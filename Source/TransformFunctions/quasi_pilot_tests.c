#include <math.h>
#include <stdint.h>
#include <stddef.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#include "gnss_codes.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define R2D (180.0 / M_PI)
#define D2R (M_PI / 180.0)

#define NUM_CHIPS 1023
#define SPC 4                     // samples per chip
//#define NSAMP (NUM_CHIPS * SPC)   // total samples in a code period
enum { NSAMP = (NUM_CHIPS * SPC) };
#define CHIP_RATE 1.023e6                  // Hz
#define FS (CHIP_RATE * SPC)               // sample rate (Hz)

// -------------------------------
// External replicas (provide these)
// E1B_Code[PRN][4092] in {+1,-1}, PRN index 1..50
extern int8_t E1B_Code[51][4092];

// -------------------------------
// Utility complex helpers
//typedef struct { float r, i; } c32;
typedef struct { double r, i; } c64;

//typedef struct { float re, im; } c32;
static inline c32 c_add(c32 a, c32 b){ return (c32){a.r+b.r, a.i+b.i}; }
static inline c32 c_mul(c32 a, c32 b){ return (c32){a.r*b.r - a.i*b.i, a.r*b.i + a.i*b.r}; }
static inline c32 c_conj(c32 a){ return (c32){a.r, -a.i}; }
static inline float c_abs2(c32 a){ return a.r*a.r + a.i*a.i; }
static c32 c_add64(c32 a, c32 b) { return (c32) { a.r + b.r, a.i + b.i }; }
static c32 c_scale(c32 a, float k) { return (c32) { a.r* k, a.i* k }; }
extern c32 c_mul_conj(c32 a, c32 b) { // a * conj(b)
  return (c32) { a.r* b.r + a.i * b.i, a.i* b.r - a.r * b.i };
}


// Fast sincos
static inline void sincosf_fast(float x, float *s, float *c){
    *s = sinf(x); *c = cosf(x);
}

// -------------------------------
// E1-B CBOC composite subcarrier sign for a fractional chip position u in [0,1).
// We model the subcarrier signs at the E1-B CBOC(6,1, 1/11) composite:
// v = sqrt(10/11)*BOC(1,1)  - sqrt(1/11)*BOC(6,1)
// For an efficient replica and mixer, we only need the sign at the current fractional chip.
// Here we use square-wave BOC signs with MSK-like chip-centered transitions.

static inline int boc_sign(int m, float u){
    if (m <= 0) return +1; // BPSK
    int k = (int)floorf(2.0f * m * u); // toggles 2*m times per chip
    return (k & 1) ? -1 : +1;
}

static inline float cboc_e1b_weight(float u){
    const float w1 = 0.9534625892455923f; // sqrt(10/11)
    const float w6 = 0.3015113445777636f; // sqrt(1/11)
    int s11 = boc_sign(1, u);
    int s61 = boc_sign(6, u);
    // E1-B uses minus sign for the 6,1 branch
    return w1 * (float)s11 - w6 * (float)s61;
}

// -------------------------------
// Code chipM_PIcker: returns +/-1 code value for PRN at fractional chip phase (wrapped).
// code_phase_chips in [0, L), L=4092
static inline int code_chip_at(const int8_t *code, float code_phase_chips, int L){
    int idx = (int)floorf(code_phase_chips);
    if (idx >= L) idx -= L;
    else if (idx < 0) idx += L;
    return code[idx];
}


// Per-PRN Doppler mixer:
// - Each PRN has its own Doppler (Hz) and initial carrier phase (deg).
// - The PRN chip value (+1/-1) scales that PRN's complex phasor.
// - Outputs NSAMP complex samples.
extern void mix_two_prns_oversampled_per_prn(const int32_t* prn_a,
  const int32_t* prn_b,
  double doppler_a_hz, double doppler_b_hz,
  double phase_a_deg, double phase_b_deg,
  c32* out_iandq, int size, float samp_rate, float sigma, int sign)
{
  // PRN A phasor increment and initial phase
  const double dth_a = 2.0 *M_PI * (double)doppler_a_hz / (double)samp_rate;
  const double ca_inc = (double)cos(dth_a);
  const double sa_inc = (double)sin(dth_a);
  double pca = (double)cos(phase_a_deg * (double)M_PI / 180.0f);
  double psa = (double)sin(phase_a_deg * (double)M_PI / 180.0f);

  // PRN B phasor increment and initial phase
  const double dth_b = 2.0 *M_PI * (double)doppler_b_hz / (double)samp_rate;
  const double cb_inc = (double)cos(dth_b);
  const double sb_inc = (double)sin(dth_b);
  double pcb = (double)cos(phase_b_deg * (double)M_PI / 180.0f);
  double psb = (double)sin(phase_b_deg * (double)M_PI / 180.0f);

  for (int samp = 0; samp < size; ++samp) {
    double a = (double)prn_a[samp];  // +/- 1
    double b = (double)prn_b[samp];  // +/- 1

    // s[n] = a*e^{j theta_a[n]} + b*e^{j theta_b[n]}
    double ia = a * pca, qa = a * psa;
    double ib = b * pcb, qb = b * psb;
    int quant = 0;
    out_iandq[samp].r = quant ? quantize_pm13(ia + ib + noise(sigma))* sign : ((ia + ib) + noise(sigma))* sign;
    out_iandq[samp].i = quant ? quantize_pm13(qa + qb + noise(sigma))* sign : ((qa + qb) + noise(sigma))* sign;

    // advance both phasors
    double npca = pca * ca_inc - psa * sa_inc;
    double npsa = pca * sa_inc + psa * ca_inc;
    pca = npca; psa = npsa;

    double npcb = pcb * cb_inc - psb * sb_inc;
    double npsb = pcb * sb_inc + psb * cb_inc;
    pcb = npcb; psb = npsb;

    // renormalize occasionally to limit float drift (not necessary for 1 ms)
    double na = 1.0 / sqrt(pca * pca + psa * psa);
    pca *= na; psa *= na;
    double nb = 1.0 / sqrt(pcb * pcb + psb * psb);
    pcb *= nb; psb *= nb;

  }

}

extern void rotate_fwd(int array[], int size, int offset) {
  int* temp = (int*)malloc(sizeof(int) * size);
  if (temp == NULL) { printf("rotate_fwd's malloc failed");  return; }
  for (int i = 0; i < size; i++) {
    if (i - offset < 0) {
      temp[i] = array[size + (i - offset)];
    }
    else {
      temp[i] = array[i - offset];
    }
  }
  memcpy(array, temp, sizeof(int) * size);
  free(temp);
}

extern void make_replica(const int32_t* prn_a, c32* out_iandq,float doppler, int size, float samp_freq)
{
  const double phase_deg = 0.0001;
  // PRN phasor increment and initial phase
  const double dth_a = 2.0 * M_PI * (double)doppler / (double)samp_freq;
  const double ca_inc = (double)cos(dth_a);
  const double sa_inc = (double)sin(dth_a);
  double pca = (double)cos(phase_deg * (double)M_PI / 180.0f);
  double psa = (double)sin(phase_deg * (double)M_PI / 180.0f);

  for (int samp = 0; samp < size; ++samp) {
    double a = (double)prn_a[samp];  // +/- 1
    
    // s[n] = a*e^{j theta_a[n]} + b*e^{j theta_b[n]}
    double ia = a * pca, qa = a * psa;
    bool quant = false;
    out_iandq[samp].r = quant ? (float) (quantize_pm13(ia ) + noise(1)) : (float)(ia );
    out_iandq[samp].i = quant ? (float) (quantize_pm13(qa ) + noise(1)) : (float)(qa );

    // advance both phasors
    double npca = pca * ca_inc - psa * sa_inc;
    double npsa = pca * sa_inc + psa * ca_inc;
    pca = npca; psa = npsa;

    // renormalize occasionally to limit float drift (not necessary for 1 ms)
    double na = 1.0 / sqrt(pca * pca + psa * psa);
    pca *= na; psa *= na;
  }
}

// -------------------------------
// Signal synthesis for two PRNs (PRN2 and PRN4) into one complex baseband stream.
// Inputs:
//   fs_hz            : sample rate
//   N                : number of samples
//   fIF_hz           : baseband center (usually 0)
//   prn2/prn4        : PRN numbers (2 and 4)
//   doppler2/4       : carrier Doppler (Hz) per PRN
//   phi2/4_rad       : initial carrier phase (rad) at sample 0
//   code_rate_cps    : nominal code chip rate (1.023e6 for E1), set per-PRN Doppler via scaling if desired
//   code_doppler_scale: optional fractional change of code rate due to Doppler (e.g. 1 + fd/c * f_code/f_carr). Use 1.0 if unknown.
//   code_phase2/4    : initial code phase (chips) in [0,4092)
//   out              : output buffer of c32 length N
// Notes:
//   - E1-B data bit is assumed +1 for this short synthesis. If needed, pass a per-sample sign array and multiply.
extern void synth_e1b_two_prns(
    float fs_hz, size_t N, float fIF_hz,
    int prn_a, float doppler_a, float phia_rad, float code_rate_cps_a, float code_doppler_scalea, float code_phase_a,
    int prn_b, float doppler_b, float phib_rad, float code_rate_cps_b, float code_doppler_scaleb, float code_phase_b,
    c32 *out
){
    const int L = 4092;
    const int8_t *ca = E1B_Code[prn_a];
    const int8_t *cb = E1B_Code[prn_b];

    float dchipsa = (code_rate_cps_a * code_doppler_scalea) / fs_hz;
    float dchipsb = (code_rate_cps_b * code_doppler_scaleb) / fs_hz;

    float dphia = 2.0f * (float)M_PI * (fIF_hz + doppler_a) / fs_hz; // phase increment per sample
    float dphib = 2.0f * (float)M_PI * (fIF_hz + doppler_b) / fs_hz;

    float chipsa = code_phase_a;
    float chipsb = code_phase_b;
    float phia = phia_rad;
    float phib = phib_rad;

    for (size_t n = 0; n < N; ++n) {
        // PRN a
        float fraca = chipsa - floorf(chipsa);
        int code_a = code_chip_at(ca, chipsa, L);
        float sca = cboc_e1b_weight(fraca);
        float ampa = (float)code_a * sca; // data assumed +1

        float sa, caph; sincosf_fast(phia, &sa, &caph);
        c32 xa = { ampa * caph, ampa * sa };

        // PRN b
        float fracb = chipsb - floorf(chipsb);
        int codeb = code_chip_at(cb, chipsb, L);
        float scb = cboc_e1b_weight(fracb);
        float ampb = (float)codeb * scb;

        float sb, cbph; sincosf_fast(phib, &sb, &cbph);
        c32 xb = { ampb * cbph, ampb * sb };

        out[n].r = quantize_pm13(xa.r + xb.r + noise(4.01)); // x2.r + x4.r;
        out[n].i = quantize_pm13(xa.i + xb.i + noise(4.01)); // x2.i + x4.i;

        // advance
        chipsa += dchipsa; 
        if (chipsa >= L) { chipsa -= L; } 
        else if (chipsa < 0) { chipsa += L; }

        chipsb += dchipsb; 
        
        if (chipsb >= L) { chipsb -= L; } 
        else if (chipsb < 0) { chipsb += L;}

        phia += dphia; 
        if (phia > 1e9f || phia < -1e9f) { phia = fmodf(phia, 2.0f * (float)M_PI); }
        phib += dphib; 
        if (phib > 1e9f || phib < -1e9f) { phib = fmodf(phib, 2.0f * (float)M_PI); }
    }
}

// -------------------------------
// Acquisition/tracking aid for PRN 2: estimate code phase (chips), carrier freq (Hz), and carrier phase (rad).
// Strategy:
//   1) Code-phase search (coarse): slide the prompt correlator over a code-phase grid (step ~ 0.5 chips or finer).
//      For each hypothesis, wipe code*subcarrier and sum complex samples over a short coherent T (e.g., 4 ms), but without carrier wipe.
//      Then at each hypothesis, we do a small frequency bin search around f_search_center ± f_search_span to find best |sum|.
//   2)M_PIck the max across code-phase and frequency. Report:
//        - code_phase_chips (chips in [0,4092)),
//        - carrier_freq_hz (offset from fIF, i.e., Doppler estimate),
//        - carrier_phase_rad at midpoint of coherent span.
// Notes:
//   - This is a compact reference, not a production-optimized acquisition.
//   - You can refine with finer code steps, more ms, and 2D FFT approaches (circular correlation + FFT over freq).
//
// Inputs:
//   buf[]: complex baseband samples
//   fs_hz: sampling rate
//   fIF_hz: baseband center (Hz) (usually 0)
//   T_coh_ms: coherent integration time in milliseconds (e.g., 4)
//   code: PRN 2 code array (±1, length 4092)
//   code_rate_cps: nominal code rate (1.023e6)
//   code_guess guess within +/- 10 chips of try code
//   code_step_chips: step for code-phase grid (e.g., 0.5 or 0.25)
//   f_search_center_hz: expected Doppler center (Hz) relative to fIF
//   f_search_span_hz: half-span (search from center-span to center+span)
//   f_search_step_hz: frequency bin step (e.g., 250 Hz)
// Outputs:
//   out_code_phase_chips
//   out_carrier_freq_hz (relative to fIF, i.e., Doppler)
//   out_carrier_phase_rad (referenced to midpoint of coherent span)
extern void estimate_prn_code_and_carrier(
    const c32 *buf, size_t N, float fs_hz, float c_search_span,
    int8_t const *code, float code_rate_cps, float code_guess,
    float T_coh_ms, float code_step_chips,
    float f_search_center_hz, float f_search_span_hz, float f_search_step_hz,
    float *out_code_phase_chips, float *out_carrier_freq_hz, float *out_carrier_phase_rad, c32* c32_sum
){
    const int L = 4092;
    // number of samples for coherent sum
    size_t Ncoh = (size_t) llroundf((T_coh_ms * 1e-3f) * fs_hz);
    if (Ncoh > N) Ncoh = N;
    if (Ncoh < 100) { // too short
        *out_code_phase_chips = 0.0f;
        *out_carrier_freq_hz = 0.0f;
        *out_carrier_phase_rad = 0.0f;
        return;
    }

    float dchips = code_rate_cps / fs_hz;
    float best_re = 0.0, best_im = 0.0;
    float best_pow = -1.0f;
    float best_code = 0.0f;
    float best_f = 0.0f;
    c32 best_sum = {0,0};

    float code_low  = code_guess - c_search_span;
    float code_high = code_guess + c_search_span;
    uint8_t wrap = 0;
    if (code_low < 0) { code_low += L; wrap = 1; }
    if (code_high > L) { code_high -= L;  wrap = 1; }
    //printf("wrap %d low %f high %f \n", wrap, code_low, code_high);
    
    // Code-phase grid
    for (float code0 = 0.0f; code0 < (float)L; code0 += code_step_chips) {
        // Precompute code*subcarrier wave over Ncoh samples for this starting phase
        // We only need a per-sample real weight w[n] = code_chip * cboc_weight(frac) (± amplitude)
        // Then we multiply buf[n] by w[n] (code wipe) and test frequency bins.
        if (wrap == 0 && (code0 < code_low || code0 > code_high)) { continue; }
        if (wrap == 1 && (code0 > code_high && code0 < code_low)) { continue; }

        // To keep it memory-light, generate on the fly inside freq loop.
        for (float f = f_search_center_hz - f_search_span_hz; f <= f_search_center_hz + f_search_span_hz + 1e-3f; f += f_search_step_hz) {
            float dphi = 2.0f * (float)M_PI * f / fs_hz;
            float phi = 0.0f; // phase start arbitrary; we only care about sum angle and magnitude

            float chips = code0;
            c32 S = {0,0};
            for (size_t n = 0; n < Ncoh; ++n) {
                float frac = chips - floorf(chips);
                int c = code_chip_at(code, chips, L);
                float w = (float)c * cboc_e1b_weight(frac); // code*subcarrier

                // code wipeoff: multiply received by w (real). Then rotate by -carrier to accumulate at DC.
                float sph, cph; sincosf_fast(-phi, &sph, &cph); // rotator for carrier wipe-off
                c32 x = buf[n];
                // Multiply by w (real) first
                x.r *= w; x.i *= w;
                // Rotate by -phi
                float xr = x.r*cph - x.i*sph;
                float xi = x.r*sph + x.i*cph;
                S.r += xr; S.i += xi;

                chips += dchips; 
                if (chips >= L) { chips -= L; }
                phi += dphi;
            }
            float pow = c_abs2(S);
            if (pow > best_pow) {
                best_pow = pow; best_code = code0; best_f = f; best_sum = S;
                best_re = S.r; best_im = S.i;
            }
        }
    }

    c32_sum->r = (double) best_re;
    c32_sum->i = (double) best_im;
    // Estimate carrier phase at coherent midpoint (t_mid ~ Ncoh/(2 fs))
    float phase = atan2f(best_sum.i, best_sum.r);
    *out_code_phase_chips = fmodf(best_code, (float)L);
    if (*out_code_phase_chips < 0) { *out_code_phase_chips += (float)L; }
    *out_carrier_freq_hz = best_f;
    *out_carrier_phase_rad = phase;
}


// Recover PRN A’s phase with per-PRN Doppler: wipe PRN A Doppler, then code-correlate.
// Returns residual phase (deg) and best 1/4-chip offset.
extern double recover_gal_prn_phase_deg_with_doppler(const c32* iandq,
  int SIZE,
  const int prn, 
  int  code_offset,
  double doppler_a_hz,
  int* best_offset_out,
  double* doppler_out,
  double* pwr,
  double* pwr_ratio,
  c32* c32_sum)
{
  double best_pow = -1.0, second_best_pow = -1, best_re = 0.0, best_im = 0.0, best_doppler = 0.0;
  int best_off = 0, second_off = -100;

  float code_search_span = 0; //chip span was 2
  float code_low = code_offset - code_search_span;
  float code_high = code_offset + code_search_span;
  uint8_t wrap = 0;
  if (code_low < 0) { code_low += SIZE; wrap = 1; }
  if (code_high > SIZE) { code_high -= SIZE;  wrap = 1; }

  float f_search_step_hz = 1.0; // Hz
  float f_search_span_hz = 0; // Hz span was 2
  for (float f = -f_search_span_hz; f <= +f_search_span_hz + 1e-3f; f += f_search_step_hz) {
    // Frequency wipe for PRN A: multiply by e^{-j 2π f_a n / fs}
    double dtheta = 2.0 * M_PI * (doppler_a_hz + f) / (double)FS;

    const float c_inc = (float)cos(dtheta);
    const float s_inc = (float)sin(dtheta);

    float pc = 1.0f, ps = 0.0f;  // start at zero phase; PRN A's initial phase stays in residual
    static float iw[1023 * 4], qw[1023 * 4]; // put replica with Doppler; code will be handled in next loop
    for (int n = 0; n < SIZE; ++n) {
      float I = iandq[n].r, Q = iandq[n].i;
      // (I + jQ) * e^{-jθ} = (I*pc + Q*ps) + j(-I*ps + Q*pc)
      iw[n] = +I * pc + Q * ps;
      qw[n] = -I * ps + Q * pc;
      // rotate by dtheta
      float npc = pc * c_inc - ps * s_inc;
      float nps = pc * s_inc + ps * c_inc;
      // cache for next round
      pc = npc;
      ps = nps;
    }

    // 1/4-chip code-phase search across one code period
    for (int off = 0; off < SIZE; ++off) {
      if (wrap == 0 && (off < code_low || off > code_high)) { continue; } // handle wrap around of code search limits
      if (wrap == 1 && (off > code_high && off < code_low)) { continue; }
      double re = 0.0, im = 0.0;
      for (int n = 0; n < SIZE; ++n) {
        int idx = n - off;
        if (idx < 0) { idx += SIZE; } // handle code rollovers
        re += iw[n] * E1B_Code[prn][idx];
        im += qw[n] * E1B_Code[prn][idx];
      }
      double p = sqrt(re * re + im * im);
      if (p > best_pow) {
        best_pow = p;
        best_doppler = (double)f;
        best_re = re;
        best_im = im;
        best_off = off;
      }
      //printf("%d, %.1f\n",off,  p);
      if (p > second_best_pow && fabs(best_off - off) > 20 && p < best_pow) {
        second_best_pow = p;
        second_off = off;
      }
    } // for off 
  } // end for f_search_span_hz
  *doppler_out = best_doppler;
  *best_offset_out = best_off;
  *pwr = best_pow;
  *pwr_ratio = (best_pow / second_best_pow);
  c32_sum->r = best_re;
  c32_sum->i = best_im;
  //printf(" best %.3f sec %.3f %d %d %d \n", best_pow, second_best_pow, best_off, second_off, (int)fabs(best_off - second_off));
  double phase_deg = atan2(best_im, best_re) * R2D;
  return phase_deg;
}

/* 


Synthesizes complex baseband I/Q for two Galileo E1‑B signals (PRN 2 and PRN 4).

Each PRN has its own carrier Doppler (Hz) and initial carrier phase (rad).
Uses provided code replicas: arrays of ±1 chips for each PRN, shaped by the E1 CBOC combination in a simple, correlator-friendly way.
Lets you set code rate (with Doppler), sampling rate, and initial code phases (in chips).
Provides a function to estimate PRN 2’s code phase and carrier phase from a buffer:

Performs code wipeoff using the known replica for PRN 2 over a search of code phase (via early/prompt/late or simple grid).
Performs a short coherent carrier frequency/phase search (1D FFT-like single-bin grid) around a provided search band to recover residual Doppler and phase at the coherent midpoint.
Returns estimated code phase (chips), carrier frequency (Hz), and carrier phase (rad).
Assumptions

You already have the E1‑B primary code replicas in ±1 (int8_t) arrays for PRN 2 and PRN 4 (length 4092 chips), e.g., e1b_code[prn][4092]. You can load them using the tiny loader we built earlier.
For E1‑B, modulation is CBOC on E1 with weights sqrt(10/11) for BOC(1,1) and −sqrt(1/11) for BOC(6,1). For correlation, many receivers use separate BOC(1,1) and BOC(6,1) branches and combine. Here, to keep it compact, we generate a single composite time-domain replica that matches the ICD’s E1‑B sign (negative on the 6,1 branch).
Navigation data on E1‑B is 250 sps; for short windows we assume the data bit is constant (or known). You can pass a per-sample data sign if desired; here we default to +1.


Notes and practical tips

Data bits: For E1‑B, the 250 sps data can flip sign; for short coherents (≤4 ms) you’re within a single symbol and safe. For longer coherents across unknown bit edges, use either noncoherent combining or differential techniques, or wipe known bits if you demodulate I/NAV/SSP first.

Code Doppler: The function accepts a code_doppler_scale to slightly stretch the code rate based on your Doppler hypothesis. A simple first-order model is scale ≈ 1 + fd/c · f_code/f_carr, or just keep 1.0 and let the correlator tolerate small misalignment over a few ms.

Search performance: This brute-force loop is clear but not optimal. For speed, use:

FFT-based circular correlation over code phase for a bank of frequency bins.
Separate BOC(1,1) and BOC(6,1) correlators and combine outputs by weights sqrt(10/11) and −sqrt(1/11).
Coherent combine multiple 4 ms blocks after compensating the best Doppler bin.
Phase reference: The returned carrier phase is at the coherent midpoint. If you track over time, account for phase evolution with the estimated Doppler to maintain continuity.

Can you provide with an exmaple of how to implement early/prompt/late correlators?
How would I implement FFT-based circular correlation for faster aquisition?
Can you show me how to combine multiple 4 ms blocks non-coherently?

*/

/*
extern void test_galileo(float true_phase_deg, float true_code, float true_doppler) {
  //static int8_t e1b[E1B_MAX_PRN + 1][E1B_CODE_LEN];
  int num_prn = load_e1b_primary_codes("C:/work/Baseband/HEX_E1B.txt", E1B_Code);
  if (num_prn < 0) { fprintf(stderr, "Failed to open file\n"); return 1; }
  fprintf(stdout, "Loaded %2d PRNs\n", num_prn);

  float fs = 4.092e6f; // sampling rate
  size_t N = (size_t)(0.05f * fs); // 5 ms buffer
  float fIF = 0.0f;
  enum EPOCHS { epoc1 = 60 };
  enum BITS { bits1 = (epoc1 / 20)};
  // Allocate buffer
  c32* buf = (c32*)malloc(N * sizeof(c32));
  int signum = -1;
  c32 segment[bits1] = { 0 }; // segment accumulator
  c32 sum_iq[epoc1] = { 0 };
  int sgn[bits1] = { 1 };
  int seg_cnt = 0;
  c32 c32_sum = { 0, 0 };
  int PRN_A = 4;
  int PRN_B = 34;
  for (int i = 0; i < epoc1; i++) {

    // Synthesize two PRNs
    synth_e1b_two_prns(
      fs, N, fIF,
      PRN_A, true_doppler, true_phase_deg * D2R, 1.023e6f, 1.0f, true_code,   // PRN to be searched
      PRN_B, 4444, -55 * D2R, 1.023e6f, 1.0f, 2222,  // other PRN 
      buf
    );

    for (int k = 0; k < N; k++) { buf[k].i *= signum; buf[k].r *= signum; }

    //printf("prn mix done \n");
    // Estimate PRN 2 code and carrier
    float code_phase, fd, ph, span_hz = 5.0f, span_chips = 10.0f;
    estimate_prn_code_and_carrier(
      buf, N, fs, span_chips,
      E1B_Code[PRN_A], 1.023e6f, true_code,
      4.0f, 1.0f,             // 4 ms coherent, 0.5-chip code step
      true_doppler, span_hz, 1.0f,  // Doppler search ±3 kHz, 0.5 Hz step
      &code_phase, &fd, &ph, &c32_sum);

    sum_iq[i] =  c32_sum;
    segment[seg_cnt] = c_add64(segment[seg_cnt], c32_sum);
    if (i % 20 == 0) {
      signum *= -1;
      segment[seg_cnt++] = c32_sum;
      c32_sum.i = 0; c32_sum.r = 0;
    }
    // get the sums of I and Q to check the bit flips
    printf("PRN%d[%2d]: phase=%7.3f [deg], code_phase=%7.3f [chips], Doppler=%7.1f [Hz] sign %d\n",
      PRN_A, i, ph * R2D, code_phase, fd, signum);
  }
  // work on revering correct phases:
  for (int k = 0; k < bits1; ++k) { sgn[k] = +1; }
  for (int k = 1; k < bits1; ++k) {
    c32 ref = segment[0];
    c32 prod_same = c_mul_conj(segment[k], ref);
    c32 prod_flip = c_mul_conj(c_scale(segment[k], -1.0f), ref);
    float score_same = prod_same.r;
    float score_flip = prod_flip.r;
    //printf("score_same %f score_flip %f \n", score_same, score_flip);
    sgn[k] = (score_flip > score_same) ? -1 : +1;
  }
  for (int i = 0; i < bits1; i++) {
    printf("Segment %d, sgn %d ang %f \n", i, sgn[i], atan2(segment[i].i, segment[i].r) * 180 /M_PI);
  }

  c32_sum.i = c32_sum.r = 0;
  for (int i = 0; i < epoc1; i++) {
    int idx = i / 20;
    sum_iq[i].r *= sgn[idx];
    sum_iq[i].i *= sgn[idx];
    c32_sum = c_add64(sum_iq[i], c32_sum);
  }
  printf("Average ang %f \n",  atan2(c32_sum.i, c32_sum.r) * 180 /M_PI);

  free(buf);
}*/

