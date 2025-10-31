
#include "up_sample.h"
#include "cplx_types.h"

#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif



static inline double sinc_pi(double x) {
  double ax = fabs(x);
  if (ax < 1e-12) return 1.0;
  return sin(M_PI * x) / (M_PI * x);
}

// Blackman window, symmetric for n in [-half .. +half-1]
static inline double blackman_sym(int n, int half) {
  int N = 2 * half;
  int k = n + half; // 0..N-1
  double a0 = 0.42, a1 = 0.5, a2 = 0.08;
  double phi = (2.0 * M_PI * k) / (double)(N - 1);
  return a0 - a1 * cos(phi) + a2 * cos(2.0 * phi);
}

iq_resamp_zp_t* iq_resamp_zp_init(double fin_hz, double fout_hz, int ntaps, double cutoff) {
  if (ntaps < 4 || (ntaps & 1)) return NULL; // even, >=4
  if (!(cutoff > 0.0 && cutoff <= 0.5)) return NULL;

  iq_resamp_zp_t* st = (iq_resamp_zp_t*)calloc(1, sizeof(iq_resamp_zp_t));
  if (!st) return NULL;
  st->fin = fin_hz;
  st->fout = fout_hz;
  st->dt = fin_hz / fout_hz;
  st->ntaps = ntaps;
  st->half = ntaps / 2;
  st->cutoff = cutoff;
  return st;
}

void iq_resamp_zp_free(iq_resamp_zp_t* st) {
  free(st);
}

// Reflect index for padding: maps idx outside [0, N-1] back via reflection.
static inline int reflect_index(int idx, int N) {
  // Handle small N gracefully
  if (N <= 1) return 0;
  while (idx < 0 || idx >= N) {
    if (idx < 0) idx = -idx - 1;       // reflect at -0.5
    if (idx >= N) idx = 2 * N - 1 - idx; // reflect at N-0.5
  }
  return idx;
}

// One-pass fractional-delay resampling (forward only), with reflection padding at boundaries.
// - in:  input array length Nin (offline block)
// - out: output array length Nout = floor(Nin * Fout/Fin + 0.5) or as you choose
// - st->dt = Fin/Fout
static void iq_resamp_forward(const iq_resamp_zp_t* st,
  const c32* in, size_t Nin,
  c32* out, size_t Nout)
{
  const double dt = st->dt;
  const int half = st->half;
  const int ntaps = st->ntaps;
  const double c = st->cutoff; // <= 0.5

  // Output times in input-sample units: t_k = k * dt
  // We'll compute each out[k] from samples around floor(t_k) with reflection at edges.
  for (size_t k = 0; k < Nout; ++k) {
    double t = (double)k * dt;
    int ti = (int)floor(t);
    double frac = t - (double)ti;

    double acc_i = 0.0, acc_q = 0.0;
    double wsum = 0.0;

    // taps m in [-half .. +half-1]
    for (int m = -half; m < half; ++m) {
      double x = (double)m - frac; // fractional offset
      // Lowpass, normalized to input Nyquist:
      // h(x) = (2*c) * sinc(2*c*x), then window
      double h = (2.0 * c) * sinc_pi(2.0 * c * x);
      h *= blackman_sym(m, half);

      int idx = ti + m;
      idx = reflect_index(idx, (int)Nin);

      acc_i += h * (double)in[idx].r;
      acc_q += h * (double)in[idx].i;
      wsum += h;
    }

    double inv = (fabs(wsum) > 1e-15) ? (1.0 / wsum) : 1.0;
    out[k].r = (float)(acc_i * inv);
    out[k].i = (float)(acc_q * inv);
  }
}

// Zero-phase resampling via forward-backward (filtfilt):
// 1) Forward resample Fin -> Fout to y_fwd
// 2) Reverse y_fwd, apply same kernel with dt = 1.0 (since input and output rates are equal in the reversed domain),
//    but we must still sample on the same time grid; simpler and more consistent: run the same forward routine
//    treating the reversed signal as "input at Fout" and producing "output at Fout" with dt_back = 1.0.
// 3) Reverse back to get y_zp.
//
// To keep code compact, we reuse iq_resamp_forward with a temporary state whose dt=1.0 and same taps/cutoff.
// Nout is chosen as round(Nin * Fout/Fin). You can pass any Nout you need; commonly round-to-nearest.
size_t iq_resamp_process_zero_phase(iq_resamp_zp_t* st,
  const c32* in, size_t Nin,
  c32* out, size_t Nout_cap)
{
  if (!st || !in || !out || Nin == 0) return 0;

  // Determine output length (nearest integer ratio)
  double ratio = st->fout / st->fin;
  size_t Nout = (size_t)llround((double)Nin * ratio);
  if (Nout > Nout_cap) Nout = Nout_cap;
  if (Nout == 0) return 0;

  // 1) Forward resample: x (Fin) -> y (Fout)
  c32* y = (c32*)malloc(Nout * sizeof(c32));
  if (!y) return 0;
  iq_resamp_forward(st, in, Nin, y, Nout);

  // 2) Backward (filtfilt-style) on y:
  // Reverse y -> yr
  c32* yr = (c32*)malloc(Nout * sizeof(c32));
  if (!yr) { free(y); return 0; }
  for (size_t k = 0; k < Nout; ++k) yr[k] = y[Nout - 1 - k];

  // Create a temporary state for identity-rate filtering on the reversed sequence:
  // dt_back = 1.0 (input and output both at Fout)
  iq_resamp_zp_t st_back = *st;
  st_back.fin = st->fout;
  st_back.fout = st->fout;
  st_back.dt = 1.0;

  // Run forward kernel on reversed signal with dt=1 to smooth phase (equivalent to applying same lowpass)
  c32* yr2 = (c32*)malloc(Nout * sizeof(c32));
  if (!yr2) { free(y); free(yr); return 0; }
  iq_resamp_forward(&st_back, yr, Nout, yr2, Nout);

  // Reverse back to get zero-phase result
  for (size_t k = 0; k < Nout; ++k) out[k] = yr2[Nout - 1 - k];

  free(y);
  free(yr);
  free(yr2);
  return Nout;
}

/* Example usage:

#include <stdio.h>

int main() {
    const double Fin = 10230e3;
    const double Fout = 16384e3;
    const int    NTAPS = 12;      // try 8, 12, or 16
    const double CUTOFF = 0.45;

    iq_resamp_zp_t *rs = iq_resamp_zp_init(Fin, Fout, NTAPS, CUTOFF);
    if (!rs) { fprintf(stderr, "init failed\n"); return 1; }

    size_t Nin = 10230; // 1 ms at 10.23 Msps
    cf32 *in = (cf32*)calloc(Nin, sizeof(cf32));
    // fill 'in' ...

    size_t Nout_cap = (size_t)llround(Nin * (Fout/Fin)) + 8;
    cf32 *out = (cf32*)malloc(Nout_cap * sizeof(cf32));

    size_t Nout = iq_resamp_process_zero_phase(rs, in, Nin, out, Nout_cap);
    printf("Produced %zu samples (expected ~16384)\n", Nout);

    free(in);
    free(out);
    iq_resamp_zp_free(rs);
    return 0;
}
*/

/*

// Modified Bessel I0 for Kaiser window 
static float i0f_approx(float x) {
  // Polynomial/rational approximation of I0(x)
  // Reference: Blair approximation
  float ax = fabsf(x);
  if (ax < 3.75f) {
    float y = (x / 3.75f) * (x / 3.75f);
    float r = 1.0f + y * (3.5156229f + y * (3.0899424f + y * (1.2067492f
      + y * (0.2659732f + y * (0.0360768f + y * 0.0045813f)))));
    return r;
  }
  else {
    float y = 3.75f / ax;
    float r = (expf(ax) / sqrtf(ax)) * (0.39894228f + y * (0.01328592f
      + y * (0.00225319f + y * (-0.00157565f + y * (0.00916281f
        + y * (-0.02057706f + y * (0.02635537f + y * (-0.01647633f + y * 0.00392377f))))))));
    return r;
  }
}

static float kaiser_w(float n, float N, float beta) {
  // n in [0..N], symmetric, N = taps-1
  float x = 2.0f * n / (N)-1.0f; // map to [-1,1]
  float arg = beta * sqrtf(1.0f - x * x);
  return i0f_approx(arg) / i0f_approx(beta);
}

// Build time-reversed taps for each fractional delay phase.
//   taps[phase * P + k] corresponds to sample offset k centered at 0 with fractional delay f in [0,1).
//   We store time-reversed so convolution becomes dot-product over ascending k with history arranged newest-first. 
  const int P = rs->taps_per_phase;
  const int H = rs->hlen;
  const int Nph = rs->nphases;
  float* ph = rs->phases;

  // Precompute window (same for all phases)
  float win[2 * RESAMP_HLEN + 1];
  int Nw = 2 * H;
  for (int k = -H, idx = 0; k <= H; ++k, ++idx) {
    win[idx] = kaiser_w((float)(k + H), (float)Nw, RESAMP_KAISER_BETA);
  }

  // For each phase, fractional delay f = phase / Nph
  for (int p = 0; p < Nph; ++p) {
    float f = (float)p / (float)Nph; // fractional delay in [0,1)
    // Build taps centered at t=0 with delay f: h[k] = sinc((k - f)*fc_norm) * window
    // Normalized cutoff relative to input sample rate Nyquist: RESAMP_CUTOFF in (0, 0.5]
    // Use sinc(x) = sin(pi*x*2*fc)/(pi*x) with fc = RESAMP_CUTOFF (relative to Fs/2 -> absolute Fc = RESAMP_CUTOFF * Fs/2)
    // Equivalent normalized digital frequency multiplier = 2*RESAMP_CUTOFF
    const float W = 2.0f * RESAMP_CUTOFF; // <= 1.0
    float sum = 0.0f;
    for (int k = -H, idx = 0; k <= H; ++k, ++idx) {
      float x = (float)k - f;
      float sinc;
      if (fabsf(x) < 1e-6f) {
        sinc = W;
      }
      else {
        float pix = (float)M_PI * x;
        sinc = sinf(pix * W) / (pix);
      }
      float h = sinc * win[idx];
      ph[p * P + (H - k)] = h; // time-reversed: index (H - k) maps k=-H..H to 0..P-1
      sum += h;
    }
    // Normalize DC gain to 1.0
    if (sum != 0.0f) {
      float inv = 1.0f / sum;
      for (int k = 0; k < P; ++k) {
        ph[p * P + k] *= inv;
      }
    }
  }
}

resamp_cpx_t* resamp_cpx_10230_to_16384_init(void) {
  resamp_cpx_t* rs = (resamp_cpx_t*)calloc(1, sizeof(resamp_cpx_t));
  if (!rs) return NULL;
  rs->hlen = RESAMP_HLEN;
  rs->taps_per_phase = 2 * RESAMP_HLEN + 1;
  rs->nphases = RESAMP_NPHASES;
  rs->phases = (float*)malloc(sizeof(float) * rs->nphases * rs->taps_per_phase);
  rs->hist_len = rs->taps_per_phase + 64; // small margin to simplify wrap
  rs->hist = (c32*)calloc(rs->hist_len, sizeof(c32));
  rs->hist_pos = 0;
  rs->t = 0.0;
  rs->dt = (double)RESAMP_IN_RATE / (double)RESAMP_OUT_RATE; // input samples per output sample
  if (!rs->phases || !rs->hist) {
    free(rs->phases); free(rs->hist); free(rs);
    return NULL;
  }
  build_polyphase(rs);
  return rs;
}

void resamp_cpx_free(resamp_cpx_t* rs) {
  if (!rs) return;
  free(rs->phases);
  free(rs->hist);
  free(rs);
}

// Process input block:
//   - in:  input samples (length Nin)
//   - out: output buffer (capacity at least ceil(Nin * Fout/Fin) + a few)
//   - returns number of output samples written
//   You can call repeatedly; state persists. 
size_t resamp_cpx_10230_to_16384(resamp_cpx_t* rs, const c32* in, size_t Nin, c32* out, size_t Nout_max)
{
  if (!rs) return 0;
  size_t out_count = 0;

  const int P = rs->taps_per_phase;
  const int H = rs->hlen;
  const int Nph = rs->nphases;
  float* ph = rs->phases;

  // For efficient dot-product, we maintain a circular history buffer where hist_pos is the next write index.
  for (size_t n = 0; n < Nin; ++n) {
    // Push input sample
    rs->hist[rs->hist_pos] = in[n];
    rs->hist_pos = (rs->hist_pos + 1) % rs->hist_len;

    // After inserting each input, emit as many outputs as fall before the next input time.
    // We track fractional time t measured in input-sample units relative to current write position.
    // We align convolution so that the center tap corresponds to time t = 0 at the latest written sample boundary.
    // When rs->t < 1.0, an output is due using fractional delay f = rs->t.
    while (rs->t < 1.0) {
      if (out_count >= Nout_max) return out_count;

      // Phase index
      double f = rs->t; // [0,1)
      int p = (int)floor(f * Nph);
      if (p >= Nph) p = Nph - 1;

      // Gather dot-product over P taps centered around the newest sample.
      // We need samples at offsets k = -H..+H from the "center" position (one before hist_pos).
      int center = (rs->hist_pos - 1 + rs->hist_len) % rs->hist_len;

      float acc_i = 0.0f, acc_q = 0.0f;
      const float* taps = &ph[p * P];

      // Walk history backward for k = -H..+H using time-reversed taps
      int idx = center + H; // start at center + H, move backward with wrap
      for (int k = 0; k < P; ++k) {
        if (--idx < 0) idx += rs->hist_len; // pre-decrement to map k=0 to center+H-1
        const c32 s = rs->hist[idx];
        float h = taps[k];
        acc_i += h * s.r;
        acc_q += h * s.i;
      }

      out[out_count].r = acc_i;
      out[out_count].i = acc_q;
      out_count++;

      rs->t += rs->dt; // advance output time by input-sample units per output
    }

    // Consume one input-sample time
    rs->t -= 1.0;
    if (rs->t < 0.0) rs->t = 0.0; // numeric safety
  }

  return out_count;
}

// Helper to estimate required output capacity for a block of Nin inputs 
static size_t resamp_out_capacity(size_t Nin) {
  double ratio = RESAMP_OUT_RATE / RESAMP_IN_RATE; // ≈ 1.601955
  return (size_t)ceil(Nin * ratio) + 8;
}

*/