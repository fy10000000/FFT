/** 
 * CMSIS.cpp : This file contains the 'main' function for testing modifications made
 * to CMSIS to handle up to 16K FFT size for the CFFT functions.
 * 
 * Note that only float and q15_t types are supported and tested here
 * For some strange reason the q31_t types do not compile
 * 
 * This derivative work is proprietary and confidential to Baseband Technologies, Inc.
 * 
 */

#include <iostream>
#include <chrono> // for profiling
#include "transform_functions.h"
#include "gnss_codes.h"
#include "up_sample.h"
#include "msb_funcs.h"

#include <complex> // not to be confused with complex.h


#define Q31_MIN   ((int32_t)0x80000000)  // -2147483648
#define Q31_MAX   ((int32_t)0x7FFFFFFF)  //  2147483647
#define Q15_MIN   (-32768)
#define Q15_MAX   ( 32767)
#define R2D (180.0f / 3.14159265358979f)
#define D2R (3.14159265358979f / 180.0f)
#define FS (4.092e6f) // sample rate

typedef struct {
  int num_errors;
  int num_trials;
} results_s;

typedef struct {
  float  num;
  float  data[20];
  int    idx;
  int    window;
} stat_s;

void stat_init(stat_s* s) {
  s->num = 0.0f;
  s->idx = 0;
  s->window = 3;
  memset(s->data, 0, sizeof(s->data));
} 

void stat_add(stat_s* s, float x) {
  if (s->num < s->window) {
    s->num++;
  }
  s->data[s->idx] = x;
  s->idx++;
  s->idx %= s->window;
}

float stat_var(stat_s* s) {
  if (s->num < 3) { return 0.0f; }
  float sum = 0.0f, sum_sq = 0.0f;
  for (int i = 0; i < s->num; i++) {
    sum += s->data[i];
    sum_sq += s->data[i] * s->data[i];
  }
  return ((sum_sq - sum * sum / s->num) / (s->num - 1));
}

float stat_mean(stat_s* s) {
  float sum = 0.0f;
  for (int i = 0; i < s->num; i++) {
    sum += s->data[i];
  }
  return  sum / s->num;
}

static inline int16_t q15_mul(int16_t a, int16_t b) {
  // Promote to 32-bit for the product
  int32_t prod = (int32_t)a * (int32_t)b; // range: [0x4000_0000, 0x4000_0000]
  // Round to nearest: add 0.5 ulp before shifting.
  // For signed values, use "round away from zero" approach:
  // add 0x4000 for positive, subtract 0x4000 for negative.
  int32_t rounded = prod >= 0 ? (prod + (1 << 14)) : (prod - (1 << 14));
  int32_t q15 = rounded >> 15;
  // Saturate to Q15 range
  if (q15 > Q15_MAX) { 
    return Q15_MAX; 
  }
  if (q15 < Q15_MIN) {
    return Q15_MIN;
  }
  return (int16_t)q15;
}

static inline int16_t q15_add(int16_t a, int16_t b) {
  int32_t sum = (int32_t)a + (int32_t)b;
  if (sum > Q15_MAX) {
    return Q15_MAX;
  }
  if (sum < Q15_MIN) {
    return Q15_MIN;
  }
  return (int16_t)sum;
}

static inline int16_t q15_sub(int16_t a, int16_t b) {
  int32_t diff = (int32_t)a - (int32_t)b;
  if (diff > Q15_MAX) {
    return Q15_MAX;
  }
  if (diff < Q15_MIN) {
    return Q15_MIN;
  }
  return (int16_t)diff;
}

static inline int32_t q31_add(int32_t a, int32_t b) {
  int64_t sum = (int64_t)a + (int64_t)b;
  if (sum > Q31_MAX) {
    return Q31_MAX;
  }
  if (sum < Q31_MIN) {
    return Q31_MIN;
  }
  return (int32_t)sum;
}

static inline int32_t q31_sub(int32_t a, int32_t b) {
  int64_t diff = (int64_t)a - (int64_t)b;
  if (diff > Q31_MAX) {
    return Q31_MAX;
  }
  if (diff < Q31_MIN) {
    return Q31_MIN;
  }
  return (int32_t)diff;
}

static inline int32_t q31_mul(int32_t a, int32_t b) {
  // 64-bit product in Q62
  int64_t prod = (int64_t)a * (int64_t)b;
  // Round to nearest, symmetric: add/subtract 0.5 LSB before shifting.
  // 0.5 LSB at Q31 corresponds to 1<<30 in Q62 domain.
  int64_t rounded = prod >= 0 ? (prod + (int64_t)1 << 30) : (prod - ((int64_t)1 << 30));
  // Shift back to Q31
  int64_t q31 = rounded >> 31;
  // Saturate
  if (q31 > Q31_MAX) { 
    return Q31_MAX; 
  }
  if (q31 < Q31_MIN) { 
    return Q31_MIN; 
  }
  return (int32_t)q31;
}

static inline float q31tof(int32_t x) { return (float)x / 2147483648.0f; } // 2^31
static inline int32_t f2q31(float x) {
  if (x >= 0.9999999995343387f) {
    return Q31_MAX; // (Q31_MAX / 2^31)
  }
  if (x <= -1.0f) {
    return Q31_MIN;
  }
  double v = (double)x * 2147483648.0; // 2^31
  if (v > Q31_MAX) {
    v = Q31_MAX;
  }
  if (v < (double)Q31_MIN) {
    v = (double)Q31_MIN;
  }
  return (int32_t)v;
}


static uint16_t U2(uint8_t* p)
{
  uint16_t value = 0;
  for (int i = 0; i < 2; i++) {
    value = (value << 8) | (p[1 - i]);
  }
  return value;
}

static uint32_t U3(uint8_t* p)
{
  uint32_t value = 0;
  for (int i = 0; i < 3; i++) {
    value = (value << 8) | (p[2 - i]);
  }
  return value;
}

static uint32_t U4(uint8_t* p)
{
  uint32_t value = 0;
  for (int i = 0; i < 4; i++) {
    value = (value << 8) | (p[3 - i]);
  }
  return value;
}

void fft_cpp(int size, std::complex<float>* w, bool fwd) {
  // interleaved real, imagmalloc(size * sizeof(float);
  float* fft_data = (float*)malloc(sizeof(float) * 2 * size); 
  for (int i = 0; i < size; i++) {
    fft_data[2 * i] = w[i].real();
    fft_data[2 * i + 1] = w[i].imag();
  }
  arm_cfft_radix2_instance_f32 s;
  arm_cfft_radix2_init_f32(&s, size, fwd ? 0:1, 1);
  arm_cfft_radix2_f32(&s, fft_data);
  for (int i = 0; i < size; i++) {
    w[i] = std::complex<float>(fft_data[2 * i], fft_data[2 * i + 1]);
  }
  free(fft_data);
}

void fft_c32(int size, c32* w, bool fwd) {
  // interleaved real, imag
  float* fft_data = (float*)malloc(sizeof(float) * 2 * size);
  memset(fft_data, 0, sizeof(float) * 2 * size);
  for (int i = 0; i < size; i++) {
    fft_data[2 * i] = w[i].r; // Use the real part
    fft_data[2 * i + 1] = w[i].i; // Use the imaginary part
  }
  arm_cfft_radix2_instance_f32 s;
  arm_cfft_radix2_init_f32(&s, size, fwd ? 0 : 1, 1);
  arm_cfft_radix2_f32(&s, fft_data);
  for (int i = 0; i < size; i++) {
    w[i].r = fft_data[2 * i]; // Update the real part
    w[i].i = fft_data[2 * i + 1]; // Update the imaginary part
  }
  free(fft_data);
}

c32 get_conj(const c32 x) {
  c32 y;
  y.r = x.r;
  y.i = -x.i;
  return y;
}

c32 mult(const c32 a, const c32 b) {
  return c32{  //(a.r +j a.i) * (b.r +j b.i)
    a.r * b.r - a.i * b.i,
    a.r * b.i + a.i * b.r
  };
}

c32 add(const c32 a, const c32 b) {
  return c32{  //(a.r +j a.i) + (b.r +j b.i)
    a.r + b.r,
    a.i + b.i
  };
}

float mag(const c32 in) {
  return sqrtf(in.r * in.r + in.i * in.i);
}


double compute_gps_time(int year, int month, int day, int hour, int minute, int second)
{
  struct tm refTime, epochTime;
  time_t refTimeRaw, epochTimeRaw;
  double diff;

  refTime.tm_year = year - 1900;
  refTime.tm_mon = month - 1;
  refTime.tm_mday = day;
  refTime.tm_hour = hour;
  refTime.tm_min = minute;
  refTime.tm_sec = second;
  refTime.tm_isdst = 0;
  refTimeRaw = mktime(&refTime);

  epochTime.tm_year = 80;//start of GPS time Jan 6 1980
  epochTime.tm_mon = 1 - 1;
  epochTime.tm_mday = 6;
  epochTime.tm_hour = 0;
  epochTime.tm_min = 0;
  epochTime.tm_sec = 0;
  epochTime.tm_isdst = 0;
  epochTimeRaw = mktime(&epochTime);

  diff = difftime(refTimeRaw, epochTimeRaw);

  return diff;
}

static float q15tof(int16_t x) { return (float)x / 32768.0f; }
static int16_t f2q15(float x) {
  if (x >= 0.9999695f) {
    return Q15_MAX; // ~ (32767/32768)
  }
  if (x <= -1.0f) {
    return Q15_MIN;
  }
  int32_t v = (int32_t)(x * 32768.0f);
  if (v > Q15_MAX) {
    v = Q15_MAX;
  }
  if (v < Q15_MIN) {
    v = Q15_MIN;
  }
  return (int16_t)v;
}

int to_fixed(float f, int e) {
  float tmp = (float) pow(2.0f, e);
  int64_t b = (int64_t)llround(f * tmp);
  //int a = (int) (f * tmp);
  //int b = (int)round(aa);

  if (b > 0 && b > 2147483647) { // 2147483647 32767
    return 2147483647;
  }
  if (b < 0 && b <= -2147483648) { // -32768 
    return -2147483648;
  }
  if (b < 0) {
    // next lines turn b into it's 2's complement.
    b = abs(b);
    b = ~b;
  }
  return (int) b;
}

/**
 * @brief Function to generate the arm Bit Reversal Table
 * The main case is 16384 with log_2N is 14
 */
void armBitRevTableCalculator(void) {
#define N 16384//4096//8192//16384
#define logN2 14//12//13//14
  int a[logN2] = { 0 };
  int y[(N / 4) + 1] = { 0 };
  int i, j;
  for (int l = 1; l <= N / 4; l++) {
    for (i = 0; i < logN2; i++) {
      a[i] = l & (1 << i);
    }
    for (j = 0; j < logN2; j++) {
      if (a[j] != 0) {
        y[l] += (1 << ((logN2 - 1) - j));
      }
    }
    y[l] = y[l] >> 1;
  }

  for (i = 0; i < ((N / 4) + 1); i++) {
    printf("%#x ,", y[i]);
  }
}

double magnitude(int64_t real, int64_t imag) {
  double real_d = (double)real;
  double imag_d = (double)imag;
  return sqrt(real_d * real_d + imag_d * imag_d);
}

double recover_prn_phase_deg_with_doppler(const c32* iandq,
  int SIZE,
  const int prn_a[], // tested with 1024*4
  int  code_offset,
  double doppler_a_hz,
  double sampe_rate_hz,
  int* best_offset_out,
  float* doppler_out,
  float* pwr)
{
  double best_pow = -1.0, best_re = 0.0, best_im = 0.0, best_doppler = 0.0;
  int best_off = 0;

  float code_search_span = 10; //chip span was 2
  float code_low = code_offset - code_search_span;
  float code_high = code_offset + code_search_span;
  uint8_t wrap = 0;
  if (code_low < 0) { code_low += SIZE; wrap = 1; }
  if (code_high > SIZE) { code_high -= SIZE;  wrap = 1; }

  float f_search_step_hz = 1.0; // Hz
  float f_search_span_hz = 1; // Hz span was 2
  for (float f = -f_search_span_hz; f <= +f_search_span_hz + 1e-3f; f += f_search_step_hz) {
    // Frequency wipe for PRN A: multiply by e^{-j 2π f_a n / fs}
    double dtheta = 2.0 * PI * (doppler_a_hz + f) / (double)sampe_rate_hz;

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
        re += iw[n] * prn_a[idx];
        im += qw[n] * prn_a[idx];
      }
      double p = sqrt(re * re + im * im);
      if (p > best_pow) {
        best_pow = p;
        best_doppler = (double)f;
        best_re = re;
        best_im = im;
        best_off = off;
      }
      
    } // for off 
  } // end for f_search_span_hz
  *doppler_out = best_doppler;
  *best_offset_out = best_off;
  *pwr = best_pow;
  //printf(" best %.3f sec %.3f %d %d %d \n", best_pow, second_best_pow, best_off, second_off, (int)fabs(best_off - second_off));
  double phase_deg = atan2(best_im, best_re) * R2D;
  return phase_deg;
}


/**
 * @brief Function to generate the arm Twiddle Coefficient Table
 * The main case is 16384 with dyad = 15 for q15_t (or 31 for q31_t)
 *
 * An adjunct is the python script genBitsReversal.py (will require installing sympy.combinatorics)
 * which will generate armBitRevIndex_fixed_N (eg armBitRevIndexTable_fixed_4096 used in
 * arm_cfft_sR_f32_len4096)
 */
void twiddleCoefCalculator() {
  int fft_len = 16384;// 4096;// 8192;
  int dyad = 31; // 15 for q15_t, 31 for q31_t (nothing to do with fft_len)
  
  // the coeffs (not twiddle factors)
  //#define q31_t int64_t
  // consider using version without '\n' in printf for one line array (more readable than 16k line)
  for (int i = 0; i < (3 * fft_len / 4); i++) {
    float twiddleCoefq15Even = cos(i * 2 * PI / (float)fft_len);
    float twiddleCoefq15Odd = sin(i * 2 * PI / (float)fft_len);
    //printf("(q15_t) %#04hx, (q15_t) %#04hx,", (q15_t)to_fixed(twiddleCoefq15Even, dyad), (q15_t)to_fixed(twiddleCoefq15Odd, dyad));
    printf("(q31_t) %d, (q31_t) %d, \n", (int)to_fixed(twiddleCoefq15Even, dyad), (int) to_fixed(twiddleCoefq15Odd, dyad));
  }
}

float signalf(float x, float F1, float F2, float F3, float F4) {
  float ans = 0.5 * sin(F1 * 2.0 * PI * x) + 0.5 * sin(F2 * 2.0 * PI * x) + 0.5 * sin(F3 * 2.0 * PI * x) + 0.5 * sin(F4 * 2.0 * PI * x);
  return ans;
}

void DecodeOrsIQCplx(uint8_t* data, uint32_t byteLength, c32 iqs[])
{
  /// A 2-bit look up table for decoding
  static int16_t TwoBitTable[4] = { 1,-1,3,-3 };

  for (uint32_t i = 0; i < byteLength; i++) {
    uint8_t byte = data[i];

    uint8_t bits = (byte) & 0x03;
    iqs[2 * i].r = (float)TwoBitTable[bits];

    bits = (byte >> 2) & 0x03;
    iqs[2 * i].i = (float)TwoBitTable[bits];

    bits = (byte >> 4) & 0x03;
    iqs[2 * i + 1].r = (float)TwoBitTable[bits];

    bits = (byte >> 6) & 0x03;
    iqs[2 * i + 1].i = (float)TwoBitTable[bits];
  }
}

void read_L1(char* input) {
  FILE* fp_msb = NULL;
  fopen_s(&fp_msb, input, "r");
  if (fp_msb == NULL) {
    fprintf(stderr, "Failed to open msb file %s\n", input);
    return;
  }
  fseek(fp_msb, 0L, SEEK_END);
  size_t bytes_to_read = ftell(fp_msb);
  rewind(fp_msb);
  

  FILE* fp_out = NULL; //output file
  errno_t er = fopen_s(&fp_out, "C:/Python/out5.csv", "w");
  if (er != 0 || fp_out == NULL) {
    fprintf(stderr, "Failed to open output file\n");
    return;
  }
 
  //// Dial in the prn and doppler here ////////////////
  #define SPC 1 // samples per chip
  #define SIZE 1024*SPC *1 // 1 for GPS 4 for Gal -> 16K for Galileo and 4K for GPS
  int proc_gps = 1; // 1 for GPS, 0 for Galileo
  int prn = 3;// 4;// 10;// 10;// 11;// 4;
  double doppler = -2570;// 805;// 1232;// -582;// -2263;// -912;// 67;
  /////////////////////////////////////////////////////

  c32* iandq = (c32*)malloc(SIZE * sizeof(c32));
  c32* repli = (c32*)malloc(SIZE * sizeof(c32));
  if (iandq == NULL || repli == NULL) {
    fprintf(stderr, "Memory allocation failed for q32 array.\n");
    free(iandq); free(repli);
    return;
  }

  char line[256];
  for (int i = 0; i < 0 * 1024; i++) {
    fgets(line, sizeof(line), fp_msb);
  }

  char* context = nullptr;
  // read in the csv data
  while (!feof(fp_msb)) {
    if (fgets(line, sizeof(line), fp_msb) != NULL) {
      // Process the line
      static int idx = 0;
      if (idx >= SIZE) { break; }
      char* token = strtok_s(line, ",", &context); // eat up the first ordinal
      //token = strtok_s(NULL, ",", &context); 
      if (token != NULL) {
        iandq[idx].r = (float)atof(token);
        token = strtok_s(NULL, ",", &context);
        if (token != NULL) {
          iandq[idx].i = (float)atof(token);
        }
        idx++;
      }
    }
  }
  fclose(fp_msb);

  float* actual = (float*)malloc(SIZE * 2 * sizeof(float));
  float* replica = (float*)malloc(SIZE * 2 * sizeof(float));
  float* prod = (float*)malloc(SIZE * 2 * sizeof(float));
  if (actual == NULL || replica == NULL || prod == NULL) {
    fprintf(stderr, "Memory allocation failed for 'actual'.\n");
    return; // Exit or handle the error appropriately
  }
  memset(actual, 0, sizeof(float) * SIZE * 2);
  memset(replica, 0, sizeof(float) * SIZE * 2);
  memset(prod, 0, sizeof(float) * SIZE * 2);

  if (proc_gps) {
    synth_gps_prn(prn, -doppler, SIZE, repli, SPC);
  }
  else { // Galileo
    synth_e1b_prn(prn, -doppler, SIZE, repli);
  }

  if (1) {
    arm_cfft_radix2_instance_f32 as;
    arm_cfft_radix2_instance_f32 rs;

    for (int i = 0; i < SIZE; i++) {
      actual[2 * i] = iandq[i].r * 0.25;
      actual[2 * i + 1] = iandq[i].i * 0.25;
      replica[2 * i] = repli[i].r * 0.25;
      replica[2 * i + 1] = repli[i].i * 0.25;
    }
    if (iandq != NULL) { free(iandq); }
    if (repli != NULL) { free(repli); }

    // do the float thing
    arm_cfft_radix2_init_f32(&as, SIZE, 0, 1); // Initialize the CFFT instance for 8-point FFT
    arm_cfft_radix2_f32(&as, actual);

    arm_cfft_radix2_init_f32(&rs, SIZE, 0, 1); // Initialize the CFFT instance for 8-point FFT
    arm_cfft_radix2_f32(&rs, replica);


    // conjugate the replica
    for (int i = 0; i < SIZE; i++) {
      float Ar = actual[i * 2], Ai = actual[i * 2 + 1];
      float Rr = replica[i * 2], Ri = replica[i * 2 + 1];
      // A * conj(R)
      prod[i * 2] = Ar * Rr + Ai * Ri;     // (Ar + jAi) * (Rr - jRi)
      prod[i * 2 + 1] = Ai * Rr - Ar * Ri;     // 
    }

    arm_cfft_radix2_instance_f32 conv;
    arm_cfft_radix2_init_f32(&conv, SIZE, 1, 1); // inverse FFT
    arm_cfft_radix2_f32(&conv, prod);
    float max = 0;
    int pos = 0;

    for (int i = 0; i < SIZE; i++) {
      float mag = sqrt(prod[2 * i] * prod[2 * i] + prod[2 * i + 1] * prod[2 * i + 1]);
      fprintf(fp_out, "%d, %f \n", i, mag);
      if (mag > max) { max = mag; pos = i; }
    }
    printf("max_float = %f pos=%d\n", max, pos);
  }

}

void read_ors(char* input) {
  FILE* fp_msb = NULL;
  fopen_s(&fp_msb, input, "rb");
  if (fp_msb == NULL) {
    fprintf(stderr, "Failed to open msb file %s\n", input); return;
  }
  fseek(fp_msb, 0L, SEEK_END);
  size_t bytes_to_read = ftell(fp_msb);
  rewind(fp_msb);

  uint8_t* buffer = (uint8_t*)malloc(bytes_to_read);
  size_t bytesRead;
  bytesRead = fread(buffer, 1, bytes_to_read, fp_msb);
  if (bytesRead != bytes_to_read) { printf("Error parsiing data\n"); }
  fclose(fp_msb);

  // 2 for res 1
  uint16_t hdrlen = U2(&buffer[2]); // 2 for header length
  // 6 for res2
  uint16_t yr = U2(&buffer[10]) + 1900; // 12
  uint16_t mon = buffer[12] + 1; // 13
  uint16_t day = buffer[13]; // 14
  double tods = U4(&buffer[14]) / 1000.0; // 18
  int hr = (int)floor(tods / 3600.0);
  int min = (int)floor((tods / 3600.0 - floor(tods / 3600.0)) * 60);
  int sec = (int)floor((tods / 60.0 - floor(tods / 60.0)) * 60);
  // map back to GPS time
  double dtime = compute_gps_time(yr, mon, day, hr, min, sec);
  int msec = (int)round((sec - floor(sec)) * 1000);
  uint32_t nanos = U3(&buffer[18]); // 21
  printf("%4d-%02d-%02d-%02d-%02d-%02d msec=%d nanos=%d \n", yr, mon, day, hr, min, sec, msec, nanos);
  uint16_t size = ((uint16_t)bytes_to_read) - hdrlen;
  printf("size (total bytes - hdrs) %d \n", size);

  FILE* fp_out = NULL; //output file
  errno_t er = fopen_s(&fp_out, "C:/Python/out2.csv", "w");
  if (er != 0 || fp_out == NULL) {
    fprintf(stderr, "Failed to open output file\n"); return;
  }

  //// Dial in the prn and doppler here ////////////////
#define SPC  4
#define SIZE 1024*SPC *4 // 16K for Galileo and 4K for GPS

  typedef struct {
    int prn;
    double doppler;
    int constel;
  } acq_struct;

  acq_struct prn2acq[15] = {0};
  prn2acq[0].prn = 4 ; prn2acq[0].doppler = 791;   prn2acq[0].constel = 1; // GPS
  prn2acq[1].prn = 9 ; prn2acq[1].doppler = -112;  prn2acq[1].constel = 2; // GAL
  prn2acq[2].prn = 6 ; prn2acq[2].doppler = -1243; prn2acq[2].constel = 2;
  prn2acq[3].prn = 26; prn2acq[3].doppler = -544;  prn2acq[3].constel = 1;
  prn2acq[4].prn = 3 ; prn2acq[4].doppler = -2584; prn2acq[4].constel = 1; 
  prn2acq[5].prn = 31; prn2acq[5].doppler = 1915;  prn2acq[5].constel = 2;
  prn2acq[6].prn = 5 ; prn2acq[6].doppler = 2268;  prn2acq[6].constel = 2;
  prn2acq[7].prn = 16; prn2acq[7].doppler = 2822;  prn2acq[7].constel = 1;
  prn2acq[8].prn = 9 ; prn2acq[8].doppler = 2986;  prn2acq[8].constel = 1;
  prn2acq[9].prn = 31; prn2acq[9].doppler = -2867; prn2acq[9].constel = 1;
  prn2acq[10].prn = 23; prn2acq[10].doppler = -867;  prn2acq[10].constel = 2;
  prn2acq[11].prn = 4 ; prn2acq[11].doppler = -2020; prn2acq[11].constel = 2;
  prn2acq[12].prn = 34; prn2acq[12].doppler = -1027; prn2acq[12].constel = 2;
  prn2acq[13].prn = 6 ; prn2acq[13].doppler = -882;  prn2acq[13].constel = 1;
  prn2acq[14].prn = 24; prn2acq[14].doppler = 3545;  prn2acq[14].constel = 2;
  
  /////////////////////////////////////////////////////

  c32* iandq = (c32*)malloc(SIZE * sizeof(c32));
  c32* signl = (c32*)malloc(SIZE * sizeof(c32));
  c32* repli = (c32*)malloc(SIZE * sizeof(c32));
  c32* prod  = (c32*)malloc(SIZE * sizeof(c32));
  if (iandq == NULL || repli == NULL) {
    fprintf(stderr, "Memory allocation failed for q32 array.\n");
    free(iandq); free(repli); free(prod); return;
  }

  int num_meas = 0;
  bb_meas_t meas;
  memset(&meas, 0, sizeof(bb_meas_t));  

  DecodeOrsIQCplx(&buffer[hdrlen], SIZE / 2, iandq);
  free(buffer);

  if (0) {
    FILE* fp1bt_out = NULL; //output file
    errno_t er2 = fopen_s(&fp1bt_out, "C:/Python/out-1bit.csv", "w");
    if (er2 != 0 || fp1bt_out == NULL) {
      fprintf(stderr, "Failed to open output file\n");
      return;
    }
    fprintf(fp1bt_out, "I, Q\n");
    for (int i = 0; i < 4 * 1024; i++) {
      fprintf(fp1bt_out, "%f, %f \n",0.25* iandq[i].r,0.25* iandq[i].i); 
      printf("%d: I=%f Q=%f \n", i, iandq[i].r, iandq[i].i);
    }
    fclose(fp1bt_out);
  }

  for (int loop = 0; loop < 15; loop++) {

    int size = prn2acq[loop].constel == 1 ? 4096 : 16384;
    int proc_gps = (prn2acq[loop].constel == 1) ? 1 : 0;
    if (proc_gps) {
      synth_gps_prn(prn2acq[loop].prn, -prn2acq[loop].doppler, size, repli, SPC);
    }
    else { // Galileo
      synth_e1b_prn(prn2acq[loop].prn, -prn2acq[loop].doppler, size, repli);
    }
    memcpy(signl, iandq, size * sizeof(c32));
    memset(prod, 0, size * sizeof(c32));

    fft_c32(size, repli, true); // F(repli)
    fft_c32(size, signl, true); // F(iandq) and prod=F(iandq) * conj(F(repli)) below
    for (int i = 0; i < size; i++) { prod[i] = mult(signl[i], get_conj(repli[i])); }
    fft_c32(size, prod, false); // in-place inv F(prod) 

    top2_pks peaks;
    find_top2_peaks_cplx(prod, size, 3, &peaks, fp_out);
    fclose(fp_out);
    // compute noise stats for SNR
    double BW = 0.3e3; // 10 MHz
    double cn0 = compute_snr_cplx(prod, size, peaks.val1, peaks.idx1) + 35;// +10 * log(BW);

    //printf("%s %d v1=%f idx1=%d ; v2=%f idx2=%d ratio=%f len=%d\n", (proc_gps == 1) ? "GPS" : "GAL", prn2acq[loop].prn, 
    //  peaks.val1, peaks.idx1, peaks.val2, peaks.idx2, (peaks.val1 / peaks.val2), size);
    double early = mag(prod[peaks.idx1 - 1]), prompt = peaks.val1, late = mag(prod[peaks.idx1 + 1]);
    double interp = InterpolateCodePhase(peaks.idx1, early, prompt, late);
    //printf("interpolated %f cn0 %f\n", interp, cn0);

    if ((peaks.val1 / peaks.val2) > 1.3) {
      meas.sats[meas.num_sat].prn = prn2acq[loop].prn;
      meas.sats[meas.num_sat].code_phase = float(interp) / 4096.0f;
      meas.sats[meas.num_sat].doppler = -prn2acq[loop].doppler;
      meas.sats[meas.num_sat].cno = (float)cn0;
      meas.sats[meas.num_sat].constellation = proc_gps ? SYS_GPS : SYS_GAL;
      
      float ratio = (peaks.val1 / peaks.val2);
      printf("Acquired %s %d Doppler %f Hz CodePhase %f [ms] C/N0 %f dB-Hz ratio=%f\n", (proc_gps == 1) ? "GPS" : "GAL",
        prn2acq[loop].prn, -prn2acq[loop].doppler,meas.sats[meas.num_sat].code_phase, meas.sats[meas.num_sat].cno, ratio* ratio);
      meas.num_sat++;
    }
  }
  write_msb(&meas, (char*)"C:/Python/out2.bin");
   
  free(iandq); free(repli); free(prod);
  //////////////////////////////////////////////////////
}

void sim_E5A() {
  #define FFT_SIZE 16384
  double doppler_a = 1580;
  int prn_a = 5;
  double doppler_b = -2580;
  int prn_b = 15;
  int offset = 5000;
  float up_offset = (float)offset * 16384.0f / 10230.0f; // upsampled offset
  printf("should be %f \n", up_offset);

  c32* samp_a  = (c32*)malloc(E5A_CODE_LEN * sizeof(c32));
  c32* samp_b  = (c32*)malloc(E5A_CODE_LEN * sizeof(c32));
  c32* repli_a = (c32*)malloc(E5A_CODE_LEN * sizeof(c32));
  
  synth_e5a_prn(prn_a, doppler_a, SIZE, samp_a, offset);
  synth_e5a_prn(prn_b, doppler_b, SIZE, samp_b, 0);
  synth_e5a_prn(prn_a, doppler_a + 100, SIZE, repli_a, 0);

  for (int i = 0; i < E5A_CODE_LEN; i++) {
    samp_a[i].r += samp_b[i].r; // add noise & quantize later
    samp_a[i].i += samp_b[i].i;
  }
  free(samp_b);
  // now up sample both samples and replicas

  c32* up_samp  = (c32*)malloc(FFT_SIZE * sizeof(c32));
  c32* up_repli = (c32*)malloc(FFT_SIZE * sizeof(c32));
  
  up_sample_10k_to_16k(samp_a, up_samp);
  up_sample_10k_to_16k(repli_a, up_repli);
  free(samp_a); free(repli_a);

  /*
  FILE* dbg_fp = NULL;
  fopen_s(&dbg_fp, "C:/Python/check.csv", "w");
  for (int i = 0; i < FFT_SIZE; i++) {
    fprintf(dbg_fp, "%d, %f, %f \n", i, up_samp[i].r, up_samp[i].i);
  }
  fclose(dbg_fp);
  */

  //now do the FFT based correlation
  float* actual = (float*)malloc(FFT_SIZE * 2 * sizeof(float));
  float* replica = (float*)malloc(FFT_SIZE * 2 * sizeof(float));
  float* prod = (float*)malloc(FFT_SIZE * 2 * sizeof(float));
  for (int i = 0; i < FFT_SIZE; i++) {
    actual[2 * i] = up_samp[i].r * 0.25;
    actual[2 * i + 1] = up_samp[i].i * 0.25;
    replica[2 * i] = up_repli[i].r * 0.25;
    replica[2 * i + 1] = up_repli[i].i * 0.25;
  }
  free(up_samp); free(up_repli);

  arm_cfft_radix2_instance_f32 as;
  arm_cfft_radix2_instance_f32 rs;

  // do the float thing
  arm_cfft_radix2_init_f32(&as, FFT_SIZE, 0, 1); // Initialize the CFFT instance for 8-point FFT
  arm_cfft_radix2_f32(&as, actual);

  arm_cfft_radix2_init_f32(&rs, FFT_SIZE, 0, 1); // Initialize the CFFT instance for 8-point FFT
  arm_cfft_radix2_f32(&rs, replica);

  // conjugate the replica
  for (int i = 0; i < FFT_SIZE; i++) {
    float Ar = actual[i * 2], Ai = actual[i * 2 + 1];
    float Rr = replica[i * 2], Ri = replica[i * 2 + 1];
    // A * conj(R)
    prod[i * 2] = Ar * Rr + Ai * Ri;     // (Ar + jAi) * (Rr - jRi)
    prod[i * 2 + 1] = Ai * Rr - Ar * Ri;     // 
  }
  free(actual); free(replica);
  arm_cfft_radix2_instance_f32 conv;
  arm_cfft_radix2_init_f32(&conv, FFT_SIZE, 1, 1); // inverse FFT
  arm_cfft_radix2_f32(&conv, prod);
  
  float max = 0;
  int pos = 0;

  FILE* dbg_fp = NULL;
  fopen_s(&dbg_fp, "C:/Python/out6.csv", "w");
  for (int i = 0; i < FFT_SIZE; i++) {
    float mag = sqrt(prod[2 * i] * prod[2 * i] + prod[2 * i + 1] * prod[2 * i + 1]);
    fprintf(dbg_fp, "%d, %f \n", i, mag);
    if (mag > max) { max = mag; pos = i; }
  }
  free(prod);
  fclose(dbg_fp);
  printf("max_float = %f pos=%d\n", max, pos);
}

/////////////////////////////////////////////////////////////////////////////
void read_E5A(char* input) {
  FILE* fp_1bitcsv = NULL;
  fopen_s(&fp_1bitcsv, input, "r");
  if (fp_1bitcsv == NULL) {
    fprintf(stderr, "Failed to open msb file %s\n", input);
    return;
  }
  fseek(fp_1bitcsv, 0L, SEEK_END);
  size_t bytes_to_read = ftell(fp_1bitcsv);
  rewind(fp_1bitcsv);


  FILE* fp_out = NULL; //output file
  errno_t er = fopen_s(&fp_out, "C:/Python/out3.csv", "w");
  if (er != 0 || fp_out == NULL) {
    fprintf(stderr, "Failed to open output file\n");
    return;
  }

  //// Dial in the prn and doppler here ////////////////
  int gal_proc = 1; // 1 for Galileo, 0 for GPS
#define SPC 1
#define FFT_SIZE 16384 
#define SAMP 10230 // for 1 ms at 10.23 MHz
  // only 9 and 36 with q31; 10, 6 also works with float
  int prn = 36;// 6;// 6;// 36;// 9;// 36
  double doppler = -1*(1580 + 1e6 +2500);// E36:1580 E6: 1261 ; G10:-582, G32:1232
#define DO_FLOAT // q31 ifndef
  /////////////////////////////////////////////////////
  c32* sampl = (c32*)malloc(SAMP * sizeof(c32));
  c32* repli = (c32*)malloc(SAMP * sizeof(c32));
  memset(sampl, 0, sizeof(c32) * SAMP);
  memset(repli, 0, sizeof(c32) * SAMP);

  if (sampl == NULL || repli == NULL) {
    fprintf(stderr, "Memory allocation failed for q32 array.\n");
    free(sampl); free(repli);
    return;
  }

  //up sample
  c32* up_samp  = (c32*)malloc(FFT_SIZE * sizeof(c32));
  c32* up_repli = (c32*)malloc(FFT_SIZE * sizeof(c32));
  c32* up_prod  = (c32*)malloc(FFT_SIZE * sizeof(c32));
  c32* sum_prod = (c32*)malloc(FFT_SIZE * sizeof(c32));
  float* nci_sum = (float*)malloc(FFT_SIZE * sizeof(float));
  memset(up_samp , 0, sizeof(c32) * FFT_SIZE);
  memset(up_repli, 0, sizeof(c32) * FFT_SIZE);
  memset(up_prod , 0, sizeof(c32) * FFT_SIZE); 
  memset(nci_sum , 0, sizeof(float) * FFT_SIZE);
  memset(sum_prod, 0, sizeof(c32) * FFT_SIZE);

  float* replica = (float*)malloc(FFT_SIZE * 2 * sizeof(float));
  if (replica == NULL) {
    fprintf(stderr, "Memory allocation failed for 'replica'.\n");
    return; // Exit or handle the error appropriately
  }
  memset(replica, 0, sizeof(float) * FFT_SIZE * 2); 


  if (gal_proc) {
    synth_e5a_prn(prn, -doppler, SAMP, repli, 0);
  }
  else { // GPS
    synth_L5I_prn(prn, -doppler, SAMP, repli, 0);
  }
  up_sample_10k_to_16k(repli, up_repli);
  free(repli);

  fft_c32(FFT_SIZE, up_repli, true); // forward FFT 

  char* context = nullptr;
  // read in the csv data
  char line[256];
  int LEN = SAMP;
  ///////////// main loop /////////////////////////////////////////////
  for (int loop = 0; loop < 20; loop++) {
    int idx = 0;
    while (!feof(fp_1bitcsv)) {
      if (fgets(line, sizeof(line), fp_1bitcsv) != NULL) {
        char* token = strtok_s(line, ",", &context);
        token = strtok_s(NULL, ",", &context);
        if (token != NULL) {
          sampl[idx].r = (float)atof(token);
          token = strtok_s(NULL, ",", &context);
          if (token != NULL) {
            sampl[idx].i = (float)atof(token);
          }
        }
        idx++;
        if ((idx != 0) && (idx % LEN == 0)) {
          break;
        }
      }
    }
    up_sample_10k_to_16k(sampl, up_samp);
    
    // note repli has been FFTed already
    fft_c32(FFT_SIZE, up_samp, true); // forward FFT
    for (int k = 0; k < FFT_SIZE; k++) { up_prod[k] = mult(up_samp[k], get_conj(up_repli[k])); }
    fft_c32(FFT_SIZE, up_prod, false); // IFFT
    for (int i = 0; i < FFT_SIZE; i++) { nci_sum[i] += mag(up_prod[i]); }

    printf("loop %d \n", loop);
  } // end NCI for loop


  float max = 0; float max2 = 0;
  int pos = 0; int pos2 = 0;
  for (int i = 0; i < FFT_SIZE; i++) {
    float mag = nci_sum[i];

    fprintf(fp_out, "%d, %f \n", i, mag);
    if (mag > max) { max = mag; pos = i; }
  }
  printf("max_float = %f pos=%d frac=%f\n", max, pos, (pos / 16384.0));
 
  fclose(fp_1bitcsv); fclose(fp_out);
  free(sampl); free(up_samp); free(up_repli); free(up_prod); free(nci_sum); free(sum_prod);
}

// try differential coherent integration. Insensitive to bit transitions
void test_quasi_diff_pilot() {
  srand((unsigned int)time(NULL)); // randomise seed
  // Test the quasi pilot generation

  int min_idx = 0;
  int loc_cnt = 0;
  float min_val = 1e5;

  int locations[50] = { 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280 };
  int window = 2; // 2 * window ms either side of center (window>=5 does not work with noise window > 7 does not work with no noise)
  int nci = 300;
#define SPC 4 // samples per chip
  int len = 1023 * SPC * nci; // 4 samples per chip and 100 ms
  int c_phase = 4*4096 /8 ;// 1023 * SPC / 8; // which chip to set the code phase to
  int prn1 = 5, prn2 = 15;
  float dop1 = 2000, dop2 = -3000;
  float dop_error = 250;// 10; // full 2*250 Hz error in wipeoff
  float dop_err_rate = 0.6;// 0.6;// 0.6;//Hz per ms
  float sigma = 3.4;// 3.4;// 3.5;// 3.5; // noise level
  c32* out = (c32*)malloc(len * sizeof(c32));
  if (out == NULL) {
    fprintf(stderr, "Memory allocation failed for 100 ms I&Q array.\n");
    return;
  }
  int* prn_c1 = (int*)malloc(sizeof(int) * 1023 * SPC);
  int* prn_c2 = (int*)malloc(sizeof(int) * 1023 * SPC);
  c32* replica = (c32*)malloc(sizeof(c32) * 1024 * SPC);
  memset(replica, 0, sizeof(c32) * 1024 * SPC);
  getCode(1023, SPC, prn1, prn_c1);
  getCode(1023, SPC, prn2, prn_c2);
  make_replica(prn_c1, replica, dop1 + dop_error, 1023 * SPC, 1.023e6 * SPC);
  // now advance code-phase
  rotate_fwd(prn_c1, 1023 * SPC, c_phase); // code phase 1/4 way
  int sign2 = 1; // sign applied a posteriori after finding BTT
  stat_s stat;
  stat_init(&stat); // moving average of peak values window size = 3
  for (int i = 0; i < nci; i++) {
    for (int j = 0; j < 50; j++) {
      if (locations[j] == i) { sign2 *= -1; break; } // change sign at the bit transitions
    }
    // offset doppler by 250 Hz and add a residual doppler ramp of 0.1 Hz per ms
    mix_two_prns_oversampled_per_prn(prn_c1, prn_c2, dop1 + i * dop_err_rate, dop2 - i * dop_err_rate, PI / 2, 0,
      &out[1023 * SPC * i], 1023 * SPC, 1.023e6 * SPC, sigma, sign2); // was 2.31 for -128.5 dBm 3.1 for -131.5
  }
  free(prn_c1); free(prn_c2);

  fft_c32(1024 * SPC, replica, true);
  
  // Compute circular correlation C_k(τ) = FFT^-1{ FFT[x_d,k] · conj(FFT[code]) }.
  c32* fft_data = (c32*)malloc(sizeof(c32) * 1024 * SPC);
  c32* fft_prev = (c32*)malloc(sizeof(c32) * 1024 * SPC);
  c32* diff_acc = (c32*)malloc(sizeof(c32) * 1024 * SPC);
  bool have_prev = false;
  if (fft_data == NULL ) { printf("Error allocating fft_data \n"); return; }
  for (int center = window / 2; center <= nci - window / 2; center++) {
    memset(fft_prev, 0, sizeof(c32) * 1024 * SPC);
    memset(diff_acc, 0, sizeof(c32) * 1024 * SPC);
    for (int windex = center - window / 2; windex < center + window / 2; windex++) {
      memset(fft_data, 0, sizeof(c32) * 1024 * SPC);
      for (int j = 0; j < 1023 * SPC; j++) { // xfer to cplx float array
        fft_data[j] = out[(1023 * SPC * windex) + j];
      }

      fft_c32(1024 * SPC, fft_data, true); // forward FFT
    
      for (int k = 0; k < 1024 * SPC; k++) { // pt-wise * with conj of replica
        fft_data[k] = mult(fft_data[k], get_conj(replica[k]));
      }

      fft_c32(1024 * SPC, fft_data, false); // IFFT
      
      // multiply with conjugate of previous result and accumulate
      if (have_prev) {
        for (int k = 0; k < 1024 * SPC; k++) {
          diff_acc[k] = add(diff_acc[k], mult(fft_data[k], get_conj(fft_prev[k])));
        }
      }
      
      memcpy(fft_prev, fft_data, sizeof(c32) * 1024 * SPC);
      have_prev = true;
    } // for windex 

    float max_coh = 0; int pos_coh = 0;
    for (int m = 0; m < 1024 * SPC; m++) {
      float mag_coh = mag(diff_acc[m]);
      if (mag_coh > max_coh) { max_coh = mag_coh; pos_coh = m; }
    }

    stat_add(&stat, max_coh);
    float mean2 = stat_mean(&stat);
    if (max_coh < min_val && max_coh < mean2 - 45 && max_coh < 160) { // for 4 SPC
      //if (max_coh < min_val && max_coh < mean2 - 15 && max_coh < 81) { // for 2 SPC
      //if (max_coh < min_val && max_coh < mean2 - 10 && max_coh < 42) { // sigma = 2.0 for 1 SPC (odd missed detect)
      min_val = max_coh;
      min_idx = center;
    }
    if ((center == min_idx + 1) && (max_coh > min_val)) {
      locations[loc_cnt++] = min_idx; // empirically the wider the window the earlier the bit transition appears
      min_val = 1e5;
      min_idx = 0;
    }
    printf("center=%d max=%6.1f pos=%d mean=%6.1f\n", center, max_coh, pos_coh, mean2);

    if (0) {//center == 270) {//fabs(pos_coh - c_phase) > 50) {
      printf("Warning: large code phase error at center %d pos_coh=%d c_phase=%d\n", center, pos_coh, c_phase);
      FILE* fp_out = NULL; //output file
      errno_t er = fopen_s(&fp_out, "C:/Python/diff_corr.csv", "w");
      for (int m = 0; m < 1024 * SPC; m++) {
        double magnitude = mag(diff_acc[m]);
        fprintf(fp_out, "%d, %f\n", m, magnitude);
      }
      fclose(fp_out);
    }
  } // for center

  for (int i = 0; i < loc_cnt; i++) {
    printf("Bit transition at %d ms \n", locations[i]);
  }
  printf("BTs: %d Random number: %d\n", loc_cnt, rand());
  free(out); free(fft_data);
  free(fft_prev); free(diff_acc); free(replica);
}

void test_quasi_pilot() {
  srand((unsigned int)time(NULL)); // randomise seed
  // Test the quasi pilot generation
 
  int min_idx = 0;
  int loc_cnt = 0;
  float min_val = 1e5;
  
  int locations[50] = { 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280 };
  int window = 3; // 2 * window ms either side of center (window>=6 does not work)
  int nci = 300;
#define SPC 4 // samples per chip
  int len = 1023 * SPC * nci; // 4 samples per chip and 100 ms
  int c_phase = 500; // which chip to set the code phase to
  int prn1 = 4, prn2 = 8;
  float dop1 = 2000, dop2 = -3000;
  float dop_error = 250;// 10; // full 2*250 Hz error in wipeoff
  float dop_err_rate = 0.6;// 0.6;// 0.6;//Hz per ms
  float sigma = 3.5;// 3.5; // noise level
  c32* out = (c32*)malloc(len * sizeof(c32));
  if (out == NULL) {
    fprintf(stderr, "Memory allocation failed for 100 ms I&Q array.\n");
    return;
  }
  int* prn_c1  = (int*)malloc(sizeof(int) * 1023 * SPC);
  int* prn_c2  = (int*)malloc(sizeof(int) * 1023 * SPC);
  c32* replica = (c32*)malloc(sizeof(c32) * 1024 * SPC);
  memset(replica, 0, sizeof(c32) * 1024 * SPC);
  getCode(1023, SPC, prn1, prn_c1);
  getCode(1023, SPC, prn2, prn_c2);
  make_replica(prn_c1, replica, dop1 + dop_error, 1023* SPC, 1.023e6 * SPC);
  rotate_fwd(prn_c1, 1023 * SPC, c_phase); // now advance code-phase
  int sign2 = 1; // sign applied a posteriori after finding BTT
  stat_s stat;
  stat_init(&stat); // moving average of peak values window size = 3
  for (int i = 0; i < nci; i++) {
    for (int j = 0; j < 50; j++) {
      if (locations[j] == i) { sign2 *= -1; break; } // change sign at the bit transitions
    }
    // offset doppler by 250 Hz and add a residual doppler ramp of 0.1 Hz per ms
    mix_two_prns_oversampled_per_prn(prn_c1, prn_c2,dop1 + i * dop_err_rate,dop2 - i * dop_err_rate,PI/2,0,
      &out[1023 * SPC * i],1023* SPC, 1.023e6 * SPC, sigma , sign2); // was 2.31 for -128.5 dBm 3.1 for -131.5
  }
  free(prn_c1); free(prn_c2);
  
  // use FFTs 
  fft_c32(1024 * SPC, replica, true);

  // Compute circular correlation C_k(τ) = FFT^-1{ FFT[x_d,k] · conj(FFT[code]) }.
  c32* fft_data = (c32*)malloc(sizeof(c32) * 1024 * SPC);
  c32* fft_prod = (c32*)malloc(sizeof(c32) * 1024 * SPC);
  c32* fft_sum  = (c32*)malloc(sizeof(c32) * 1024 * SPC);
  if (fft_data == NULL || fft_prod == NULL) { printf("Error allocating fft_data or fft_prod\n"); return; }
  for (int center = window/2; center <= nci - window/2; center++) {
    // use linearity to sum the ncis coherently in moving window
    memset(fft_data, 0, sizeof(c32) * 1024 * SPC);
    memset(fft_prod, 0, sizeof(c32) * 1024 * SPC);
    memset(fft_sum , 0, sizeof(c32) * 1024 * SPC);
    for (int windex = center - window /2; windex < center + window/2; windex++) {
      for (int j = 0; j < 1023 * SPC; j++) { // xfer to float array
        fft_data[j] = out[(1023 * SPC * windex) + j];
      }

      fft_c32(1024 * SPC, fft_data, true); // forward FFT

      for (int k = 0; k < 1024 * SPC; k++) { // pt-wise * with conj of replica
        fft_sum[k] = add(fft_sum[k], mult(fft_data[k], get_conj(replica[k])));
      }
    } // for windex 

    fft_c32(1024 * SPC, fft_sum, false); // IFFT 
    
    if (0) {//center == 30) {
      FILE* fp_out = NULL; //output file
      errno_t er = fopen_s(&fp_out, "C:/Python/nci_sum4.csv", "w");
      for (int m = 0; m < 1024 * SPC; m++) {
        double magn = mag(fft_sum[m]);
        fprintf(fp_out, "%d, %f\n", m, magn);
      }
      fclose(fp_out);
    }

    float max_coh = 0; int pos_coh = 0;
    for (int m = 0; m < 1024 * SPC; m++) {
      float mag_coh = mag(fft_sum[m]);
      //fprintf(fp_out, "%d, %f \n", m, mag);
      if (mag_coh > max_coh) { max_coh = mag_coh; pos_coh = m; }
    }
   
    stat_add(&stat, max_coh);
    float mean2 = stat_mean(&stat);
    if (max_coh < min_val && max_coh < mean2 - 1000 && max_coh < 2500) { // for 4 SPC
    //if (max_coh < min_val && max_coh < mean2 - 200 && max_coh < 1200) { // for 2 SPC
    //if (max_coh < min_val && max_coh < mean2 - 10 && max_coh < 600) { // sigma = 2.0 for 1 SPC (odd missed detect)
      min_val = max_coh;
      min_idx = center;
    }
    if ((center == min_idx + 1) && (max_coh > min_val) ) {
      locations[loc_cnt++] = min_idx; // empirically the wider the window the earlier the bit transition appears
      min_val = 1e5;
      min_idx = 0;
    }
    printf("center=%d max=%6.1f pos=%d mean=%6.1f\n", center, max_coh, pos_coh, mean2);
  } // for center

  for (int i = 0; i < loc_cnt; i++) {
    printf("Bit transition at %d ms \n", locations[i]);
  }
  printf("BTs: %d Random number: %d\n", loc_cnt , rand());
  free(out); free(fft_data); free(fft_prod); free(fft_sum); free(replica);
}

void find_prn_shift2( c32* prnA,  c32* prnA_Shift, const int size)
{
  int best_k = -1;
  int best_val = -1e6;
#define FFT_SZ 512 * 2
#define E5A_SIZE 330 * 2

  c32 cprnA[FFT_SZ] = { 0 };
  c32 cprnA_Shift[FFT_SZ] = { 0 };
  c32 cprod[FFT_SZ] = { 0 };
  //memcpy(cprnA, prnA, sizeof(c32) * size); // tried the padded FFT (at least twice the length) worked up to 165 of 330
  //memcpy(cprnA_Shift, prnA_Shift, sizeof(c32) * size);

  up_sample_N_to_M(prnA, E5A_SIZE, cprnA, FFT_SZ);
  up_sample_N_to_M(prnA_Shift, E5A_SIZE, cprnA_Shift, FFT_SZ);

  fft_c32(FFT_SZ, cprnA, true);
  fft_c32(FFT_SZ, cprnA_Shift, true);

  for (int i = 0; i < FFT_SZ; i++) {
    cprod[i] = mult(cprnA[i], get_conj(cprnA_Shift[i]));
  }
  fft_c32(FFT_SZ, cprod, false);

  for (int k = 0; k < FFT_SZ; ++k) {
   
    float val = mag(cprod[k]);
    //fprintf(stderr, "Shift, %d, %f\n", k, val);
    if (val > best_val) {
      best_val = val;
      best_k = k;
    }
  }
  float denom = FFT_SZ;
  float numer = E5A_SIZE;
  printf("Best2 corr2 value %d at shift %d round %d \n", best_val, best_k,(int) round(best_k *numer/ denom));
}

void test_quasi_pilot_330Up() {
  srand((unsigned int)time(NULL)); // randomise seed
  // Test the quasi pilot generation
  int min_idx = 0; int loc_cnt = 0;
  float min_val = 1e5;
#define FFT_QP_SIZE 512
  float chipping_rate = 5.115e6; // chips per sec
  int locations[50] = { -1 };// { 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280 };
  int window = 26; // 2 * window ms either side of center (window>=6 does not work)
  int nci = 300;
#define SPC 1 // samples per chip
  float up_ratio = float(E5_QP_CODE_LEN) / float(FFT_QP_SIZE);
  int len = E5_QP_CODE_LEN * SPC * nci; // 4 samples per chip and 100 ms
  int c_phase = 29;// 166;// 329; // which chip to set the code phase to
  int prn1 = 4, prn2 = 8;
  float dop1 = 2000, dop2 = 2000;
  float dop_error = 1500;// 10; // full 2*250 Hz error in wipeoff
  float dop_err_rate = 0.6;// 0.6;// 0.6;//Hz per ms
  float sigma = 11.652;// 3.5;// 3.5;// 3.5; // noise level
  int num_errors = 0;
  c32* out = (c32*)malloc(len * sizeof(c32));
  if (out == NULL) { fprintf(stderr, "Memory allocation failed for 100 ms I&Q array.\n"); return; }
  int* prn_c1 = (int*)malloc(sizeof(int) * E5_QP_CODE_LEN * SPC);
  int* prn_c2 = (int*)malloc(sizeof(int) * E5_QP_CODE_LEN * SPC);
  c32* replica = (c32*)malloc(sizeof(c32) * E5_QP_CODE_LEN * SPC);
  float* sum_mag = (float*)malloc(sizeof(float) * FFT_QP_SIZE * SPC);
  memset(prn_c1,  0, sizeof(int) * E5_QP_CODE_LEN * SPC);
  memset(prn_c2,  0, sizeof(int) * E5_QP_CODE_LEN * SPC);
  memset(replica, 0, sizeof(c32) * E5_QP_CODE_LEN * SPC);
  getE5_QPCode(E5_QP_CODE_LEN, SPC, prn1, prn_c1);
  getE5_QPCode(E5_QP_CODE_LEN, SPC, prn2, prn_c2);

  make_replica(prn_c1, replica, dop1 + dop_error, E5_QP_CODE_LEN * SPC, chipping_rate * SPC);
  rotate_fwd(prn_c1, E5_QP_CODE_LEN * SPC, c_phase); // now advance code-phase  
  int sign2 = 1; // sign applied a posteriori after finding BTT
  stat_s stat;
  stat_init(&stat); // moving average of peak values window size = 3
  for (int i = 0; i < nci; i++) {
    for (int j = 0; j < 50; j++) {
      if (locations[j] == i) { sign2 *= -1; break; } // change sign at the bit transitions
    }
    // offset doppler by 250 Hz and add a residual doppler ramp of 0.1 Hz per ms
    mix_two_prns_oversampled_per_prn(prn_c1, prn_c2, dop1 + i * dop_err_rate, dop2 - i * dop_err_rate,
      PI / 2, 0, &out[E5_QP_CODE_LEN * SPC * i], E5_QP_CODE_LEN * SPC, chipping_rate * SPC, sigma, sign2); // was 2.31 for -128.5 dBm 3.1 for -131.5
  }
  free(prn_c1); free(prn_c2);

  // Compute circular correlation C_k(τ) = FFT^-1{ FFT[signal] · conj(FFT[replica]) }.
  c32* fft_repl = (c32*)malloc(sizeof(c32) * FFT_QP_SIZE * SPC);
  c32* fft_data = (c32*)malloc(sizeof(c32) * FFT_QP_SIZE * SPC);
  c32* fft_sum  = (c32*)malloc(sizeof(c32) * FFT_QP_SIZE * SPC);
  memset(fft_repl, 0, sizeof(c32) * FFT_QP_SIZE * SPC);
  if (fft_data == NULL || fft_repl == NULL) { printf("Error allocating fft_data or fft_prod\n"); return; }

  up_sample_N_to_M(replica, E5_QP_CODE_LEN * SPC, fft_repl, FFT_QP_SIZE * SPC);
  free(replica);
  fft_c32(FFT_QP_SIZE * SPC, fft_repl, true);
  auto start = std::chrono::high_resolution_clock::now();////////////////////////////////////////
  for (int center = window / 2; center <= nci - window / 2; center++) {
    memset(fft_sum, 0, sizeof(c32) * FFT_QP_SIZE * SPC);
    memset(sum_mag, 0, sizeof(float) * FFT_QP_SIZE * SPC);
    for (int windex = center - window / 2; windex < center + window / 2; windex++) {
      memset(fft_data, 0, sizeof(c32) * FFT_QP_SIZE * SPC);
      up_sample_N_to_M(&out[E5_QP_CODE_LEN * SPC * windex], E5_QP_CODE_LEN * SPC, fft_data, FFT_QP_SIZE * SPC);
      fft_c32(FFT_QP_SIZE * SPC, fft_data, true); // forward FFT

      for (int k = 0; k < FFT_QP_SIZE * SPC; k++) { // accumulate pt-wise * with conj of replica
        fft_sum[k] = add(fft_sum[k], mult(fft_data[k], get_conj(fft_repl[k])));
      }

      fft_c32(FFT_QP_SIZE * SPC, fft_sum, false); // IFFT 
      for (int k = 0; k < FFT_QP_SIZE * SPC; k++) { sum_mag[k] += mag(fft_sum[k]); }
    } // for windex 

    //fft_c32(FFT_QP_SIZE * SPC, fft_sum, false); // IFFT 
    for (int k = 0; k < FFT_QP_SIZE * SPC; k++) { sum_mag[k] += mag(fft_sum[k]); }

    if (false) { //center == 90) {
      FILE* fp_out = NULL; //output file
      errno_t er = fopen_s(&fp_out, "C:/Python/nci_sum4.csv", "w");
      for (int m = 0; m < FFT_QP_SIZE * SPC; m++) {
        double magn = sum_mag[m];// mag(fft_sum[m]);
        fprintf(fp_out, "%d, %f\n", m, magn);
      }
      fclose(fp_out);
    }

    float max_coh = 0; int pos_coh = 0;
    // E5_QP_CODE_LEN FFT_QP_SIZE
    for (int m = 0; m < FFT_QP_SIZE * SPC; m++) {
      float mag_coh = sum_mag[m];// mag(fft_sum[m]);
      if (mag_coh > max_coh) { max_coh = mag_coh; pos_coh = m; }
    }

    if ((int)round(pos_coh * up_ratio) != c_phase) {
      num_errors++;
    }

    stat_add(&stat, max_coh);
    float mean2 = stat_mean(&stat);
    if (max_coh < min_val && max_coh < mean2 - 1000 && max_coh < 2500) { // for 4 SPC
      min_val = max_coh;
      min_idx = center;
    }
    if ((center == min_idx + 1) && (max_coh > min_val)) {
      locations[loc_cnt++] = min_idx; // empirically the wider the window the earlier the bit transition appears
      min_val = 1e5;
      min_idx = 0;
    }
    printf("center=%03d max=%6.1f pos=%d mean=%6.1f corrPos=%d\n", center, max_coh, pos_coh, mean2, (int)round(pos_coh * up_ratio));// E5_QP_CODE_LEN / FFT_QP_SIZE));
  } // for center
  auto end = std::chrono::high_resolution_clock::now();////////////////////////////////////////
  std::chrono::duration<double> duration = end - start;
  printf("Processing time for quasi pilot 330Up ms: %f seconds\n", duration.count());
  printf("num errors %d \n", num_errors);
  for (int i = 0; i < loc_cnt; i++) {
    printf("Bit transition at %d ms \n", locations[i]);
  }
  printf("BTs: %d Random number: %d\n", loc_cnt, rand());
  free(out); free(fft_data); free(fft_sum); free(fft_repl); free(sum_mag);
}


void test_quasi_pilot_330(results_s* results) {
  //srand((unsigned int)time(NULL)); // randomise seed
  auto now = std::chrono::system_clock::now();
  // Convert the time point to duration since epoch
  auto duration_since_epoch = now.time_since_epoch();
  srand((unsigned int)std::chrono::duration_cast<std::chrono::milliseconds>(duration_since_epoch).count());
  // Test the quasi pilot generation
  int min_idx = 0; int loc_cnt = 0;
  float min_val = 1e5;
  // uncomment "//shift" for type 2) method
#define FFT_QP_SIZE 512 * 2
  float chipping_rate = 5.115e6; // chips per sec
  int locations[50] = { -1 };// { 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280 };
  int window = 26; //  window/2 ms either side of center 
  int nci = 500;
#define SPC 1 // samples per chip
  int   len = E5_QP_CODE_LEN * SPC * nci; // 4 samples per chip and 100 ms
  int   c_phase = 329;// 1;// 329; // which chip to set the code phase to
  int   prn1 = 14, prn2 = 18;
  float dop1 = 2000, dop2 = 2000;
  float dop_error = 1500;// 10; // full 2*250 Hz error in wipeoff
  float dop_err_rate = 0.6;// 0.6;// 0.6;//Hz per ms
  float N0 = -174.0 + 2.5;// dBm/Hz thermal noise ktb density + noise figure
  float cn0 = -128.75 - N0; // dBm/Hz pgae 14 Galileo_OS_SIS_ICD_v2.2
  float sigma = sqrt((1.0 * chipping_rate * SPC) / (2.0f* pow(10.0,cn0/10.0)));
  //float sigma = 12.7;// 15.98;// 3.5;// 3.5; // noise level 15->6*31
  int   num_errors = 0; int num_tries = 0;
  c32*  out    = (c32*)malloc(len * sizeof(c32));
  if (out == NULL) { fprintf(stderr, "Memory allocation failed for 100 ms I&Q array.\n"); return; }
  int* prn_c1  = (int*)malloc(sizeof(int) * E5_QP_CODE_LEN * SPC);
  int* prn_c2  = (int*)malloc(sizeof(int) * E5_QP_CODE_LEN * SPC);
  int* prn_c3  = (int*)malloc(sizeof(int) * E5_QP_CODE_LEN * SPC);
  c32* replica = (c32*)malloc(sizeof(c32) * E5_QP_CODE_LEN * SPC);
  float* mag_sum = (float*)malloc(sizeof(float) * FFT_QP_SIZE * SPC);
  memset(prn_c1 , 0, sizeof(int) * E5_QP_CODE_LEN * SPC);
  memset(prn_c2 , 0, sizeof(int) * E5_QP_CODE_LEN * SPC);
  memset(replica, 0, sizeof(c32) * E5_QP_CODE_LEN * SPC);
  getE5_QPCode(E5_QP_CODE_LEN, SPC, prn1, prn_c1);
  getE5_QPCode(E5_QP_CODE_LEN, SPC, prn2, prn_c2);

  memcpy(prn_c3, prn_c1, sizeof(int) * E5_QP_CODE_LEN * SPC); // unshifted version for later
  //rotate_fwd(prn_c3, E5_QP_CODE_LEN * SPC, (c_phase > 165) ? 165 : 0);

  make_replica(prn_c3, replica, dop1 + dop_error, E5_QP_CODE_LEN * SPC, chipping_rate * SPC);
  rotate_fwd(prn_c1, E5_QP_CODE_LEN * SPC, c_phase); // now advance code-phase  
  free(prn_c3);
  int sign2 = 1; // sign applied a posteriori after finding BTT
  stat_s stat;
  stat_init(&stat); // moving average of peak values window size = 3
  for (int i = 0; i < nci; i++) {
    for (int j = 0; j < 50; j++) {
      if (locations[j] == i) { sign2 *= -1; break; } // change sign at the bit transitions
    }
    // offset doppler by 250 Hz and add a residual doppler ramp of 0.1 Hz per ms
    mix_two_prns_oversampled_per_prn(prn_c1, prn_c2, dop1 + i * dop_err_rate, dop2 - i * dop_err_rate, 
      PI / 2, 0,&out[E5_QP_CODE_LEN * SPC * i], E5_QP_CODE_LEN * SPC, chipping_rate * SPC, sigma, sign2); // was 2.31 for -128.5 dBm 3.1 for -131.5
  }
  free(prn_c1); free(prn_c2);

  // Compute circular correlation C_k(τ) = FFT^-1{ FFT[signal] · conj(FFT[replica]) }.
  c32* fft_repl = (c32*)malloc(sizeof(c32) * FFT_QP_SIZE * SPC);
  c32* fft_data = (c32*)malloc(sizeof(c32) * FFT_QP_SIZE * SPC);
  c32* fft_sum  = (c32*)malloc(sizeof(c32) * FFT_QP_SIZE * SPC);
  memset(fft_repl, 0, sizeof(c32) * FFT_QP_SIZE * SPC);
  if (fft_data == NULL || fft_repl == NULL) { printf("Error allocating fft_data or fft_prod\n"); return; }

  //up_sample_N_to_M(replica, E5_QP_CODE_LEN * SPC, fft_repl, FFT_QP_SIZE * SPC);
  memcpy(fft_repl, replica, sizeof(c32) * E5_QP_CODE_LEN * SPC);
  // this is not good memcpy(&fft_repl[E5_QP_CODE_LEN * SPC], replica, sizeof(c32) * (FFT_QP_SIZE - E5_QP_CODE_LEN) * SPC);
  free(replica);
  fft_c32(FFT_QP_SIZE * SPC, fft_repl, true);
  auto start = std::chrono::high_resolution_clock::now();////////////////////////////////////////
  //for (int center = window / 2; center <= nci - window / 2; center++) {
  for (int center = window / 2; center <= nci - window; center++) {
    memset(fft_sum, 0, sizeof(c32) * FFT_QP_SIZE * SPC);
    memset(mag_sum, 0, sizeof(float) * FFT_QP_SIZE * SPC);
    for (int windex = center - window / 2; windex < center + window / 2; windex +=1) {
      memset(fft_data, 0, sizeof(c32) * FFT_QP_SIZE * SPC);
      memcpy(fft_data, &out[E5_QP_CODE_LEN * SPC * windex], sizeof(c32) * SPC * (E5_QP_CODE_LEN));
      memcpy(&fft_data[E5_QP_CODE_LEN * SPC], &out[E5_QP_CODE_LEN * SPC * windex], sizeof(c32) * SPC * (FFT_QP_SIZE - E5_QP_CODE_LEN)); 
      //up_sample_N_to_M(&out[E5_QP_CODE_LEN * SPC * windex], E5_QP_CODE_LEN * SPC, fft_data, FFT_QP_SIZE * SPC);
      fft_c32(FFT_QP_SIZE * SPC, fft_data, true); // forward FFT

      for (int k = 0; k < FFT_QP_SIZE * SPC; k++) { // accumulate pt-wise * with conj of replica
        fft_sum[k] = add(fft_sum[k], mult(fft_data[k], get_conj(fft_repl[k])));
        //fft_sum[k] = mult(fft_data[k], get_conj(fft_repl[k]));
      }

      //fft_c32(FFT_QP_SIZE * SPC, fft_sum, false); // IFFT
      //for (int k = 0; k < FFT_QP_SIZE * SPC; k++) { mag_sum[k] += mag(fft_sum[k]); }
    } // for windex 

    // used to have the IFFT here
    fft_c32(FFT_QP_SIZE * SPC, fft_sum, false); // IFFT
    for (int k = 0; k < FFT_QP_SIZE * SPC; k++) { mag_sum[k] += mag(fft_sum[k]); }

    if (false) {//center == 17) {
      FILE* fp_out = NULL; //output file
      errno_t er = fopen_s(&fp_out, "C:/Python/nci_sum4.csv", "w");
      for (int m = 0; m < FFT_QP_SIZE * SPC; m++) {
        fprintf(fp_out, "%d, %f\n", m, mag_sum[m]);
      }
      fclose(fp_out);
    }

    float max_coh = 0; int pos_coh = 0;
    int idx1 = -1, idx2 = -1;
    float v1 = -1e6, v2 = -1e6;
    // E5_QP_CODE_LEN FFT_QP_SIZE
    for (int m = 0; m < E5_QP_CODE_LEN * SPC; m++) {
      float mag_coh = mag_sum[m];// mag(fft_sum[m]);
      if (mag_coh > max_coh) { max_coh = mag_coh; pos_coh = m; }

      if (mag_coh >= v1) {
        // Promote current best to second, insert new best
        if (idx1 != -1) { v2 = v1; idx2 = idx1;}
        v1 = mag_coh; idx1 = m;
      }
      else if (mag_coh > v2) {
        // Update second best if it doesn't collide with best 
        if (mag_coh != v1) { v2 = mag_coh;  idx2 = m; }
      }
    }

    num_tries++;
    //pos_coh = (c_phase > 165) ? pos_coh + 165 : pos_coh;
    if (pos_coh != c_phase) {
      num_errors++;
    }

    stat_add(&stat, max_coh);
    float mean2 = stat_mean(&stat);
    if (max_coh < min_val && max_coh < mean2 - 1000 && max_coh < 2500) { // for 4 SPC
      min_val = max_coh;
      min_idx = center;
    }
    if ((center == min_idx + 1) && (max_coh > min_val)) {
      locations[loc_cnt++] = min_idx; // empirically the wider the window the earlier the bit transition appears
      min_val = 1e5;
      min_idx = 0;
    }
    //printf("center=%03d max=%6.1f pos=%d mean=%6.1f ratio=%3.1f sep=%d\n", center, max_coh, pos_coh, mean2, (v1 / v2), (int)abs(idx1 - idx2));// E5_QP_CODE_LEN / FFT_QP_SIZE));
  } // for center
  auto end = std::chrono::high_resolution_clock::now();////////////////////////////////////////
  std::chrono::duration<double> duration = end - start;
  printf("Processing time for quasi pilot 330 ms: %f seconds\n", duration.count());
  printf("cn0=%4.2f sigma=%6.3f num errors pilot_330: %d out of %d\n",cn0,sigma, num_errors, num_tries);
  results->num_errors += num_errors;
  results->num_trials += num_tries;
  for (int i = 0; i < loc_cnt; i++) {
    printf("Bit transition at %d ms \n", locations[i]);
  }
  //printf("BTs: %d Random number: %d\n", loc_cnt, rand());
  free(out); free(fft_data); free(fft_sum); free(fft_repl); free(mag_sum);
}

void test_quasi_pilot2() {
  srand((unsigned int)time(NULL)); // randomise seed
  // Test the quasi pilot generation

  int min_idx = 0;
  int loc_cnt = 0;
  float min_val = 1e5;

  int locations[50] = { 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280 };
  int window = 2; // 2 * window ms either side [window works from 2 - 37]
  int nci = 300;
#define SPC 4 // samples per chip
  int len = 1023 * SPC * nci; // 4 samples per chip and 100 ms
  int c_phase = 3333; // which chip to set the code phase to
  int prn1 = 4, prn2 = 8;
  float dop1 = 2000, dop2 = -3000;
  float dop_error = 250;// 10; // full 250 Hz error in wipeoff
  float dop_err_rate = 0.6;// 0.6;//Hz per ms
  float sigma = 3.5;// 3.5; // noise level
  c32* out = (c32*)malloc(len * sizeof(c32));
  if (out == NULL) {
    fprintf(stderr, "Memory allocation failed for 100 ms I&Q array.\n");
    return;
  }
  int prn_c1[1023 * SPC], prn_c2[1023 * SPC];
  c32* replica = (c32*)malloc(sizeof(c32) * 1024 * SPC);
  memset(replica, 0, sizeof(c32) * 1024 * SPC);
  getCode(1023, SPC, prn1, prn_c1);
  getCode(1023, SPC, prn2, prn_c2);
  make_replica(prn_c1, replica, dop1 + dop_error, 1023 * SPC, 1.023e6 * SPC);
  // now advance code-phase
  rotate_fwd(prn_c1, 1023 * SPC, c_phase); // code phase 1/4 way
  int sign2 = 1; // sign applied a posteriori after finding BTT
  stat_s stat;
  stat_init(&stat); // moving average of peak values window size = 3
  for (int i = 0; i < nci; i++) {
    for (int j = 0; j < 50; j++) {
      if (locations[j] == i) { sign2 *= -1; break; } // change sign at the bit transitions
    }
    // offset doppler by 250 Hz and add a residual doppler ramp of 0.1 Hz per ms
    mix_two_prns_oversampled_per_prn(prn_c1, prn_c2, dop1 + dop_error + i * dop_err_rate, dop2 - i * dop_err_rate, PI / 2, 0,
      &out[1023 * SPC * i], 1023 * SPC, 1.023e6 * SPC, sigma, sign2); // was 2.31 for -128.5 dBm 3.1 for -131.5
  }

  // use FFTs 
  fft_c32(1024 * SPC, replica, true);

  // Compute circular correlation C_k(τ) = FFT^-1{ FFT[x_d,k] · conj(FFT[code]) }.
  c32* fft_data = (c32*)malloc(sizeof(c32) * 1024 * SPC);
  c32* fft_sum = (c32*)malloc(sizeof(c32) * 1024 * SPC);
  if (fft_data == NULL || fft_sum == NULL) { printf("Error allocating fft_data or fft_prod\n"); return; }
  for (int center = window / 2; center <= nci - window / 2; center++) {
    // use linearity to sum the ncis coherently in moving window
    memset(fft_data, 0, sizeof(c32) * 1024 * SPC);
    memset(fft_sum, 0, sizeof(c32) * 1024 * SPC);
    for (int windex = center - window / 2; windex < center + window / 2; windex++) {
      for (int j = 0; j < 1023 * SPC; j++) { // xfer to float array
        // additive coherent sum (would it work with real data?)
        fft_data[j] = add(fft_data[j], out[(1023 * SPC * windex) + j]);
      }
    } // for windex 

    fft_c32(1024 * SPC, fft_data, true); // forward FFT

    for (int k = 0; k < 1024 * SPC; k++) { // pt-wise * with conj of replica
      fft_sum[k] = mult(fft_data[k], get_conj(replica[k]));
    }
    
    fft_c32(1024 * SPC, fft_sum, false); // IFFT

    if (0) {//center == 170) {
      FILE* fp_out = NULL; //output file
      errno_t er = fopen_s(&fp_out, "C:/Python/nci_sum3.csv", "w");
      for (int k = 0; k < 1024 * SPC; k++) {
        float mag1 = mag(fft_sum[k]);
        fprintf(fp_out, "%d, %f \n", k, mag1);
      }
      fclose(fp_out);
    }

    float max_nci = 0; int pos_nci = 0;
    for (int m = 0; m < 1024 * SPC; m++) {
      float mag_nci = mag(fft_sum[m]);
      if (mag_nci > max_nci) { max_nci = mag_nci; pos_nci = m; }
    }
    
    //printf("nci pos =%d max=%5.0f \n", pos_nci, max_nci);
    stat_add(&stat, (max_nci));

    float mean2 = stat_mean(&stat);
    //if (max_nci < min_val && max_nci < mean2 - 20) {
    if (max_nci < min_val && max_nci < mean2 - 900) {
      min_val = max_nci;
      min_idx = center;
    }
    if (loc_cnt > 49) { loc_cnt = 49; }
    if ((center == min_idx + 1) && (max_nci > min_val)) {
      locations[loc_cnt++] = min_idx; // empirically the wider the window the earlier the bit transition appears
      min_val = 1e5;
      min_idx = 0;
    }
    printf("center=%d max=%5.0f pos=%d mean=%2.6f\n", center, max_nci, pos_nci, mean2);
  } // for center

  for (int i = 0; i < loc_cnt; i++) {
    printf("Bit transition at %d ms \n", locations[i]);
  }
  printf("Random number: %d\n", rand());
  free(out); free(fft_data); free(fft_sum); free(replica);
}

void test_quasi_pilot3() {
  srand((unsigned int)time(NULL)); // randomise seed
  // Test the quasi pilot generation
  int nci = 300;
  int len = 1023 * 2 * nci; // 2 samples per chip and 100 ms
  int min_idx = 0;
  int loc_cnt = 19;
  float min_val = 1e5;
  c32* out = (c32*)malloc(len * sizeof(c32));
  if (out == NULL) {
    fprintf(stderr, "Memory allocation failed for 100 ms I&Q array.\n");
    return;
  }
  // was {-1}
  int locations[50] =  { 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 
    190, 200, 210, 220, 230, 240, 250, 260, 270, 280, 290 };
  int window = 3; // 3 ms either side
  
  int code_phase = 555; // which chip to set the code phase to
  float carrier_phase = 90; // deg
  int prn1 = 4, prn2 = 8;
  float dop1 = 2000, dop2 = -3000;
  float dop_error = 20; // full 250 Hz error in wipeoff
  float sigma = 2.00;// 3.4; // noise level
  int prn_c1[1023 * 2], prn_c2[1023 * 2];
  
  getCode(1023, 2, prn1, prn_c1);
  getCode(1023, 2, prn2, prn_c2);
  rotate_fwd(prn_c1, 1023 * 2, code_phase); // code phase 1/4 way

  int sign2 = 1; // sign applied a posteriori after finding BTT
  stat_s stat;
  stat_init(&stat); // moving average of peak values window size = 3
  int num_flips = 0;
  for (int i = 0; i < nci; i++) {
    for (int j = 0; j < 50; j++) {
      if (locations[j] == i) { sign2 *= -1; num_flips++;  break; } // change sign at the bit transitions
    }
    // offset doppler by 250 Hz and add a residual doppler ramp of 0.1 Hz per ms
    mix_two_prns_oversampled_per_prn(prn_c1, prn_c2, dop1 + dop_error + i * 0.0, dop2 - i * 0.0, carrier_phase, 0,
      &out[1023 * 2 * i], 1023 * 2, 1.023e6 * 2, sigma, sign2); // was 2.31 for -128.5 dBm 3.1 for -131.5
  }

  getCode(1023, 2, prn1, prn_c1);
  float dop_out, pwr_out;
  int code_out;
  for (int windex = 0; windex < nci; windex++) {
    float phase = recover_prn_phase_deg_with_doppler(&out[(1023 * 2 * windex)], 1023 * 2, 
      prn_c1, code_phase - 9, dop1, 2.046e6, &code_out, &dop_out, &pwr_out);
    printf("window=%d best_off=%d dop_out=%f pwr_out=%f phase=%f\n", windex, code_out, dop_out, pwr_out, phase);
  } // for center
  
  for (int i = 0; i < num_flips -1; i++) {
    printf("Bit transition at %d ms \n", locations[i]);
  }
  printf("Random number: %d\n", rand());
  free(out);
}

/**
 * Main for testing and developing under Visual Studio 2022
 */
int main(int argc,char* argv[])
{
  float T = (float) (1.0 / 500.0); // Period for the test sine wave (1/500 seconds = 2ms period)
  const int multK = 1; // multiple of 1K size can be used for 1,2,4,8,or 16 (has not been tested for less than 1K)
  float F1 = 50.0f; // Frequency of the first sine wave
  float F2 = 80.0f;// 80.0f; // Frequency of the second sine wave
  float F3 = 155.0f;// 155.0f; // Frequency of the third sine wave (not used in this example, but can be added for more complexity)
  float F4 = 230.0f; // Frequency of the third sine wave (not used in this example, but can be added for more complexity)
  
// File lengths are too long to do float and radix4 and radix2 serially: use compile flags
//#define DO_FLOAT
//#define DO_Q15_RADIX2
//#define DO_Q15_RADIX4
//#define DO_Q31_RADIX2 
//#define DO_Q31_RADIX4

  if (0) {
    // set to 1 above if need to recalculate some of the twiddleCoefs
    // e.g. twiddleCoef_8192_q15 in arm_common_tables.c
    // minor code changes would be necessary to support q31_t once compile issues are fixed
    twiddleCoefCalculator();
    return 0;
  }
  if (0) {
    // set to 1 above if need to recalculate some of the armBitRevTables
    // e.g. armBitRevTable[1024 * 4] in arm_common_tables.c which is on one big line!
    armBitRevTableCalculator();
    return 0;
  }

  if (0) {
    sim_E5A();
    return 0;
  }

  if (1) {
    read_E5A((char*)"C:/work/Baseband/TestData/E5/t14/G_2024_10_21_22_29_43.047.csv");
    //read_E5A((char*)"C:/work/Baseband/TestData/E5/t14/G_2024_10_21_22_32_30.500_resampled_16368Hz.csv");
    //read_E5A((char*)"C:/work/Baseband/TestData/E5/t14/G_2024_10_21_22_29_43.047_resampled_16368Hz.csv");
    return 0;
  }

  if (0) {
    read_L1((char*)"C:/Python/out-1bit_1spc_1bit.csv");// "C:/work/Baseband/TestData/E5/t14/G_2024_10_21_22_29_43.047.csv");
  }

  if (0) {
    read_ors((char*)"C:/work/Baseband/TestData/100ms/bw25/G_2025_09_03_23_04_45.ors");
    //read_ors((char*)"C:/work/Baseband/TestData/100ms/bw25/G_2025_09_03_23_22_41.ors");
    //read_ors((char*)"C:/work/Baseband/TestData/100ms/bw25/G_2025_09_03_23_04_45.ors");
    //read_ors((char*)"C:/work/Baseband/TestData/G_2025_06_05_22_11_26.ors");
    return 0;
  }

  if (0) {
    //test_quasi_pilot_330Up();
    results_s results = {0};
    for (int i = 0; i < 50; i++) {
      test_quasi_pilot_330(&results);
    }
    printf("Total tries: %d, Total errors: %d, Avg errors per try: %f\n", results.num_trials, results.num_errors,
      (float)results.num_errors / (float)results.num_trials);
    //test_quasi_pilot2();
    return 0;
  }

#ifdef DO_FLOAT
  float test[1024 * multK * 2];

  memset(test, 0, sizeof(test));
  int code[1024 * 1] = { 1 };
  getCode(1024, 1, 23, code);
  for (int i = 0; i < 1024 * multK; i++) {
    double x = i * T;
    test[i * 2 + 0] =  cos(F1 * 2.0 * PI * x) + code[i * 2 + 0] * cos(F2 * 2.0 * PI * x);// 0.5 * signalf(x, F1, F2, F3, F4);
    test[i * 2 + 1] =  sin(F1 * 2.0 * PI * x) + code[i * 2 + 1] * sin(F2 * 2.0 * PI * x);// cos(F1 * 2.0 * PI * x); // set the imaginary part to 0
  }
  
  arm_cfft_radix2_instance_f32 s;
  std::cout << "************ float32 FFT ************************ \n";
 
  arm_cfft_radix2_init_f32(&s, 1024 * multK, 0, 1); // Initialize the CFFT instance for 8-point FFT
  arm_cfft_radix2_f32(&s, test);
  for (int i = 0; i < 1024 * multK / 2; i++) {
    std::cout << "float [" << i << "] , " << (double(i)/double(T* multK*1024)) << " , " << sqrt(test[2 * i] * test[2 * i] + test[2 * i + 1] * test[2 * i + 1]) << "\n";
  }
#endif // DO_FLOAT

  
  q15_t src_q15[1024 * 2 * multK];
  memset(src_q15, 0, sizeof(src_q15));

#ifdef DO_Q15_RADIX4
  // try the Q15 for radix4 
  arm_cfft_radix4_instance_q15 s_q15_radix4;
  std::cout << "************ FFT q15_t Radix4 ************************ \n";
  for (int i = 0; i < 1024 * multK; i++) {
    double x = i * T;
    src_q15[i * 2] = (q15_t)(f2q15(signalf(x, F1, F2, F3, F4)));
    src_q15[i * 2 + 1] = 0; // set the imaginary part to 0
  }
  arm_cfft_radix4_init_q15(&s_q15_radix4, 1024 * multK, 0, 1); // Initialize the CFFT instance for 8-point FFT
  arm_cfft_radix4_q15(&s_q15_radix4, src_q15); 
  for (int i = 0; i < 1024 * multK / 2; i++) {
    std::cout << "rad4_q15[" << i << "] , " << (double(i)/double(T* multK*1024)) << " , " << sqrt(src_q15[2 * i] * src_q15[2 * i] + src_q15[2 * i + 1] * src_q15[2 * i + 1]) << "\n";
  }
#endif // DO_Q15_RADIX4

#ifdef DO_Q15_RADIX2 
  arm_cfft_radix2_instance_q15 s_q15;  
  std::cout << "************ FFT q15_t Radix2 ************************ \n";
  for (int i = 0; i < 1024 * multK; i++) {
    double x = i * T;
    src_q15[i*2] = (q15_t)f2q15(signalf(x, F1, F2, F3, F4));
    src_q15[i*2+1] = 0; // set the imaginary part to 0
  }
 
  arm_cfft_radix2_init_q15(&s_q15, 1024 * multK, 0, 1); // Initialize the CFFT instance 
  arm_cfft_radix2_q15(&s_q15, src_q15); // fft done in place result is symmetric about 1/2 point

  for (int i = 0; i < 1024 * multK / 2; i++) {
    std::cout << "rad2_q15[" << i << "] , " << (double(i)/double(T* multK*1024)) << " , " << sqrt(src_q15[2 * i] * src_q15[2 * i] + src_q15[2 * i + 1] * src_q15[2 * i + 1]) << "\n";
  }
#endif // DO_Q15_RADIX2

#if defined(DO_Q31_RADIX2) || defined(DO_Q31_RADIX4)
  q31_t src_q31[1024 * 2 * multK];
  memset(src_q31, 0, sizeof(src_q31));
#endif

#ifdef DO_Q31_RADIX2
  // now try the Q31 version warning does not compile right
  arm_cfft_radix2_instance_q31 s_q31;

  std::cout << "************ FFT q31_t ************************ \n";
  for (int i = 0; i < 1024 * multK; i++) {
    double x = i * T;
    src_q31[i*2] = (q31_t)f2q31(signalf(x, F1, F2, F3, F4));
    src_q31[i * 2 + 1] = 0;
  }
  arm_cfft_radix2_init_q31(&s_q31, 1024 * multK, 0, 1); // Initialize the CFFT instance for 8-point FFT
  arm_cfft_radix2_q31(&s_q31, src_q31);
  for (int i = 0; i < 1024 * multK / 2; i++) {
    std::cout << "rad2_q31[" << i << "] , " << (double(i) / double(T * multK * 1024)) << " , " << magnitude(src_q31[2 * i], src_q31[2 * i + 1]) << "\n";
  }
#endif //DO_Q31_RADIX2

#ifdef DO_Q31_RADIX4
  // now try the Q31 version warning does not compile right
  arm_cfft_radix4_instance_q31 s_q31_r4;
 
  std::cout << "************ FFT q31_t rad4 ************************ \n";
  for (int i = 0; i < 1024 * multK; i++) {
    double x = i * T;
    src_q31[i * 2] = (q31_t)f2q31(signalf(x, F1, F2, F3, F4));
    src_q31[i * 2 + 1] = 0;
  }
  arm_cfft_radix4_init_q31(&s_q31_r4, 1024 * multK, 0, 1); // Initialize the CFFT instance for 8-point FFT
  arm_cfft_radix4_q31(&s_q31_r4, src_q31);
  for (int i = 0; i < 1024 * multK / 2; i++) {
    std::cout << "rad4_q31[" << i << "] , " << (double(i) / double(T * multK * 1024)) << " , " << magnitude(src_q31[2 * i], src_q31[2 * i + 1]) << "\n";
  }
#endif //DO_Q31_RADIX2

  return 0;
}

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file


