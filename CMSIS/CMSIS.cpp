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
#include <string>
#include <chrono> // for profiling
#include "transform_functions.h"
#include "gnss_codes.h"
#include "up_sample.h"
#include "msb_funcs.h"

#include "string_helpers.h"

#include <complex> // not to be confused with complex.h

#define Q31_MIN   ((int32_t)0x80000000)  // -2147483648
#define Q31_MAX   ((int32_t)0x7FFFFFFF)  //  2147483647
#define Q15_MIN   (-32768)
#define Q15_MAX   ( 32767)
#define R2D (180.0f / 3.14159265358979f)
#define D2R (3.14159265358979f / 180.0f)
#define FS (4.092e6f) // sample rate

FILE* fpResults = NULL;

typedef enum
{
  OPTION_ERROR = 0x00,
  OPTION_L1 = 0x01,
  OPTION_E1 = 0x02,
  OPTION_L5 = 0x10,
  OPTION_E5a = 0x20,
  OPTION_QP = 0x40,
  OPTION_QX = 0x80 // use only qp satellites, but not as qp
} optionSigEnum;

static int optionSig = OPTION_ERROR;
static int optionMsCoh = -1;
static int optionMsps = -1;
static int optionDUnc = -1; // doppler uncertainty in Hz

typedef struct {
  int num_errors;
  int num_trials;
  int btt_missed_detects;
  int btt_false_alarms;
  int differences[100];
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

c32 get_minus_conj(const c32 x) {
  c32 y;
  y.r = -x.r;
  y.i = +x.i;
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


double compute_gps_time(int year, int month, int day, int hour, int minute, double second)
{
  struct tm refTime, epochTime;
  time_t refTimeRaw, epochTimeRaw;
  double diff;

  refTime.tm_year = year - 1900;
  refTime.tm_mon = month - 1;
  refTime.tm_mday = day;
  refTime.tm_hour = hour;
  refTime.tm_min = minute;
  refTime.tm_sec = (int)second;
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
  diff += second - (int)second; // include fractional seconds

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

bool findMax(double input) {
  const int MAX_VALUES = 21;
  static double values[MAX_VALUES];
  static int num_values = 0;
  if (num_values < MAX_VALUES) {
    values[num_values++] = input;
  }
  else {
    for (int i = 0; i < MAX_VALUES - 1; i++) {
      values[i] = values[i + 1]; // shift values left
    }
    values[MAX_VALUES - 1] = input;
  }
  bool primed = (num_values == MAX_VALUES);
  bool is_max = false;
  if (primed) {
    //for (int i = 0; i < MAX_VALUES; i++) {  printf("%f,", values[i]); } printf("\n");
    if (values[0] < values[5] && values[5] < values[10] && values[10] > values[15] && values[15] > values[20]) {
      if (values[10] > 10000) { //was 16000
        is_max = true;
      }
    }
  }
  return is_max;
}

bool checkBT(int center, int* locations,int num_bt) {
  for (int i = 0; i < num_bt; i++) {
    if (locations[i] == center) {
      return true;
    }
  }
  return false;
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

  float code_search_span = 25; //chip span was 2
  float code_low = code_offset - code_search_span;
  float code_high = code_offset + code_search_span;
  uint8_t wrap = 0;
  if (code_low < 0) { code_low += SIZE; wrap = 1; }
  if (code_high > SIZE) { code_high -= SIZE;  wrap = 1; }

  float f_search_step_hz = 1.0; // Hz
  float f_search_span_hz = 0; // Hz span was 2
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

void read_L1_CSV(char* input) {
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
  errno_t er = fopen_s(&fp_out, "out5.csv", "w");
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
    synth_gps_prn(prn, -doppler, SIZE, repli, SPC, 0);
  }
  else { // Galileo
    synth_e1b_prn(prn, -doppler, SIZE, repli, SPC, 0);
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

// copy doppler info from .ast file corresp to input .ors
int read_assist(char* input, acq_struct* prn2acq) {
  char filename[256];
  strcpy_s(filename, input);
  *strrchr(filename, '.') = 0;
  strcat_s(filename, ".ast");
  FILE* fpAssist = NULL;
  fopen_s(&fpAssist, filename, "r");
  if (fpAssist == NULL) {
    fprintf(stderr, "Failed to open assitance file %s\n", filename);
    return 0;
  }
  int cnt = 0;
  char line[256];
  while (fgets(line, sizeof(line), fpAssist) != NULL && cnt < 40)
  {
    char c = line[0];
    if (c == 'G')
      prn2acq[cnt].constel = 1;
    else if (c == 'E')
      prn2acq[cnt].constel = 2;
    else
      continue;
    sscanf_s(line + 1, "%d %lf", &prn2acq[cnt].prn, &prn2acq[cnt].doppler);
    cnt++;
  }
  fclose(fpAssist);
  return cnt;
}


void read_ors(char* directory, char* input) {
  char filename[_MAX_PATH];
  sprintf_s(filename, "%s/%s", directory, input);

  FILE* fp_msb = NULL;
  fopen_s(&fp_msb, filename, "rb");
  if (fp_msb == NULL) {
    fprintf(stderr, "Failed to open ors file %s\n", filename); return;
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
  dtime -= 18; // leap seconds
  printf("week %d tow %f \n", (int)floor(dtime / 604800.0), dtime - floor(dtime / 604800.0) * 604800.0);
  int msec = (int)round((sec - floor(sec)) * 1000);
  uint32_t nanos = U3(&buffer[18]); // 21
  printf("%4d-%02d-%02d-%02d-%02d-%02d msec=%d nanos=%d \n", yr, mon, day, hr, min, sec, msec, nanos);
  int payload_size = bytes_to_read - hdrlen;
  printf("size (total bytes - hdrs) %d \n", payload_size);
  printf("Working on %s\n", input);

  FILE* fp_out = NULL; //output file
  char str[128];
  top2_pks peaks;

  //// Dial in the prn and doppler here ////////////////
#define SPC  4
#define SIZE 1024*SPC *4 // 16K for Galileo and 4K for GPS

  acq_struct prn2acq[30] = {0};
  //sprintf_s(filename, "%s/replay_bw25.log", directory);
  //int cnt = parse_log((char*)filename, input, prn2acq);
  //qsort(prn2acq, cnt, sizeof(acq_struct), compareConstelPRN); // sort by constel & prn for expected order
  int cnt = read_assist(filename, prn2acq);

  c32* iandq = (c32*)malloc(payload_size * sizeof(c32));
  c32* signl = (c32*)malloc(SIZE * sizeof(c32));
  c32* repli = (c32*)malloc(SIZE * sizeof(c32));
  c32* prod  = (c32*)malloc(SIZE * sizeof(c32));
  //c32* sums  = (c32*)malloc(SIZE * sizeof(c32));
  float* sums = (float*)malloc(SIZE * sizeof(c32));

  if (iandq == NULL || prod == NULL) {
    fprintf(stderr, "Memory allocation failed for q32 array.\n");
    free(iandq); free(repli); free(prod); free(signl); return;
  }

  bb_meas_t meas;
    
  bool gal16k = true;

  // use iandq as the main buffer (don't modify)
  DecodeOrsIQCplx(&buffer[hdrlen], payload_size / 2, iandq);
  free(buffer);
  static int first_pass[30] = {0}, first_pass_cnt = 0;
  memset(first_pass, 0, sizeof(int) * 30);
  // on second pass cyclically advance code phase to before pos in 1st pass
  // NB: do not cycle past the code phase otherswise errors will occur!
  for (int pass = 0; pass < 2; pass++) {
    memset(&meas, 0, sizeof(bb_meas_t)); first_pass_cnt = 0;
    for (int loop = 0; loop < cnt; loop++) {
      bool is_gps = (prn2acq[loop].constel == 1) ? true : false;
      int rep_offset = first_pass[first_pass_cnt];
      int size = is_gps ? 4092 : 16368; // the samples per epoch
      int fft_size = is_gps ? 4096 : 16384;
      memset(repli, 0, SIZE * sizeof(c32));
      memset(sums, 0, SIZE * sizeof(float));
      if (rep_offset > 200) {
        rep_offset -= 150;
      }
      int spc = SPC;
      if (!is_gps && !gal16k)
      {
        size = 4092;
        fft_size = 4096;
        spc = 1;
      }

      if (is_gps) {
        synth_gps_prn(prn2acq[loop].prn, -prn2acq[loop].doppler, size, repli, spc, rep_offset);
      }
      else { // Galileo
        synth_e1b_prn(prn2acq[loop].prn, -prn2acq[loop].doppler, size, repli, spc, rep_offset);
      }
      
      fft_c32(fft_size, repli, true); // F(repli)
      for (int j = 0; j < 4; j++) {
        memset(prod, 0, SIZE * sizeof(c32));
        memset(signl, 0, SIZE * sizeof(c32));
        if (is_gps || gal16k)
        {
          memcpy(signl, &iandq[j * size], size * sizeof(c32));
        }
        else
        {
          // decimate by 4
          c32* s = &iandq[j * size * 4];
          for (int k = 0; k < size; k++)
            signl[k] = s[k * 4];
        }
       
        fft_c32(fft_size, signl, true); // F(iandq) and prod=F(iandq) * conj(F(repli)) below
        for (int i = 0; i < fft_size; i++) { prod[i] = mult(signl[i], get_conj(repli[i])); }
        fft_c32(fft_size, prod, false); // in-place inv F(prod) 
        for (int i = 0; i < fft_size; i++) {
          sums[i] += mag(prod[i]);// add(sums[i], prod[i]);
        }
      }
      snprintf(str, sizeof(str), "%s/out%s%02d.csv", directory, (is_gps) ? "G" : "E", prn2acq[loop].prn);
      //errno_t err= fopen_s(&fp_out, str, "w");
      find_top2_peaks_real(sums, size, 3, &peaks, fp_out); if (fp_out) { fclose(fp_out); }
      if (false) {//is_gps == false) {
        c32 sum = { 0 };
        float c_span = 3;
        int8_t e1_code[16368];
        float code_out = 0, freq_out = 0, phase_out = 0;
        GetE1BCode(prn2acq[loop].prn, SPC, e1_code);
        estimate_prn_code_and_carrier(&iandq[0], size, SPC * 1.023e6f, c_span, e1_code, 1.023e6f,
          peaks.idx1, 4.f, 1, -prn2acq[loop].doppler, 3, 1, &code_out, &freq_out, &phase_out, &sum);
        printf("code_out=%f, freq_out=%f, phase_out=%f \n", code_out, freq_out, phase_out);
      }
      // compute noise stats for SNR
      double BW = 3.1623e3; // Hz
      double cn0 = compute_snr_real(sums, size, peaks) + 35;// +10 * log(BW)
      double interp = interpCodePhaseFloat(sums, size, &peaks);
      float thresh = (pass == 0) ? 1.001 : 1.3;
      float ratio = (peaks.val1 / peaks.val2);
      if (ratio > thresh) {
        meas.sats[meas.num_sat].prn = prn2acq[loop].prn;
        meas.sats[meas.num_sat].code_phase = float(interp + rep_offset) / (spc * 1023.0f);
        meas.sats[meas.num_sat].doppler = prn2acq[loop].doppler;
        meas.sats[meas.num_sat].cno = (float)cn0;
        meas.sats[meas.num_sat].constellation = is_gps ? SYS_GPS : SYS_GAL;
        printf("Acquired %s %d Doppler %f Hz Code %f [ms] Chips %f  C/N0 %f dB-Hz ratio=%f\n", is_gps ? "GPS" : "GAL",
          prn2acq[loop].prn, -prn2acq[loop].doppler, meas.sats[meas.num_sat].code_phase, meas.sats[meas.num_sat].code_phase * 4092.0, meas.sats[meas.num_sat].cno, ratio * ratio);
        meas.num_sat++;
      }
      first_pass[first_pass_cnt] = (int)interp;
      first_pass_cnt++; // must increment here for things to stay in syn
    }
    printf("Done inner loop\n");
  }

  sprintf_s(filename, "%s/%s", directory,input);
  *strrchr(filename, '.') = 0;
  strcat_s(filename,_MAX_PATH, ".msb");
  printf("writing to %s \n", filename);
  int num_bytes = write_msb(&meas, (char*)filename);
  bb_meas_t check = { 0 };
  FILE* test = NULL;
  errno_t er2 = fopen_s(&test, (char*)filename, "rb");
  uint8_t tbuff[128] = { 0 };
  fread(tbuff,1, num_bytes,test);
  read_bb_msb(tbuff, num_bytes, &check);
  fclose(test);
   
  free(iandq); free(repli); free(prod); free(signl); free(sums);
  //////////////////////////////////////////////////////
}

void sim_E5A() {
  #define FFT_SIZE 16384
  int    prn_a = 5, prn_b = 15;
  double doppler_a = 1580, doppler_b = -2580;
  int    offset = 300;// 10229;
  float  up_offset = (float)offset * 16384.0f / 10230.0f; // upsampled offset
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
  c32* up_prod  = (c32*)malloc(FFT_SIZE * sizeof(c32));
  
  up_sample_10k_to_16k(samp_a, up_samp);
  up_sample_10k_to_16k(repli_a, up_repli);
  free(samp_a); free(repli_a);

  fft_c32(FFT_SIZE, up_samp, true); 
  fft_c32(FFT_SIZE, up_repli, true);
  for (int i = 0; i < FFT_SIZE; i++) { up_prod[i] = mult(up_samp[i], get_conj(up_repli[i])); }
  fft_c32(FFT_SIZE, up_prod, false); // in-place inv F(prod)

  FILE* dbg_fp = NULL;
  fopen_s(&dbg_fp, "C:/Python/out6.csv", "w");
  top2_pks peaks = { 0 };
  find_top2_peaks_cplx(up_prod, FFT_SIZE, 3, &peaks, dbg_fp);
  double early = mag(up_prod[peaks.idx1 - 1]), prompt = peaks.val1, late = mag(up_prod[peaks.idx1 + 1]);
  double interp = InterpolateCodePhase(round(peaks.idx1 * 1.0), early * early, prompt * prompt, late * late);

  fclose(dbg_fp);
  free(up_prod); free(up_samp); free(up_repli);
  printf("max_float=%f offset=%d pos=%d interp=%f check=%f\n", peaks.val1, offset, peaks.idx1, interp, interp * (10230.0/16384.0));
}

void sim_L1() {
#define SPC  4
  bool   do_gps = false;
  float  scale = do_gps ? 1 : 4;
  int    fft_size = 1024 * scale;
  int    prn_a = 2, prn_b = 8;
  double doppler_a = 1580, doppler_b = -2580;
  int    offset = 4020;// 1022;//8000 bad is threshold for error
  int    rep_cylce = 0;// 1024 / 2;// 512;
  top2_pks peaks = { 0 };
  //printf("should be %d \n", offset);
  c32* samp_a = (c32*)malloc(fft_size * SPC * sizeof(c32));
  c32* samp_b = (c32*)malloc(fft_size * SPC * sizeof(c32));
  c32* repl_a = (c32*)malloc(fft_size * SPC * sizeof(c32));
  c32* prod   = (c32*)malloc(fft_size * SPC * sizeof(c32));

  for (int offset = 0; offset < 1024 * SPC * scale; offset += 20) {
    if (offset > 1024 * SPC * scale / 2) { rep_cylce = 1024 * SPC * scale / 2; }
    memset(samp_a, 0, fft_size * SPC * sizeof(c32));
    memset(samp_b, 0, fft_size * SPC * sizeof(c32));
    memset(repl_a, 0, fft_size * SPC * sizeof(c32));
    memset(prod, 0, fft_size * SPC * sizeof(c32));

    if (do_gps) {
      synth_gps_prn(prn_a, doppler_a, 1023 * SPC, samp_a, SPC, offset);
      synth_gps_prn(prn_b, doppler_b, 1023 * SPC, samp_b, SPC, 300);
      synth_gps_prn(prn_a, doppler_a + 10, 1023 * SPC, repl_a, SPC, 0);
    }
    else {
      synth_e1b_prn(prn_a, doppler_a, 4092 * SPC, samp_a, SPC, offset);
      synth_e1b_prn(prn_b, doppler_b, 4092 * SPC, samp_b, SPC, 300);
      synth_e1b_prn(prn_a, doppler_a + 10, 4092 * SPC, repl_a, SPC, 0);
    }
    rotate_fwd_c32(repl_a,(do_gps? 1023 : 4092) * SPC, rep_cylce);// fft_size* SPC - rep_cylce);

    for (int i = 0; i < fft_size * SPC; i++) {
      //samp_a[i].r += samp_b[i].r; // add noise & quantize later
      //samp_a[i].i += samp_b[i].i;
    }

    fft_c32(fft_size * SPC, samp_a, true);
    fft_c32(fft_size * SPC, repl_a, true);
    for (int i = 0; i < fft_size * SPC; i++) { prod[i] = mult(samp_a[i], get_conj(repl_a[i])); }
    fft_c32(fft_size * SPC, prod, false); // in-place inv F(prod)
    FILE* dbg_fp = NULL; //fopen_s(&dbg_fp, "C:/Python/out6.csv", "w");
    find_top2_peaks_cplx(prod, fft_size * SPC, 3, &peaks, dbg_fp); if (dbg_fp) { fclose(dbg_fp); }
    double interp = interpCodePhase(prod, fft_size * SPC, &peaks);
    printf("max_float,%f offset,%f pos,%d interp,%f error,%f\n", peaks.val1, double(offset), peaks.idx1, interp, SPEED_LIGHT * (offset - (rep_cylce + interp)) / (1023.0 * scale * SPC));
  }
  
  free(prod); free(samp_a); free(repl_a); free(samp_b);
}

#define USE_FFT 1
void TimeDomainCorrelate(int size, c32* samples, c32* code, float* power)
{
  uint32_t i, j;

  //memset(power, 0, sizeof(*power) * size);

  for (i = 0; i < size; i++)
  {
    double sumI = 0.0;
    double sumQ = 0.0;
    for (j = 0; j < size; j++)
    {
      sumI += samples[(i + j) % size].r * code[j].r;
      sumQ += samples[(i + j) % size].i * code[j].i;
    }
    power[i] += (float)sqrt(sumI * sumI + sumQ * sumQ);

  }

}
/////////////////////////////////////////////////////////////////////////////
void read_E5A(char* input) {
  FILE* fp_1bitcsv = NULL;
  fopen_s(&fp_1bitcsv, input, "r");
  if (fp_1bitcsv == NULL) {
    fprintf(stderr, "Failed to open sample file %s\n", input);
    return;
  }
  fseek(fp_1bitcsv, 0L, SEEK_END);
  size_t bytes_to_read = ftell(fp_1bitcsv);
  rewind(fp_1bitcsv);
  FILE* fp_out = NULL; //output file
 
#define SPC 1
#define FFT_SIZE 16384 
#define SAMP 10230 // for 1 ms at 10.23 MHz
  //#define SAMP 16384
  //#define SAMP 15345
    // only 9 and 36 with q31; 10, 6 also works with float
    //int prn = 36;// 6;// 6;// 36;// 9;// 36
    //double doppler = -1*(1580 + 1e6 +2500);// E36:1580 E6: 1261 ; G10:-582, G32:1232

  char filename[256];
  strcpy_s(filename, input);
  *strrchr(filename, '.') = 0;
  strcat_s(filename, ".ast");
  FILE* fpAssist = NULL;
  fopen_s(&fpAssist, filename, "r");
  if (fpAssist == NULL) {
    fprintf(stderr, "Failed to open assitance file %s\n", filename);
    return;
  }

  acq_struct prn2acq[40] = { 0 };
  int cnt = 0;
  char line[256];
  /**/
  while (fgets(line, sizeof(line), fpAssist) != NULL)
  {
    char c = line[0];
    if (c == 'G')
      prn2acq[cnt].constel = 1;
    else if (c == 'E')
      prn2acq[cnt].constel = 2;
    else
      continue;
    sscanf_s(line + 1, "%d %lf", &prn2acq[cnt].prn, &prn2acq[cnt].doppler);
    if (c == 'G' && !HasGPSL5(prn2acq[cnt].prn)) continue;
    prn2acq[cnt].doppler *= 1176.45 / 1575.42; // adjust to L5

    cnt++;
  }
  fclose(fpAssist);

  /////////////////////////////////////////////////////
  int corrLength = USE_FFT ? FFT_SIZE : SAMP;

  c32* sampl = (c32*)malloc(SAMP * sizeof(c32));
  c32* repli = (c32*)malloc(SAMP * sizeof(c32));
  c32* up_samp = (c32*)malloc(corrLength * sizeof(c32));
  c32* up_repli = (c32*)malloc(corrLength * sizeof(c32));
  c32* up_prod = (c32*)malloc(corrLength * sizeof(c32));
  c32* sum_prod = (c32*)malloc(corrLength * sizeof(c32));
  float* nci_sum = (float*)malloc(corrLength * sizeof(float));

  bb_meas_t meas;
  memset(&meas, 0, sizeof(bb_meas_sat_t));

  char* pc = strrchr(input, '/') + 1;
  strcpy_s(filename, pc);
  *strrchr(filename, '.') = 0;
  strcat_s(filename, "-pwr.txt");
  errno_t er = fopen_s(&fp_out, filename, "w");
  if (er != 0 || fp_out == NULL) { fprintf(stderr, "Failed to open output file\n"); return; }


  ///////////////////// main prn loop ////////////////////////////////
  for (int prn_loop = 0; prn_loop < cnt; prn_loop++) {
    int prn = prn2acq[prn_loop].prn;
    int gal_proc = (prn2acq[prn_loop].constel == 2) ? 1 : 0;
    //printf("Processing PRN %d Doppler %f constel %d \n", prn, doppler, prn2acq[prn_loop].constel);

    memset(sampl, 0, sizeof(c32) * SAMP);
    memset(repli, 0, sizeof(c32) * SAMP);
    memset(up_samp, 0, sizeof(c32) * corrLength);
    memset(up_repli, 0, sizeof(c32) * corrLength);
    memset(up_prod, 0, sizeof(c32) * corrLength);
    memset(sum_prod, 0, sizeof(c32) * corrLength);

    int dopplerOffset;
    double ratio = 0.0;
    double bestRatio = 0.0;
    for (dopplerOffset = -600; dopplerOffset <= 600; dopplerOffset += 100)
    {
      double doppler = -1 * (prn2acq[prn_loop].doppler + 1e6 + 4100 + dopplerOffset);
      //double doppler = -1 * (prn2acq[prn_loop].doppler + 4100 + dopplerOffset);

      int codeoff = 0;
      for (codeoff = 0; codeoff < 10230; codeoff += 3410)
      {
        if (gal_proc) {
          synth_e5a_prn(prn, -doppler, SAMP, repli, 0);
        }
        else { // GPS
          synth_L5I_prn(prn, -doppler, SAMP, repli, 0);
        }
        if (USE_FFT)
        {
          up_sample_10k_to_16k(repli, up_repli);
          //memcpy(repli, up_repli, sizeof(c32)* SAMP);

          fft_c32(FFT_SIZE, up_repli, true); // forward FFT 
        }

        char* context = nullptr;
        // read in the csv data
        //char line[256];
        int LEN = SAMP;
        int NCI = 20;

        //if (0)
        if (codeoff > 0)
        {
          NCI = 19;
          for (int skip = 0; skip < codeoff; skip++)
            fgets(line, sizeof(line), fp_1bitcsv);
        }

        ///////////// NCI loop /////////////////////////////////////////////
        memset(nci_sum, 0, sizeof(float) * corrLength);
        for (int loop = 0; loop < NCI; loop++) {
          int idx = 0;
          while (!feof(fp_1bitcsv)) {
            if (fgets(line, sizeof(line), fp_1bitcsv) != NULL) {
              char* token = strtok_s(line, ",", &context);
              token = strtok_s(NULL, ",", &context);
              if (token != NULL) {
                sampl[idx].r = (float)atof(token);
                token = strtok_s(NULL, ",", &context);
                if (token != NULL) { sampl[idx].i = (float)atof(token); }
              }
              idx++;
              if ((idx != 0) && (idx % LEN == 0)) { break; }
            }
          }
          if (USE_FFT)
          {
            //up_sample_10k_to_16k(sampl, up_samp);
            up_sample_N_to_M(sampl, SAMP, up_samp, FFT_SIZE);
            //memcpy(sampl, up_samp, sizeof(c32) * SAMP);

            // note repli has been FFTed already
            fft_c32(FFT_SIZE, up_samp, true); // forward FFT
            for (int k = 0; k < FFT_SIZE; k++) { up_prod[k] = mult(up_samp[k], get_conj(up_repli[k])); }
            fft_c32(FFT_SIZE, up_prod, false); // IFFT
            for (int i = 0; i < FFT_SIZE; i++) { nci_sum[i] += mag(up_prod[i]); }
          }
          else
          {
            TimeDomainCorrelate(SAMP, sampl, repli, nci_sum);
          }

          //printf("loop %d \n", loop);
        } // end NCI for loop

        top2_pks peaks;
        //fprintf(fp_out, "%c%02d %6.0f %5d %5d\n", (gal_proc == 1) ? 'E' : 'G', prn2acq[prn_loop].prn, prn2acq[prn_loop].doppler + dopplerOffset, dopplerOffset, codeoff);
        find_top2_peaks_real(nci_sum, corrLength, 3, &peaks, fp_out);
        double cn0 = compute_snr_real(nci_sum, corrLength, peaks) + 35.0;
        double early = nci_sum[peaks.idx1 - 1], prompt = nci_sum[peaks.idx1], late = nci_sum[peaks.idx1 + 1];
        double interp = InterpolateCodePhase(peaks.idx1, early * early, prompt * prompt, late * late) * SAMP / corrLength;


        interp += codeoff;
        if (interp > SAMP)
          interp -= SAMP;

        ratio = peaks.ratio;
        printf("Avail PRN %c%02d doppler %6.0f ratio %5.2f loc %5d interp %10.4f CN0 %5.2f %5d %5d %c\n", (gal_proc == 1) ? 'E' : 'G', prn2acq[prn_loop].prn, prn2acq[prn_loop].doppler + dopplerOffset, ratio, (int)peaks.idx1, interp, cn0, dopplerOffset, codeoff, (ratio > bestRatio) ? '*' : ' ');

        if (ratio > bestRatio) {
          meas.sats[meas.num_sat].prn = prn2acq[prn_loop].prn;
          meas.sats[meas.num_sat].code_phase = float(interp / SAMP); // ie (interp / 16384) * (16384 /16368) fixme!
          meas.sats[meas.num_sat].doppler = -(prn2acq[prn_loop].doppler + dopplerOffset);
          meas.sats[meas.num_sat].cno = (float)cn0;
          meas.sats[meas.num_sat].constellation = gal_proc ? SYS_GAL : SYS_GPS;

          bestRatio = ratio;
          //printf("Acquired %s %d Doppler %f Hz CodePhase %f [ms] C/N0 %f dB-Hz ratio=%f\n", (gal_proc == 1) ? "GAL" : "GPS",
          //  prn2acq[prn_loop].prn, -prn2acq[prn_loop].doppler, meas.sats[meas.num_sat].code_phase, meas.sats[meas.num_sat].cno, ratio * ratio);
        }
        rewind(fp_1bitcsv);
      }
      rewind(fp_1bitcsv);
    }
    if (bestRatio > 1.29)
      meas.num_sat++;
  } // end for prn_loop
  fclose(fp_out);

  pc = strrchr(input, '/') + 1;
  strcpy_s(filename, pc);
  *strrchr(filename, '.') = 0;
  strcat_s(filename, ".msb");
  int num_bytes = write_msb(&meas, (char*)filename);
  bb_meas_t check = { 0 };
  FILE* test = NULL;
  errno_t er2 = fopen_s(&test, filename, "rb");
  uint8_t tbuff[128] = { 0 };
  fread(tbuff, 1, num_bytes, test);
  read_bb_msb(tbuff, num_bytes, &check);
  fclose(test);


  fclose(fp_1bitcsv);
  free(sampl); free(repli); free(up_samp); free(up_repli); free(up_prod); free(nci_sum); free(sum_prod);
}

/////////////////////////////////////////////////////////////////////////////
void read_L5E5AE5QP(char* input) {
  FILE* fp_1bitcsv = NULL;
  fopen_s(&fp_1bitcsv, input, "r");
  if (fp_1bitcsv == NULL) {
    fprintf(stderr, "Failed to open sample file %s\n", input);
    return;
  }
  fseek(fp_1bitcsv, 0L, SEEK_END);
  size_t bytes_to_read = ftell(fp_1bitcsv);
  rewind(fp_1bitcsv);
  FILE* fp_out = NULL; //output file

#define SPC 2
#define FFT_SIZE 16384 
//#define SAMP 10230 // for 1 ms at 10.23 MHz
//#define SAMP 5115 // for 1 ms at 10.23 MHz
//#define SAMP 16384
//#define SAMP 15345
  int msps = 10230;
  if (optionMsps > 0)
    msps = optionMsps;
  // only 9 and 36 with q31; 10, 6 also works with float
  //int prn = 36;// 6;// 6;// 36;// 9;// 36
  //double doppler = -1*(1580 + 1e6 +2500);// E36:1580 E6: 1261 ; G10:-582, G32:1232

  char filename[256];
  strcpy_s(filename, input);
  *strrchr(filename, '.') = 0;
  strcat_s(filename, ".ast");
  FILE* fpAssist = NULL;
  fopen_s(&fpAssist, filename, "r");
  if (fpAssist == NULL) {
    fprintf(stderr, "Failed to open assitance file %s\n", filename);
    return;
  }

  if (msps < 10230 && !(optionSig & OPTION_QP))
  {
    fprintf(stderr, "msps too low\n");
    return;
  }

  acq_struct prn2acq[40] = { 0 }; 
  int cnt = 0;
  char line[256];
  /**/
  bool hasQp = false;

  while (fgets(line, sizeof(line), fpAssist) != NULL)
  {
    char c = line[0];
    if (c == 'G')
      prn2acq[cnt].constel = 1;
    else if (c == 'E')
      prn2acq[cnt].constel = 2;
    else
      continue;
    sscanf_s(line + 1, "%d %lf", &prn2acq[cnt].prn, &prn2acq[cnt].doppler);
    if (c == 'G' && !HasGPSL5(prn2acq[cnt].prn)) continue; // won't be able to get these
    hasQp = (optionSig & OPTION_QP) && (c == 'E') && HasE5QP(prn2acq[cnt].prn);
    if (!hasQp && msps < 10230) continue; // won't work
    //if (!hasQp) continue; // only do QP

    prn2acq[cnt].doppler *= 1176.45 / 1575.42; // adjust to L5

    cnt++;
  }
  fclose(fpAssist);

  /////////////////////////////////////////////////////
  int maxCorrLength = USE_FFT ? FFT_SIZE : msps;

  c32* sampl     = (c32*)malloc(msps * sizeof(c32));
  c32* repli     = (c32*)malloc(msps * sizeof(c32));
  c32* up_samp   = (c32*)malloc(maxCorrLength * sizeof(c32));
  c32* up_repli  = (c32*)malloc(maxCorrLength * sizeof(c32));
  c32* up_prod   = (c32*)malloc(maxCorrLength * sizeof(c32));
  c32* sum_prod  = (c32*)malloc(maxCorrLength * sizeof(c32));
  float* nci_sum = (float*)malloc(maxCorrLength * sizeof(float));

  bb_meas_t meas;
  memset(&meas, 0, sizeof(bb_meas_sat_t)); 

  char* dateTimeStr = strchr(input,'G'); // sample filename
  int year, month, day, hour, minute;
  double sec;
  sscanf_s(dateTimeStr, "G_%d_%d_%d_%d_%d_%lf", &year, &month, &day, &hour, &minute, &sec);
  double gpstime = compute_gps_time(year, month, day, hour, minute, sec);
  meas.time = (uint64_t)gpstime;
  meas.sec = gpstime - (uint64_t)gpstime;
  int week = (int)(gpstime / 604800.0);
  double tow = fmod(gpstime, 604800.0);
  
  char* pc = 0;
  if (0)
  {
    pc = strrchr(input, '/') + 1;
    strcpy_s(filename, pc);
    *strrchr(filename, '.') = 0;
    strcat_s(filename, "-pwr.txt");
    errno_t er = fopen_s(&fp_out, filename, "w");
    if (er != 0 || fp_out == NULL) { fprintf(stderr, "Failed to open output file\n"); return; }
  }

  ///////////////////// main prn loop ////////////////////////////////
  for (int prn_loop = 0; prn_loop < cnt; prn_loop++) {
    int prn = prn2acq[prn_loop].prn;
    int gal_proc = (prn2acq[prn_loop].constel == 2) ? 1 : 0;
    hasQp = (optionSig & OPTION_QP) && gal_proc && HasE5QP(prn);
    //EDif (gal_proc == 0 || hasQp == false) { continue; }
    //EDif ( (prn2acq[prn_loop].prn == 15 || prn2acq[prn_loop].prn == 34) == false) { continue; }
    //printf("Processing PRN %d Doppler %f constel %d \n", prn, doppler, prn2acq[prn_loop].constel);
    printf("Processing %c%02d%c around %.1lf Hz\n", gal_proc ? 'E' : 'G', prn, gal_proc ? (hasQp ? 'q' : 'a') : ' ', prn2acq[prn_loop].doppler);
    int codeLength = (hasQp) ? E5_QP_CODE_LEN : E5A_CODE_LEN;
    int codeRate = (hasQp) ? 5115 : 10230;
    int spc = msps / codeRate;

    int mixcount = 0;
    int fftcount = 0;
    int ifftcount = 0;
    int interpcount = 0;

    int sampLength = codeLength * spc;
    int corrLength = sampLength;
    if (USE_FFT)
    {
      // find the next power of 2
      for (corrLength = 128; corrLength < sampLength; corrLength *= 2)
        ;
      //printf("fftSize = %d\n", corrLength);
    }

    int LEN = sampLength;
    int msSkip = 0; // skip this many ms worth of samples
    int msUse = 20; // use this many ms of samples
    int tc = 1; // number of code periods to use for coherent integration (same as number of ms for L5/E5a)

    int dopplerStep = 500 / tc;
    if (hasQp)
    {
      int msCoh = 2;
      if (optionMsCoh > 0)
        msCoh = optionMsCoh;
      tc = msCoh * 31 / 2; // n*31/2 for n ms
      dopplerStep = 500 / (tc * 2.0/ 31);
    }
    //dopplerStep = 200;
    

    memset(sampl   , 0, sizeof(c32) * sampLength);
    memset(repli   , 0, sizeof(c32) * sampLength);
    memset(up_samp, 0, sizeof(c32)* corrLength);
    memset(up_repli, 0, sizeof(c32) * corrLength);
    memset(up_prod , 0, sizeof(c32) * corrLength);
    memset(sum_prod, 0, sizeof(c32) * corrLength);
    
    int dopplerOffset;
    int dopplerUncertainty = 3 * dopplerStep;
    if (optionDUnc >= 0) dopplerUncertainty = optionDUnc;
    double ratio = 0.0;
    double bestRatio = 0.0;
    int bestCodeOffset = 0.0;
    int bestDoppOffset = 0.0;
    int codeoff = 0;
    double div = 2;
    if (hasQp) div = 1;

    // Generate code replica

    if (gal_proc) {
      if (hasQp)
      {
        //getE5_QPCode(E5_QP_CODE_LEN, 1, prn2acq[cnt].prn, prn_code);
        //make_replica(prn_code, replica, doppler, E5_QP_CODE_LEN * spc, chipping_rate * SPC);
        //memcpy(fft_repl, replica, sizeof(c32) * E5_QP_CODE_LEN);
        //synth_e5qp_prn(prn, -doppler, msps, repli, 0);
        //synth_e5qp_prn(prn, 0.0, msps, repli, 0);
        int* prn_code = (int*)malloc(E5_QP_CODE_LEN * spc * sizeof(int));
        getE5_QPCode(E5_QP_CODE_LEN, spc, prn, prn_code);
        memset(repli, 0, sizeof(c32) * sampLength);
        for (int sampi = 0; sampi < sampLength; sampi++)
          repli[sampi].r = prn_code[sampi];
        free(prn_code);

      }
      else
      {
        //synth_e5a_prn(prn, -doppler, sampLength, repli, 0);
        synth_e5a_prn(prn, 0.0, msps, repli, 0);
      }
    }
    else { // GPS
      synth_L5I_prn(prn, 0.0, msps, repli, 0);
    }
    if (USE_FFT)
    {
      //up_sample_10k_to_16k(repli, up_repli);
      //memcpy(repli, up_repli, sizeof(c32)* SAMP);
      memset(up_repli, 0, sizeof(c32) * corrLength);
      //memcpy(up_repli, repli, sampLength);
      //up_sample_N_to_M(repli, sampLength, up_repli, corrLength);
      upsample_linear_c32(repli, sampLength, up_repli, corrLength, false);
      interpcount++;
      fft_c32(corrLength, up_repli, true); // forward FFT 
      fftcount++;
    }


    //for (codeoff = bestCodeOffset - (sampLength/div); codeoff <= (sampLength/div); codeoff += (sampLength/(2*div))) // try different code offsets
    for (codeoff = 0; codeoff < sampLength; codeoff += sampLength / div) // try different code offsets
    {
      for (dopplerOffset = -dopplerUncertainty; dopplerOffset <= dopplerUncertainty; dopplerOffset += dopplerStep)
      {
        //double doppler = -1 * (prn2acq[prn_loop].doppler + 1e6 + 4100 + dopplerOffset);
        double doppler = -1 * (prn2acq[prn_loop].doppler + 4100 + dopplerOffset);


        char* context = nullptr;
        // read in the csv data
        //char line[256];
        int NCI = msUse / tc;
        if (hasQp)
        {
          NCI = msUse / (tc * 2.0 / 31);
        }

        if (msSkip > 0)
        {
          // skip samples at the beginning of the file
          for (int skip = 0; skip < msSkip*SAMP; skip++)
            fgets(line, sizeof(line), fp_1bitcsv);
        }
        if (codeoff > 0)
        {
          NCI--; // we'll run out of space if there's a code offset because we are skipping samples to do this
          for (int skip = 0; skip < codeoff; skip++)
            fgets(line, sizeof(line), fp_1bitcsv);
        }

        ///////////// NCI loop /////////////////////////////////////////////
        memset(nci_sum, 0, sizeof(float) * corrLength);
        //double upsampSecs = 1.0 / (msps * 1000.0 * corrLength / sampLength);
        //double dphi = 2.0 * PI_F64 * doppler * upsampSecs;
        double sampSecs = 1.0 / (msps * 1000.0);
        double dphi = 2.0 * PI_F64 * doppler * sampSecs;
        if (dphi < 0) dphi += 2.0 * PI_F64;
        uint32_t totalSampleCount = 0;
        for (int nci = 0; nci < NCI; nci++)
        {
          memset(sampl, 0, sizeof(c32) * sampLength);
          memset(up_samp, 0, sizeof(c32) * corrLength);
          memset(up_prod, 0, sizeof(c32) * corrLength);
          memset(sum_prod, 0, sizeof(c32) * corrLength);
          for (int tci = 0; tci < tc; tci++)
          {
            int idx = 0;
            /*
            while (!feof(fp_1bitcsv)) {
              if (fgets(line, sizeof(line), fp_1bitcsv) != NULL) {
                char* token = strtok_s(line, ",", &context);
                token = strtok_s(NULL, ",", &context);
                if (token != NULL) {
                  sampl[idx].r = (float)atof(token);
                  token = strtok_s(NULL, ",", &context);
                  if (token != NULL) { sampl[idx].i = (float)atof(token); }
                }
                idx++;
                if ((idx != 0) && (idx % LEN == 0)) { break; }
              }
            }
            // upsample
            c32* temp_up_samp = (c32*)malloc(maxCorrLength * sizeof(c32));
            memset(temp_up_samp, 0, sizeof(c32) * corrLength);
            //memcpy(temp_up_samp, sampl, sampLength);
            //up_sample_N_to_M(sampl, sampLength, temp_up_samp, corrLength);
            upsample_linear_c32(sampl, sampLength, temp_up_samp, corrLength, false);

            // wipe Doppler
            for (int sampi = 0; sampi < corrLength; sampi++)
            {
              c32 phase;
              double phi = fmod(dphi * totalSampleCount,2*PI_F64);
              phase.r = cos(phi);
              phase.i = sin(phi);
              c32 tempsamp1 = mult(temp_up_samp[sampi], phase);
              c32 tempsamp2 = add(up_samp[sampi], tempsamp1);
              up_samp[sampi] = tempsamp2;
              totalSampleCount++;
            }
            free(temp_up_samp);
            */
            while (!feof(fp_1bitcsv)) {
              if (fgets(line, sizeof(line), fp_1bitcsv) != NULL) {
                char* token = strtok_s(line, ",", &context);
                token = strtok_s(NULL, ",", &context);
                if (token != NULL) {
                  sampl[idx].r = (float)atof(token);
                  token = strtok_s(NULL, ",", &context);
                  if (token != NULL) { sampl[idx].i = (float)atof(token); }
                }
                c32 phase;
                double phi = fmod(dphi * totalSampleCount, 2 * PI_F64);
                phase.r = cos(phi);
                phase.i = sin(phi);
                c32 tempsamp1 = mult(sampl[idx], phase);
                c32 tempsamp2 = add(up_samp[idx], tempsamp1);
                up_samp[idx] = tempsamp2;
                totalSampleCount++;

                idx++;

                if ((idx != 0) && (idx % LEN == 0)) { break; }
              }
            }
            mixcount++;

          }
          if (USE_FFT)
          {
            /*
            // note repli has been FFTed already
            fft_c32(corrLength, up_samp, true); // forward FFT
            for (int k = 0; k < corrLength; k++) { up_prod[k] = mult(up_samp[k], get_conj(up_repli[k])); }
            fft_c32(corrLength, up_prod, false); // IFFT
            for (int i = 0; i < corrLength; i++) { nci_sum[i] += mag(up_prod[i]); }
            */

            // upsample
            c32* temp_up_samp = (c32*)malloc(maxCorrLength * sizeof(c32));
            memset(temp_up_samp, 0, sizeof(c32)* corrLength);
            //memcpy(temp_up_samp, up_samp, sampLength);
            //up_sample_N_to_M(up_samp, sampLength, temp_up_samp, corrLength);
            upsample_linear_c32(up_samp, sampLength, temp_up_samp, corrLength, false);
            interpcount++;
            fftcount++;
            ifftcount++;

            // note repli has been FFTed already
            fft_c32(corrLength, temp_up_samp, true); // forward FFT
            for (int k = 0; k < corrLength; k++) { up_prod[k] = mult(temp_up_samp[k], get_conj(up_repli[k])); }
            fft_c32(corrLength, up_prod, false); // IFFT
            for (int i = 0; i < corrLength; i++) { nci_sum[i] += mag(up_prod[i]); }
            free(temp_up_samp);

          }
          else
          {
            //TimeDomainCorrelate(sampLength, sampl, repli, nci_sum);
          }
          if (0)
          {
            top2_pks peaks;
            find_top2_peaks_real(nci_sum, corrLength, 10, &peaks, fp_out);
            double cn0 = compute_snr_real(nci_sum, corrLength, peaks) + 35.0;
            double early = nci_sum[peaks.idx1 - 1], prompt = nci_sum[peaks.idx1], late = nci_sum[peaks.idx1 + 1];
            double interp = InterpolateCodePhase(peaks.idx1, early * early, prompt * prompt, late * late) * sampLength / corrLength;
            interp += codeoff;
            if (interp > sampLength)
              interp -= sampLength;
            ratio = peaks.ratio;
            printf("NC%03d PRN %c%02d doppler %6.0f ratio %5.2f loc %5d interp %10.4f CN0 %5.2f %5d %5d %c\n", nci, (gal_proc == 1) ? 'E' : 'G', prn2acq[prn_loop].prn, prn2acq[prn_loop].doppler + dopplerOffset, ratio, (int)peaks.idx1, interp, cn0, dopplerOffset, codeoff, (ratio > bestRatio) ? '*' : ' ');
          }
        }

        top2_pks peaks;
        //fprintf(fp_out, "%c%02d %6.0f %5d %5d\n", (gal_proc == 1) ? 'E' : 'G', prn2acq[prn_loop].prn, prn2acq[prn_loop].doppler + dopplerOffset, dopplerOffset, codeoff);
        find_top2_peaks_real(nci_sum, corrLength, 10, &peaks, fp_out);
        double cn0 = compute_snr_real(nci_sum, corrLength, peaks) + 35.0;
        double early = nci_sum[peaks.idx1 - 1], prompt = nci_sum[peaks.idx1], late = nci_sum[peaks.idx1 + 1];
        double interp = InterpolateCodePhase(peaks.idx1, early * early, prompt * prompt, late * late) * sampLength / corrLength;


        interp += codeoff;
        if (interp > sampLength)
          interp -= sampLength;

        ratio = peaks.ratio;
        printf("Avail PRN %c%02d doppler %6.0f ratio %5.2f loc %5d interp %10.4f CN0 %5.2f %5d %5d %c\n", 
          (gal_proc == 1) ? 'E' : 'G', prn2acq[prn_loop].prn, prn2acq[prn_loop].doppler + dopplerOffset, ratio, (int)peaks.idx1, interp, cn0, dopplerOffset, codeoff, (ratio > bestRatio) ? '*' : ' ');

        if (ratio > bestRatio) {
          meas.sats[meas.num_sat].prn = prn2acq[prn_loop].prn;
          meas.sats[meas.num_sat].code_phase = float(interp / msps); // [ms] QP will have some offset of 2/31 segments
          meas.sats[meas.num_sat].doppler = -(prn2acq[prn_loop].doppler + dopplerOffset);
          meas.sats[meas.num_sat].cno = (float)cn0;
          meas.sats[meas.num_sat].ratio = (float)ratio;
          meas.sats[meas.num_sat].constellation = gal_proc ? SYS_GAL : SYS_GPS;
          meas.sats[meas.num_sat].signal = gal_proc ? (hasQp ? SIGNAL_E5aQP : SIGNAL_E5a) : SIGNAL_L5;
          if ((optionSig & OPTION_QX) && gal_proc && HasE5QP(prn))
            meas.sats[meas.num_sat].signal = SIGNAL_E5aQP;
          bestDoppOffset = dopplerOffset;
          bestCodeOffset = codeoff;

          if (fp_out != NULL)
          {
            fprintf(fp_out, "PRN %c%02d, doppler %6.0f, ratio %5.2f, loc %5d, interp %10.4f, CN0 %5.2f, %5d, %5d\n", (gal_proc == 1) ? 'E' : 'G', prn2acq[prn_loop].prn, prn2acq[prn_loop].doppler + dopplerOffset, ratio, (int)peaks.idx1, interp, cn0, dopplerOffset, codeoff);
            for (int pwri = 0; pwri < corrLength; pwri++) fprintf(fp_out, "%03d, %10.4lf\n", pwri, nci_sum[pwri]);
            fflush(fp_out);
          }
            
          bestRatio = ratio;
          //printf("Acquired %s %d Doppler %f Hz CodePhase %f [ms] C/N0 %f dB-Hz ratio=%f\n", (gal_proc == 1) ? "GAL" : "GPS",
          //  prn2acq[prn_loop].prn, -prn2acq[prn_loop].doppler, meas.sats[meas.num_sat].code_phase, meas.sats[meas.num_sat].cno, ratio * ratio);
        } 
        rewind(fp_1bitcsv);
      } // for codeoff
      //div *= 2;
      rewind(fp_1bitcsv);
    } // for doppler
#ifdef COUNT_OPERATIONS
    fprintf(fpResults, "%c%02d, %d, %3d, %2.0lf, %3d, %3d, %3d, %3d\n",
      meas.sats[meas.num_sat].constellation == SYS_GPS ? 'G' : 'E',
      meas.sats[meas.num_sat].prn,
      meas.sats[meas.num_sat].signal,
      dopplerUncertainty,
      div,
      mixcount,
      interpcount,
      fftcount,
      ifftcount);
#endif

    printf("Avail PRN %c%02d doppler %6.0f ratio %5.2f loc %5d CN0 %5.2f  %5d %5d \n", (gal_proc == 1) ? 'E' : 'G', 
      meas.sats[meas.num_sat].prn, meas.sats[meas.num_sat].doppler, bestRatio,(int) (meas.sats[meas.num_sat].code_phase *msps), meas.sats[meas.num_sat].cno, bestDoppOffset, bestCodeOffset);
    if (bestRatio > 1.29) {
      meas.num_sat++;
    }
  } // end for prn_loop
  if (fp_out != NULL) fclose(fp_out);

#ifndef COUNT_OPERATIONS
  if (fpResults != NULL)
  {
    for (int imeas = 0; imeas < meas.num_sat; imeas++)
      fprintf(fpResults, "%c%02d, %d, %04d, %9.3lf, %6.0lf, %5.2f, %.6lf, %10.3lf\n",
        meas.sats[imeas].constellation == SYS_GPS ? 'G' : 'E',
        meas.sats[imeas].prn,
        meas.sats[imeas].signal,
        week,
        tow,
        meas.sats[imeas].doppler,
        meas.sats[imeas].ratio,
        meas.sats[imeas].code_phase,
        meas.sats[imeas].code_phase * SPEED_LIGHT * 0.001
      );
  }
#endif

  pc = strrchr(input, '/')+1;
  strcpy_s(filename, pc);
  *strrchr(filename, '.') = 0;
  strcat_s(filename, ".msb");
  printf("Writing %s\n", filename);
  int num_bytes = write_msb(&meas, (char*)filename);
  bb_meas_t check = { 0 };
  FILE* test = NULL;
  errno_t er2 = fopen_s(&test, filename, "rb");
  uint8_t tbuff[128] = { 0 };
  fread(tbuff, 1, num_bytes, test);
  read_bb_msb(tbuff, num_bytes, &check);
  fclose(test);

 
  fclose(fp_1bitcsv); 
  free(sampl); free(repli); free(up_samp); free(up_repli); free(up_prod); free(nci_sum); free(sum_prod);
}
/////////////////////////////////////////////////////////////////////////////
void read_L1E1(char* input) {
  FILE* fp_1bitcsv = NULL;
  fopen_s(&fp_1bitcsv, input, "r");
  if (fp_1bitcsv == NULL) {
    fprintf(stderr, "Failed to open sample file %s\n", input);
    return;
  }
  fseek(fp_1bitcsv, 0L, SEEK_END);
  size_t bytes_to_read = ftell(fp_1bitcsv);
  rewind(fp_1bitcsv);
  FILE* fp_out = NULL; //output file

#define SPC 2
#define FFT_SIZE 16384 
  //#define SAMP 10230 // for 1 ms at 10.23 MHz
//#define SAMP (4092*4) // for 1 ms at 10.23 MHz
//#define SAMP 16384
//#define SAMP 15345
  int msps = 4092;
  if (optionMsps > 0)
    msps = optionMsps;
  // only 9 and 36 with q31; 10, 6 also works with float
  //int prn = 36;// 6;// 6;// 36;// 9;// 36
  //double doppler = -1*(1580 + 1e6 +2500);// E36:1580 E6: 1261 ; G10:-582, G32:1232

  char filename[256];
  strcpy_s(filename, input);
  *strrchr(filename, '.') = 0;
  strcat_s(filename, ".ast");
  FILE* fpAssist = NULL;
  fopen_s(&fpAssist, filename, "r");
  if (fpAssist == NULL) {
    fprintf(stderr, "Failed to open assitance file %s\n", filename);
    return;
  }

  acq_struct prn2acq[40] = { 0 };
  int cnt = 0;
  char line[256];
  /**/
  bool hasQp = false;

  while (fgets(line, sizeof(line), fpAssist) != NULL)
  {
    char c = line[0];
    if (c == 'G')
      prn2acq[cnt].constel = 1;
    else if (c == 'E')
      prn2acq[cnt].constel = 2;
    else
      continue;
    sscanf_s(line + 1, "%d %lf", &prn2acq[cnt].prn, &prn2acq[cnt].doppler);
    if (c == 'G' && !HasGPSL5(prn2acq[cnt].prn)) continue; // keep the same sats as qp testing
    hasQp = (c == 'E') && HasE5QP(prn2acq[cnt].prn);
    //if (!hasQp && SAMP < 10230) continue; // won't work
    //if (!hasQp) continue; // only do QP

    //prn2acq[cnt].doppler *= 1176.45 / 1575.42; // adjust to L5

    cnt++;
  }
  fclose(fpAssist);

  /////////////////////////////////////////////////////
  int maxSampleLength = msps * 4; // allow for E1
  int maxCorrLength = USE_FFT ? FFT_SIZE : maxSampleLength;


  c32* sampl = (c32*)malloc(maxSampleLength * sizeof(c32));
  c32* repli = (c32*)malloc(maxSampleLength * sizeof(c32));
  c32* up_samp = (c32*)malloc(maxCorrLength * sizeof(c32));
  c32* up_repli = (c32*)malloc(maxCorrLength * sizeof(c32));
  c32* up_prod = (c32*)malloc(maxCorrLength * sizeof(c32));
  c32* sum_prod = (c32*)malloc(maxCorrLength * sizeof(c32));
  float* nci_sum = (float*)malloc(maxCorrLength * sizeof(float));

  bb_meas_t meas;
  memset(&meas, 0, sizeof(bb_meas_sat_t));

  char* dateTimeStr = strchr(input, 'G'); // sample filename
  int year, month, day, hour, minute;
  double sec;
  sscanf_s(dateTimeStr, "G_%d_%d_%d_%d_%d_%lf", &year, &month, &day, &hour, &minute, &sec);
  double gpstime = compute_gps_time(year, month, day, hour, minute, sec);
  meas.time = (uint64_t)gpstime;
  meas.sec = gpstime - (uint64_t)gpstime;
  int week = (int)(gpstime / 604800.0);
  double tow = fmod(gpstime, 604800.0);

  char* pc = 0;
  if (0)
  {
    pc = strrchr(input, '/') + 1;
    strcpy_s(filename, pc);
    *strrchr(filename, '.') = 0;
    strcat_s(filename, "-pwr.txt");
    errno_t er = fopen_s(&fp_out, filename, "w");
    if (er != 0 || fp_out == NULL) { fprintf(stderr, "Failed to open output file\n"); return; }
  }

  ///////////////////// main prn loop ////////////////////////////////
  for (int prn_loop = 0; prn_loop < cnt; prn_loop++) {
    int prn = prn2acq[prn_loop].prn;
    int gal_proc = (prn2acq[prn_loop].constel == 2) ? 1 : 0;
    hasQp = gal_proc && HasE5QP(prn);
    //EDif (gal_proc == 0 || hasQp == false) { continue; }
    //EDif ( (prn2acq[prn_loop].prn == 15 || prn2acq[prn_loop].prn == 34) == false) { continue; }
    //printf("Processing PRN %d Doppler %f constel %d \n", prn, doppler, prn2acq[prn_loop].constel);
    printf("Processing %c%02d%c around %.1lf Hz\n", gal_proc ? 'E' : 'G', prn, gal_proc ? (hasQp ? 'q' : 'a') : ' ', prn2acq[prn_loop].doppler);
    int msps = 4092;
    int codeRate = 1023;
    int spc = msps / codeRate;
    int codeLengthMs = (gal_proc) ? 4 : 1;
    int codeLength = codeRate * codeLengthMs;

    int mixcount = 0;
    int fftcount = 0;
    int ifftcount = 0;
    int interpcount = 0;

    int sampLength = codeLength * spc;
    int corrLength = sampLength;
    if (USE_FFT)
    {
      // find the next power of 2
      for (corrLength = 128; corrLength < sampLength; corrLength *= 2)
        ;
      //printf("fftSize = %d\n", corrLength);
    }


    int LEN = sampLength;
    int msSkip = 0; // skip this many ms worth of samples
    int msUse = 20; // use this many ms of samples
    int msCoherent = 4; // number of ms to use for coherent integration
    if (optionMsCoh > 0)
      msCoherent = optionMsCoh;

    int dopplerStep = 500 / msCoherent;
    //dopplerStep = 200;
    int tc = msCoherent / codeLengthMs; // number of code lengths for coherent integration

    memset(sampl, 0, sizeof(c32) * sampLength);
    memset(repli, 0, sizeof(c32) * sampLength);
    memset(up_samp, 0, sizeof(c32) * corrLength);
    memset(up_repli, 0, sizeof(c32) * corrLength);
    memset(up_prod, 0, sizeof(c32) * corrLength);
    memset(sum_prod, 0, sizeof(c32) * corrLength);

    int dopplerOffset;
    int dopplerUncertainty = 3 * dopplerStep;
    if (optionDUnc >= 0) dopplerUncertainty = optionDUnc;
    double ratio = 0.0;
    double bestRatio = 0.0;
    int bestCodeOffset = 0.0;
    int bestDoppOffset = 0.0;
    int codeoff = 0;
    double div = 1;
    if (gal_proc) div = 4;

    // Generate code replica

    if (gal_proc) {
      //synth_e1b_prn(prn, 0, sampLength, repli, spc, 0);
      int8_t* prn_code = (int8_t*)malloc(codeLength * sizeof(int8_t));
      GetE1BCode(prn, 1, prn_code);
      for (int sampi = 0; sampi < sampLength; sampi++)
        repli[sampi].r = prn_code[sampi / spc];
      free(prn_code);
    }
    else {
      //synth_gps_prn(prn, 0, sampLength, repli, spc, 0);
      int* prn_code = (int*)malloc(codeLength * sizeof(int));
      getCode(codeLength, 1, prn, prn_code);
      for (int sampi = 0; sampi < sampLength; sampi++)
        repli[sampi].r = prn_code[sampi / spc];
      free(prn_code);

    }

    if (USE_FFT)
    {
      //up_sample_10k_to_16k(repli, up_repli);
      //memcpy(repli, up_repli, sizeof(c32)* sampLength);
      memset(up_repli, 0, sizeof(c32) * corrLength);
      //memcpy(up_repli, repli, sampLength);
      //up_sample_N_to_M(repli, sampLength, up_repli, corrLength);
      upsample_linear_c32(repli, sampLength, up_repli, corrLength, false);
      interpcount++;

      fft_c32(corrLength, up_repli, true); // forward FFT 
      fftcount++;
    }


    //for (codeoff = bestCodeOffset - (sampLength/div); codeoff <= (sampLength/div); codeoff += (sampLength/(2*div))) // try different code offsets
    for (codeoff = 0; codeoff < sampLength; codeoff += sampLength / div) // try different code offsets
    {
      for (dopplerOffset = -dopplerUncertainty; dopplerOffset <= dopplerUncertainty; dopplerOffset += dopplerStep)
      {
        //double doppler = -1 * (prn2acq[prn_loop].doppler + 1e6 + 4100 + dopplerOffset);
        double doppler = -1 * (prn2acq[prn_loop].doppler + 4100 * 1575.42 / 1176.45 + dopplerOffset);


        char* context = nullptr;
        // read in the csv data
        //char line[256];
        int NCI = msUse / msCoherent;

        if (msSkip > 0)
        {
          // skip samples at the beginning of the file
          for (int skip = 0; skip < msSkip * msps; skip++)
            fgets(line, sizeof(line), fp_1bitcsv);
        }
        if (codeoff > 0)
        {
          NCI--; // we'll run out of space if there's a code offset because we are skipping samples to do this
          for (int skip = 0; skip < codeoff; skip++)
            fgets(line, sizeof(line), fp_1bitcsv);
        }

        ///////////// NCI loop /////////////////////////////////////////////
        memset(nci_sum, 0, sizeof(float) * corrLength);
        //double upsampSecs = 1.0 / (msps * 1000.0 * corrLength / sampLength);
        //double dphi = 2.0 * PI_F64 * doppler * upsampSecs;
        double sampSecs = 1.0 / (msps * 1000.0);
        double dphi = 2.0 * PI_F64 * doppler * sampSecs;
        if (dphi < 0) dphi += 2.0 * PI_F64;
        uint32_t totalSampleCount = 0;
        for (int nci = 0; nci < NCI; nci++)
        {
          memset(sampl, 0, sizeof(c32) * sampLength);
          memset(up_samp, 0, sizeof(c32) * corrLength);
          memset(up_prod, 0, sizeof(c32) * corrLength);
          memset(sum_prod, 0, sizeof(c32) * corrLength);
          for (int tci = 0; tci < tc; tci++)
          {
            int idx = 0;
            /*
            while (!feof(fp_1bitcsv)) {
              if (fgets(line, sizeof(line), fp_1bitcsv) != NULL) {
                char* token = strtok_s(line, ",", &context);
                token = strtok_s(NULL, ",", &context);
                if (token != NULL) {
                  sampl[idx].r = (float)atof(token);
                  token = strtok_s(NULL, ",", &context);
                  if (token != NULL) { sampl[idx].i = (float)atof(token); }
                }
                idx++;
                if ((idx != 0) && (idx % LEN == 0)) { break; }
              }
            }
            // upsample
            c32* temp_up_samp = (c32*)malloc(maxCorrLength * sizeof(c32));
            memset(temp_up_samp, 0, sizeof(c32) * corrLength);
            //memcpy(temp_up_samp, sampl, sampLength);
            //up_sample_N_to_M(sampl, sampLength, temp_up_samp, corrLength);
            upsample_linear_c32(sampl, sampLength, temp_up_samp, corrLength, false);

            // wipe Doppler
            for (int sampi = 0; sampi < corrLength; sampi++)
            {
              c32 phase;
              double phi = fmod(dphi * totalSampleCount,2*PI_F64);
              phase.r = cos(phi);
              phase.i = sin(phi);
              c32 tempsamp1 = mult(temp_up_samp[sampi], phase);
              c32 tempsamp2 = add(up_samp[sampi], tempsamp1);
              up_samp[sampi] = tempsamp2;
              totalSampleCount++;
            }
            free(temp_up_samp);
            */
            while (!feof(fp_1bitcsv)) {
              if (fgets(line, sizeof(line), fp_1bitcsv) != NULL) {
                char* token = strtok_s(line, ",", &context);
                token = strtok_s(NULL, ",", &context);
                if (token != NULL) {
                  sampl[idx].r = (float)atof(token);
                  token = strtok_s(NULL, ",", &context);
                  if (token != NULL) { sampl[idx].i = (float)atof(token); }
                }
                c32 phase;
                double phi = fmod(dphi * totalSampleCount, 2 * PI_F64);
                phase.r = cos(phi);
                phase.i = sin(phi);
                c32 tempsamp1 = mult(sampl[idx], phase);
                c32 tempsamp2 = add(up_samp[idx], tempsamp1);
                up_samp[idx] = tempsamp2;
                totalSampleCount++;

                idx++;

                if ((idx != 0) && (idx % LEN == 0)) { break; }
              }
            }
            mixcount++;

          }
          if (USE_FFT)
          {
            /*
            // note repli has been FFTed already
            fft_c32(corrLength, up_samp, true); // forward FFT
            for (int k = 0; k < corrLength; k++) { up_prod[k] = mult(up_samp[k], get_conj(up_repli[k])); }
            fft_c32(corrLength, up_prod, false); // IFFT
            for (int i = 0; i < corrLength; i++) { nci_sum[i] += mag(up_prod[i]); }
            */

            // upsample
            c32* temp_up_samp = (c32*)malloc(maxCorrLength * sizeof(c32));
            memset(temp_up_samp, 0, sizeof(c32) * corrLength);
            //memcpy(temp_up_samp, up_samp, sampLength);
            //up_sample_N_to_M(up_samp, sampLength, temp_up_samp, corrLength);
            upsample_linear_c32(up_samp, sampLength, temp_up_samp, corrLength, false);
            interpcount++;
            fftcount++;
            ifftcount++;

            // note repli has been FFTed already
            fft_c32(corrLength, temp_up_samp, true); // forward FFT
            for (int k = 0; k < corrLength; k++) { up_prod[k] = mult(temp_up_samp[k], get_conj(up_repli[k])); }
            fft_c32(corrLength, up_prod, false); // IFFT
            for (int i = 0; i < corrLength; i++) { nci_sum[i] += mag(up_prod[i]); }
            free(temp_up_samp);

          }
          else
          {
            //TimeDomainCorrelate(sampLength, sampl, repli, nci_sum);
          }
          if (0)
          {
            top2_pks peaks;
            find_top2_peaks_real(nci_sum, corrLength, 3, &peaks, fp_out);
            double cn0 = compute_snr_real(nci_sum, corrLength, peaks) + 35.0;
            double early = nci_sum[peaks.idx1 - 1], prompt = nci_sum[peaks.idx1], late = nci_sum[peaks.idx1 + 1];
            double interp = InterpolateCodePhase(peaks.idx1, early * early, prompt * prompt, late * late) * sampLength / corrLength;
            interp += codeoff;
            if (interp > sampLength)
              interp -= sampLength;
            ratio = peaks.ratio;
            printf("NC%03d PRN %c%02d doppler %6.0f ratio %5.2f loc %5d interp %10.4f CN0 %5.2f %5d %5d %c\n", nci, (gal_proc == 1) ? 'E' : 'G', prn2acq[prn_loop].prn, prn2acq[prn_loop].doppler + dopplerOffset, ratio, (int)peaks.idx1, interp, cn0, dopplerOffset, codeoff, (ratio > bestRatio) ? '*' : ' ');
          }
        }

        top2_pks peaks;
        //fprintf(fp_out, "%c%02d %6.0f %5d %5d\n", (gal_proc == 1) ? 'E' : 'G', prn2acq[prn_loop].prn, prn2acq[prn_loop].doppler + dopplerOffset, dopplerOffset, codeoff);
        find_top2_peaks_real(nci_sum, corrLength, 10, &peaks, fp_out);
        double cn0 = compute_snr_real(nci_sum, corrLength, peaks) + 35.0;
        double early = nci_sum[peaks.idx1 - 1], prompt = nci_sum[peaks.idx1], late = nci_sum[peaks.idx1 + 1];
        double interp = InterpolateCodePhase(peaks.idx1, early * early, prompt * prompt, late * late) * sampLength / corrLength;


        interp += codeoff;
        if (interp > sampLength)
          interp -= sampLength;

        ratio = peaks.ratio;
        printf("Avail PRN %c%02d doppler %6.0f ratio %5.2f loc %5d interp %10.4f CN0 %5.2f %5d %5d %c\n",
          (gal_proc == 1) ? 'E' : 'G', prn2acq[prn_loop].prn, prn2acq[prn_loop].doppler + dopplerOffset, ratio, (int)peaks.idx1, interp, cn0, dopplerOffset, codeoff, (ratio > bestRatio) ? '*' : ' ');

        if (ratio > bestRatio) {
          meas.sats[meas.num_sat].prn = prn2acq[prn_loop].prn;
          meas.sats[meas.num_sat].code_phase = float(interp / msps); // [ms] QP will have some offset of 2/31 segments
          meas.sats[meas.num_sat].doppler = -(prn2acq[prn_loop].doppler + dopplerOffset);
          meas.sats[meas.num_sat].cno = (float)cn0;
          meas.sats[meas.num_sat].ratio = (float)ratio;
          meas.sats[meas.num_sat].constellation = gal_proc ? SYS_GAL : SYS_GPS;
          meas.sats[meas.num_sat].signal = gal_proc ? SIGNAL_E1 : SIGNAL_L1;
          if ((optionSig & OPTION_QX) && gal_proc && HasE5QP(prn))
            meas.sats[meas.num_sat].signal = SIGNAL_E5aQP;
          bestDoppOffset = dopplerOffset;
          bestCodeOffset = codeoff;

          if (fp_out != NULL)
          {
            fprintf(fp_out, "PRN %c%02d, doppler %6.0f, ratio %5.2f, loc %5d, interp %10.4f, CN0 %5.2f, %5d, %5d\n", (gal_proc == 1) ? 'E' : 'G', prn2acq[prn_loop].prn, prn2acq[prn_loop].doppler + dopplerOffset, ratio, (int)peaks.idx1, interp, cn0, dopplerOffset, codeoff);
            for (int pwri = 0; pwri < corrLength; pwri++) fprintf(fp_out, "%03d, %10.4lf\n", pwri, nci_sum[pwri]);
            fflush(fp_out);
          }

          bestRatio = ratio;
          //printf("Acquired %s %d Doppler %f Hz CodePhase %f [ms] C/N0 %f dB-Hz ratio=%f\n", (gal_proc == 1) ? "GAL" : "GPS",
          //  prn2acq[prn_loop].prn, -prn2acq[prn_loop].doppler, meas.sats[meas.num_sat].code_phase, meas.sats[meas.num_sat].cno, ratio * ratio);
        }
        rewind(fp_1bitcsv);
      } // for codeoff
      //div *= 2;
      rewind(fp_1bitcsv);
    } // for doppler
#ifdef COUNT_OPERATIONS
    fprintf(fpResults, "%c%02d, %d, %3d, %2.0lf, %3d, %3d, %3d, %3d\n",
      meas.sats[meas.num_sat].constellation == SYS_GPS ? 'G' : 'E',
      meas.sats[meas.num_sat].prn,
      meas.sats[meas.num_sat].signal,
      dopplerUncertainty,
      div,
      mixcount,
      interpcount,
      fftcount,
      ifftcount);
#endif

    printf("Avail PRN %c%02d doppler %6.0f ratio %5.2f loc %5d CN0 %5.2f  %5d %5d \n", (gal_proc == 1) ? 'E' : 'G',
      meas.sats[meas.num_sat].prn, meas.sats[meas.num_sat].doppler, bestRatio, (int)(meas.sats[meas.num_sat].code_phase * msps), meas.sats[meas.num_sat].cno, bestDoppOffset, bestCodeOffset);
    if (bestRatio > 1.29) {
      meas.num_sat++;
    }
  } // end for prn_loop
  if (fp_out != NULL) fclose(fp_out);

#ifndef COUNT_OPERATIONS
  if (fpResults != NULL)
  {
    for (int imeas = 0; imeas < meas.num_sat; imeas++)
      fprintf(fpResults, "%c%02d, %d, %04d, %9.3lf, %6.0lf, %5.2f, %.6lf, %10.3lf\n",
        meas.sats[imeas].constellation == SYS_GPS ? 'G' : 'E',
        meas.sats[imeas].prn,
        meas.sats[imeas].signal,
        week,
        tow,
        meas.sats[imeas].doppler,
        meas.sats[imeas].ratio,
        meas.sats[imeas].code_phase,
        meas.sats[imeas].code_phase * SPEED_LIGHT * 0.001
      );
  }
#endif

  pc = strrchr(input, '/') + 1;
  strcpy_s(filename, pc);
  *strrchr(filename, '.') = 0;
  strcat_s(filename, ".msb");
  printf("Writing %s\n", filename);
  int num_bytes = write_msb(&meas, (char*)filename);
  bb_meas_t check = { 0 };
  FILE* test = NULL;
  errno_t er2 = fopen_s(&test, filename, "rb");
  uint8_t tbuff[128] = { 0 };
  fread(tbuff, 1, num_bytes, test);
  read_bb_msb(tbuff, num_bytes, &check);
  fclose(test);


  fclose(fp_1bitcsv);
  free(up_samp);
  free(sampl);
  free(repli);
  free(up_repli);
  free(up_prod);
  free(nci_sum);
  free(sum_prod);
}

void read_L1(char* input, TestCase test_case) {
  FILE* fp_qp = NULL;
  fopen_s(&fp_qp, input, "r");
  if (fp_qp == NULL) {
    fprintf(stderr, "Failed to open msb file %s\n", input);
    return;
  }
  fseek(fp_qp, 0L, SEEK_END);
  size_t bytes_to_read = ftell(fp_qp);
  rewind(fp_qp);

  uint8_t* buffer = (uint8_t*)malloc(bytes_to_read);
  size_t bytesRead;
  bytesRead = fread(buffer, 1, bytes_to_read, fp_qp);
  if (bytesRead != bytes_to_read) { printf("Error parsiing data\n"); }
  fclose(fp_qp);

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
  dtime -= 18; // leap seconds
  printf("week %d tow %f \n", (int)floor(dtime / 604800.0), dtime - floor(dtime / 604800.0) * 604800.0);
  int msec = (int)round((sec - floor(sec)) * 1000);
  uint32_t nanos = U3(&buffer[18]); // 21
  printf("%4d-%02d-%02d-%02d-%02d-%02d msec=%d nanos=%d \n", yr, mon, day, hr, min, sec, msec, nanos);
  int payload_size = bytes_to_read - hdrlen;
  printf("size (total bytes - hdrs) %d \n", payload_size);
  printf("Working on %s\n", input);


  c32* samples = (c32*)malloc(payload_size * sizeof(c32));
  DecodeOrsIQCplx(&buffer[hdrlen], payload_size / 2, samples);
  free(buffer);

  // Test the quasi pilot generation
  int min_idx = 0; int loc_cnt = 0; int missed = 0;
  float min_val = 1e5;
//#define FFT_QP_SIZE 4096 * 1 // was 2
  int fft_size = test_case.fft_size;
  float chipping_rate = 1.023e6; // chips per sec
  int window = 4; //  window/2 ms either side of center 
//#define SPC 4 // samples per chip was 3
  int spc = test_case.spc;
  int nci = payload_size / (L1C_CODE_LEN * spc);
  double IF_OFFSET = test_case.IF;// 5400;// 5400;// 1e6 + 4100; // March 31 540 April 1st 4100
  printf("Using %d len FFT and %d SPC and window size %d (with %d NCI avail) IF_OFFSET=%d\n", fft_size, spc, window, nci, (int)IF_OFFSET);
  acq_struct prn2acq[30] = { 0 }; int cnt = 0;
  int num_prns = read_assist(input, prn2acq);


  ///////////////////// main prn loop ////////////////////////////////

  // Compute circular correlation C_k(τ) = FFT^-1{ FFT[signal] · conj(FFT[replica]) }.
  c32* replica  = (c32*)malloc(sizeof(c32) * fft_size);
  c32* fft_repl = (c32*)malloc(sizeof(c32) * fft_size);
  c32* fft_data = (c32*)malloc(sizeof(c32) * fft_size);
  c32* fft_sum  = (c32*)malloc(sizeof(c32) * fft_size);
  c32* fft_prod = (c32*)malloc(sizeof(c32) * fft_size);
  float* nci_sum = (float*)malloc(sizeof(float) * fft_size);
  int* prn_code = (int*)malloc(sizeof(int) * L1C_CODE_LEN * spc);
  for (int cnt = 0; cnt < num_prns; cnt++) {
    // only GAL=2 right PRNs
    if (prn2acq[cnt].constel == 2) { continue; }
    //if (prn2acq[cnt].constel == 1 || (prn2acq[cnt].prn != 13 && prn2acq[cnt].prn != 23)) { continue; }
    //if (prn2acq[cnt].constel == 1 || (prn2acq[cnt].prn != 25 && prn2acq[cnt].prn != 36)) { continue; }
    memset(fft_repl, 0, sizeof(c32) * fft_size);
    memset(replica, 0, sizeof(c32) * fft_size);
    
    getCode(L1C_CODE_LEN, spc, prn2acq[cnt].prn, prn_code);
    double doppler = IF_OFFSET + (prn2acq[cnt].doppler);
    make_replica(prn_code, replica, doppler, L1C_CODE_LEN * spc, chipping_rate * spc);
    memcpy(fft_repl, replica, sizeof(c32) * L1C_CODE_LEN);

    fft_c32(fft_size, fft_repl, true);

    int found[50] = { 0 };
    int last_location = 0;
    for (int center = window / 2; center <= nci - window / 2 - 2; center++) {
      memset(fft_sum, 0, sizeof(c32) * fft_size);
      memset(nci_sum, 0, sizeof(float) * fft_size);
      for (int windex = center - window / 2; windex < center + window / 2; windex += 1) {
        memset(fft_data, 0, sizeof(c32) * fft_size);
        memcpy(fft_data, &samples[L1C_CODE_LEN * spc * windex], sizeof(c32) * spc * (L1C_CODE_LEN));
        memcpy(&fft_data[L1C_CODE_LEN * spc], &samples[L1C_CODE_LEN * spc * windex], sizeof(c32) * (fft_size - spc * L1C_CODE_LEN));
        fft_c32(fft_size, fft_data, true); // forward FFT

        for (int k = 0; k < fft_size; k++) { // accumulate pt-wise * with conj of replica
          if (windex < center) { // cheaper method     //get_minus_conj
            fft_sum[k] = add(fft_sum[k], mult(fft_data[k], get_conj(fft_repl[k])));
          } else {
            fft_sum[k] = add(fft_sum[k], mult(fft_data[k], get_conj(fft_repl[k])));
          }
          fft_prod[k] = mult(fft_data[k], get_conj(fft_repl[k]));
        }
        // parallel track prod and IFFT then square and sum
        fft_c32(fft_size, fft_prod, false);
        for (int k = 0; k < fft_size; k++) { nci_sum[k] += mag(fft_prod[k]); }

      } // for windex 

      // used to have the IFFT here
      fft_c32(fft_size, fft_sum, false); // IFFT // cheaper method

      FILE* fp_out = NULL;
      top2_pks peaks;
      find_top2_peaks_cplx(fft_sum, L1C_CODE_LEN * spc, 10, &peaks, NULL);
      top2_pks peaks2;
      find_top2_peaks_real(nci_sum, L1C_CODE_LEN * spc, 10, &peaks2, fp_out);
      bool isMax = findMax(peaks.val1);
      double ang = atan2(fft_sum[peaks.idx1].i, fft_sum[peaks.idx1].r) * 57.2957795;
      float ratio = peaks.val1 / peaks.val2;
      float nci_ratio = peaks2.val1 / peaks2.val2;
      // nix if last max was less than 20 points ago 
      if (isMax && fabs(last_location - center + 8) > 20) {
        found[loc_cnt++] = center - 8;
        last_location = center - 8;
      }
      char constel = (prn2acq[cnt].constel == 1) ? 'G' : 'E';
      printf("prn, %c%02d, center,%03d, max,%6.1f,p1:,%4d,p2:,%4d,r1,%4.2f,r2,%4.2f,ang,%4.2f,dop,%6.1f\n",
        constel,prn2acq[cnt].prn, center, peaks.val1, peaks.idx1, peaks2.idx1, ratio, nci_ratio, ang, doppler);
    } // for center
  }

  //for (int i = 0; i < loc_cnt; i++) { printf("Bit transition at %d ms \n", found[i]); }
  //printf("BTs: %d \n", loc_cnt);
  free(samples); free(replica); free(nci_sum); free(prn_code);
  free(fft_data); free(fft_sum); free(fft_repl); free(fft_prod);
}

// make sure to edit the line if (prn2acq[cnt].constel == 1 to allow the
// right prns are used or use all
void read_QP(char* input, TestCase test_case) {
  FILE* fp_qp = NULL;
  fopen_s(&fp_qp, input, "r");
  if (fp_qp == NULL) {
    fprintf(stderr, "Failed to open msb file %s\n", input);
    return;
  }
  fseek(fp_qp, 0L, SEEK_END);
  size_t bytes_to_read = ftell(fp_qp);
  rewind(fp_qp);

  uint8_t* buffer = (uint8_t*)malloc(bytes_to_read);
  size_t bytesRead;
  bytesRead = fread(buffer, 1, bytes_to_read, fp_qp);
  if (bytesRead != bytes_to_read) { printf("Error parsiing data\n"); }
  fclose(fp_qp);

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
  dtime -= 18; // leap seconds
  printf("week %d tow %f \n", (int)floor(dtime / 604800.0), dtime - floor(dtime / 604800.0) * 604800.0);
  int msec = (int)round((sec - floor(sec)) * 1000);
  uint32_t nanos = U3(&buffer[18]); // 21
  printf("%4d-%02d-%02d-%02d-%02d-%02d msec=%d nanos=%d \n", yr, mon, day, hr, min, sec, msec, nanos);
  int payload_size = bytes_to_read - hdrlen;
  printf("size (total bytes - hdrs) %d \n", payload_size);
  printf("Working on %s\n", input);
  

  c32* samples = (c32*)malloc(payload_size * sizeof(c32));
  DecodeOrsIQCplx(&buffer[hdrlen], payload_size / 2, samples);
  free(buffer);

  // Test the quasi pilot generation
  int min_idx = 0; int loc_cnt = 0; int missed = 0;
  float min_val = 1e5;
#define FFT_QP_SIZE 512 * 2 // was 2
  float chipping_rate = 5.115e6; // chips per sec
  int window = 20; //  window/2 ms either side of center 
  int spc = test_case.spc; // samples per chip was 3
  int nci = payload_size / (E5_QP_CODE_LEN * spc);
  double IF_OFFSET = test_case.IF;// +5400;// 1e6 + 4100;
  printf("Using %d len FFT and %d SPC and window size %d (with %d NCI avail) IF_OFFSET=%d spc=%d\n", 
    FFT_QP_SIZE, SPC, window, nci, (int)IF_OFFSET, spc);
  acq_struct prn2acq[30] = { 0 }; int cnt = 0;
  int num_prns = read_assist(input, prn2acq);
  
  
  ///////////////////// main prn loop ////////////////////////////////

  // Compute circular correlation C_k(τ) = FFT^-1{ FFT[signal] · conj(FFT[replica]) }.
  c32* replica  = (c32*)malloc(sizeof(c32) * FFT_QP_SIZE);
  c32* fft_repl = (c32*)malloc(sizeof(c32) * FFT_QP_SIZE);
  c32* fft_data = (c32*)malloc(sizeof(c32) * FFT_QP_SIZE);
  c32* fft_sum  = (c32*)malloc(sizeof(c32) * FFT_QP_SIZE);
  c32* fft_prod = (c32*)malloc(sizeof(c32) * FFT_QP_SIZE);
  float* nci_sum = (float*)malloc(sizeof(float) * FFT_QP_SIZE);
  int* prn_code = (int*)malloc(sizeof(int) * E5_QP_CODE_LEN * spc);
  double fac = 1176.45 / 1575.42;// values were quoted at L1, need at L5
  for (int cnt = 0; cnt < num_prns; cnt++) {
    // only GAL=2 right PRNs
    if (prn2acq[cnt].constel == 1) { continue; } // only do Gal
    bool proceed = false;
    for (int i = 0; i < 9; i++) {
      if (prn2acq[cnt].prn == test_case.good_svs[i]) { proceed = true; break; }
    }
    if (!proceed) { continue; }
    memset(fft_repl, 0, sizeof(c32) * FFT_QP_SIZE);
    memset(replica, 0, sizeof(c32) * FFT_QP_SIZE);
    
    getE5_QPCode(E5_QP_CODE_LEN, spc, prn2acq[cnt].prn, prn_code);
    double doppler = IF_OFFSET + (prn2acq[cnt].doppler) * fac;
    make_replica(prn_code, replica, doppler, E5_QP_CODE_LEN * spc, chipping_rate * SPC);
    memcpy(fft_repl, replica, sizeof(c32) * E5_QP_CODE_LEN);
    
    fft_c32(FFT_QP_SIZE, fft_repl, true);

    int found[50] = { 0 };
    int last_location = 0;
    for (int center = window / 2; center <= nci - window / 2 - 2; center++) {
      memset(fft_sum, 0, sizeof(c32) * FFT_QP_SIZE);
      memset(nci_sum, 0, sizeof(float) * FFT_QP_SIZE);
      for (int windex = center - window / 2; windex < center + window / 2; windex += 1) {
        memset(fft_data, 0, sizeof(c32) * FFT_QP_SIZE);
        memcpy(fft_data, &samples[E5_QP_CODE_LEN * spc * windex], sizeof(c32) * spc * (E5_QP_CODE_LEN));
        memcpy(&fft_data[E5_QP_CODE_LEN * spc], &samples[E5_QP_CODE_LEN * spc * windex], sizeof(c32) * (FFT_QP_SIZE - spc * E5_QP_CODE_LEN));
        fft_c32(FFT_QP_SIZE, fft_data, true); // forward FFT

        for (int k = 0; k < FFT_QP_SIZE; k++) { // accumulate pt-wise * with conj of replica
          if (windex < center) { // cheaper method     //get_minus_conj
            fft_sum[k] = add(fft_sum[k], mult(fft_data[k], get_conj(fft_repl[k])));
          } else {
            fft_sum[k] = add(fft_sum[k], mult(fft_data[k], get_conj(fft_repl[k])));
          }
          fft_prod[k] = mult(fft_data[k], get_conj(fft_repl[k]));
        }
        // parallel track prod and IFFT then square and sum
        fft_c32(FFT_QP_SIZE, fft_prod, false);
        for (int k = 0; k < FFT_QP_SIZE; k++) { nci_sum[k] += mag(fft_prod[k]); }

      } // for windex 

      // used to have the IFFT here
      fft_c32(FFT_QP_SIZE, fft_sum, false); // IFFT // cheaper method

      FILE* fp_out = NULL;
      if (false) {} //center == 50) {  errno_t er = fopen_s(&fp_out, "C:/Python/nci_sum4.csv", "w"); }
      top2_pks peaks;
      find_top2_peaks_cplx(fft_sum, E5_QP_CODE_LEN * spc, 10, &peaks, fp_out);
      top2_pks peaks2;
      find_top2_peaks_real(nci_sum, E5_QP_CODE_LEN * spc, 10, &peaks2, fp_out);
      bool isMax = findMax(peaks.val1);
      double ang = atan2(fft_sum[peaks.idx1].i, fft_sum[peaks.idx1].r) * 57.2957795;
      float ratio = peaks.val1 / peaks.val2;
      float ratio2 = peaks2.val1 / peaks2.val2;
      // nix if last max was less than 20 points ago 
      if (isMax && fabs(last_location - center + 8) > 20) {
        found[loc_cnt++] = center - 8;
        last_location = center - 8;
      }
      char constel = (prn2acq[cnt].constel == 1) ? 'G' : 'E';
      printf("prn,%c%02d, center, %03d,max,%5.1f, p1,%4d,p2,%4d,r1,%4.2f,r2,%4.2f,ang,%4.2f,dop,%6.1f\n", 
        constel, prn2acq[cnt].prn, center, peaks.val1, peaks.idx1, peaks2.idx1, ratio, ratio2, ang, doppler);
    } // for center
  }

  //for (int i = 0; i < loc_cnt; i++) { printf("Bit transition at %d ms \n", found[i]); }
  //printf("BTs: %d \n", loc_cnt);
  free(samples); free(replica); free(nci_sum); free(prn_code);
  free(fft_data); free(fft_sum); free(fft_repl); free(fft_prod);
}

// try differential coherent integration. Insensitive to bit transitions
void test_quasi_diff_pilot() {
  srand((unsigned int)time(NULL)); // randomise seed
  // Test the quasi pilot generation

  int min_idx = 0;
  int loc_cnt = 0;
  float min_val = 1e5;

  int locations[50] = { 50, 100, 150, 200, 250, 300, 350 };
  int window = 20; // 2 * window ms either side of center (window>=5 does not work with noise window > 7 does not work with no noise)
  int nci = 400;
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
      errno_t er = fopen_s(&fp_out, "diff_corr.csv", "w");
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
      errno_t er = fopen_s(&fp_out, "nci_sum4.csv", "w");
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
      errno_t er = fopen_s(&fp_out, "nci_sum4.csv", "w");
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

// shift left
void shiftDown(c32* array,int gap, int size) {
  for (int i = 0; i < size; i++) {
    array[i] = array[i + gap];
  }
}

void test_quasi_pilot_330(results_s* results) {
  // todo adapt for live signals eg skip samples
  // for 1 sample/chip (but 
  // buffer fft (short ones) for moving window to find BT
  // refine and double check the exact location of the 
  auto now = std::chrono::system_clock::now();
  // Convert the time point to duration since epoch
  auto duration_since_epoch = now.time_since_epoch();
  srand((unsigned int)std::chrono::duration_cast<std::chrono::milliseconds>(duration_since_epoch).count());
  // Test the quasi pilot generation
  int min_idx = 0; int loc_cnt = 0; int missed = 0;
  float min_val = 1e5;
#define FFT_QP_SIZE 512 * 2 // was 2
  float chipping_rate = 5.115e6; // chips per sec
  // note can do code-phase error stats or BTT detection stats but not both at once!!!! use {-1} for no BTTs
  int locations[50] = { 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000 };// { 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280 };
  int found[50] = { 0 };
  int window = 70; //  window/2 ms either side of center 
  int nci = 1100;
#define SPC 2 // samples per chip was 3
  int   len = E5_QP_CODE_LEN * SPC * nci; // 4 samples per chip and 100 ms
  int   c_phase = 300;// 1;// 329; // which chip to set the code phase to
  int   prn1 = 14, prn2 = 18;
  float dop1 = 2000, dop2 = -2000;
  float dop_error = 500;// 10; // full 2*250 Hz error in wipeoff
  float dop_err_rate = 0.6;// 0.6;// 0.6;//Hz per ms
  float pwr = -128.5;// -129.5;// -128.5;// -128.75;// dBm signal power
  float N0 = -174.0 + 2.5;// dBm/Hz thermal noise ktb density + noise figure
  float cn0 = pwr - N0; // dBm/Hz pgae 14 Galileo_OS_SIS_ICD_v2.2
  float sigma = sqrt((1.0 * chipping_rate * SPC) / (2.0f* pow(10.0,cn0/10.0)));
  //printf("sigma %f \n", sigma);
  //sigma = 0.1;// 12.7;// 15.98;// 3.5;// 3.5; // noise level 15->6*31
  int   num_errors = 0; int num_tries = 0;
  c32*  out    = (c32*)malloc(len * sizeof(c32));
  c32*  cache  = (c32*)malloc(FFT_QP_SIZE * window * sizeof(c32));
  memset(cache, 0, sizeof(c32) * FFT_QP_SIZE * window);
  if (out == NULL) { fprintf(stderr, "Memory allocation failed for 100 ms I&Q array.\n"); return; }
  int* prn_c1  = (int*)malloc(sizeof(int) * E5_QP_CODE_LEN * SPC);
  int* prn_c2  = (int*)malloc(sizeof(int) * E5_QP_CODE_LEN * SPC);
  int* prn_c3  = (int*)malloc(sizeof(int) * E5_QP_CODE_LEN * SPC);
  c32* replica = (c32*)malloc(sizeof(c32) * E5_QP_CODE_LEN * SPC);
  //float* nci_sum = (float*)malloc(sizeof(float) * FFT_QP_SIZE);
  memset(prn_c1 , 0, sizeof(int) * E5_QP_CODE_LEN * SPC);
  memset(prn_c2 , 0, sizeof(int) * E5_QP_CODE_LEN * SPC);
  memset(replica, 0, sizeof(c32) * E5_QP_CODE_LEN * SPC);
  getE5_QPCode(E5_QP_CODE_LEN, SPC, prn1, prn_c1);
  getE5_QPCode(E5_QP_CODE_LEN, SPC, prn2, prn_c2);

  memcpy(prn_c3, prn_c1, sizeof(int) * E5_QP_CODE_LEN * SPC); // unshifted version for later

  make_replica(prn_c3, replica, dop1 + dop_error, E5_QP_CODE_LEN * SPC, chipping_rate * SPC);
  rotate_fwd(prn_c1, E5_QP_CODE_LEN * SPC, c_phase); // now advance code-phase  
  
  int sign2 = 1; // sign applied a posteriori after finding BTT
  int num_btts = -1;
  for (int i = 0; i < nci; i++) {
    for (int j = 0; j < 50; j++) {
      if (locations[j] == i) { sign2 *= -1; num_btts++;  break; } // change sign at the bit transitions
    }
    // offset doppler by 250 Hz and add a residual doppler ramp of 0.1 Hz per ms
    mix_two_prns_oversampled_per_prn(prn_c1, prn_c2, dop1 + i * dop_err_rate, dop2 - i * dop_err_rate, 
      PI / 2, 0,&out[E5_QP_CODE_LEN * SPC * i], E5_QP_CODE_LEN * SPC, chipping_rate * SPC, sigma, sign2); // was 2.31 for -128.5 dBm 3.1 for -131.5
  }
  free(prn_c1); free(prn_c2);

  // Compute circular correlation C_k(τ) = FFT^-1{ FFT[signal] · conj(FFT[replica]) }.
  c32* fft_repl = (c32*)malloc(sizeof(c32) * FFT_QP_SIZE);
  c32* fft_data = (c32*)malloc(sizeof(c32) * FFT_QP_SIZE);
  c32* fft_sum  = (c32*)malloc(sizeof(c32) * FFT_QP_SIZE);
  c32* fft_prod = (c32*)malloc(sizeof(c32) * FFT_QP_SIZE);
  memset(fft_repl, 0, sizeof(c32) * FFT_QP_SIZE);
  if (fft_data == NULL || fft_repl == NULL) { printf("Error allocating fft_data or fft_prod\n"); return; }

  memcpy(fft_repl, replica, sizeof(c32) * E5_QP_CODE_LEN);
  free(replica);
  fft_c32(FFT_QP_SIZE, fft_repl, true);
  auto start = std::chrono::high_resolution_clock::now();////////////////////////////////////////
  int last_location = 0;
  int cptr = 0;
  for (int center = window / 2; center <= nci - window / 2 - 2; center++) {
    if (center == window / 2) { // priming cache speeds things up 4x
      for (int windex = center - window / 2; windex < center + window / 2; windex += 1) {
        memcpy(&cache[cptr * FFT_QP_SIZE], &out[E5_QP_CODE_LEN * SPC * windex], sizeof(c32) * SPC * (E5_QP_CODE_LEN));
        memcpy(&cache[cptr * FFT_QP_SIZE + E5_QP_CODE_LEN * SPC], &out[E5_QP_CODE_LEN * SPC * windex], sizeof(c32) * (FFT_QP_SIZE - SPC * E5_QP_CODE_LEN));
        fft_c32(FFT_QP_SIZE, &cache[cptr * FFT_QP_SIZE], true); // forward FFT
        cptr = (cptr+1) % window;
      }
    } else {
      int indx = window - 1;// center + window / 2 - 1; // index of new window member
      int indx2 = center + (window / 2) - 1; // index of old window member to be removed
      //memcpy(&cache[cptr * FFT_QP_SIZE], &out[E5_QP_CODE_LEN * SPC * indx], sizeof(c32) * SPC * (E5_QP_CODE_LEN));
      //memcpy(&cache[cptr * FFT_QP_SIZE + E5_QP_CODE_LEN * SPC], &out[E5_QP_CODE_LEN * SPC * indx], sizeof(c32) * SPC * (FFT_QP_SIZE - E5_QP_CODE_LEN));
      //fft_c32(FFT_QP_SIZE, &cache[cptr], true); // forward FFT
      //cptr = (cptr + 1) % window;

      //printf("B %f %f %f %f %f %f %f %f \n", cache[0 * FFT_QP_SIZE].r, cache[1 * FFT_QP_SIZE].r, cache[64 * FFT_QP_SIZE].r, cache[65 * FFT_QP_SIZE].r, cache[66 * FFT_QP_SIZE].r, cache[67 * FFT_QP_SIZE].r, cache[68 * FFT_QP_SIZE].r, cache[69 * FFT_QP_SIZE].r);
      shiftDown(cache, FFT_QP_SIZE, FFT_QP_SIZE * window);
      //printf("A %f %f %f %f %f %f %f %f \n", cache[0 * FFT_QP_SIZE].r, cache[1 * FFT_QP_SIZE].r, cache[64 * FFT_QP_SIZE].r, cache[65 * FFT_QP_SIZE].r, cache[66 * FFT_QP_SIZE].r, cache[67 * FFT_QP_SIZE].r, cache[68 * FFT_QP_SIZE].r, cache[69 * FFT_QP_SIZE].r);
      memcpy(&cache[indx * FFT_QP_SIZE], &out[E5_QP_CODE_LEN * SPC * indx2], sizeof(c32) * SPC * (E5_QP_CODE_LEN));
      memcpy(&cache[indx * FFT_QP_SIZE + E5_QP_CODE_LEN * SPC], &out[E5_QP_CODE_LEN * SPC * indx2], sizeof(c32) * (FFT_QP_SIZE - SPC * E5_QP_CODE_LEN));
      fft_c32(FFT_QP_SIZE, &cache[indx * FFT_QP_SIZE], true); // forward FFT
    }
 
    //memset(nci_sum, 0, sizeof(float) * FFT_QP_SIZE);
    memset(fft_sum, 0, sizeof(c32) * FFT_QP_SIZE);
    int tmp_idx = 0;
    for (int windex = center - window / 2; windex < center + window / 2; windex +=1) {
      //memset(fft_data, 0, sizeof(c32) * FFT_QP_SIZE);
      //memcpy(fft_data, &out[E5_QP_CODE_LEN * SPC * windex], sizeof(c32) * SPC * (E5_QP_CODE_LEN));
      //memcpy(&fft_data[E5_QP_CODE_LEN * SPC], &out[E5_QP_CODE_LEN * SPC * windex], sizeof(c32) * (FFT_QP_SIZE - SPC * E5_QP_CODE_LEN)); 
      //fft_c32(FFT_QP_SIZE, fft_data, true); // forward FFT

      //for (int lp = 0; lp < FFT_QP_SIZE; lp++) { 
      //  if (fft_data[lp].i != cache[tmp_idx * FFT_QP_SIZE + lp].i) {
      //    printf("cache and fftp disagree at %d \n", lp);
      //  }
      //}
     
      for (int k = 0; k < FFT_QP_SIZE; k++) { // accumulate pt-wise * with conj of replica
        c32 fft_data_ins = cache[tmp_idx * FFT_QP_SIZE + k];// cache[tmp_idx * FFT_QP_SIZE + k]; fft_data[k];
        if (windex < center) { // cheaper method           //get_minus_conj
          fft_sum[k] = add(fft_sum[k], mult(fft_data_ins, get_minus_conj(fft_repl[k])));
        } else {
          fft_sum[k] = add(fft_sum[k], mult(fft_data_ins, get_conj(fft_repl[k])));
        }
        //fft_prod[k] = mult(fft_data_ins, get_conj(fft_repl[k]));
      }

      // parallel track prod and IFFT then square and sum
      //fft_c32(FFT_QP_SIZE, fft_prod, false);
      //for (int k = 0; k < FFT_QP_SIZE; k++) { nci_sum[k] += mag(fft_prod[k]); }

      tmp_idx = (tmp_idx + 1) % window;
    } // for windex 

    // used to have the IFFT here
    fft_c32(FFT_QP_SIZE, fft_sum, false); // IFFT // cheaper method

    FILE* fp_out = NULL;
    if (false) {} //center == 50) {  errno_t er = fopen_s(&fp_out, "C:/Python/nci_sum4.csv", "w"); }
    top2_pks peaks;
    find_top2_peaks_cplx(fft_sum, E5_QP_CODE_LEN * SPC, 4, &peaks, fp_out);
    top2_pks peaks2 = {0,0,1,0};
    //find_top2_peaks_real(nci_sum, E5_QP_CODE_LEN * SPC, 4, &peaks2, fp_out);

 
    bool isMax = findMax(peaks.val1);
    int best_code = 0; float best_pwr = 0; float best_dop = 0; double phase_deg = 0;
    if (false) {
      phase_deg = recover_prn_phase_deg_with_doppler(&out[E5_QP_CODE_LEN * SPC * (center )], E5_QP_CODE_LEN * SPC,
        prn_c3, peaks.idx1, dop1, chipping_rate * SPC, &best_code, &best_dop, &best_pwr);
    }
    bool isBT = checkBT(center, locations, num_btts);
    if (peaks.idx1 != c_phase && isBT) { num_errors++; printf("code %d \n", peaks.idx1); }
    //if (peaks2.idx1 != c_phase && isBT) { num_errors++; printf("code2 %d \n", peaks2.idx1); }
    
    double ang = atan2(fft_sum[peaks.idx1].i, fft_sum[peaks.idx1].r) * 57.2957795;
    
    float ratio = peaks.val1 / peaks.val2;
    float ratio2 = peaks2.val1 / peaks2.val2;
    // nix if last max was less than 20 points ago 
    if (isMax && fabs(last_location - center + 8) > 20) {
      found[loc_cnt++] = center - 8;
      last_location = center - 8;
    }
    if (center > last_location + 2*100) {
      last_location = locations[loc_cnt++];
      missed++;
      //printf("padd %d\n", last_location);
    }
    //printf("center,%03d, max,%6.1f,p1,%d,p2,%d,bestPwr,%4.2f,ratio,%4.2f,ratio2,%4.2f,ang,%4.2f,cnt,%d,phase,%f\n", 
    //  center, peaks.val1, peaks.idx1, peaks2.idx1, best_pwr, ratio, ratio2, ang, isMax, phase_deg);
  } // for center
  num_tries++;
  printf("Found %d btt there are %d btt \n", num_btts, num_btts - missed);
  for (int i = 0; i < num_btts; i++) {
    printf("Diff %d %d %d \n", i, locations[i] , found[i]);
    int diff = abs(locations[i] - found[i]);
    if (diff < 100  && found[i] != 0) {
      results->differences[diff]++;
    }
  }
  
  results->num_errors += num_errors;
  results->num_trials += num_btts;// num_tries;
  if ((loc_cnt - num_btts) > 0) { results->btt_false_alarms+= (loc_cnt - num_btts); }
  if (missed) { results->btt_missed_detects += missed; }
  
  free(prn_c3); free(cache); //free(nci_sum);
  free(out); free(fft_data); free(fft_sum); free(fft_repl); free(fft_prod);// free(fft_prev); //free(mag_sum);
}

void test_quasi_pilot2() {
  srand((unsigned int)time(NULL)); // randomise seed
  // Test the quasi pilot generation

  int min_idx = 0;
  int loc_cnt = 0;
  float min_val = 1e5;

  int locations[50] = { 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280 };
  int window = 4; // 2 * window ms either side [window works from 2 - 37]
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

  for (int argi = 1; argi < argc; argi++)
  {
    if (strcmp(argv[argi], "L1") == 0) optionSig |= OPTION_L1;
    if (strcmp(argv[argi], "E1") == 0) optionSig |= OPTION_E1;
    if (strcmp(argv[argi], "L5") == 0) optionSig |= OPTION_L5;
    if (strcmp(argv[argi], "E5") == 0) optionSig |= OPTION_E5a;
    if (strcmp(argv[argi], "QP") == 0) optionSig |= OPTION_QP;
    if (strcmp(argv[argi], "QX") == 0) optionSig |= OPTION_QX;

    if (strcmp(argv[argi], "mscoh") == 0)
    {
      argi++; 
      optionMsCoh = atoi(argv[argi]);
    }

    if (strcmp(argv[argi], "msps") == 0)
    {
      argi++;
      optionMsps = atoi(argv[argi]);
    }

    if (strcmp(argv[argi], "dunc") == 0)
    {
      argi++;
      optionDUnc = atoi(argv[argi]);
    }
  }


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

  if (0) {
     read_E5A((char*)"C:/work/support/eric/E5/t14/G_2024_10_21_22_29_43.047.csv");
    //read_E5A((char*)"C:/work/Baseband/TestData/E5/t14/G_2024_10_21_22_32_30.500_resampled_16368Hz.csv");
    //read_E5A((char*)"C:/work/Baseband/TestData/E5/t14/G_2024_10_21_22_29_43.047_resampled_16368Hz.csv");
    return 0;
  }

  if (1) {

    // G_2025_09_03_23_04_45.ors G_2025_09_03_23_04_56.ors G_2025_09_03_23_05_33.ors G_2025_09_03_23_04_56.ors G_2025_09_03_23_12_45.ors G_2025_09_03_23_19_10.ors
    char line[256] = { 0 };
    char path[256] = { 0 };
    //const char* folder = "C:/work/support/eric/E5/t14";
    //const char* folder = "C:/work/support/esa/qp/data/t11/2026-04-01/L5-1_16";
    //const char* folder = "C:/work/support/esa/qp/data/t11/2026-04-01/L5-1_10";
    //const char* folder = "C:/work/baseband/utilities/2026-04-07/L5-1_10"; //old
    const char* folder = "C:/work/Baseband/Utilities/2026-04-20/L5_10_87";  
    //const char* folder = "C:/work/support/esa/qp/data/t12/2026-04-07/L5-5_05_42";
    //const char* folder = "C:/work/support/esa/qp/data/t12/2026-04-07/L5_05_42";
    //const char* folder = "C:/work/support/esa/qp/data/t12/2026-04-07/L5_05_87";
    //const char* folder = "C:/work/support/esa/qp/data/t12/2026-04-07/L5-1_05_87";
    //const char* folder = ".";
    sprintf_s(path, 256, "%s/results.csv", folder);
    errno_t er =  fopen_s(&fpResults, path, "w");
    sprintf_s(path, 256, "%s/csvFiles.txt", folder);
    FILE* fp_in = NULL; //output file
    er = fopen_s(&fp_in, path, "r");
    while (fgets(line, sizeof(line), fp_in)) {
      *strchr(line, '\n') = 0;
      sprintf_s(path, 256, "%s/%s", folder, line);
      printf("Processing %s\n", path);
      //read_E5A(path);
      if (optionSig & (OPTION_L1 | OPTION_E1))
        read_L1E1(path);
      else
        read_L5E5AE5QP(path);
      if (fpResults != NULL) fflush(fpResults);
      break;
    }
    if (fp_in != NULL) { fclose(fp_in); }
    if (fpResults != NULL) { fclose(fpResults); }

    //read_ors((char*)"C:/work/Baseband/TestData/100ms/bw25/G_2025_09_03_23_05_33.ors");
    //read_ors((char*)"C:/work/Baseband/TestData/100ms/bw25/G_2025_09_03_23_04_56.ors");
    //read_ors((char*)"C:/work/Baseband/TestData/100ms/bw25/G_2025_09_03_23_04_45.ors");
    //read_ors((char*)"C:/work/Baseband/TestData/G_2025_06_05_22_11_26.ors");
    return 0;
  }
  if (false) { // read_L1 data .ors with .ast 
    TestName day = APRIL1;
    // choices are march31_cases april1_cases april7_cases
    for (int j = 0; j < l1_test_wrappers[day].num_folders; j++) {
      for (int i = 0; i < l1_test_wrappers[day].tests[j].num_files; i++) {
        char buff[256] = "C:/work/Baseband/Utilities/"; // path to folder containing folders containing the .ors and .ast files
        strcat_s(buff, l1_test_wrappers[day].tests[j].name);
        strcat_s(buff, l1_test_wrappers[day].tests[j].files[i]);
        read_L1(buff, l1_test_wrappers[day].tests[j]);
      }
      printf("==========================================================\n");
    }
    return 0;
  }

  if (false) { // read_QP data .ors with .ast
    TestName day = APRIL7;
    // choices are march31_cases april1_cases april7_cases
    for (int j = 0; j < test_wrappers[day].num_folders; j++) {
      for (int i = 0; i < test_wrappers[day].tests[j].num_files; i++) {
        char buff[256] = "C:/work/Baseband/Utilities/"; // path to folder containing folders containing the .ors and .ast files
        strcat_s(buff, test_wrappers[day].tests[j].name);
        strcat_s(buff, test_wrappers[day].tests[j].files[i]);
        read_QP(buff, test_wrappers[day].tests[j]);
      }
      printf("==========================================================\n");
    }
    return 0;
  }

  if (0) {
    read_E5A((char*)"C:/work/Baseband/TestData/E5/t14/G_2024_10_21_22_29_43.047.csv");
    //read_E5A((char*)"C:/work/Baseband/TestData/E5/t14/G_2024_10_21_22_32_30.500_resampled_16368Hz.csv");
    //read_E5A((char*)"C:/work/Baseband/TestData/E5/t14/G_2024_10_21_22_29_43.047_resampled_16368Hz.csv");
    return 0;
  }

  if (0) {
    read_L1_CSV((char*)"C:/Python/out-1bit_1spc_1bit.csv");// "C:/work/Baseband/TestData/E5/t14/G_2024_10_21_22_29_43.047.csv");
    return 0;
  }

  if (0) {
    sim_L1();
    return 0;
  }

  if (0) {
    // G_2025_09_03_23_04_45.ors G_2025_09_03_23_04_56.ors G_2025_09_03_23_05_33.ors G_2025_09_03_23_04_56.ors G_2025_09_03_23_12_45.ors G_2025_09_03_23_19_10.ors
    FILE* fp_in = NULL; //output file
    //char directory[_MAX_PATH] = "C:/work/support/eric/data/100ms/bw25";
    char directory[_MAX_PATH] = "C:/work/support/esa/qp/data/t11/2026-04-01/L1_04";
    char filename[_MAX_PATH];
    sprintf_s(filename, "%s/orsfiles.txt", directory);

    errno_t er = fopen_s(&fp_in, filename, "r");
    char line[256] = {0};
    while (fgets(line, sizeof(line), fp_in)) {
      for (int i = 0; i < strlen(line); i++) {
        if (line[i] == '\\') { line[i] = '/'; } // swap delimiter
        if (line[i] == '\n') { line[i] = NULL; } // swap terminator 
      }
      //printf("about to process %s\n", line);
      read_ors(directory,line);
    }
    if (fp_in != NULL) { fclose(fp_in); }
    //read_ors((char*)"C:/work/Baseband/TestData/100ms/bw25/G_2025_09_03_23_05_33.ors");
    //read_ors((char*)"C:/work/Baseband/TestData/100ms/bw25/G_2025_09_03_23_04_56.ors");
    //read_ors((char*)"C:/work/Baseband/TestData/100ms/bw25/G_2025_09_03_23_04_45.ors");
    //read_ors((char*)"C:/work/Baseband/TestData/G_2025_06_05_22_11_26.ors");
    return 0;
  }

  if (0) {
    test_quasi_pilot();
    return 0;
  }

  if (0) {
    //test_quasi_diff_pilot();
    //test_quasi_pilot_330Up();
    
    results_s results = {0};
    auto start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < 100; i++) {
      test_quasi_pilot_330(&results);
    }
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    printf("Time taken: %d \n",(long)duration.count());
    printf("Total tries: %d, Total errors: %d, code errors: %d #fa %d #md %d\n", results.num_trials, (results.btt_false_alarms + results.btt_missed_detects),  results.num_errors,
      results.btt_false_alarms, results.btt_missed_detects);
    for (int i = 0; i < 100; i++) {
      if (results.differences[i] == 0) { continue; }
      printf("Num with diff %d is %d \n", i, results.differences[i]);
    }
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


