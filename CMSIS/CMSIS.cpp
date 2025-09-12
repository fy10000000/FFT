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
#include "transform_functions.h"
#include "gnss_codes.h"


#define Q31_MIN   ((int32_t)0x80000000)  // -2147483648
#define Q31_MAX   ((int32_t)0x7FFFFFFF)  //  2147483647
#define Q15_MIN   (-32768)
#define Q15_MAX   ( 32767)

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


void read_ors(char* input) {
  FILE* fp_msb = NULL;
  fopen_s(&fp_msb, input, "rb");
  if (fp_msb == NULL) {
    fprintf(stderr, "Failed to open msb file %s\n", input);
    return;
  }
  fseek(fp_msb, 0L, SEEK_END);
  size_t bytes_to_read = ftell(fp_msb);
  rewind(fp_msb);

  uint8_t* buffer = (uint8_t*)malloc(bytes_to_read);
  size_t bytesRead;
  bytesRead = fread(buffer, 1, bytes_to_read, fp_msb);
  if (bytesRead != bytes_to_read) {
    // Process the read data
    printf("Error parsiing data\n");
    for (size_t i = 0; i < 50; i++) {
      printf("%02X ", buffer[i]);
    }
    printf("\n");
  }
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
    fprintf(stderr, "Failed to open output file\n");
    return;
  }
  
  //// Dial in the prn and doppler here ////////////////
#define SIZE 1024*4 *4 // 16K for Galileo and 4K for GPS
  int proc_gps = 0; // 1 for GPS, 0 for Galileo
  int prn = 5;// 15;// 24;// 4;// 23;// 5;// 31;// 6;// 9;// 34;// 36;// 24;// 34;// 4;// 15;// 36;// 9;// 23;// 9;// 6;// 26;// 9;// 29;// 11;// 6;// 7;// 31;// 3;// 26;// 4;// 9;// 7;// 11;// 6;// 31;// 9;// 3;// 26;// 4;// 10;// 23;// 13;// 27;// 18;// 15;// 10;//15
  double doppler = 2277.32;// 1557.36;// 3553.67;// -2011.07;// -854.19;// 2277.32;// 1927.59;// -1231.94;// -100.53;// -1013.74;// -2730.37;// 3553.67;// -1013.74;// -2011.07;// 1557.36;// -2730.37;// -100.53;// -854.19;// -1927.59;// -100.53;// -1231.94;// -1279.33;// 2651.04;// -82.61;// 646.60;// -1487.68;// 3213.83;// -3087.45;// -3052.57;// -1279.33;// 272.41;// 2651.04;// 3341.25;// 1312.40;// -903.17;// -2885.45;// 2967.03;// -2603.85;// -565.13;// 770.87;// 2023;// -530.0;// -3314.0;// 154.0;// -2034.0;// -2490;// 2030; // -2490
  
  /////////////////////////////////////////////////////

  c32* iandq = (c32*)malloc(SIZE * sizeof(c32));
  c32* repli = (c32*)malloc(SIZE * sizeof(c32));
  if (iandq == NULL || repli == NULL) {
    fprintf(stderr, "Memory allocation failed for q32 array.\n");
    free(iandq); free(repli);
    return;
  }

  DecodeOrsIQCplx(&buffer[hdrlen], SIZE/2, iandq);
  free(buffer);

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
    synth_gps_prn(prn, -doppler, SIZE, repli);
  }
  else { // Galileo
    synth_e1b_prn(prn, -doppler, SIZE, repli);
  }

  if (1) {
    arm_cfft_radix2_instance_f32 as;
    arm_cfft_radix2_instance_f32 rs;

    for (int i = 0; i < SIZE; i++) {
      actual[2 * i]      = iandq[i].r * 0.25;
      actual[2 * i + 1]  = iandq[i].i * 0.25;
      replica[2 * i]     = repli[i].r * 0.25;
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
      float Ar = actual[i * 2],  Ai = actual[i * 2 + 1];
      float Rr = replica[i * 2], Ri = replica[i * 2 + 1];
      // A * conj(R)
      prod[i * 2]     = Ar * Rr + Ai * Ri;     // (Ar + jAi) * (Rr - jRi)
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

  if (0) { // q31 
    arm_cfft_radix4_instance_q31 aq;
    arm_cfft_radix4_instance_q31 rq;
    q31_t  actual[SIZE * 2] = { 0 };
    q31_t replica[SIZE * 2] = { 0 };
    for (int i = 0; i < SIZE; i++) {
      actual[2 * i]      = f2q31(iandq[i].r * 0.00001);
      actual[2 * i + 1]  = f2q31(iandq[i].i * 0.00001);
      replica[2 * i]     = f2q31(repli[i].r * 0.00001);
      replica[2 * i + 1] = f2q31(repli[i].i * 0.00001);
    }
    if (iandq != NULL) { free(iandq); }
    if (repli != NULL) { free(repli); }

    arm_cfft_radix4_init_q31(&aq, SIZE, 0, 1); // Initialize the CFFT instance for 8-point FFT
    arm_cfft_radix4_q31(&aq, actual);

    arm_cfft_radix4_init_q31(&rq, SIZE, 0, 1); // Initialize the CFFT instance for 8-point FFT
    arm_cfft_radix4_q31(&rq, replica);

    q31_t prod[SIZE * 2] = { 0 };
    // Mult Actual with conjugate of Replica
    for (int i = 0; i < SIZE; i++) {
      float Ar = actual[i * 2], Ai = actual[i * 2 + 1];
      float Rr = replica[i * 2], Ri = replica[i * 2 + 1];
      // A * conj(R)
      prod[i * 2] = q31_add(q31_mul(Ar, Rr), q31_mul(Ai, Ri));     // (Ar + jAi) * (Rr - jRi)
      prod[i * 2 + 1] = q31_sub(q31_mul(Ai, Rr), q31_mul(Ar, Ri));     // 
    }

    arm_cfft_radix4_instance_q31 conv;
    arm_cfft_radix4_init_q31(&conv, SIZE, 1, 1); // inverse FFT
    arm_cfft_radix4_q31(&conv, prod);
    float max = 0;
    int pos = 0;
    for (int i = 0; i < SIZE; i++) {
      float mag = sqrt((float)prod[2 * i] * (float)prod[2 * i] + (float)prod[2 * i + 1] * (float)prod[2 * i + 1]);
      if (i <= 4 || i >= SIZE - 4) { mag = 0; }
      fprintf(fp_out, "%d, %f \n", i, mag);
      if (mag > max) { max = mag; pos = i; }
    }
    printf("max_q31 = %f pos=%d\n", max, pos);
  }

  if (0) { // q15 

    arm_cfft_radix4_instance_q15 aq;
    arm_cfft_radix4_instance_q15 rq;
    q15_t  actual[SIZE * 2] = { 0 };
    q15_t replica[SIZE * 2] = { 0 };
    for (int i = 0; i < SIZE; i++) {
      actual[2 * i] = f2q15(iandq[i].r * 0.2);
      actual[2 * i + 1] = f2q15(iandq[i].i * 0.2);
      replica[2 * i] = f2q15(repli[i].r * 0.2);
      replica[2 * i + 1] = f2q15(repli[i].i * 0.2);
    }
    if (iandq != NULL) { free(iandq); }
    if (repli != NULL) { free(repli); }

    arm_cfft_radix4_init_q15(&aq, SIZE, 0, 1); // Initialize the CFFT instance for 8-point FFT
    arm_cfft_radix4_q15(&aq, actual);

    arm_cfft_radix4_init_q15(&rq, SIZE, 0, 1); // Initialize the CFFT instance for 8-point FFT
    arm_cfft_radix4_q15(&rq, replica);

    q15_t prod[SIZE * 2] = { 0 };
    // Mult Actual with conjugate of Replica
    for (int i = 0; i < SIZE; i++) {
      float Ar = actual[i * 2], Ai = actual[i * 2 + 1];
      float Rr = replica[i * 2], Ri = replica[i * 2 + 1];
      // A * conj(R)
      prod[i * 2] = q15_add(q15_mul(Ar, Rr), q15_mul(Ai, Ri));     // (Ar + jAi) * (Rr - jRi)
      prod[i * 2 + 1] = q15_sub(q15_mul(Ai, Rr), q15_mul(Ar, Ri));     // 
    }

    arm_cfft_radix4_instance_q15 conv;
    arm_cfft_radix4_init_q15(&conv, SIZE, 1, 1); // inverse FFT
    arm_cfft_radix4_q15(&conv, prod);
    float max = 0;
    int pos = 0;
    for (int i = 0; i < SIZE; i++) {
      float mag = sqrt((float)prod[2 * i] * (float)prod[2 * i] + (float)prod[2 * i + 1] * (float)prod[2 * i + 1]);
      fprintf(fp_out, "%d, %f \n", i, mag);
      if (mag > max) { max = mag; pos = i; }
    }
    printf("max_q15 = %f pos=%d\n", max, pos);
  }

  fclose(fp_out);
  free(actual); free(replica); free(prod);
  //////////////////////////////////////////////////////
}

void test_quasi_pilot() {
  srand((unsigned int)time(NULL)); // randomise seed
  // Test the quasi pilot generation
  int len = 1023 * 2 * 240; // 2 samples per chip and 100 ms
  int min_idx = 0;
  int loc_cnt = 0;
  float min_val = 1e5;
  c32* out = (c32*)malloc(len * sizeof(c32));
  if (out == NULL) {
    fprintf(stderr, "Memory allocation failed for 100 ms I&Q array.\n");
    return;
  }
  int locations[50] = { -1 };// { 19, 39, 60, 79, 99, 119, 140, 159, 180 };// { -1 };// { 19, 39, 59, 79, 99, 119, 139, 159, 179 };// index into locations of bit transitions
  int window = 6; // 5 ms either side
  int nci = 240;
  int c_phase = 888; // which chip to set the code phase to
  int prn1 = 4, prn2 = 8;
  float dop1 = 2000, dop2 = -3000;
  int prn_c1[1023 * 2], prn_c2[1023 * 2];
  c32 replica[1023 * 2] = { 0 };
  getCode(1023, 2, prn1, prn_c1);
  getCode(1023, 2, prn2, prn_c2);
  make_replica(prn_c1, replica, dop1, 1023*2, 1.023e6 * 2);
  // now advance code-phase
  rotate_fwd(prn_c1, 1023 * 2, c_phase); // code phase 1/4 way
  int sign = 1;
  int sign2 = 1;
  stat_s stat;
  stat_init(&stat);
  for (int i = 0; i < nci; i++) {
    
    sign *= ((i +1) % 10) == 0 ? -1 : 1; // change sign every 20 ms
    for (int j = 0; j < 50; j++) {
      if (locations[j] == i) { sign2 *= -1; break; } // change sign at the bit transitions
    }
    // offset doppler by 250 Hz and add a residual doppler ramp of 0.1 Hz per ms
    mix_two_prns_oversampled_per_prn(prn_c1, prn_c2,dop1 + 250.0 + i * 0.1 ,dop2 - i * 0.1,0,0,
      &out[1023 * 2 * i],1023*2, 1.023e6 * 2 , 2.26, sign * sign2); // was 2.26 for -118 dBm
    printf("%d sign %d \n", i + 1, sign * sign2);
  }
  // use FFTs 
  float fft_replica[1024 * 2 * 2] = { 0 };
  arm_cfft_radix2_instance_f32 s;
  memset(fft_replica, 0, sizeof(fft_replica));
  arm_cfft_radix2_init_f32(&s, 1024 * 2, 0, 1);
  // xfer to float array
  for (int i = 0; i < 1023 * 2; i++) {
    fft_replica[i * 2 + 0] = replica[i].r * 0.25;
    fft_replica[i * 2 + 1] = replica[i].i * 0.25;
  }
  arm_cfft_radix2_f32(&s, fft_replica);
  
  float fft_data[1024 * 2 * 2] = { 0 };
  for (int center = window/2; center <= nci - window/2; center++) {
    // use linearity to sum the ncis coherently in moving window
    memset(fft_data, 0, sizeof(fft_data));
    for (int windex = center - window /2; windex < center + window/5; windex++) {
      for (int j = 0; j < 1023 * 2; j++) { // xfer to float array
        fft_data[j * 2 + 0] += out[(1023 * 2 * windex) + j].r * 0.25;
        fft_data[j * 2 + 1] += out[(1023 * 2 * windex) + j].i * 0.25;
      }
    } // for coherent integrations
    // need a better normalization strategy if using fixed point to avoid overflow
    for (int j = 0; j < 1023 * 2 * 2; j++) { fft_data[j] /= (float)window; }
    arm_cfft_radix2_init_f32(&s, 1024 * 2, 0, 1);
    arm_cfft_radix2_f32(&s, fft_data);
    // pt-wise multiply with conj of replica
    float fft_prod[1024 * 2 * 2] = { 0 };
    for (int k = 0; k < 1024 * 2; k++) {
      float Ar = fft_data[k * 2 + 0], Ai = fft_data[k * 2 + 1];
      float Rr = fft_replica[k * 2 + 0], Ri = fft_replica[k * 2 + 1]; // conj
      // A * conj(R) and add this product coherently
      fft_prod[k * 2 + 0] += (Ar * Rr + Ai * Ri);     // (Ar + jAi) * (Rr - jRi)
      fft_prod[k * 2 + 1] += (Ai * Rr - Ar * Ri);     // 
    }

    // inverse FFT
    arm_cfft_radix2_init_f32(&s, 1024 * 2, 1, 1);
    arm_cfft_radix2_f32(&s, fft_prod);

    float max = 0; int pos = 0;
    for (int m = 0; m < 1024 * 2; m++) {
      float mag = sqrt(fft_prod[m * 2] * fft_prod[m * 2] + fft_prod[m * 2 + 1] * fft_prod[m * 2 + 1]);
      //fprintf(fp_out, "%d, %f \n", m, mag);
      if (mag > max) { max = mag; pos = m; }
    }
    stat_add(&stat, max);
  
    float mean2 = stat_mean(&stat);
    if (max < min_val && max < mean2 - 0.5 && max < 20) {
      min_val = max;
      min_idx = center;
    }
    if ((center == min_idx + 1) && (max > min_val) ) {
      locations[loc_cnt++] = min_idx; // empirically the wider the window the earlier the bit transition appears
      min_val = 1e5;
      min_idx = 0;
    }
    printf("center=%d max=%f pos=%d mean=%f\n", center, max, pos, mean2);
  } // for center

  for (int i = 0; i < loc_cnt; i++) {
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
  const int multK = 16; // multiple of 1K size can be used for 1,2,4,8,or 16 (has not been tested for less than 1K)
  float F1 = 50.0f; // Frequency of the first sine wave
  float F2 = 0;// 80.0f; // Frequency of the second sine wave
  float F3 = 0;// 155.0f; // Frequency of the third sine wave (not used in this example, but can be added for more complexity)
  float F4 = 230.0f; // Frequency of the third sine wave (not used in this example, but can be added for more complexity)
  
// File lengths are too long to do float and radix4 and radix2 serially: use compile flags
//#define DO_FLOAT
#define DO_Q15_RADIX2
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
    read_ors((char*)"C:/work/Baseband/TestData/100ms/bw25/G_2025_09_03_23_04_45.ors");
    //read_ors((char*)"C:/work/Baseband/TestData/100ms/bw25/G_2025_09_03_23_22_41.ors");
    //read_ors((char*)"C:/work/Baseband/TestData/100ms/bw25/G_2025_09_03_23_04_45.ors");
    //read_ors((char*)"C:/work/Baseband/TestData/G_2025_06_05_22_11_26.ors");
    return 0;
  }

  if (1) {
    test_quasi_pilot();
    return 0;
  }

#ifdef DO_FLOAT
  float test[1024 * multK * 2];

  memset(test, 0, sizeof(test));
  for (int i = 0; i < 1024 * multK; i++) {
    double x = i * T;
    test[i*2] = 0.5 * signalf(x,F1,F2,F3,F4);
    test[i * 2 + 1] = 0; // set the imaginary part to 0
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


