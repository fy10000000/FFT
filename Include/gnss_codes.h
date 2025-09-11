/**
 * expose the codes in gal_codes.c
 *
 */

#pragma once
#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "ctype.h"
#include <stdint.h> // Include this for int8_t




#define E1B_CODE_LEN 4092
#define E1B_HEX_LEN  1023
#define E1B_MAX_PRN  50

#ifndef PI
#define PI               3.14159265358979f
#endif

typedef struct { float  r, i; } c32; // complex float
int8_t E1B_Code[E1B_MAX_PRN + 1][E1B_CODE_LEN]; // one based indexing

int load_e1b_primary_codes(char* path, int8_t out[E1B_MAX_PRN + 1][E1B_CODE_LEN]);



#ifdef __cplusplus
extern "C" {
#endif

void getCode(int num, int samplesPerChip, const int prn, int* out);

// this loads Gal prns under the covers
void synth_e1b_prn(
  int prn, // one based indexing
  float doppler,
  size_t N,
  c32* out
);

void synth_gps_prn(int prn, float doppler, size_t size, c32* replica);

void mix_two_prns_oversampled_per_prn(const int32_t* prn_a, const int32_t* prn_b,
  double doppler_a_hz, double doppler_b_hz,
  double phase_a_deg, double phase_b_deg,
  c32* out_iandq, int size, float samp_rate, float sigma, int sign);

void make_replica(const int32_t* prn_a, c32* out_iandq, float doppler, int size, float samp_freq);
void rotate_fwd(int array[], int size, int offset);

#ifdef __cplusplus
}
#endif