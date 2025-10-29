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

#include "cplx_types.h"


#define E1B_CODE_LEN 4092
#define E5A_CODE_LEN 10230

#define E1B_HEX_LEN  1023
#define E5A_HEX_LEN  2558

#define E1B_MAX_PRN  50
#define E5A_MAX_PRN  50

#ifndef PI
#define PI               3.14159265358979f
#endif

int8_t E1B_Code[E1B_MAX_PRN + 1][E1B_CODE_LEN]; // one based indexing

//int8_t E5AICodes[E5A_MAX_PRN][E5A_CODE_LEN];// zero based indexing


int load_e1b_primary_codes(char* path, int8_t out[E1B_MAX_PRN + 1][E1B_CODE_LEN]);

int load_e5a_primary_codes(char* path, uint8_t out[E5A_CODE_LEN], int prn);

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

void synth_e5a_prn(
  int prn, // one based indexing
  float doppler,
  size_t N,
  c32* out,
  int rotate_offset
);

void up_sample_10k_to_16k(c32* in, c32* out);

void getGalCode(int prn, int* out, int size);

void synth_gps_prn(int prn, float doppler, size_t size, c32* replica, int spc);

void mix_two_prns_oversampled_per_prn(const int32_t* prn_a, const int32_t* prn_b,
  double doppler_a_hz, double doppler_b_hz,
  double phase_a_deg, double phase_b_deg,
  c32* out_iandq, int size, float samp_rate, float sigma, int sign);

void make_replica(const int32_t* prn_a, c32* out_iandq, float doppler, int size, float samp_freq);
void rotate_fwd(int array[], int size, int offset);

int8_t quantize_pm13(double x);
int8_t quantize_pm1(double x);

double noise(double sigma);

#ifdef __cplusplus
}
#endif