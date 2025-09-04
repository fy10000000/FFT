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

// this loads Gal prns under the covers
void synth_e1b_prn(
  int prn, // one based indexing
  float doppler,
  float phi_rad,
  float code_phase,
  float fs_hz, // typically 4.092e6f
  float code_rate_cps, // typically 1.023e6f
  size_t N,
  c32* out
);

void synth_gps_prn(int prn, float doppler, size_t size, c32* replica);

#ifdef __cplusplus
}
#endif