/**
 * read write functions for the msb binary data
 *
 */

#pragma once
#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "ctype.h"
#include <stdint.h> // Include this for int8_t
#include <stdbool.h>


#define SYS_GPS     0x01                /* navigation system: GPS */
#define SYS_GAL     0x02                /* navigation system: Galileo */

typedef struct {
  uint8_t* buff;
  int      bit_pos;
} buff_str;


typedef struct {
  double time_tag;   /* time (s) expressed by standard time_t */
  float pseudorange; /* pseudorange (m) */
  float code_phase;  /* code phase (ms) */
  float doppler;     /* doppler frequency (Hz) */
  float carr_phase;  /* cycles at L1 90.5,0.5] */
  float cno;         /* carrier-to-noise density ratio (dB-Hz) */
  uint16_t prn;      /* satellite PRN number */
  uint16_t constellation; /* constellation type (GPS, GLONASS, Galileo, etc.) */
} bb_meas_sat_t;

typedef struct {
  uint64_t time;          /* time (s) expressed by standard time_t */
  double sec;             /* fraction of second under 1 s */
  int num_sat;            /* number of satellites collected */
  bb_meas_sat_t sats[20]; /* satellite measurements */
} bb_meas_t;

#ifdef __cplusplus
extern "C" {
#endif

  extern int write_bb_msb(const bb_meas_t* measurements, uint8_t* bin_buff, int buff_size);

  extern void read_bb_msb(uint8_t* bin_buff, int bin_size, bb_meas_t* measurements);

#ifdef __cplusplus
}
#endif