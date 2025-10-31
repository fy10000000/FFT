#pragma once

#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <string.h>

#include "cplx_types.h"


typedef struct {
  double fin;      // input sample rate (Hz)
  double fout;     // output sample rate (Hz)
  double dt;       // input samples per output sample = Fin/Fout
  int    ntaps;    // total taps (even), e.g., 8, 12, or 16
  int    half;     // ntaps/2
  double cutoff;   // normalized to input Nyquist (0 < cutoff <= 0.5), e.g., 0.45
} iq_resamp_zp_t;

iq_resamp_zp_t* iq_resamp_zp_init(double fin_hz, double fout_hz, int ntaps, double cutoff);

size_t iq_resamp_process_zero_phase(iq_resamp_zp_t* st, const c32* in, size_t Nin, c32* out, size_t Nout_cap);

void iq_resamp_zp_free(iq_resamp_zp_t* st);