
#include "gnss_codes.h"

/* Convert one hex digit to value 0..15, or -1 on error. */
static int hex_nibble(int c) {
  if (c >= '0' && c <= '9') { return c - '0'; }
  if (c >= 'A' && c <= 'F') { return 10 + (c - 'A'); }
  if (c >= 'a' && c <= 'f') { return 10 + (c - 'a'); }
  return -1;
}

/* Map logic bit to signal chip per Table 11:
   logic 1 -> signal -1, logic 0 -> signal +1. */
static inline int8_t logic_to_chip_icd(int bit) {
  return bit ? (int8_t)-1 : (int8_t)+1;
}

static inline int boc_sign(int m, float u) {
  if (m <= 0) return +1; // BPSK
  int k = (int)floorf(2.0f * m * u); // toggles 2*m times per chip
  return (k & 1) ? -1 : +1;
}

static inline float cboc_e1b_weight(float u) {
  const float w1 = 0.9534625892455923f; // sqrt(10/11)
  const float w6 = 0.3015113445777636f; // sqrt(1/11)
  int s11 = boc_sign(1, u);
  int s61 = boc_sign(6, u);
  // E1-B uses minus sign for the 6,1 branch
  return w1 * (float)s11 - w6 * (float)s61;
}

// -------------------------------
// Code chip picker: returns +/-1 code value for PRN at fractional chip phase (wrapped).
// code_phase_chips in [0, L), L=4092
static inline int code_chip_at(const int8_t* code, float code_phase_chips, int L) {
  int idx = (int)floorf(code_phase_chips);
  if (idx >= L) idx -= L;
  else if (idx < 0) idx += L;
  return code[idx];
}


/* Convert an in-memory hex stream (any separators are allowed; only [0-9A-Fa-f] are consumed)
   into E1B_CODE_LEN chips in {+1, -1}, MSB-first within each nibble (Annex C.2).
   Returns 0 on success, <0 on error (e.g., insufficient hex). */
static int hexstream_to_e1b_chips(const char* hex_stream, int8_t chips[E1B_CODE_LEN]) {
  int chip_idx = 0;
  for (const char* p = hex_stream; *p && chip_idx < E1B_CODE_LEN; ++p) {
    int v = hex_nibble((unsigned char)*p);
    if (v < 0) continue; // skip non-hex
    // MSB-first within nibble: bits 3,2,1,0
    for (int b = 3; b >= 0 && chip_idx < E1B_CODE_LEN; --b) {
      int bit = (v >> b) & 1;
      chips[chip_idx++] = logic_to_chip_icd(bit);
    }
  }
  return (chip_idx == E1B_CODE_LEN) ? 0 : -1;
}

/* Trim leading/trailing whitespace in-place; returns pointer to first nonspace. */
static char* trim_inplace(char* s) {
  if (!s) { return s; }
  while (isspace((unsigned char)*s)) { ++s; }
  if (*s == 0) { return s; }
  char* e = s + strlen(s) - 1;
  while (e >= s && isspace((unsigned char)*e)) { *e-- = '\0'; }
  return s;
}

/* Try to parse an integer PRN at the start of line and locate start of hex field.
   Accepts separators: comma, semicolon, space, tab.
   On success: *out_prn set 1..E1B_MAX_PRN, returns pointer into the line where hex likely starts.
   On failure (no PRN at start): returns NULL. */
static const char* parse_prnnd_hex_start(const char* line, int* out_prn) {
  const char* p = line;
  // skip leading spaces
  while (isspace((unsigned char)*p)) { ++p; }
  // must start with digit to consider PRN
  if (!isdigit((unsigned char)*p)) { return NULL; }

  char* endptr = NULL;
  long prn = strtol(p, &endptr, 10);
  if (endptr == p || prn < 1 || prn > E1B_MAX_PRN) { return NULL; }

  // skip separators to reach hex
  const char* q = endptr;
  while (*q == ',' || *q == ';' || isspace((unsigned char)*q)) { ++q; }
  *out_prn = (int)prn;
  return q;
}

/* Load E1-B codes from a file.
   Formats supported:
     - CSV/TXT: "<PRN><sep><1023-hex>" per line (sep can be ',', ';', space, or tab)
     - HEX-only: "<1023-hex>" per line, PRNs assigned by line number (1..50)
   Returns number of PRNs successfully loaded (0..50), or <0 on file open/read error. */
int load_e1b_primary_codes(char* path, int8_t out[E1B_MAX_PRN + 1][E1B_CODE_LEN]) {
  FILE* f = fopen(path, "rb");
  if (!f) { return -1; }

  // zero out to mark "not loaded"
  for (int prn = 0; prn <= E1B_MAX_PRN; ++prn) {
    for (int i = 0; i < E1B_CODE_LEN; ++i) {
      out[prn][i] = 0;
    }
  }

  char line[16384];
  int loaded = 0;
  int next_prn_seq = 1;

  while (fgets(line, (int)sizeof(line), f)) {
    char* raw = trim_inplace(line);
    if (*raw == '\0') continue;           // empty
    if (*raw == '#')  continue;           // comment

    // If the line is suspiciously short in hex, skip early
    // (still allow separators; we’ll count hex later)
    int hex_chars = 0;
    for (char* c = raw; *c; ++c) {
      if (isxdigit((unsigned char)*c)) { ++hex_chars; }
    }
    if (hex_chars == 0) {
      continue;
    }

    int prn = 0;
    const char* hex_start = parse_prnnd_hex_start(raw, &prn);

    if (!hex_start) {
      // No explicit PRN; treat as hex-only line
      if (next_prn_seq > E1B_MAX_PRN) continue; // extra lines ignored
      prn = next_prn_seq++;
      hex_start = raw;
    }

    // Copy only hex digits into a temporary buffer to ensure exactly 1023 nibbles are considered.
    char hexbuf[E1B_HEX_LEN + 1];
    int h = 0;
    for (const char* p = hex_start; *p && h < E1B_HEX_LEN; ++p) {
      if (isxdigit((unsigned char)*p)) hexbuf[h++] = *p;
    }
    hexbuf[h] = '\0';

    if (h != E1B_HEX_LEN) {
      // Try to keep reading more hex from the remainder of the current line (already done) – not enough.
      // Some files may wrap the hex across multiple lines; this loader assumes one line per PRN.
      // You can extend here to accumulate across lines if needed.
      // Skip malformed line
      continue;
    }

    if (hexstream_to_e1b_chips(hexbuf, out[prn]) == 0) {
      loaded++;
    }
  }

  fclose(f);
  return loaded;
}

void sincosf_fast(float phia, float* sine, float* cosine) {
  *sine = sin(phia);
  *cosine = cos(phia);
}

// -------------------------------
// Signal synthesis for two PRNs (PRN2 and PRN4) into one complex baseband stream.
// Inputs:

//   prn              : PRN number
//   doppler          : carrier Doppler (Hz) per PRN
//   phi   _rad       : initial carrier phase (rad) at sample 0
//   code_phase       : initial code phase (chips) in [0,4092)
//   fs_hz            : sample rate
//   code_rate_cps    : nominal code chip rate (1.023e6 for E1), set per-PRN Doppler via scaling if desired
//   N                : number of samples
//   out              : output buffer of c32 length N
// Notes:
//   - E1-B data bit is assumed +1 for this short synthesis. If needed, pass a per-sample sign array and multiply.
extern void synth_e1b_prn(
  int prn, // one based indexing
  float doppler, 
  float phi_rad, 
  float code_phase,
  float fs_hz, // typically 4.092e6f
  float code_rate_cps, // typically 1.023e6f
  size_t N,
  c32* out
) {
  int ret = load_e1b_primary_codes((char*)"C:/work/Baseband/HEX_E1B.txt", E1B_Code);
  if (ret < 0) { printf("Error loading Galileo codes; check path in synth_e1b_prn()\n"); return; }
  const int L = 4092;
  const int8_t* ca = E1B_Code[prn];

  float dchips = (code_rate_cps) / fs_hz;
  float dphia = 2.0f * (float) PI * (doppler) / fs_hz; // phase increment per sample
  float chips = code_phase;
  float phia = phi_rad;

  for (size_t n = 0; n < N; ++n) {
    // PRN a
    float frac = chips - floorf(chips);
    int code = code_chip_at(ca, chips, L);
    float sca = cboc_e1b_weight(frac);
    float ampa = (float)code * sca; // data assumed +1

    float sa, caph; 
    sincosf_fast(phia, &sa, &caph);
    c32 xa = { ampa * caph, ampa * sa };
 
    out[n].r = xa.r; 
    out[n].i = xa.i; 

    // advance
    chips += dchips;
    if (chips >= L) { chips -= L; }
    else if (chips < 0) { chips += L; }

    phia += dphia;
    if (phia > 1e9f || phia < -1e9f) { 
      phia = fmodf(phia, 2.0f * PI);
    }
    
  }
}

////////////////////////////////////////////////////////////////

void getCode(int num, int samplesPerChip, const int prn, int* out)
{
  //Feedback taps as defined in GPS spec
  int g1tap[] = { 2,9 };
  int g2tap[] = { 1,2,5,7,8,9 };
  int g1tap_len = 2;
  int g2tap_len = 6;
  int sats[37][2] = { {1, 5}, {2, 6}, {3, 7}, {4, 8}, {0, 8}, {1, 9}, {0, 7}, {1, 8}, {2, 9}, {1, 2},
                      {2, 3}, {4, 5}, {5, 6}, {6, 7}, {7, 8}, {8, 9}, {0, 3}, {1, 4}, {2, 5}, {3, 6},
                      {4, 7}, {5, 8}, {0, 2}, {3, 5}, {4, 6}, {5, 7}, {6, 8}, {7, 9}, {0, 5}, {1, 6},
                      {2, 7}, {3, 8}, {4, 9}, {3, 9}, {0, 6}, {1, 7}, {3, 9} };

  int g[1023 * 10]; // max C/A code length (1023) * max oversample (10)
  int g1[10], g2[10];
  int i, j, k;
  for (i = 0; i < 10; ++i) {
    g1[i] = 1;
    g2[i] = 1;
  }

  for (i = 0; i < num; ++i) {
    int val = (g1[9] + g2[sats[prn - 1][0]] + g2[sats[prn - 1][1]]) % 2;
    g[i] = val;

    // Shift g1
    int g1_fb = 0;
    for (j = 0; j < g1tap_len; ++j) {
      g1_fb += g1[g1tap[j]];
    }
    g1_fb %= 2;
    for (j = 9; j > 0; --j) {
      g1[j] = g1[j - 1];
    }
    g1[0] = g1_fb;

    int g2_fb = 0;
    for (j = 0; j < g2tap_len; ++j) {
      g2_fb += g2[g2tap[j]];
    }
    g2_fb %= 2;
    for (j = 9; j > 0; --j) {
      g2[j] = g2[j - 1];
    }
    g2[0] = g2_fb;

  }

  for (i = 0; i < num; ++i) {
    if (g[i] == 0) { g[i] = -1; }
  }

  if (samplesPerChip > 1) {
    for (i = 0, k = 0; i < num; ++i) {
      for (j = 0; j < samplesPerChip; ++j) {
        out[k++] = g[i];
      }
    }
  }
  else {
    memcpy(out, g, num * sizeof(int));
  }
}

void mix_prn(const int32_t prn_a[1023 * 4],
  const double doppler_a_hz, const double phase_a_deg, c32* signal, int size)
{
  double FS = 1.023e6 * 4;
  // PRN A phasor increment and initial phase
  const double dth_a = 2.0 * PI * (double)doppler_a_hz / (double)FS;
  const double ca_inc = (double)cos(dth_a);
  const double sa_inc = (double)sin(dth_a);

  double pca = cos(phase_a_deg * PI / 180.0);
  double psa = sin(phase_a_deg * PI / 180.0);

  int n = 0;
  for (int samp = 0; samp < size; ++samp) {
    double a = (double)prn_a[samp];  // +/- 1

    // s[n] = a*e^{j theta_a[n]} 
    double ia = a * pca, qa = a * psa;

    signal[samp].r = (float)(ia);
    signal[samp].i = (float)(qa);

    // advance both phasors
    double npca = pca * ca_inc - psa * sa_inc;
    double npsa = pca * sa_inc + psa * ca_inc;
    pca = npca; psa = npsa;

    // renormalize occasionally to limit float drift (not necessary for 1 ms)
    double na = 1.0 / sqrt(pca * pca + psa * psa);
    pca *= na; psa *= na;
  }
}


void synth_gps_prn(int prn, float doppler, size_t size, c32*  replica) {
  int* code = (int*)malloc(sizeof(int) * size);
  getCode(size/4, 4, prn, code);
  mix_prn(code,doppler,0, replica, size);
  free(code);
}

