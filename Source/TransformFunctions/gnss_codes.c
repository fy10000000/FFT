
#include "gnss_codes.h"
#include "up_sample.h"
#include <stdbool.h>
#include <complex.h>

//#include "transform_functions.h"
//void fft_simple2(int size, complex float* w, bool fwd) {
//  // interleaved real, imagmalloc(size * sizeof(float);
//  float* fft_data = (float*)malloc(sizeof(float) * 2 * size);
//  for (int i = 0; i < size; i++) {
//    fft_data[2 * i] = crealf(w[i]);
//    fft_data[2 * i + 1] = cimagf(w[i]);
//  }
//  arm_cfft_radix2_instance_f32 s;
//  arm_cfft_radix2_init_f32(&s, size, fwd ? 0 : 1, 1);
//  arm_cfft_radix2_f32(&s, fft_data);
//  for (int i = 0; i < size; i++) {
//    w[i] = fft_data[2 * i] + I * fft_data[2 * i + 1];
//  }
//  free(fft_data);
//}




//  GPS L5 I-channel spreading code generator.
//
//  References:
//    - IS-GPS-705 (L5), publicly available ICD.
//    - Chip rate: 10.23 MHz, code length: 10,230 chips (1 ms).
//    - Two 13-stage LFSRs, G1 and G2, with all-ones initialization.
//    - Output chip c[n] = g1[n] XOR g2[s1[n]] XOR g2[s2[n]] where s1,s2 are PRN-specific taps.
//      Here we implement as XOR of two selected G2 tap outputs (Gold-like with phase selectors).
//    - This implementation returns chips as +1/-1 (bipolar) or 0/1 depending on flag.
//
//  Conventions:
//    - LFSR state bits: bit 0 = stage 1 (rightmost), bit 12 = stage 13 (leftmost).
//    - At each step, output is the bit of stage 13 (leftmost) OR as per ICD definition; the
//      XOR tap positions below are defined accordingly so that generated sequence matches ICD.
//    - Initialization: all stages = 1.
//
//  WARNING:
//    Different sources label stages left-to-right differently. The taps below are defined for
//    the shift/update implemented here. If you use a different shifting convention, adjust.



// Choose output mapping: bipolar (+1/-1) or binary (0/1) 
//typedef enum { L5_OUT_BINARY = 0, L5_OUT_BIPOLAR = 1 } l5_out_t;

// LFSR polynomials (tap masks) for 13-stage registers.
//   We implement right-shift with new bit injected at bit 12 (MSB), output taken from bit 0 (LSB).
//   Feedback bit = XOR of tapped bits of current state.
//   Taps chosen to realize: c.f. p 14 IS-GPS-705J.pdf
//   XA: 1+ x9 + x10 + x12 + x13
//   XBIi or XBQi: 1 + x + x3 + x4 + x6 + x7 + x8 + x12 + x13
//     G1: x^13 + x^12 + x^10 + x^9 + 1
//     G2: x^13 + x^12 + x^8 + x^7 + x^6 + x^4 + x^3 + x^1 + 1
//
//   With LSB as stage 1, MSB as stage 13, the feedback mask marks stages that are XORed.
//   Mask bit i corresponds to stage (i+1).
static inline uint16_t g1_feedback(uint16_t s) {
  // taps at stages: 13,12,10,9,1 => bits 12,11,9,8,0 
  unsigned b = ((s >> 12) ^ (s >> 11) ^ (s >> 9) ^ (s >> 8) ^ (s >> 0)) & 1u;
  return (uint16_t)b;
}
static inline uint16_t g2_feedback(uint16_t s) {
  // taps at stages: 13,12,8,7,6,4,3,1 => bits 12,11,7,6,5,3,2,0 
  unsigned b = ((s >> 12) ^ (s >> 11) ^ (s >> 7) ^ (s >> 6) ^ (s >> 5) ^ (s >> 3) ^ (s >> 2) ^ (s >> 0)) & 1u;
  return (uint16_t)b;
}

// One LFSR step: right shift, inject new bit at MSB (bit 12), output old LSB (bit 0) 
static inline unsigned lfsr_step(uint16_t* state, unsigned (*fb)(uint16_t)) {
  unsigned out = *state & 1u;                 // output = stage 1 (LSB) 
  unsigned f = fb(*state);
  *state = (uint16_t)((*state >> 1) | (f << 12));
  return out;
}

// PRN-specific G2 tap selector pairs for L5 I-channel.
//   Each PRN uses XOR of two delayed taps from G2. Values are 1..13 (stage numbers).
//   Note: Below is an example table for a subset of PRNs. Fill as per ICD annex for full set.
//   For demonstration, we include PRNs 1..10. Extend to all satellites as needed.
//
typedef struct { uint8_t tap1; uint8_t tap2; } L5TapPair;

// Example subset; replace/extend with full ICD mapping for operational use 
static const L5TapPair L5I_TAPS[] = {
  // index 0 unused to make PRN index 1-based 
  {0,0}, // 1-37
  {1, 5}, {2, 6}, {3, 7}, {4, 8}, {0, 8}, {1, 9}, {0, 7}, {1, 8}, {2, 9}, {1, 2},
  {2, 3}, {4, 5}, {5, 6}, {6, 7}, {7, 8}, {8, 9}, {0, 3}, {1, 4}, {2, 5}, {3, 6},
  {4, 7}, {5, 8}, {0, 2}, {3, 5}, {4, 6}, {5, 7}, {6, 8}, {7, 9}, {0, 5}, {1, 6},
  {2, 7}, {3, 8}, {4, 9}, {3, 9}, {0, 6}, {1, 7}, {3, 9},
  // ... fill out to PRN as needed 
};

// Number of PRNs implemented in table (max index) 
#define L5I_MAX_PRN ((int)(sizeof(L5I_TAPS)/sizeof(L5I_TAPS[0])) - 1)

// Extract stage k (1..13) bit from LFSR state with LSB as stage 1 
static inline unsigned stage_bit(uint16_t s, unsigned k) {
  return (s >> (k - 1u)) & 1u;
}

// Generate N chips of L5 I-channel code for a given PRN.
//   - prn: 1..L5I_MAX_PRN (extend table for full constellation)
//   - out: destination buffer of length N
//   - start_at_epoch: if true, start from code phase 0 (all ones initial state)
//   - init_state optional: if start_at_epoch is false and init_state provided, use it to resume.
//
//   Returns false if PRN unsupported.
typedef struct {
  uint16_t g1; // 13-bit state, LSB=stage1, initialize to 0x1FFF (13 ones) 
  uint16_t g2;
  uint32_t chip_index; // modulo 10230 if desired 
} L5State;

static void l5_reset(L5State* st) {
  st->g1 = 0x1FFFu; // 13 ones 
  st->g2 = 0x1FFFu;
  st->chip_index = 0;
}

bool l5_generate_I(int prn, int8_t* out, size_t N, bool start_at_epoch, L5State* state_io)
{
  if (prn <= 0 || prn > L5I_MAX_PRN) return false;

  L5State st_local;
  L5State* st = &st_local;

  if (!start_at_epoch && state_io && state_io->g1 && state_io->g2) {
    *st = *state_io;
  }
  else {
    l5_reset(st);
  }

  const L5TapPair tp = L5I_TAPS[prn];

  for (size_t i = 0; i < N; ++i) {
    // Output taps from current states BEFORE stepping 
    unsigned g1_out = stage_bit(st->g1, 13); // Some conventions use stage 13; adjust if needed 
    // For robustness, use the same output definition for both: 
    // Alternatively use LSB as output; then adjust taps accordingly. 

    // G2 selector outputs: tap1 XOR tap2 (stage indices 1..13) 
    unsigned g2_t1 = stage_bit(st->g2, tp.tap1);
    unsigned g2_t2 = stage_bit(st->g2, tp.tap2);
    unsigned g2_sel = g2_t1 ^ g2_t2;

    unsigned chip = g1_out ^ g2_sel;

    if (1) {
      out[i] = chip ? -1 : +1; // map 0->+1, 1->-1 
    }
    else {
      out[i] = (int8_t)chip;
    }

    // Advance both LFSRs one chip 
    (void)lfsr_step(&st->g1, g1_feedback);
    (void)lfsr_step(&st->g2, g2_feedback);

    st->chip_index++;
    if (st->chip_index == L5_CODE_LEN) {
      st->chip_index = 0;
      // Optionally, reinitialize at each ms if desired: l5_reset(st); 
    }
  }

  if (state_io) {
    *state_io = *st;
  }
  return true;
}

// Helper to generate one full 1 ms epoch for a PRN (10,230 chips) 
bool l5_generate_I_epoch(int prn, int8_t* out)
{
  L5State st; l5_reset(&st);
  return l5_generate_I(prn, out, L5_CODE_LEN, true, &st);
}

/* Demo main: generate first 64 chips of PRN 1, print as +1/-1 */
#ifdef L5_DEMO
int main(void) {
  const int PRN = 1;
  int8_t chips[64];

  //if (!l5_generate_I(PRN, chips, 64, L5_OUT_BIPOLAR, true, NULL)) {
  if (!l5_generate_I(PRN, chips, 64, true, NULL)) {
    fprintf(stderr, "Unsupported PRN\n");
    return 1;
  }
  for (int i = 0; i < 64; ++i) {
    printf("%d%s", (int)chips[i], (i % 32 == 31) ? "\n" : " ");
  }
  return 0;
}
#endif

///--------------------------------------------------------

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

/* Unpack MSB-first packed bits into +/-1 chips.
   Input:
     packed: pointer to packed bits, MSB-first in each byte (bit 7 is first)
     nbits:  number of valid bits; must be >= 10230
   Output:
     code: length 10230, values +1 or -1 (int8_t)
   Mapping: bit 0 -> +1, bit 1 -> -1 (change if your convention differs)
*/
static int hexstream_to_e5a_chips(const char* hex_stream, int8_t chips[E5A_CODE_LEN]) {
  int chip_idx = 0;
  for (const char* p = hex_stream; *p && chip_idx < E5A_CODE_LEN; ++p) {
    int v = hex_nibble((unsigned char)*p);
    if (v < 0) continue; // skip non-hex
    // MSB-first within nibble: bits 3,2,1,0
    for (int b = 3; b >= 0 && chip_idx < E5A_CODE_LEN; --b) {
      int bit = (v >> b) & 1;
      chips[chip_idx++] = logic_to_chip_icd(bit);
    }
  }
  return (chip_idx == E5A_CODE_LEN) ? 0 : -1;
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

/* Load E5-A codes from a file.
   Formats supported:
     - CSV/TXT: "<PRN><sep><1023-hex>" per line (sep can be ',', ';', space, or tab)
     - HEX-only: "<1023-hex>" per line, PRNs assigned by line number (1..50)
   Returns number of PRNs successfully loaded (0..50), or <0 on file open/read error. */
int load_e5a_primary_codes(char* path, int8_t out[E5A_CODE_LEN], int target_prn) {
  FILE* f = fopen(path, "rb");
  if (!f) { return -1; }  

  char line[E5A_HEX_LEN + 512];
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

    if (prn != target_prn) {
      continue;
    }

    if (!hex_start) {
      // No explicit PRN; treat as hex-only line
      if (next_prn_seq > E5A_MAX_PRN) continue; // extra lines ignored
      prn = next_prn_seq++;
      hex_start = raw;
    }

    // Copy only hex digits into a temporary buffer to ensure exactly 1023 nibbles are considered.
    char hexbuf[E5A_HEX_LEN + 1];
    int h = 0;
    for (const char* p = hex_start; *p && h < E5A_HEX_LEN; ++p) {
      if (isxdigit((unsigned char)*p)) hexbuf[h++] = *p;
    }
    hexbuf[h] = '\0';

    if (h != E5A_HEX_LEN) {
      // Try to keep reading more hex from the remainder of the current line (already done) – not enough.
      // Some files may wrap the hex across multiple lines; this loader assumes one line per PRN.
      // You can extend here to accumulate across lines if needed.
      // Skip malformed line
      continue;
    }

    if (hexstream_to_e5a_chips(hexbuf, out) == 0) {
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

extern void getGalCode(int prn, int* out, int size) {
  for (int i = 0; i < size; i++) {
    out[i] = E1B_Code[prn][i];
  }
}

static size_t gcd_size_t(size_t a, size_t b) {
  while (b) {
    size_t t = a % b;
    a = b;
    b = t;
  }
  return a;
}

// Rotate array a[0..len-1] right by N positions (forward).
void rotate_right_i8_cycles(int8_t* a, size_t len, long long N) {
  if (!a || len == 0) return;

  long long kll = N % (long long)len;
  if (kll < 0) kll += (long long)len; // handle negative
  size_t k = (size_t)kll;
  if (k == 0) return;

  size_t g = gcd_size_t(len, k);
  for (size_t start = 0; start < g; ++start) {
    int8_t tmp = a[start];
    size_t i = start;

    while (1) {
      size_t next = i + k;
      if (next >= len) next -= len;
      if (next == start) break;

      a[i] = a[next];
      i = next;
    }
    a[i] = tmp;
  }
}

//   prn              : PRN number
//   doppler          : carrier Doppler (Hz) per PRN
//   phi   _rad       : initial carrier phase (rad) at sample 0
//   code_phase       : initial code phase (chips) in [0,4092)
//   fs_hz            : sample rate
//   code_rate_cps    : nominal code chip rate (1.023e6 for E1), set per-PRN Doppler via scaling if desired
//   N                : number of samples
//   out              : output buffer of c32 length N
extern void synth_e5a_prn(
  int prn, // one based indexing
  float doppler,
  size_t N,
  c32* out,
  int rotate_offset
  ) {
  float phi_rad = 0.0;
  float code_phase = 0.0;
  const float fs_hz = 10.23e6;// Mspc
  const float code_rate_cps = 10.23e6f;

  int8_t code_loc[E5A_CODE_LEN] = { 0 };
  int ret = load_e5a_primary_codes((char*)"C:/work/Baseband/HEX_E5AI.txt", code_loc, prn);
  if (ret < 0) { printf("Error loading Galileo codes; check path in synth_e5b_prn()\n"); return; }

  if (rotate_offset) {
    //rotate_right_i8_cycles(code_loc, E5A_CODE_LEN, rotate_offset); //still a bug in this
    int size = E5A_CODE_LEN;
    int8_t* temp = (int8_t*)malloc(sizeof(int8_t) * size);
    if (temp == NULL) { printf("rotate_fwd's malloc failed");  return; }
    for (int i = 0; i < size; i++) {
      if (i - rotate_offset < 0) { temp[i] = code_loc[size + (i - rotate_offset)]; }
      else { temp[i] = code_loc[i - rotate_offset]; }
    }
    memcpy(code_loc, temp, sizeof(int8_t) * size);
    free(temp);
  }

  //FILE* fp_out = NULL; //output file
  //errno_t er = fopen_s(&fp_out, "C:/Python/prn36.csv", "w");
  //for (int i = 0; i < E5A_CODE_LEN; i++) {
  //  fprintf(fp_out, "%d, %d\n",i, -code_loc[i]);
  //}
  //fclose(fp_out);

  
  const int L = E5A_CODE_LEN;
  const int8_t* ca = code_loc;

  float dchips = (code_rate_cps) / fs_hz;
  float dphia = 2.0f * (float)PI * (doppler) / fs_hz; // phase increment per sample
  float chips = code_phase;
  float phia = phi_rad;
  c32 in[E5A_CODE_LEN];
  int8_t no_quant = 1; // 0 for quantization
  for (size_t n = 0; n < L; ++n) {
    // PRN a
    float frac = chips - floorf(chips);
    int code = code_chip_at(ca, chips, L);
    float sca = 1.0;// cboc_e1b_weight(frac);
    float ampa = (float)code * sca; // data assumed +1

    float sa, ca;
    sincosf_fast(phia, &sa, &ca);
    c32 xa = { ampa * ca, ampa * sa };

    out[n].r = no_quant ? xa.r : quantize_pm1(xa.r + noise(1.0));
    out[n].i = no_quant ? xa.i : quantize_pm1(xa.i + noise(1.0));

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

void readL5Icode(int prn, int8_t* out) {
  FILE* f = NULL;
  errno_t er = fopen_s(&f, "C:/work/Baseband/codes_L5I.csv", "r");
  if (er != 0 || f == NULL) {
    printf("Error opening L5I code file\n");
    return;
  }
  char line[128];
  // prn N is in the N-1 th column
  int line_num = 0; 
  char* token = NULL;
  char* context = NULL;
  while (fgets(line, sizeof(line), f)) {
    for (int i = 0; i < 32; i++) {
      token = (i == 0) ? strtok_s(line, ",", &context) : strtok_s(NULL,",", &context);
      if (i == (prn - 1)) {
        out[line_num] = atoi(token);
        break;
      }
    }
    line_num++;
  }
  if (line_num != L5_CODE_LEN) {
    printf("Error reading L5I code for PRN %d\n", prn);
  }
  fclose(f);
}

extern void synth_L5I_prn(
  int prn, // one based indexing
  float doppler,
  size_t N,
  c32* out,
  int rotate_offset
) {
  float phi_rad = 0.0;
  float code_phase = 0.0;
  const float fs_hz = 10.23e6;// Mspc
  const float code_rate_cps = 10.23e6f;

  int8_t code_loc[L5_CODE_LEN] = { 0 };

  readL5Icode(prn, code_loc);
  //if (!l5_generate_I(prn, code_loc, L5_CODE_LEN, true, NULL)) {
  //  fprintf(stderr, "Unsupported PRN\n");
  //  return 1;
  //}
 

  if (rotate_offset) {
    int size = L5_CODE_LEN;
    int8_t* temp = (int8_t*)malloc(sizeof(int8_t) * size);
    if (temp == NULL) { printf("rotate_fwd's malloc failed");  return; }
    for (int i = 0; i < size; i++) {
      if (i - rotate_offset < 0) { temp[i] = code_loc[size + (i - rotate_offset)]; }
      else { temp[i] = code_loc[i - rotate_offset]; }
    }
    memcpy(code_loc, temp, sizeof(int8_t) * size);
    free(temp);
  }

 
  const int L = L5_CODE_LEN;
  const int8_t* ca = code_loc;

  float dchips = (code_rate_cps) / fs_hz;
  float dphia = 2.0f * (float)PI * (doppler) / fs_hz; // phase increment per sample
  float chips = code_phase;
  float phia = phi_rad;
  c32 in[L5_CODE_LEN];
  int8_t no_quant = 1; // 0 for quantization
  for (size_t n = 0; n < L; ++n) {
    // PRN a
    float frac = chips - floorf(chips);
    int code = code_chip_at(ca, chips, L);
    float sca = 1.0;
    float ampa = (float)code * sca; // data assumed +1

    float sa, ca;
    sincosf_fast(phia, &sa, &ca);
    c32 xa = { ampa * ca, ampa * sa };

    out[n].r = no_quant ? xa.r : quantize_pm1(xa.r + noise(1.0));
    out[n].i = no_quant ? xa.i : quantize_pm1(xa.i + noise(1.0));

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

void up_sample_10k_to_16k(c32* in , c32* out) {
  // wrapper around the function doing this
  const double Fin = 10230e3;
  const double Fout = 16384e3;
  const int    NTAPS = 12;      // try 8, 12, or 16
  const double CUTOFF = 0.45;

  iq_resamp_zp_t* rs = iq_resamp_zp_init(Fin, Fout, NTAPS, CUTOFF);
  if (!rs) { fprintf(stderr, "init failed\n"); return ; }

  size_t Nin = 10230; // 1 ms at 10.23 Msps
  //c32* in = (c32*)calloc(Nin, sizeof(c32));
  // fill 'in' ...

  size_t Nout_cap = (size_t)llround(Nin * (Fout / Fin)) + 8;
  //c32* out = (c32*)malloc(Nout_cap * sizeof(c32));

  size_t Nout = iq_resamp_process_zero_phase(rs, in, Nin, out, Nout_cap);
  printf("Produced %zu samples (expected ~16384)\n", Nout);

  iq_resamp_zp_free(rs);
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
  size_t N,
  c32* out
) {
  float phi_rad = 0.0;
  float code_phase = 0.0;
  const float fs_hz = 4.092e6f;
  const float code_rate_cps = 1.023e6f;
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

extern void getCode(int num, int samplesPerChip, const int prn, int* out)
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

void mix_prn(const int32_t* prn_a,
  const double doppler_a_hz, const double phase_a_deg, c32* signal, int size, int spc)
{
  double FS = 1.023e6 * spc;
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


void synth_gps_prn(int prn, float doppler, size_t size, c32*  replica, int spc) {
  int* code = (int*)malloc(sizeof(int) * size);
  getCode(size/ spc, spc, prn, code);
  mix_prn(code,doppler,0, replica, size, spc);
  free(code);
}

#define Q13_THRESHOLD 0.568 //0.5
int8_t quantize_pm13(double x) {
  double thr = Q13_THRESHOLD;
  if (thr < 0.0) thr = 0.0;
  if (thr > 1.0) thr = 1.0;
  int8_t mag = (fabs(x) >= thr) ? 3 : 1;
  return (x < 0.0) ? (int8_t)(-mag) : mag;
}

int8_t quantize_pm1(double x) {
  // Symmetric tie-break: zero maps to +1
  return (x >= 0.0f) ? +1 : -1;
}

double noise(double sigma) {
  double u1 = 0.0;
  do { u1 = (double)rand() / RAND_MAX; } while (u1 == 0.0); // avoid log(0)
  double u2 = (double)rand() / RAND_MAX;
  double ans = sigma * sqrt(-2.0 * log(u1)) * cos(2.0 * PI * u2);
  //printf("%f ", ans);
  return ans;
}