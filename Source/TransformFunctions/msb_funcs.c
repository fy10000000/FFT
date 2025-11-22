
#include <msb_funcs.h>

int cur_bit = 0; //not mod 8 but the bit length of msg
int buff_size = 0; // size of the buffer in bytes

void reset_bit_reader(int buffer_size) {
  cur_bit = 0;
  buff_size = buffer_size;
}


// primitive to read the actual bits
int64_t read_bit(uint8_t* buff) {
  if (cur_bit >= buff_size * 8) {
    fprintf(stderr, "Bit reader out of bounds\n");
    return -1; // Error
  }

  int byte_index_loc = cur_bit / 8;
  int bit_index_loc = cur_bit % 8;
  int64_t current_byte = buff[byte_index_loc];
  int bit = (current_byte >> (7 - bit_index_loc)) & 0x01;
  cur_bit++;
  return bit;
}

int64_t bit_reader(uint8_t* buff, unsigned num_bit, int is_signed) {

  int64_t ans = 0;
  for (unsigned i = 0; i < num_bit; i++) {
    ans = ans | (read_bit(buff) << (num_bit - (i + 1))); // read from LSB bit
  }

  // test num_bit is high for negative
  int64_t test = (0x1LL << (num_bit - 1)) & ans;
  if (test && is_signed) {
    //printf("is negative\n");
    int tmp = 0xFFFFFFFFFFFFFFFF << num_bit;
    ans = ans | tmp;
  }

  return ans;
}


void bit_writer(buff_str* bs, uint64_t val, int num_bits)
{
  if (!bs || !bs->buff || num_bits <= 0) {
    return;
  }
  if (num_bits > 64) {
    num_bits = 64; // clamp to 64 bits to avoid undefined behavior 
  }

  int bit_index = bs->bit_pos;  // absolute bit index from start of buffer 

  // Write from the most significant of the requested bits down to LSB 
  for (int i = num_bits - 1; i >= 0; --i) {
    // Extract the next bit to write (MSB-first of the num_bits window) 
    uint8_t bit = (uint8_t)((val >> i) & 0x1u);

    int byte_idx = bit_index >> 3;               // which byte 
    int bit_in_byte_from_msb = 7 - (bit_index & 7); // position in that byte, MSB-first 

    if (bit) {
      bs->buff[byte_idx] |= (uint8_t)(1u << bit_in_byte_from_msb);
    }
    else {
      bs->buff[byte_idx] &= (uint8_t)~(1u << bit_in_byte_from_msb);
    }

    ++bit_index;
  }

  bs->bit_pos = bit_index;
}

extern int write_msb(const bb_meas_t* measurements, char* file_path) {
  uint8_t buff[2048];
  memset(buff, 0, sizeof(buff));
  int num_written = write_bb_msb(measurements, buff, sizeof(buff));
  FILE* fp_out2 = NULL; //output file
  errno_t er2 = fopen_s(&fp_out2, file_path, "wb");
  if (er2 != 0 || fp_out2 == NULL) { fprintf(stderr, "Failed to open output file\n"); return; }
  fwrite(buff, 1, num_written, fp_out2);
  fclose(fp_out2); 
  return num_written; 
}

// returns number of bytes written
extern int write_bb_msb(const bb_meas_t* measurements, uint8_t* bin_buff, int buff_size) {
  buff_str buff;
  buff.bit_pos = 0;
  buff.buff = bin_buff;
  memset(bin_buff, 0, buff_size);
  // write the prn mask
  int64_t prn_mask = 0;
  for (int i = 0; i < measurements->num_sat; i++) {
    if (measurements->sats[i].constellation == SYS_GPS) {
      prn_mask |= (0x1LL << (measurements->sats[i].prn - 1));
    }
  }
  bit_writer(&buff, prn_mask, 32);
  // write per prn data
  for (int i = 0; i < measurements->num_sat; i++) {
    if (measurements->sats[i].constellation == SYS_GPS) {
      int64_t code_phase = (int64_t)(measurements->sats[i].code_phase / 9.54e-7f);
      int64_t doppler = (int64_t)(measurements->sats[i].doppler / 2.0f);
      int64_t cno = (int64_t)((measurements->sats[i].cno - 28) / 4);
      bit_writer(&buff, code_phase, 20);
      bit_writer(&buff, doppler, 13);
      bit_writer(&buff, cno, 4);
    }
  }
  // write Galileo presence
  int gal_present = 0;
  for (int i = 0; i < measurements->num_sat; i++) {
    if (measurements->sats[i].constellation == SYS_GAL) {
      gal_present = 1;
      break;
    }
  }
  bit_writer(&buff, gal_present, 3);
  if (gal_present) {
    int64_t gal_prn_mask = 0;
    for (int i = 0; i < measurements->num_sat; i++) {
      if (measurements->sats[i].constellation == SYS_GAL) {
        gal_prn_mask |= (0x1LL << (measurements->sats[i].prn - 1));
      }
    }
    bit_writer(&buff, gal_prn_mask, 50);
    // write per prn data
    for (int i = 0; i < measurements->num_sat; i++) {
      if (measurements->sats[i].constellation == SYS_GAL) {
        int64_t code_phase = (int64_t)(measurements->sats[i].code_phase / 9.54e-7f);
        int64_t doppler = (int64_t)(measurements->sats[i].doppler / 2.0f);
        int64_t cno = (int64_t)((measurements->sats[i].cno - 28) / 4);
        bit_writer(&buff, code_phase, 22);
        bit_writer(&buff, doppler, 13);
        bit_writer(&buff, cno, 4);
      }
    }
  }
  return buff.bit_pos / 8;
}

/// <summary>
/// A debug utility to read the msb file and print the contents. Input is the path to the msb file.
/// </summary>
/// <param name="input"></param>
extern void read_bb_msb(uint8_t* bin_buff, int bin_size, bb_meas_t* measurements) {

  printf("Size of buffer %d \n", (int)bin_size);

  // read header for time and then the binary data after converting from hex to binary

  reset_bit_reader((int)bin_size);
  ////////////////////

  int64_t prn_mask = bit_reader(bin_buff, 32, 0);

  //memset(measurements, 0, sizeof(bb_meas_t)); //assume previously inited
  int64_t ans = prn_mask;
  int num_svs = 0;
  int used_gps_svs[30] = { 0 };
  int used_gal_svs[30] = { 0 };
  int cnt = 0;
  for (int i = 0; i < 32; i++) {
    if (ans & 0x1) {
      num_svs++;
      //printf("PRN %d is present \n", i + 1);
      used_gps_svs[cnt++] = i + 1;
    }
    ans >>= 1;
  }
  printf("num GPS svs %d \n", num_svs);
  measurements->num_sat = num_svs;


  for (int i = 0; i < num_svs; i++) {
    int64_t code_phase = bit_reader(bin_buff, 20, 0);
    int64_t doppler = bit_reader(bin_buff, 13, 1);
    int64_t cno = bit_reader(bin_buff, 4, 0);
    measurements->sats[i].time_tag = measurements->time + measurements->sec;
    measurements->sats[i].pseudorange = 0;
    measurements->sats[i].code_phase = ((float)code_phase) * 9.54e-7f;
    measurements->sats[i].cno = (float)(cno * 4 + 28);
    measurements->sats[i].doppler = (float)(2 * doppler);
    measurements->sats[i].prn = used_gps_svs[i];
    measurements->sats[i].constellation = SYS_GPS; // GPS constellation
    printf("GPS %d Code phase %f, Doppler %.1f, CNO %.1f \n",
      measurements->sats[i].prn, measurements->sats[i].code_phase * 4 * 1023, measurements->sats[i].doppler, measurements->sats[i].cno);
  }
  int64_t gal = bit_reader(bin_buff, 3, 0);
  if (gal == 1) {
    printf("Galileo is present \n");
    int64_t gal_prns = bit_reader(bin_buff, 50, 0);
    printf("Galileo PRNs %lld \n", gal_prns);
    int64_t ans = gal_prns;
    int num_svs = 0;
    int cnt = 0;
    for (int i = 0; i < 50; i++) {
      if (ans & 0x1) {
        num_svs++;
        used_gal_svs[cnt++] = i + 1;
      }
      ans >>= 1;
    }
    printf("num GAL svs %d \n", num_svs);
    int num_gps = measurements->num_sat;
    measurements->num_sat += num_svs;

    for (int i = 0; i < num_svs; i++) {
      int64_t code_phase = bit_reader(bin_buff, 22, 0);
      int64_t doppler = bit_reader(bin_buff, 13, 1);
      int64_t cno = bit_reader(bin_buff, 4, 0);
      measurements->sats[num_gps + i].time_tag = measurements->time + measurements->sec;
      measurements->sats[num_gps + i].pseudorange = 0;
      measurements->sats[num_gps + i].code_phase = ((float)code_phase) * 9.54e-7f;
      measurements->sats[num_gps + i].cno = (float)(cno * 4 + 28);
      measurements->sats[num_gps + i].doppler = (float)(2 * doppler);
      measurements->sats[num_gps + i].prn = used_gal_svs[i];
      measurements->sats[num_gps + i].constellation = SYS_GAL; // Galileo constellation
      printf("GAL %d Code phase %f, Doppler %.1f, CNO %.1f \n",
        measurements->sats[num_gps + i].prn, measurements->sats[num_gps + i].code_phase, measurements->sats[num_gps + i].doppler, measurements->sats[num_gps + i].cno);
    }
  }
  else {
    printf("Galileo is not present \n");
  }
  printf("Read in %d bits \n", cur_bit);
}
