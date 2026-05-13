/* USER CODE BEGIN Header */
/**
  ******************************************************************************
  * @file           : main.c
  * @brief          : Main program body
  ******************************************************************************
  * @attention
  *
  * Copyright (c) 2025 STMicroelectronics.
  * All rights reserved.
  *
  * This software is licensed under terms that can be found in the LICENSE file
  * in the root directory of this software component.
  * If no LICENSE file comes with this software, it is provided AS-IS.
  *
  ******************************************************************************
  */
/* USER CODE END Header */
/* Includes ------------------------------------------------------------------*/
#include "main.h"

/* Private includes ----------------------------------------------------------*/
/* USER CODE BEGIN Includes */
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include "dsp/gnss_codes.h"
#include "dsp/cplx_types.h"
#include "dsp/msb_funcs.h"
#include "dsp/controller_functions.h"
#include "dsp/transform_functions.h"
//#define PI 3.14159265358979323846
/* USER CODE END Includes */

/* Private typedef -----------------------------------------------------------*/
/* USER CODE BEGIN PTD */

/* USER CODE END PTD */

/* Private define ------------------------------------------------------------*/
/* USER CODE BEGIN PD */

/* USER CODE END PD */

/* Private macro -------------------------------------------------------------*/
/* USER CODE BEGIN PM */

/* USER CODE END PM */

/* Private variables ---------------------------------------------------------*/

/* USER CODE BEGIN PV */

/* USER CODE END PV */

/* Private function prototypes -----------------------------------------------*/
void SystemClock_Config(void);
/* USER CODE BEGIN PFP */

/* USER CODE END PFP */

/* Private user code ---------------------------------------------------------*/
/* USER CODE BEGIN 0 */


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

void convert_crlf_to_lf(char *str) {
    char *src = str;
    char *dst = str;
    while (*src) {
        if (*src == '\r' && *(src + 1) == '\n') {
            // Skip the '\r'
            src++;
        }
        *dst++ = *src++;
    }
    *dst = '\0';
}

void read_L5E5AE5QP(char* input) {
  convert_crlf_to_lf(input);
  FILE* fp_1bitcsv = fopen(input, "r");
  if (fp_1bitcsv == NULL) {
    fprintf(stderr, "Failed to open sample file %s\n", input);
    return;
  }
  fseek(fp_1bitcsv, 0L, SEEK_END);
  //size_t bytes_to_read = ftell(fp_1bitcsv);
  rewind(fp_1bitcsv);
  FILE* fp_out = NULL; //output file

#define USE_FFT 1
#define SPC 1
//#define FFT_SIZE 16384
#define FFT_SIZE 1024 // temp for qp only
#define SAMP 10230 // for 1 ms at 10.23 MHz
//#define SAMP 5115 // for 1 ms at 10.23 MHz
//#define SAMP 16384
//#define SAMP 15345
  // only 9 and 36 with q31; 10, 6 also works with float
  //int prn = 36;// 6;// 6;// 36;// 9;// 36
  //double doppler = -1*(1580 + 1e6 +2500);// E36:1580 E6: 1261 ; G10:-582, G32:1232

  char filename[256];
  strcpy(filename, input);
  *strrchr(filename, '.') = 0;
  strcat(filename, ".ast");
  FILE* fpAssist = fopen(filename, "r");
  if (fpAssist == NULL) {
    fprintf(stderr, "Failed to open assistance file %s\n", filename);
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
    sscanf(line + 1, "%d %lf", &prn2acq[cnt].prn, &prn2acq[cnt].doppler);
    if (c == 'G' && !HasGPSL5(prn2acq[cnt].prn)) continue; // won't be able to get these
    hasQp = (c == 'E') && HasE5QP(prn2acq[cnt].prn);
    if (!hasQp && SAMP < 10230) continue; // won't work
    if (!hasQp) continue; // only do QP

    prn2acq[cnt].doppler *= 1176.45 / 1575.42; // adjust to L5

    cnt++;
  }
  fclose(fpAssist);

  /////////////////////////////////////////////////////
  int maxCorrLength = USE_FFT ? FFT_SIZE : SAMP;

  c32* sampl     = (c32*)malloc(SAMP * sizeof(c32));
  c32* repli     = (c32*)malloc(SAMP * sizeof(c32));
  c32* up_samp   = (c32*)malloc(maxCorrLength * sizeof(c32));
  c32* up_repli  = (c32*)malloc(maxCorrLength * sizeof(c32));
  c32* up_prod   = (c32*)malloc(maxCorrLength * sizeof(c32));
  c32* sum_prod  = (c32*)malloc(maxCorrLength * sizeof(c32));
  float* nci_sum = (float*)malloc(maxCorrLength * sizeof(float));

  bb_meas_t meas;
  memset(&meas, 0, sizeof(bb_meas_sat_t));

  char* pc = 0;


  ///////////////////// main prn loop ////////////////////////////////
  for (int prn_loop = 0; prn_loop < cnt; prn_loop++) {
    int prn = prn2acq[prn_loop].prn;
    int gal_proc = (prn2acq[prn_loop].constel == 2) ? 1 : 0;
    hasQp = gal_proc && HasE5QP(prn);
    if (!hasQp) {
    	continue;
    }
    //EDif (gal_proc == 0 || hasQp == false) { continue; }
    //EDif ( (prn2acq[prn_loop].prn == 15 || prn2acq[prn_loop].prn == 34) == false) { continue; }
    //printf("Processing PRN %d Doppler %f constel %d \n", prn, doppler, prn2acq[prn_loop].constel);
    printf("Processing %c%02d%c around %.1lf Hz\n", gal_proc ? 'E' : 'G', prn, gal_proc ? (hasQp ? 'q' : 'a') : ' ', prn2acq[prn_loop].doppler);
    int msps = SAMP;
    int codeLength = (hasQp) ? E5_QP_CODE_LEN : E5A_CODE_LEN;
    int codeRate = (hasQp) ? 5115 : 10230;
    int spc = msps / codeRate;


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
      tc = 2 * 31 / 2; // n*31/2 for n ms
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
    int dopplerUncertainty = 0;//3 * dopplerStep;
    double ratio = 0.0;
    double bestRatio = 0.0;
    for (dopplerOffset = -dopplerUncertainty; dopplerOffset <= dopplerUncertainty; dopplerOffset += dopplerStep)
    {
      //double doppler = -1 * (prn2acq[prn_loop].doppler + 1e6 + 4100 + dopplerOffset);
      double doppler = -1 * (prn2acq[prn_loop].doppler + 1e6 + 4100 + dopplerOffset);

      int codeoff=0;
      //for (codeoff = 0; codeoff < sampLength; codeoff += sampLength/2) // try different code offsets
      {
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
            memset(repli, 0, sizeof(c32)* sampLength);
            for (int sampi = 0; sampi < sampLength; sampi++) {
              repli[sampi].r = prn_code[sampi];
            }
            free(prn_code);

          }
          else
          {
            synth_e5a_prn(prn, -doppler, sampLength, repli, 0);
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
          up_sample_N_to_M(repli, sampLength, up_repli, corrLength);

          fft_c32(corrLength, up_repli, true); // forward FFT
        }

        //char* context = NULL;
        // read in the csv data
        //char line[256];
        int NCI = msUse / tc;
        if (hasQp)
        {
          NCI = msUse / (tc * 2.0 / 31);
        }
        NCI = 1;//fixme

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
        double upsampSecs = 1.0 / (msps * 1000.0 * corrLength / sampLength);
        double dphi = 360.0 * doppler * upsampSecs;//2.0 * PI_F64 * doppler* upsampSecs; // degrees for arm_sin_cos_f32
        if (dphi < 0) dphi += 360.0;//2.0 * PI_F64;
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
            while (!feof(fp_1bitcsv)) {
              if (fgets(line, sizeof(line), fp_1bitcsv) != NULL) {
                char* token = strtok(line, ",");
                token = strtok(NULL, ",");
                if (token != NULL) {
                  sampl[idx].r = (float)atof(token);
                  token = strtok(NULL, ",");
                  if (token != NULL) { sampl[idx].i = (float)atof(token); }
                }
                idx++;
                if ((idx != 0) && (idx % LEN == 0)) { break; }
              }
            }
            // upsample
            //c32* temp_up_samp = (c32*)malloc(maxCorrLength * sizeof(c32));
            //memset(temp_up_samp, 0, sizeof(c32) * corrLength);
            up_sample_N_to_M(sampl, sampLength, up_samp, corrLength);
            // wipe Doppler
            for (int sampi = 0; sampi < corrLength; sampi++)
            {
              c32 phase;
              double phi = fmod(dphi * totalSampleCount,360.0);//2*PI_F64); // degrees for arm_sin_cos_f32
              //phase.r = cos(phi);
              //phase.i = sin(phi);
              arm_sin_cos_f32((float)phi,&phase.i,&phase.r); // phi is degrees
              c32 tempsamp1 = mult(up_samp[sampi], phase);
              c32 tempsamp2 = add(up_samp[sampi], tempsamp1);
              up_samp[sampi] = tempsamp2;
              totalSampleCount++;
            }
            //free(temp_up_samp);
          }
          if (USE_FFT)
          {

            // note repli has been FFTed already
            fft_c32(corrLength, up_samp, true); // forward FFT
            for (int k = 0; k < corrLength; k++) { up_prod[k] = mult(up_samp[k], get_conj(up_repli[k])); }
            fft_c32(corrLength, up_prod, false); // IFFT
            for (int i = 0; i < corrLength; i++) { nci_sum[i] += mag(up_prod[i]); }
          }
          else
          {
            //TimeDomainCorrelate(sampLength, sampl, repli, nci_sum);
          }
          if (1)
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
        find_top2_peaks_real(nci_sum, corrLength, 3, &peaks, fp_out);
        double cn0 = compute_snr_real(nci_sum, corrLength, peaks) + 35.0;
        double early = nci_sum[peaks.idx1 - 1], prompt = nci_sum[peaks.idx1], late = nci_sum[peaks.idx1 + 1];
        double interp = InterpolateCodePhase(peaks.idx1, early * early, prompt * prompt, late * late) * sampLength / corrLength;


        interp += codeoff;
        if (interp > sampLength)
          interp -= sampLength;

        ratio = peaks.ratio;
        printf("Avail PRN %c%02d doppler %6.0f ratio %5.2f loc %5d interp %10.4f CN0 %5.2f %5d %5d %c\n", (gal_proc == 1) ? 'E' : 'G', prn2acq[prn_loop].prn, prn2acq[prn_loop].doppler + dopplerOffset, ratio, (int)peaks.idx1, interp, cn0, dopplerOffset, codeoff, (ratio > bestRatio) ? '*' : ' ');

        if (ratio > bestRatio) {
          meas.sats[meas.num_sat].prn = prn2acq[prn_loop].prn;
          meas.sats[meas.num_sat].code_phase = (interp / msps); // [ms] QP will have some offset of 2/31 segments
          meas.sats[meas.num_sat].doppler = -(prn2acq[prn_loop].doppler + dopplerOffset);
          meas.sats[meas.num_sat].cno = (float)cn0;
          meas.sats[meas.num_sat].constellation = gal_proc ? SYS_GAL : SYS_GPS;

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
      }
      rewind(fp_1bitcsv);
    }
    if (bestRatio > 1.29)
      meas.num_sat++;
  } // end for prn_loop
  if (fp_out != NULL) fclose(fp_out);

  pc = strrchr(input, '/')+1;
  strcpy(filename, pc);
  *strrchr(filename, '.') = 0;
  strcat(filename, ".msb");
  printf("Writing %s\n", filename);
  int num_bytes = write_msb(&meas, (char*)filename);
  bb_meas_t check = { 0 };
  FILE* test = fopen(filename, "rb");
  uint8_t tbuff[128] = { 0 };
  fread(tbuff, 1, num_bytes, test);
  read_bb_msb(tbuff, num_bytes, &check);
  fclose(test);


  fclose(fp_1bitcsv);
  free(sampl); free(repli); free(up_samp); free(up_repli); free(up_prod); free(nci_sum); free(sum_prod);
}






// Provided by rdimon/newlib for semihosting
extern void initialise_monitor_handles(void);

static void copy_file_semihosting(const char* in_path, const char* out_path)
{

    FILE* fo = fopen(out_path, "wb");
    if (!fo) {
        printf("Create failed: %s\n", out_path);
        return;
    }
    //printf("I am here ############################### \n");
    //fprintf(stdout, "I am here !!!!!!!!!!!!!!!!!!! \n");
    //fprintf(stderr, "I am here %%%%%%%%%%%%%%%%%% \n");

    size_t n, total = 0;
    char buf2[] = "trying to debug this darned interface. Mary had a little lamb, and the name was\n";
    n = sizeof(buf2);
    printf("here %d ############################### \n", n);

    setvbuf(fo, NULL, _IONBF, 0);// no buffer
    int m = fwrite(buf2, 1,n, fo);
    int ret = fflush(fo);
    int ret2 = fclose(fo);
    printf("wrote %d %d %d ############################### \n", m, ret, ret2);


    FILE* fi = fopen(in_path, "r");
	if (!fi) {
		printf("Open failed: %s\n", in_path);
		return;
	}
    char buf[1024] = {0};
    while ((n = fread(buf, 1, sizeof(buf), fi)) > 0) {
        total += n;
    }
    printf("Done. Copied %lu bytes.\n", (unsigned long)total);
    fclose(fi);
}

int16_t f2q15(float value) {
  // Clamp the value to the Q15 range [-1.0, 1.0)
  if (value > 0.999969) value = 0.999969; // Slightly less than 1.0
  if (value < -1.0) value = -1.0;

  // Convert to Q15 format
  return (int16_t)round(value * 32768.0f);
}

int32_t f2q31(float value) {
  // Clamp the value to the Q15 range [-1.0, 1.0)
  if (value > 0.999969) value = 0.999969; // Slightly less than 1.0
  if (value < -1.0) value = -1.0;

  // Convert to Q15 format
  return (int32_t)round(value * 65536.0f);
}
/* USER CODE END 0 */

/**
  * @brief  The application entry point.
  * @retval int
  */
int main(void)
{

  /* USER CODE BEGIN 1 */

  /* USER CODE END 1 */

  /* MCU Configuration--------------------------------------------------------*/

  /* Reset of all peripherals, Initializes the Flash interface and the Systick. */
  HAL_Init();

  /* USER CODE BEGIN Init */

  /* USER CODE END Init */

  /* Configure the system clock */
  SystemClock_Config();

  /* USER CODE BEGIN SysInit */
  initialise_monitor_handles(); // for writing /reading to files
  setvbuf(stdout, NULL, _IONBF, 0);

  //printf("STM32 semihosting file copy start \n\r");
  //fprintf(stderr,"trying with fprintf**************** \n\r");
  //puts("trying with puts**************** \n\r");
  //write(1, "test\n", 5);

  // Use double backslashes in C strings for Windows paths
  const char* in_path  = "C:/work/temp.txt";
  const char* out_path = "result.txt";
  copy_file_semihosting(in_path, out_path); // was for testing file ops


  /* USER CODE END SysInit */

  /* Initialize all configured peripherals */
  /* USER CODE BEGIN 2 */
  char line[256] = { 0 };
  char path[512] = { 0 };
#if IS_ERIC
  const char* folder = "C:/work/baseband/utilities/2026-04-07/L5-1_10"; // system does not like the /r at end of the lines
#else
  const char* folder = "C:/work/support/esa/qp/data/t13/2026-04-20/L5_10_87"; // system does not like the /r at end of the lines
  //const char* folder = ".";
#endif
  sprintf(path, "%s/csvFiles.txt", folder);
  FILE* fp_in = fopen(path, "r");
  while (fgets(line, sizeof(line), fp_in)) {
	line[strcspn(line, "\r\n")] = 0; // handle either windows or linux line endings
    sprintf(path, "%s/%s", folder, line);
    printf("Processing %s\n", path);
    //read_E5A(path);
    // please fix path to 2026-04-07
    //strcpy(path,"C:/work/baseband/utilities/2026-04-07/L5-1_10/G_2026_04_07_17_07_51.758.smp");
    read_L5E5AE5QP(path);
  }
  if (fp_in != NULL) { fclose(fp_in); }

  //float T = 1.0 / 500.0; // Period for the test sine wave (1/500 seconds = 2ms period)
  //const int multK = 1; // multiple of 1K size can be used for 1,2,4,8,or 16 (has not been tested for less than 1K)
  //float F1 = 50.0; // Frequency of the first sine wave
  //float F2 = 80.0; // Frequency of the second sine wave
//#define DO_FLOAT
//#define DO_Q15_RADIX2
//#define DO_Q15_RADIX4
//#define DO_Q31_RADIX2
//#define DO_Q31_RADIX4

  CoreDebug->DEMCR |= CoreDebug_DEMCR_TRCENA_Msk;
  DWT->CYCCNT = 0;
  DWT->CTRL |= DWT_CTRL_CYCCNTENA_Msk;
  //unsigned long t1,t2, diff;
  // or could use

  printf("******************Begin the profiling of FFTs ********************************\n");
  /*
  FILE* fp = fopen("/myfile.txt","w");
  if (fp == NULL) {
	perror("Error opening file descriptor");
	return 1;
  }
  fprintf(fp, "This is testing for fprintf...\n");
  fputs("This is testing for fputs...\n", fp);
  fclose(fp);
  */


#ifdef DO_FLOAT
  float test[1024 * multK * 2];

  memset(test, 0, sizeof(test));
  for (int i = 0; i < 1024 * multK; i++) {
    double x = i * T;
    test[i*2] = 2000 * sin(F1 * 2.0 * PI * x) + 2000 * sin(F2 * 2.0 * PI * x);
    test[i * 2 + 1] = 0; // set the imaginary part to 0
  }

  arm_cfft_radix2_instance_f32 s;
  printf("************ float32 FFT ************************ \n");

  //FILE* fp = fopen("test.txt",'w');
  fprintf(stdout,"trying this");
  //fclose(fp);

  char buffer[50];
  sprintf(buffer, "Temp: %.2f C\r\n", 25.5);
  //HAL_UART_Transmit(&huart2, (uint8_t*)buffer, strlen(buffer), 100);


  arm_cfft_radix2_init_f32(&s, 1024 * multK, 0, 1); // Initialize the CFFT instance for 8-point FFT
  t1 = DWT->CYCCNT;
  arm_cfft_radix2_f32(&s, test);
  t2 = DWT->CYCCNT;
  diff = t2 - t1;
  printf("Timer %d\n", (int)diff); // divide by clk freq 80 MHz to get time
  //for (int i = 0; i < 1024 * multK / 2; i++) {
	//printf("float [%d] , freq , %f , val. %f \n", i, ((double)i/(double)(T* multK*1024)) , sqrt(test[2 * i] * test[2 * i] + test[2 * i + 1] * test[2 * i + 1]) );
  //}
#endif // DO_FLOAT

#if defined(DO_Q15_RADIX2) || defined(DO_Q15_RADIX4)
  q15_t src_q15[1024 * 2 * multK];
  memset(src_q15, 0, sizeof(src_q15));
#endif

#ifdef DO_Q15_RADIX4
  // try the Q15 for radix4
  arm_cfft_radix4_instance_q15 s_q15_radix4;
  printf("************ FFT q15_t Radix4 ************************ \n");
  for (int i = 0; i < 1024 * multK; i++) {
    double x = i * T;
    src_q15[i * 2] = (q15_t)(f2q15 (0.5 * sin(F1 * 2.0 * PI * x) + 0.5 * sin(F2 * 2.0 * PI * x)));
    src_q15[i * 2 + 1] = 0; // set the imaginary part to 0
  }
  arm_cfft_radix4_init_q15(&s_q15_radix4, 1024 * multK, 0, 1); // Initialize the CFFT instance for 8-point FFT

  t1 = DWT->CYCCNT;
  arm_cfft_radix4_q15(&s_q15_radix4, src_q15);
  t2 = DWT->CYCCNT;
  diff = t2 - t1;
  printf("Timer %d \n", (int)diff);
  //for (int i = 0; i < 1024 * multK / 2; i++) {
	//printf("rad4 [%d] , freq , %f , val. %f \n", i, ((double)i/(double)(T* multK*1024)) , sqrt(src_q15[2 * i] * src_q15[2 * i] + src_q15[2 * i + 1] * src_q15[2 * i + 1]) );
  //}
#endif // DO_Q15_RADIX4

#ifdef DO_Q15_RADIX2
  arm_cfft_radix2_instance_q15 s_q15;
  printf("************ FFT q15_t Radix2 ************************ \n");
  for (int i = 0; i < 1024 * multK; i++) {
    double x = i * T;
    src_q15[i*2] = (q15_t)(f2q15(0.5*sin(F1 * 2.0 * PI * x) + 0.5*sin(F2 * 2.0 * PI * x)));
    src_q15[i*2+1] = 0; // set the imaginary part to 0
  }
  arm_cfft_radix2_init_q15(&s_q15, 1024 * multK, 0, 1); // Initialize the CFFT instance

  t1 = DWT->CYCCNT;
  arm_cfft_radix2_q15(&s_q15, src_q15); // fft done in place result is symmetric about 1/2 point
  t2 = DWT->CYCCNT;
  diff = t2 - t1;
  printf("Timer %d \n", (int)diff);

  //for (int i = 0; i < 1024 * multK / 2; i++) {
	//printf("rad2 [%d] , freq , %f , val. %f \n", i, ((double)i/(double)(T* multK*1024)) , sqrt(src_q15[2 * i] * src_q15[2 * i] + src_q15[2 * i + 1] * src_q15[2 * i + 1]) );
  //}
#endif // DO_Q15_RADIX2

#if defined(DO_Q31_RADIX2) || defined(DO_Q31_RADIX4)
  q31_t src_q31[1024 * 2 * multK];
  memset(src_q31, 0, sizeof(src_q31));
#endif

#ifdef DO_Q31_RADIX2
  // now try the Q31 version warning does not compile right
  arm_cfft_radix2_instance_q31 s_q31;

  printf("************ FFT q31_t radix 2 ************************ \n");
  for (int i = 0; i < 1024 * multK; i++) {
    double x = i * T;
    src_q31[i*2] = (q31_t)f2q31(0.5 * sin(F1 * 2.0 * PI * x) + 0.5 * sin(F2 * 2.0 * PI * x));
    src_q31[i * 2 + 1] = 0;
  }
  arm_cfft_radix2_init_q31(&s_q31, 1024 * multK, 0, 1); // Initialize the CFFT instance for 8-point FFT
  t1 = DWT->CYCCNT;
  arm_cfft_radix2_q31(&s_q31, src_q31);
  t2 = DWT->CYCCNT;
  diff = t2 - t1;
  printf("Timer %d \n", (int)diff);

  //for (int i = 0; i < 1024 * multK / 2; i++) {
  //  printf("rad2 [%d] , freq , %f , val. %f \n", i, ((double)i/(double)(T* multK*1024)) , sqrt(src_q31[2 * i] * src_q31[2 * i] + src_q31[2 * i + 1] * src_q31[2 * i + 1]) );
  //}
#endif //DO_Q31_RADIX2

#ifdef DO_Q31_RADIX4
  // now try the Q31 version warning does not compile right
  arm_cfft_radix4_instance_q31 s_q31_r4;

  printf("************ FFT q31_t radix 2 ************************ \n");
  for (int i = 0; i < 1024 * multK; i++) {
    double x = i * T;
    src_q31[i*2] = (q31_t)f2q31(0.5 * sin(F1 * 2.0 * PI * x) + 0.5 * sin(F2 * 2.0 * PI * x));
    src_q31[i * 2 + 1] = 0;
  }
  arm_cfft_radix4_init_q31(&s_q31_r4, 1024 * multK, 0, 1); // Initialize the CFFT instance for 8-point FFT
  t1 = DWT->CYCCNT;
  arm_cfft_radix4_q31(&s_q31_r4, src_q31);
  t2 = DWT->CYCCNT;
  diff = t2 - t1;
  printf("Timer %d \n", (int)diff);

  //for (int i = 0; i < 1024 * multK / 2; i++) {
  //  printf("rad4 [%d] , freq , %f , val. %f \n", i, ((double)i/(double)(T* multK*1024)) , sqrt(src_q31[2 * i] * src_q31[2 * i] + src_q31[2 * i + 1] * src_q31[2 * i + 1]) );
  //}
#endif //DO_Q31_RADIX2

  printf("Done profiling\n");
  /* USER CODE END 2 */

  /* Infinite loop */
  /* USER CODE BEGIN WHILE */
  //while (1)
  //{
    /* USER CODE END WHILE */

    /* USER CODE BEGIN 3 */
  //}
  /* USER CODE END 3 */
  return 0;
}

/**
  * @brief System Clock Configuration
  * @retval None
  */
void SystemClock_Config(void)
{
  RCC_OscInitTypeDef RCC_OscInitStruct = {0};
  RCC_ClkInitTypeDef RCC_ClkInitStruct = {0};

  /** Configure the main internal regulator output voltage
  */
  if (HAL_PWREx_ControlVoltageScaling(PWR_REGULATOR_VOLTAGE_SCALE1) != HAL_OK)
  {
    Error_Handler();
  }

  /** Initializes the RCC Oscillators according to the specified parameters
  * in the RCC_OscInitTypeDef structure.
  */
  RCC_OscInitStruct.OscillatorType = RCC_OSCILLATORTYPE_MSI;
  RCC_OscInitStruct.MSIState = RCC_MSI_ON;
  RCC_OscInitStruct.MSICalibrationValue = 0;
  RCC_OscInitStruct.MSIClockRange = RCC_MSIRANGE_6;
  RCC_OscInitStruct.PLL.PLLState = RCC_PLL_NONE;
  if (HAL_RCC_OscConfig(&RCC_OscInitStruct) != HAL_OK)
  {
    Error_Handler();
  }

  /** Initializes the CPU, AHB and APB buses clocks
  */
  RCC_ClkInitStruct.ClockType = RCC_CLOCKTYPE_HCLK|RCC_CLOCKTYPE_SYSCLK
                              |RCC_CLOCKTYPE_PCLK1|RCC_CLOCKTYPE_PCLK2;
  RCC_ClkInitStruct.SYSCLKSource = RCC_SYSCLKSOURCE_MSI;
  RCC_ClkInitStruct.AHBCLKDivider = RCC_SYSCLK_DIV1;
  RCC_ClkInitStruct.APB1CLKDivider = RCC_HCLK_DIV1;
  RCC_ClkInitStruct.APB2CLKDivider = RCC_HCLK_DIV1;

  if (HAL_RCC_ClockConfig(&RCC_ClkInitStruct, FLASH_LATENCY_0) != HAL_OK)
  {
    Error_Handler();
  }
}

/* USER CODE BEGIN 4 */
/*
int _write(int file, char *ptr, int len)
{
  (void)file;
  int DataIdx;

  for (DataIdx = 0; DataIdx < len; DataIdx++)
  {
    ITM_SendChar(*ptr++);
  }
  return len;
}
*/
/* USER CODE END 4 */

/**
  * @brief  This function is executed in case of error occurrence.
  * @retval None
  */
void Error_Handler(void)
{
  /* USER CODE BEGIN Error_Handler_Debug */
  /* User can add his own implementation to report the HAL error return state */
  __disable_irq();
  while (1)
  {
  }
  /* USER CODE END Error_Handler_Debug */
}
#ifdef USE_FULL_ASSERT
/**
  * @brief  Reports the name of the source file and the source line number
  *         where the assert_param error has occurred.
  * @param  file: pointer to the source file name
  * @param  line: assert_param error line source number
  * @retval None
  */
void assert_failed(uint8_t *file, uint32_t line)
{
  /* USER CODE BEGIN 6 */
  /* User can add his own implementation to report the file name and line number,
     ex: printf("Wrong parameters value: file %s on line %d\r\n", file, line) */
  /* USER CODE END 6 */
}
#endif /* USE_FULL_ASSERT */
