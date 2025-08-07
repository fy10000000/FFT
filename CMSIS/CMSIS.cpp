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

int32_t to_fixed(float f, int e) {
  int a = f * pow(2.0, e);
  int b = (int)round(a);

  if (b > 0 && b > 32767) {
    return 32767;
  }
  if (b < 0 && b <= -32768) {
    return -32768;
  }
  if (b < 0) {
    // next lines turn b into it's 2's complement.
    b = abs(b);
    b = ~b;
  }
  return b;
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
  int dyad = 15; // 15 for q15_t, 31 for q31_t (nothing to do with fft_len)
  float factor = pow(2, dyad);

  // the coeffs (not twiddle factors)
#define q15_t int16_t
// consider using version without '\n' in printf for one line array (more readable than 16k line)
  for (int i = 0; i < (3 * fft_len / 4); i++) {
    float twiddleCoefq15Even = cos(i * 2 * PI / (float)fft_len);
    float twiddleCoefq15Odd = sin(i * 2 * PI / (float)fft_len);
    //printf("(q15_t) %#04hx, (q15_t) %#04hx,", (q15_t)to_fixed(twiddleCoefq15Even, dyad), (q15_t)to_fixed(twiddleCoefq15Odd, dyad));
    printf("(q15_t) %#04hx, (q15_t) %#04hx, \n",
      (q15_t)to_fixed(twiddleCoefq15Even, dyad), (q15_t)to_fixed(twiddleCoefq15Odd, dyad));
  }
}



/**
 * Main for testing and developing under Visual Studio 2022
 */
int main(int argc,char* argv[])
{
  float T = 1.0 / 500.0; // Period for the test sine wave (1/500 seconds = 2ms period)
  const int multK = 16; // multiple of 1K size can be used for 1,2,4,8,or 16 (has not been tested for less than 1K)
  float F1 = 0.0; // Frequency of the first sine wave
  float F2 = 200.0; // Frequency of the second sine wave
  float F3 = 0.0; // Frequency of the third sine wave

// File lengths are too long to do float and radix4 and radix2 serially: use compile flags
//#define DO_FLOAT
#define DO_Q15_RADIX2
//#define DO_Q15_RADIX4
//#define DO_Q31_RADIX2 need to fix compilation issues first!!!

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
  }

#ifdef DO_FLOAT
  float test[1024 * multK * 2];

  memset(test, 0, sizeof(test));
  for (int i = 0; i < 1024 * multK; i++) {
    double x = i * T;
    test[i*2] = 2000 * sin(F1 * 2.0 * PI * x) + 2000 * sin(F2 * 2.0 * PI * x);
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
    src_q15[i * 2] = (q15_t)(f2q15 (0.5 * sin(F1 * 2.0 * PI * x) + 0.5 * sin(F2 * 2.0 * PI * x)));
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
    src_q15[i*2] = (q15_t)(f2q15(0.5*sin(F1 * 2.0 * PI * x) + 0.5*sin(F2 * 2.0 * PI * x) + 0.5 * sin(F3 * 2.0 * PI * x)));
    src_q15[i*2+1] = 0; // set the imaginary part to 0
  }
 
  arm_cfft_radix2_init_q15(&s_q15, 1024 * multK, 0, 1); // Initialize the CFFT instance 
  arm_cfft_radix2_q15(&s_q15, src_q15); // fft done in place result is symmetric about 1/2 point

  for (int i = 0; i < 1024 * multK / 2; i++) {
    std::cout << "rad2_q15[" << i << "] , " << (double(i)/double(T* multK*1024)) << " , " << sqrt(src_q15[2 * i] * src_q15[2 * i] + src_q15[2 * i + 1] * src_q15[2 * i + 1]) << "\n";
  }
#endif // DO_Q15_RADIX2

#ifdef DO_Q31_RADIX2
  // now try the Q31 version warning does not compile right
  arm_cfft_radix2_instance_q31 s_q31;
  q31_t src_q31[1024 * 2 * multK];
  std::cout << "************ FFT q31_t ************************ \n";
  for (int i = 0; i < 1024 * multK; i++) {
    double x = i * T;
    src_q31[i*2] = (q31_t)f2q31(0.5 * sin(F1 * 2.0 * PI * x) + 0.5 * sin(F2 * 2.0 * PI * x));
    src_q31[i * 2 + 1] = 0;
  }
  arm_cfft_radix2_init_q31(&s_q31, 1024 * multK, 0, 1); // Initialize the CFFT instance for 8-point FFT
  arm_cfft_radix2_q31(&s_q31, src_q31);
  for (int i = 0; i < 1024 * multK / 2; i++) {
    std::cout << "rad2_q31[" << i << "] , " << (double(i)/double(T* multK*1024)) << " , " << sqrt(src_q31[2 * i] * src_q31[2 * i] + src_q31[2 * i + 1] * src_q31[2 * i + 1]) << "\n";
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


