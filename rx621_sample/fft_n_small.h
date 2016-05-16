#ifndef __FFT_N_SMALL__
#define __FFT_N_SMALL__

#include    "etc.h"

void fft_n_small_init(void);
void fft_real_n2(   double din_re[],        // npoint=2 with interva. size = 2*intrvl
                    int    intrv,           // sampling interval
                    double out_re[],        // 0, 1
                    double out_im[]);       // 0, 1
void fft_real_n3(   double din_re[],        // npoint=3 with interval
                    int    intrv,           // sampling interval
                    double out_re[],        // 0, 1
                    double out_im[]);       // 0, 1
void fft_real_n4(   double din_re[],        // npoint=4 with interva. size = 4*intrvl
                    int    intrv,           // sampling interval
                    double out_re[],        // 0, 1, 2
                    double out_im[]);       // 0, 1, 2
void fft_real_n5(   double din_re[],        // npoint=5 with interva. size = 4*intrvl
                    int    intrv,           // sampling interval
                    double out_re[],        // 0, 1, 2
                    double out_im[]);       // 0, 1, 2
void fft_real_n8(   double din_re[],        // npoint=8 with interva. size = 2*intrvl
                    int    intrv,           // sampling interval
                    double out_re[],        // 0, 1, 2, 3, 4
                    double out_im[]);       // 0, 1, 2, 3, 4
void fft_real_n16(  double din_re[],        // npoint=16 with interval. size = 16*d
                    int    d,               // sampling interval
                    double out_re[],        // 0, 1, 2, 3, ..., 8 
                    double out_im[]);       // 0, 1, 2, 3, ..., 8
void fft_mul(   int     fi_max,
                double  k,
                double  re[],
                double  im[]);
#endif
