//
// FFT real n-general functions
//

#ifndef __FFT_REAL_N_RECURSIVE__
#define __FFT_REAL_N_RECURSIVE__

#include    "etc.h"

extern  int g_fft_status;
extern  int g_fft_init_status;

void fft_real_n_init(int npoint);
void fft_real_init_by_preset(int no);
int  fft_real_get_npoint(void);
int  fft_real_default_n_available(int npoint);

void fft_real_entry(double din_re[], double out_re[], double out_im[]);
void fft_real_recursive(    double din_re[],    // size=npoint*dintrv
                            int    npoint,      // 128, 120, ... 32, ..., 8, etc.
                            int    dintrv,      // interval for valid data
                            double out_re[],    // range={0, 1, 2, ..., fNYQ},  fNYQ=npoint/2
                            double out_im[]);   // same above
//--- for debug
void    fft_real_mon(       double      din_re[],   // size=npoint*interval
                            const int   npoint,     // e.g. 128, 120, 32, 8, 2, etc.
                            const int   interval,   // interval for valid data
                            double      out_re[],   // range={0, 1, 2, ..., fNYQ},  fNYQ=npoint/2
                            double      out_im[]);  // same above
void    fft_real_radix_tbl_mon(void);
void    fft_real_radix_mon(int cr_sw);  // 0: CR off, 1: CR on
void    fft_real_radix_preset_tbl_mon(
            int no,     // >=0: show the number, -1: show all
            int cr_sw); // 0: CR off, 1: CR on
void    fft_address_show(void);

#endif
