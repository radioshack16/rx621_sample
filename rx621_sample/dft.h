//
// DFT functions.
//

#ifndef __DFT__
#define __DFT__

#include    "etc.h"

extern double        g_cos_tbl[];
extern double *g_minus_sin_tbl;

void    dft_init(int npoint);
void    dft_test(   int     method,     // 0: DFT-raw, 1: DFT, 2: FFT
                    int     mon_sw);    // 0: mon off, 1: mon on

void    dft_input_real_raw(
                        int npoint,
                        int interval,       // for interval samplng (mainly for debug purpose)
                        double din_re[],    // npoint
                        double out_re[],    // 0..floor(npoint/2)
                        double out_im[]);   // 0..floor(npoint/2)

void    dft_input_real( int npoint,
                        int interval,       // for interval samplng (mainly for debug purpose)
                        double din_re[],    // npoint
                        double out_re[],    // 0..floor(npoint/2)
                        double out_im[]);   // 0..floor(npoint/2)

void    dft_make_conjugate_part(    int     npoint,
                                    int     fi_max,
                                    double  inout_re[],     // npoint. Take care the size.
                                    double  inout_im[]);    // npoint. Take care the size.

void    dft_address_show(void);

#endif
