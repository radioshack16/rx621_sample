//
// FFT
//  N=16, 8, 5, 4, 3, 2.
//
// Kazuyuki Hirooka
//
// 20160420-20160515
//
#include    <math.h>

#include    "dft.h"     // get_cos(), get_minus_sin()
#include    "sci.h"     // SCIprintf()
#include    "stack_var_adrs.h"

//--------------------------------------------------------------------------------
// init
//--------------------------------------------------------------------------------
// for fft_real_n3
static  double  w_n3_re;    // cos(-2*M_PI/3)
static  double  w_n3_im;    // sin(-2*M_PI/3)

// for fft_real_n5
static  double  w_n5_1_re;  // cos(  -2*M_PI/5)
static  double  w_n5_1_im;  // sin(  -2*M_PI/5)
static  double  w_n5_2_re;  // cos(-2*2*M_PI/5)
static  double  w_n5_2_im;  // sin(-2*2*M_PI/5)

// for fft_real_n8
static  double  w_n8_re;    // cos(-2*M_PI/8)
static  double  w_n8_im;    // sin(-2*M_PI/8)

// for fft_real_n16
static  double    w_n16_a;
static  double    w_n16_b;
static  double    w_n16_c;
static  double  m_w_n16_a;
static  double  m_w_n16_b;
static  double  m_w_n16_c;

//--------------------------------------------------------------------------------
// init entry
//--------------------------------------------------------------------------------
void    fft_n_small_init(void)
{
double th;

    th = -2.0*M_PI/3.0;
    w_n3_re = cos(th);
    w_n3_im = sin(th);

    th = -2.0*M_PI/5.0;
    w_n5_1_re = cos(th);
    w_n5_1_im = sin(th);
    w_n5_2_re = cos(2.0*th);
    w_n5_2_im = sin(2.0*th);

    w_n8_re = cos(-M_PI/4.0);
    w_n8_im = -w_n8_re;

      w_n16_a = cos(M_PI/4.0);
      w_n16_b = cos(M_PI/8.0);
      w_n16_c = sin(M_PI/8.0);
    m_w_n16_a = -w_n16_a;
    m_w_n16_b = -w_n16_b;
    m_w_n16_c = -w_n16_c;
}

//--------------------------------------------------------------------------------
// FFT, npoint=2
//
//  input:  2-point real
//  output: frequencies: {0, 1}, full.
//--------------------------------------------------------------------------------
void fft_real_n2(   double din_re[],        // npoint=2 with interval. size = 2*d
                    int    d,               // sampling interval
                    double out_re[],        // 0, 1
                    double out_im[])        // 0, 1
{
const double x0 = din_re[0*d];
const double x1 = din_re[1*d];

    out_re[0] = x0 + x1;
    out_im[0] = 0.0;

    out_re[1] = x0 - x1;
    out_im[1] = 0.0;

    stack_var_adrs_min_update();
}

//--------------------------------------------------------------------------------
// FFT, npoint=3
//
//  input:  3-point real
//  output: frequencies: {0, 1}
//--------------------------------------------------------------------------------
void fft_real_n3(   double din_re[],        // npoint=3 with interval. size = 3*d
                    int    d,               // sampling interval
                    double out_re[],        // 0, 1
                    double out_im[])        // 0, 1
{
const double x0 = din_re[0*d];
const double x1 = din_re[1*d];
const double x2 = din_re[2*d];

const double x1_pls_x2 = x1 + x2;

    out_re[0] = x0 + x1_pls_x2;
    out_im[0] = 0.0;

    out_re[1] = x0 + w_n3_re* x1_pls_x2;
    out_im[1] =      w_n3_im*(x1  -  x2);

    stack_var_adrs_min_update();
}

//--------------------------------------------------------------------------------
// FFT, npoint=4
//
//  input:  4-point real
//  output: frequencies: {0, 1, 2}          // 2=Nyquist frequency
//--------------------------------------------------------------------------------
void fft_real_n4(   double din_re[],        // npoint=4 with interval. size = 4*d
                    int    d,               // sampling interval
                    double out_re[],        // 0, 1, 2
                    double out_im[])        // 0, 1, 2
{
const double x0 = din_re[0*d];
const double x1 = din_re[1*d];
const double x2 = din_re[2*d];
const double x3 = din_re[3*d];

const double x0_p_x2 = x0 + x2;
const double x1_p_x3 = x1 + x3;

    out_re[0] = x0_p_x2 + x1_p_x3;
    out_im[0] = 0.0;

    out_re[1] = x0 - x2;
    out_im[1] = x3 - x1;

    out_re[2] = x0_p_x2 - x1_p_x3;
    out_im[2] = 0.0;

    stack_var_adrs_min_update();
}

//--------------------------------------------------------------------------------
// FFT, npoint=5
//
//  input:  5-point real
//  output: frequencies: {0, 1, 2}
//--------------------------------------------------------------------------------
void fft_real_n5(   double din_re[],        // npoint=5 with interval. size = 5*d
                    int    d,               // sampling interval
                    double out_re[],        // 0, 1, 2
                    double out_im[])        // 0, 1, 2
{
const double x0 = din_re[0*d];
const double x1 = din_re[1*d];
const double x2 = din_re[2*d];
const double x3 = din_re[3*d];
const double x4 = din_re[4*d];

const double x1_p_x4 = x1 + x4;
const double x2_p_x3 = x2 + x3;

const double x1_m_x4 = x1 - x4;
const double x2_m_x3 = x2 - x3;

    out_re[0] = x0 + x1_p_x4 + x2_p_x3;
    out_im[0] = 0.0;

    out_re[1] = x0 + (w_n5_1_re * x1_p_x4) + (w_n5_2_re * x2_p_x3);
    out_im[1] =      (w_n5_1_im * x1_m_x4) + (w_n5_2_im * x2_m_x3);

    out_re[2] = x0 + (w_n5_2_re * x1_p_x4) + (w_n5_1_re * x2_p_x3);
    out_im[2] =      (w_n5_2_im * x1_m_x4) - (w_n5_1_im * x2_m_x3);

    stack_var_adrs_min_update();
}

//--------------------------------------------------------------------------------
// FFT, npoint=8
//
//  input:  8-point real
//  output: frequencies: {0, 1, 2, 3, 4}    // 4=Nyqyist frequency
//--------------------------------------------------------------------------------
void fft_real_n8(   double din_re[],        // npoint=8 with interval. size = 8*d
                    int    d,           // sampling interval
                    double out_re[],        // 0, 1, 2, 3, 4
                    double out_im[])        // 0, 1, 2, 3, 4
{
const double x0 = din_re[0*d];
const double x1 = din_re[1*d];
const double x2 = din_re[2*d];
const double x3 = din_re[3*d];
const double x4 = din_re[4*d];
const double x5 = din_re[5*d];
const double x6 = din_re[6*d];
const double x7 = din_re[7*d];

const double x0_p_x4 = x0 + x4;
const double x1_p_x5 = x1 + x5;
const double x2_p_x6 = x2 + x6;
const double x3_p_x7 = x3 + x7;

const double x0_m_x4 = x0 - x4;
const double x1_m_x5 = x1 - x5;
const double x2_m_x6 = x2 - x6;
const double x3_m_x7 = x3 - x7;

const double x04_p_x26 = x0_p_x4 + x2_p_x6;
const double x15_p_x37 = x1_p_x5 + x3_p_x7;

const double tmp0 = w_n8_re*(x1_m_x5 - x3_m_x7);
const double tmp1 = w_n8_im*(x1_m_x5 + x3_m_x7);

    out_im[0] = 0.0;
    out_im[4] = 0.0;

    out_re[0] = x04_p_x26 + x15_p_x37;
    out_re[4] = x04_p_x26 - x15_p_x37;

    out_re[2] = x0_p_x4 - x2_p_x6;
    out_im[2] = x3_p_x7 - x1_p_x5;

    out_re[1] = x0_m_x4 + tmp0;
    out_re[3] = x0_m_x4 - tmp0;

    out_im[1] = tmp1 - x2_m_x6;
    out_im[3] = tmp1 + x2_m_x6;

    stack_var_adrs_min_update();
}

//--------------------------------------------------------------------------------
// FFT, npoint=16
//
//  input:  16-point real
//  output: frequencies: {0, 1, 2, 3, ..., 8}    // 8=Nyqyist frequency
//--------------------------------------------------------------------------------
void fft_real_n16(  double din_re[],        // npoint=16 with interval. size = 16*d
                    int    d,               // sampling interval
                    double out_re[],        // 0, 1, 2, 3, ..., 8
                    double out_im[])        // 0, 1, 2, 3, ..., 8
{
const double x_00p08 = din_re[0*d] + din_re[ 8*d];
const double x_00m08 = din_re[0*d] - din_re[ 8*d];

const double x_02p10 = din_re[2*d] + din_re[10*d];
const double x_02m10 = din_re[2*d] - din_re[10*d];

const double x_04p12 = din_re[4*d] + din_re[12*d];
const double x_04m12 = din_re[4*d] - din_re[12*d];

const double x_06p14 = din_re[6*d] + din_re[14*d];
const double x_06m14 = din_re[6*d] - din_re[14*d];

const double x_01p09 = din_re[1*d] + din_re[ 9*d];
const double x_01m09 = din_re[1*d] - din_re[ 9*d];

const double x_03p11 = din_re[3*d] + din_re[11*d];
const double x_03m11 = din_re[3*d] - din_re[11*d];

const double x_05p13 = din_re[5*d] + din_re[13*d];
const double x_05m13 = din_re[5*d] - din_re[13*d];

const double x_07p15 = din_re[7*d] + din_re[15*d];
const double x_07m15 = din_re[7*d] - din_re[15*d];

const double x_01m09_m_07m15 = x_01m09 - x_07m15;
const double x_01m09_p_07m15 = x_01m09 + x_07m15;
const double x_02m10_m_06m14 = x_02m10 - x_06m14;
const double x_02m10_p_06m14 = x_02m10 + x_06m14;
const double x_03m11_m_05m13 = x_03m11 - x_05m13;
const double x_03m11_p_05m13 = x_03m11 + x_05m13;

const double x_00p08_m_04p12 = x_00p08 - x_04p12;
const double x_01p09_m_07p15 = x_01p09 - x_07p15;
const double x_02p10_m_06p14 = x_02p10 - x_06p14;
const double x_03p11_m_05p13 = x_03p11 - x_05p13;

const double x_00p08_p_04p12 = x_00p08 + x_04p12;
const double x_01p09_p_07p15 = x_01p09 + x_07p15;
const double x_02p10_p_06p14 = x_02p10 + x_06p14;
const double x_03p11_p_05p13 = x_03p11 + x_05p13;

const double sum_even = x_00p08_p_04p12 + x_02p10_p_06p14;
const double sum_odd  = x_01p09_p_07p15 + x_03p11_p_05p13;

const double m_tmp_asum = m_w_n16_a*(x_03p11_m_05p13 + x_01p09_m_07p15);
const double   tmp_adif =   w_n16_a*(x_03p11_p_05p13 - x_01p09_p_07p15);

const double   tmp_am =     w_n16_a*x_02m10_m_06m14;
const double   tmp_ap =     w_n16_a*x_02m10_p_06m14;

const double   tmp_b1m =    w_n16_b*x_01m09_m_07m15;
const double   tmp_b1p =    w_n16_b*x_01m09_p_07m15;
const double   tmp_b3m =    w_n16_b*x_03m11_m_05m13;
const double m_tmp_b3p =  m_w_n16_b*x_03m11_p_05m13;

const double   tmp_c1m =    w_n16_c*x_01m09_m_07m15;
const double m_tmp_c1p =  m_w_n16_c*x_01m09_p_07m15;
const double   tmp_c3m =    w_n16_c*x_03m11_m_05m13;
const double   tmp_c3p =    w_n16_c*x_03m11_p_05m13;

const double   b3m_m_c1m =   tmp_b3m - tmp_c1m;
const double m_b3p_p_c1p = m_tmp_b3p + m_tmp_c1p;

const double c3m_p_b1m =  tmp_c3m + tmp_b1m;
const double c3p_m_b1p =  tmp_c3p - tmp_b1p;

const double tmp_am_sum = x_00m08 + tmp_am;
const double tmp_am_dif = x_00m08 - tmp_am;

const double tmp_ap_sum = x_04m12 + tmp_ap;
const double tmp_ap_dif = x_04m12 - tmp_ap;
// double * 52 = 8*52 = 416byte.

    out_im[0] = 0.0;
    out_im[8] = 0.0;

    out_re[0] = sum_even + sum_odd;
    out_re[8] = sum_even - sum_odd;
    out_re[4] = x_00p08_p_04p12 - x_02p10_p_06p14;
    out_im[4] = x_03p11_m_05p13 - x_01p09_m_07p15;
    out_im[6] = m_tmp_asum + x_02p10_m_06p14;
    out_im[2] = m_tmp_asum - x_02p10_m_06p14;
    out_re[6] = x_00p08_m_04p12 + tmp_adif;
    out_re[2] = x_00p08_m_04p12 - tmp_adif;
    out_re[1] = tmp_am_sum + c3m_p_b1m;
    out_re[7] = tmp_am_sum - c3m_p_b1m;
    out_re[5] = tmp_am_dif + b3m_m_c1m;
    out_re[3] = tmp_am_dif - b3m_m_c1m;
    out_im[1] = m_b3p_p_c1p - tmp_ap_sum;
    out_im[7] = m_b3p_p_c1p + tmp_ap_sum;
    out_im[5] =   c3p_m_b1p - tmp_ap_dif;
    out_im[3] =   c3p_m_b1p + tmp_ap_dif;

    stack_var_adrs_min_update();
}

void fft_mul(   int     fi_max,
                double  k,
                double  re[],
                double  im[])
{
int fi;

    for (fi=0;fi<=fi_max;fi++) {
        re[fi] *= k;
        im[fi] *= k;
    }
}
