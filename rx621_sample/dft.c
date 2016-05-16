//
// DFT
//
// Kazuyuki Hirooka
//
// 20160409
//
#include    "iodefine.h"
#include    <math.h>

#include    "dft_fft_param.h"       // DFT/FFT size max

#include    "dft.h"
#include    "util.h"                // ms_abs_count_show()
#include    "sci.h"                 // SCIprintf()
#include    "fft_real_n_recur.h"

#define     FI_MAX_MAX      (DFT_N_MAX/2)       // e.g. 64
#define     FOUT_SIZE_MAX   (FI_MAX_MAX+1)      // e.g. 65

//--------------------------------------------------------------------------------
// cos, -sin table for DFT, FFT
//  1.25-period of cos
#define COS_TABLE_LEN_MAX   (DFT_N_MAX*5/4)     // e.g. 128 --> 160
double        g_cos_tbl[COS_TABLE_LEN_MAX];
double *g_minus_sin_tbl                     = 0;    // to be init by dft_init().
//--------------------------------------------------------------------------------
// DFT, FFT in and out
//  in:  imaginary part is zero fixed.
//  out: upper, conjugate, part is omitted.
double  g_dft_in_re[DFT_N_MAX];         // Real part only.
double  g_dft_out_re[FOUT_SIZE_MAX];    // Real      part. {0, 1, ..., (DFT_N_MAX/2)}
double  g_dft_out_im[FOUT_SIZE_MAX];    // Imaginary part. {0, 1, ..., (DFT_N_MAX/2)}
//--------------------------------------------------------------------------------
static  int     dft_n       = 0;

//----------------------------------------
// Set test signal
//----------------------------------------
static  void dft_set_signal(    int     npoint,
                                double  din_re[],
                                int     sel,
                                int     mon_sw)
{
int     i;
double  a;
double  ofst;
double  dth, th;
double  th0;
double  fi;

    for (i=0;i<npoint;i++) {
        din_re[i] = 0.0;
    }
    switch(sel) {
        default:
        case 10:     // DC + sin(freq=fi)
            // toggle at fNyq
            for (i=0;i<npoint;i++) {
                din_re[i] += (double)((i%2)*2-1);    // {-1, +1, -1, ...}
            }
            // sine waves
            //---
            fi  = 1;
            a   = 1700.0;
            ofst= 2048.0;
            dth = 2.0*M_PI/(double)npoint;
            th0 = 0.0;
            //
            for (i=0;i<npoint;i++) {
                th = dth*((double)(i*fi)) + th0;
                din_re[i] += a*sin(th) + ofst;
            }
            //---
            fi = 3;
            a  = 60000.0;
            for (i=0;i<npoint;i++) {
                th = dth*((double)(i*fi));
                din_re[i] += a*cos(th);     // Add another cos wave
            }
            //---
            fi = 7;
            a   = 123.4*sqrt(2.0)*2.0;
            th0 = M_PI/4.0;
            for (i=0;i<npoint;i++) {
                th = dth*((double)(i*fi)) + th0;
                din_re[i] += a*sin(th);     // Add another sin wave
            }
            break;
        case 0:     // DC
            for (i=0;i<npoint;i++) {
                din_re[i] = 1.0;
            }
            break;
        case 1:     // Pulse
            for (i=0;i<npoint;i++) {
                din_re[i] = 0.0;
            }
            din_re[0] = 1.0;
            break;
        case 2:
            for (i=0;i<npoint;i++) {
                din_re[i] = 0.0;
            }
            din_re[0] = 1.0;
            din_re[1] = 2.0;
            din_re[2] = 4.0;
            din_re[3] = 8.0;
            break;
        case 3:     // toggle at fNyq
            for (i=0;i<npoint;i++) {
                din_re[i] = (double)((i%2)*2-1);    // {-1, +1, -1, ...}
            }
            break;
    }
    if (mon_sw==1) {
        // Monitor
        for (i=0;i<npoint;i++) {
            SCIprintf(" i=%3d: re = %12.3le\n", i, din_re[i]);
        }
    }
}

//
// dft_init, fft init must have been done.
//
void    dft_test(   int     method,     // 0: DFT-raw, 1: DFT, 2: FFT
                    int     mon_sw)     // 0: mon off, 1: mon on
{
int     fi, findex_max;
double  t_ms;
double  t_ms_fine;
double  t_ms_coarse;
int     f_fine;

    dft_set_signal(dft_n, g_dft_in_re, 10,     // signal select
                                        0);     // mon_sw

    ms_abs_count_mark();
    tcnt_restart();
    //----------------------------------------
    // DFT or FFT
    //----------------------------------------
    if (method==0) {
        dft_input_real_raw(
                        dft_n, 1,
                        g_dft_in_re,
                        g_dft_out_re,
                        g_dft_out_im);
        t_ms = (double)ms_abs_count_get_diff();
        SCIprintf("npoint=%d, DFT_raw done, elapsed = %6.0lf[ms]\n", dft_n, t_ms);
    } else if (method==1) {
        dft_input_real( dft_n, 1,
                        g_dft_in_re,
                        g_dft_out_re,
                        g_dft_out_im);
        t_ms_fine   = tcnt_get_in_ms();
        t_ms_coarse = (double)ms_abs_count_get_diff();
        f_fine = (t_ms_coarse<=200.0) ? 1 : 0;
        if (f_fine==1) {
            SCIprintf("npoint=%d, DFT done, elapsed = %6.2lf[ms]\n", dft_n, t_ms_fine);
        } else {
            SCIprintf("npoint=%d, DFT done, elapsed = %6.0lf[ms]\n", dft_n, t_ms_coarse);
        }
    } else {
        fft_real_entry( g_dft_in_re,
                        g_dft_out_re,
                        g_dft_out_im);
        t_ms = tcnt_get_in_ms();
        fft_real_radix_mon(0);
        SCIprintf(", FFT done, elapsed = %6.2lf[ms]", t_ms);
        SCIprintf(", g_fft_init_status=%d, g_fft_status=%d\n",
                     g_fft_init_status,    g_fft_status);
    }

    //----------------------------------------
    // Monitor
    //----------------------------------------
    if (mon_sw==1) {
        SCIprintf("----------------------------------------\n");
        SCIprintf("dft/fft_input_real output:\n");
        findex_max = dft_n/2;              // e.g. 120 --> 60
        for (fi=0;fi<=findex_max;fi++) {    // e.g. 0, 1, ..., 60. Count = 61.
            SCIprintf(" fi=%2d: re = %12.3le, im = %12.3le\n", fi, g_dft_out_re[fi], g_dft_out_im[fi]);
        }
        SCIprintf("----------------------------------------\n");
    }
}

//--------------------------------------------------------------------------------
// init cos table
//--------------------------------------------------------------------------------
void    dft_init(int npoint)
{
int     i;
int     len;
double  dth, th;
double  t_ms;

    ms_abs_count_mark();

    dft_n = npoint;
    len   = dft_n*5/4;  // +1/4

    if ((dft_n % 4)!=0) {   // NG, the -sin table, g_minus_sin_tbl[], malfuncitons.
        // Fill zero for alarm
        for (i=0;i<len;i++) {
            g_cos_tbl[i] = 0.0;
        }
        g_minus_sin_tbl = &(g_cos_tbl[0]);
        return;
    }

    dth     = 2.0*M_PI/(double)dft_n;
    for (i=0;i<len;i++) {
        th  = dth*(double)i;
        g_cos_tbl[i] = cos(th);
    }

    g_minus_sin_tbl = &(g_cos_tbl[dft_n/4]);

    t_ms = (double)ms_abs_count_get_diff();
    SCIprintf("dft_init(npoint=%d) done, elapsed = %6.0lf[ms]\n", dft_n, t_ms);
}

//--------------------------------------------------------------------------------
// DFT
//  input:  real only,  N=npoints
//  output: positive frequencies only, {0 ... (floor(npoint/2))}
//          e.g.: npoints=120 --> {0, 1, ..., 60}
//
// Memo:
//  If npoint is even, the right end corresponds to the Nyquist frequency.
//  In that case, since the positinve and the negative frequency rotators
//  get to be identical, the Nyquist component will have no imaginary part,
//  and its amplitude coefficient gets to be twice the others.
//
//  And its amplitude varies depending on the sampling phase.
//
//  Practically, for wave reconstruction or analysis, the above critical
//  frequency component is almost useless.
//--------------------------------------------------------------------------------
void    dft_input_real_raw(
                        int npoint,
                        int interval,           // for interval samplng (mainly for debug purpose)
                        double din_re[],        // size = npoint * inteval
                        double out_re[],        // 0..(floor(npoint/2))
                        double out_im[])        // 0..(floor(npoint/2))
{
int     i;
int     fi, findex_max;
double  dth, th;
double  npoint_inv;
double  sum_c, sum_s;
double  x;

    if (npoint==0) return;

    findex_max = npoint/2;      // e.g. 120 --> 60

    dth = 2.0*M_PI/(double)npoint;
    npoint_inv = 1.0/(double)npoint;

    for (fi=0;fi<=findex_max;fi++) {    // e.g. 0, 1, ..., 60. Count = 61.
        sum_c=0.0;
        sum_s=0.0;
        for (i=0;i<npoint;i++) {
            x   = din_re[interval*i];
            th  = dth*(double)(fi*i);
            sum_c += x*  cos(th);       // cos() each time
            sum_s += x*(-sin(th));      // sin() each time
        }
        out_re[fi] = npoint_inv * sum_c;
        out_im[fi] = npoint_inv * sum_s;
    }
}

//--------------------------------------------------------------------------------
// DFT
// cos, sin table version
//
// dft_init() must be done in advance.
//--------------------------------------------------------------------------------
void    dft_input_real( int npoint,
                        int interval,           // for interval samplng (mainly for debug purpose)
                        double din_re[],        // size = npoint * interval
                        double out_re[],        // 0..(floor(npoint/2))
                        double out_im[])        // 0..(floor(npoint/2))
{
int     i;
int     fi, findex_max;
int     thi;
double  npoint_inv;
double  sum_c, sum_s;
double  x;

    if (npoint==0) return;

    findex_max = npoint/2;   // e.g. 120 --> 60
    npoint_inv = 1.0/(double)npoint;

    for (fi=0;fi<=findex_max;fi++) {    // e.g. 0, 1, ..., 60. Count = 61.
        sum_c=0.0;
        sum_s=0.0;
        for (i=0;i<npoint;i++) {
            x   = din_re[interval*i];
            thi = (fi*i*interval) % dft_n;
            sum_c += x*      g_cos_tbl[thi];    // from table
            sum_s += x*g_minus_sin_tbl[thi];    // from table
        }
        out_re[fi] = npoint_inv * sum_c;
        out_im[fi] = npoint_inv * sum_s;
    }
}

//--------------------------------------------------------------------------------
// Make conjugate part
//  Fill the last half part with conjugate of the first half
//--------------------------------------------------------------------------------
void dft_make_conjugate_part(   int     npoint,
                                int     fi_max,
                                double  inout_re[],     // npoint. Take care the size.
                                double  inout_im[])     // npoint. Take care the size.
{
int fi;
int fi_end;
int fi_dst;

    fi_dst=npoint-1;
    fi_end = (npoint==fi_max*2) ? fi_max-1 : fi_max;    // Omit Nyquist freq, if it is duplicate.
    for (fi=1;fi<=fi_end;fi++, fi_dst--) {  // Omit DC, fi=0.
        inout_re[fi_dst] =      inout_re[fi];
        inout_im[fi_dst] = -1.0*inout_im[fi];
    }
}

void dft_address_show(void)
{
    SCIprintf("dft_address_show():--------------------\n");
    SCIprintf("DFT_N_MAX            = %d\n",    DFT_N_MAX);
    SCIprintf("FOUT_SIZE_MAX        = %d\n",    FOUT_SIZE_MAX);
    SCIprintf("COS_TABLE_LEN_MAX    = %d\n",    COS_TABLE_LEN_MAX);
    SCIprintf("&g_cos_tbl[]         = %08Xh - %08Xh\n", &g_cos_tbl[0],      &g_cos_tbl[COS_TABLE_LEN_MAX-1]);
    SCIprintf("&g_dft_in_re[]       = %08Xh - %08Xh\n", &g_dft_in_re[0],    &g_dft_in_re[DFT_N_MAX-1]);
    SCIprintf("&g_dft_out_re[]      = %08Xh - %08Xh\n", &g_dft_out_re[0],   &g_dft_out_re[FOUT_SIZE_MAX-1]);
    SCIprintf("&g_dft_out_im[]      = %08Xh - %08Xh\n", &g_dft_out_im[0],   &g_dft_out_im[FOUT_SIZE_MAX-1]);
}

