//
// FFT real
// N = 4*(2^p)*(3^q)*(5^r)*(7^s)
//  e.g.: N= 128, 120, 112, 100, 96, ...
//
// Recursive
//
// Kazuyuki Hirooka
//
// 20160423-20160515
//
#include    <math.h>
#include    <string.h>              // memset()

#include    "dft_fft_param.h"       // DFT/FFT size max

#include    "dft.h"                 // get_cos(), get_minus_sin()
#include    "fft_n_small.h"
#include    "fft_real_n_recur.h"
#include    "util.h"                // ms_abs_count_show()
#include    "sci.h"                 // SCIprintf()

#define     FFT_N_STG1          (FFT_N_MAX/4)           // 32 (=128/4)
#define     FFT_N_STG2          (FFT_N_STG1/4)          //  8
#define     FFT_OUT_LEN_STG1    ((FFT_N_STG1)/2+1)      // 32/2+1=17
#define     FFT_OUT_LEN_STG2    ((FFT_N_STG2)/2+1)      //  8/2+1=5

#define     FFT_DATA_HEAP_LEN_STRICT    ((FFT_OUT_LEN_STG1+FFT_OUT_LEN_STG2)*4*2)   // 4:Radix4, 2:real and imaginary
#define     FFT_DATA_HEAP_LEN           (FFT_DATA_HEAP_LEN_STRICT*12/10)            // +20%

int g_fft_init_status   = -1;
int g_fft_status        = -1;

//----------------------------------------
// data stack area
//----------------------------------------
static  double  fft_data_heap[FFT_DATA_HEAP_LEN];
static  int     fft_data_heap_ptr       = 0;
static  int     fft_data_heap_ptr_max   = 0;    // Debug purpose

    // //----------------------------------------
    // // Sample of sub FFT work memory
    // //------
    // // Stage-1: N=128, Radix-4                         // 17(=16+1) is enough, not 32.
    // double  stg1_ph0_out_re[17];  double  stg1_ph0_out_im[17];  // Half+1
    // double  stg1_ph1_out_re[17];  double  stg1_ph1_out_im[17];
    // double  stg1_ph2_out_re[17];  double  stg1_ph2_out_im[17];
    // double  stg1_ph3_out_re[17];  double  stg1_ph3_out_im[17];
    // //------
    // // Stage-2: N=32, Radix-4                          // 5(=4+1) is enough, not 8.
    // double  stg2_ph0_out_re[5];   double  stg2_ph0_out_im[5];   // Half+1
    // double  stg2_ph1_out_re[5];   double  stg2_ph1_out_im[5];
    // double  stg2_ph2_out_re[5];   double  stg2_ph2_out_im[5];
    // double  stg2_ph3_out_re[5];   double  stg2_ph3_out_im[5];
    // //------

//----------------------------------------
// radix definition for N
//----------------------------------------
#define RADIX_TBL_SIZE      (FFT_N_MAX+1)       // index can be FFT_N_MAX.
static  unsigned char radix_tbl[RADIX_TBL_SIZE];

static  int     fft_n       = 0;
static  double  w_n4_re;    // cos(-M_PI/4)

static  int     fft_radix_preset_no     = 0;

const unsigned char radix_preset[][9] = {
                                                        // *1: not work. Measuered with n16 removed from recursive()
                                                        // *2: not work. Measuered with n16, n8 removed from recursive()
                                                        // *3: not work. Measuered with n16 ,n8, n4 removed from recursive()
                                                        // *4: not work from out of memory
    // N, r0, r1, r2,r3,r4,r5,r6    // N=r0*r1*...*r<k>; 0: terminator
    // 1: exclusion mark
    {128,  8, 16,  0, 0, 0, 0, 0, 0},   // 0    16              // 12.42ms*
    {  1,  8,  8,  2, 0, 0, 0, 0, 0},   // 1    16, 2           // 19.84ms  *1
    {  1,  8,  4,  4, 0, 0, 0, 0, 0},   // 2    16, 4           // 16.90ms  *1
    {  1,  8,  4,  2, 2, 0, 0, 0, 0},   // 3    16, 4, 2        // 21.24ms  *3
    {  1,  8,  2,  8, 0, 0, 0, 0, 0},   // 4    16, 8           // 14.81ms  *2
    {  1,  8,  2,  4, 2, 0, 0, 0, 0},   // 5    16, 8, 2        // 20.69ms  *2
    {  1,  8,  2,  2, 4, 0, 0, 0, 0},   // 6    16, 8, 4        // 17.81ms  *2
    {  1,  8,  2,  2, 2, 2, 0, 0, 0},   // 7    16, 8, 4, 2     // 21.84ms  *3
    {128,  4,  8,  4, 0, 0, 0, 0, 0},   // 8    32, 4           // 16.64ms
    {  1,  4,  8,  2, 2, 0, 0, 0, 0},   // 9    32, 4, 2        // 21.09ms  *3
    {128,  4,  4,  8, 0, 0, 0, 0, 0},   // 10   32, 8           // 14.63ms
    {  1,  4,  4,  4, 2, 0, 0, 0, 0},   // 11   32, 8, 2        // 20.65ms  *2
    {  1,  4,  4,  2, 4, 0, 0, 0, 0},   // 12   32, 8, 4        // 17.77ms  *2
    {  1,  4,  4,  2, 2, 2, 0, 0, 0},   // 13   32, 8, 4, 2     // 21.84ms  *3
    {128,  4,  2, 16, 0, 0, 0, 0, 0},   // 14   32, 16          // 13.01ms
    {  1,  4,  2,  8, 2, 0, 0, 0, 0},   // 15   32, 16, 2       // 20.40ms  *2
    {  1,  4,  2,  4, 4, 0, 0, 0, 0},   // 16   32, 16, 4       // 17.45ms  *2
    {  1,  4,  2,  4, 2, 2, 0, 0, 0},   // 17   32, 16, 4, 2    // 21.77ms  *3
    {  1,  4,  2,  2, 8, 0, 0, 0, 0},   // 18   32, 16, 2       // 15.28ms  *1
    {  1,  4,  2,  2, 4, 2, 0, 0, 0},   // 19   32, 16, 8, 2    // 21.17ms  *2
    {  1,  4,  2,  2, 2, 4, 0, 0, 0},   // 20   32, 16, 8, 4    // 18.27ms  *2
    {  1,  4,  2,  2, 2, 2, 2, 0, 0},   // 21   32, 16, 8, 4, 2     //      *3, *4
    {  1,  2,  2,  2, 2, 2, 2, 2, 0},   // 22   64, 32, 16, 8, 4, 2 //      *3, *4
    //----------------------------------------
    {120,  8,  5,  3, 0, 0, 0, 0, 0},   // 23   15,  3          // 22.182ms
    {120,  8,  3,  5, 0, 0, 0, 0, 0},   // 24   15,  5          // 18.125ms*
    {120,  5,  8,  3, 0, 0, 0, 0, 0},   // 25   24,  3          // 23.117ms
    {120,  5,  4,  3, 2, 0, 0, 0, 0},   // 26   24,  6, 2       // 26.304ms
    {120,  5,  4,  2, 3, 0, 0, 0, 0},   // 27   24,  6, 3       // 24.000ms
    {120,  5,  3,  8, 0, 0, 0, 0, 0},   // 28   24,  8          // 20.774ms
    {120,  5,  2,  4, 3, 0, 0, 0, 0},   // 29   24, 12, 3       // 23.757ms
    {120,  5,  2,  3, 4, 0, 0, 0, 0},   // 30   24, 12, 4       // 23.411ms
    {120,  4,  5,  3, 2, 0, 0, 0, 0},   // 31   30,  6, 2       // 26.189ms
    {120,  4,  5,  2, 3, 0, 0, 0, 0},   // 32   30,  6, 3       // 24.064ms
    {120,  4,  3,  5, 2, 0, 0, 0, 0},   // 33   30, 10, 2       // 24.717ms
    {120,  4,  3,  2, 5, 0, 0, 0, 0},   // 34   30, 10, 5       // 19.008ms
    {120,  4,  2,  5, 3, 0, 0, 0, 0},   // 35   30, 15, 3       // 22.643ms
    {120,  4,  2,  3, 5, 0, 0, 0, 0},   // 36   30, 15, 5       // 18.573ms
    {120,  3,  8,  5, 0, 0, 0, 0, 0},   // 37   40,  5          // 18.214ms
    {120,  3,  5,  8, 0, 0, 0, 0, 0},   // 38   40,  8          // 20.570ms
    //----------------------------------------
    {112,  8,  7,  2, 0, 0, 0, 0, 0},   // 39   14, 2           // 23.67ms
    {112,  7, 16,  0, 0, 0, 0, 0, 0},   // 40   16              // 19.52ms*
    {  1,  7,  2,  8, 0, 0, 0, 0, 0},   // 41   16, 8           // 21.288ms *1
    {  1,  7,  4,  4, 0, 0, 0, 0, 0},   // 42   16, 4           // 24.000ms *2
    {  1,  4,  7,  4, 0, 0, 0, 0, 0},   // 43   28, 4           // 22.566ms *2
    //----------------------------------------
    {108,  4,  3,  3, 3, 0, 0, 0, 0},   // 44   27,  3, 3       // 18.65ms*
    {108,  3,  4,  3, 3, 0, 0, 0, 0},   // 45   36,  9, 3       // 18.73ms
    {108,  3,  3,  4, 3, 0, 0, 0, 0},   // 46   36, 12, 3       // 19.15ms
    {108,  3,  3,  3, 4, 0, 0, 0, 0},   // 47   36, 12, 4       // 18.79ms
    //----------------------------------------
    {100,  5,  5,  4, 0, 0, 0, 0, 0},   // 48   20, 4           // 21.11ms
    {100,  5,  4,  5, 0, 0, 0, 0, 0},   // 49   20, 5           // 18.07ms
    {100,  4,  5,  5, 0, 0, 0, 0, 0},   // 50   25, 5           // 17.54ms*
    //----------------------------------------
    { 96,  8,  4,  3, 0, 0, 0, 0, 0},   // 51   12, 3           // 13.491ms
    { 96,  8,  3,  4, 0, 0, 0, 0, 0},   // 52   12, 4           // 13.171ms
    { 96,  4,  8,  3, 0, 0, 0, 0, 0},   // 53   24, 3           // 13.286ms
    { 96,  4,  4,  3, 2, 0, 0, 0, 0},   // 54   24, 6, 2        // 15.834ms
    { 96,  4,  4,  2, 3, 0, 0, 0, 0},   // 55   24, 6, 3        // 14.042ms
    { 96,  4,  3,  8, 0, 0, 0, 0, 0},   // 56   24, 8           // 11.520ms
    { 96,  4,  2,  4, 3, 0, 0, 0, 0},   // 57   24, 12, 3       // 13.837ms
    { 96,  4,  2,  3, 4, 0, 0, 0, 0},   // 58   24, 12, 4       // 13.542ms
    { 96,  3,  2, 16, 0, 0, 0, 0, 0},   // 59   32, 16          // 10.368ms
    { 96,  2,  3, 16, 0, 0, 0, 0, 0},   // 60   48, 16          // 10.227ms*
                                    //    $vvv TO MODIFIED
    //----------------------------------------
    { 84,  7,  4,  3, 0, 0, 0, 0, 0},   // 61   12, 3           // 18.64ms
    { 84,  7,  3,  4, 0, 0, 0, 0, 0},   // 62   12, 4           // 18.34ms
    { 84,  4,  7,  3, 0, 0, 0, 0, 0},   // 63   21, 3           // 17.50ms*
    //----------------------------------------
    { 80,  8,  5,  2, 0, 0, 0, 0, 0},   // 64   10, 2           // 14.170ms
    { 80,  8,  2,  5, 0, 0, 0, 0, 0},   // 65   10, 5           // 10.278ms
    { 80,  5, 16,  0, 0, 0, 0, 0, 0},   // 66   16              // 10.406ms
    { 80,  4,  5,  4, 0, 0, 0, 0, 0},   // 67   20, 4           // 12.749ms
    { 80,  4,  4,  5, 0, 0, 0, 0, 0},   // 68   20, 5           // 10.227ms
    { 80,  4,  2,  5, 2, 0, 0, 0, 0},   // 69   20, 5           // 14.554ms
    { 80,  4,  2,  2, 5, 0, 0, 0, 0},   // 70   20, 5           // 10.650ms
    { 80,  2,  8,  5, 0, 0, 0, 0, 0},   // 71   40, 5           // 10.125ms*
    { 80,  2,  5,  8, 0, 0, 0, 0, 0},   // 72   40, 5           // 11.648ms
    { 80,  2,  4,  5, 2, 0, 0, 0, 0},   // 73   40, 5           // 14.490ms
    { 80,  2,  4,  2, 5, 0, 0, 0, 0},   // 74   40, 5           // 10.560ms
    { 80,  2,  2,  5, 4, 0, 0, 0, 0},   // 75   40, 5           // 13.082ms
    { 80,  2,  2,  4, 5, 0, 0, 0, 0},   // 76   40, 5           // 10.547ms
    //----------------------------------------
    { 72,  8,  3,  3, 0, 0, 0, 0, 0},   // 77    9, 3           // 10.240ms
    { 72,  4,  2,  3, 3, 0, 0, 0, 0},   // 78   18, 9, 3        // 10.560ms
    { 72,  4,  3,  2, 3, 0, 0, 0, 0},   // 79   18, 6, 3        // 11.046ms
    { 72,  4,  3,  3, 2, 0, 0, 0, 0},   // 80   18, 9, 2        // 12.390ms
    { 72,  3,  8,  3, 0, 0, 0, 0, 0},   // 81   24, 3           // 10.432ms
    { 72,  3,  4,  3, 2, 0, 0, 0, 0},   // 82   24, 6, 2        // 12.339ms
    { 72,  3,  4,  2, 3, 0, 0, 0, 0},   // 83   24, 6, 3        // 10.995ms
    { 72,  3,  3,  8, 0, 0, 0, 0, 0},   // 84   24, 8           //  9.075ms*
    { 72,  3,  2,  4, 3, 0, 0, 0, 0},   // 85   24, 12, 3       // 10.803ms
    { 72,  3,  2,  3, 4, 0, 0, 0, 0},   // 86   24, 12, 4       // 10.573ms
    //----------------------------------------
    { 64,  8,  8,  0, 0, 0, 0, 0, 0},   // 87    8              //  5.709ms
    { 64,  4, 16,  0, 0, 0, 0, 0, 0},   // 88   16              //  4.915ms*
    {  1,  4,  4,  4, 0, 0, 0, 0, 0},   // 89   16, 4           //  7.142ms *2
    {  1,  4,  2,  2, 2, 2, 0, 0, 0},   // 90   16, 8, 4, 2     //  9.408ms *3
    { 64,  2,  8,  4, 0, 0, 0, 0, 0},   // 91   32, 4           //  7.002ms
    { 64,  2,  4,  8, 0, 0, 0, 0, 0},   // 92   32, 8           //  6.016ms
    {  1,  2,  2, 16, 0, 0, 0, 0, 0},   // 93   32, 16          //  5.222ms *3
    {  1,  2,  2,  2, 2, 2, 2, 0, 0},   // 94   32, 16, 8, 4, 2 //  9.651ms *3
    //----------------------------------------
    { 60,  4,  5,  3, 0, 0, 0, 0, 0},   // 95   15, 3           //  9.882ms
    { 60,  4,  3,  5, 0, 0, 0, 0, 0},   // 96   15, 5           //  7.808ms*
    { 60,  3,  5,  4, 0, 0, 0, 0, 0},   // 97   20, 4           //  9.856ms
    { 60,  3,  4,  5, 0, 0, 0, 0, 0},   // 98   20, 5           //  7.962ms
    //----------------------------------------
    { 56,  7,  8,  0, 0, 0, 0, 0, 0},   // 99   8               //  9.267ms
    //----------------------------------------
    { 48,  8,  3,  2, 0, 0, 0, 0, 0},   // 100    6,  2         //  6.746ms
    { 48,  4,  4,  3, 0, 0, 0, 0, 0},   // 101   12,  3         //  5.747ms
    { 48,  4,  3,  4, 0, 0, 0, 0, 0},   // 102   12,  4         //  5.581ms
    { 48,  4,  2,  3, 2, 0, 0, 0, 0},   // 103   12,  6, 2      //  6.976ms
    { 48,  4,  2,  2, 3, 0, 0, 0, 0},   // 104   12,  6, 3      //  6.080ms
    { 48,  3, 16,  0, 0, 0, 0, 0, 0},   // 105   16             //  3.942ms*
    { 48,  2,  4,  3, 2, 0, 0, 0, 0},   // 106   24,  6, 2      //  6.886ms
    { 48,  2,  8,  3, 0, 0, 0, 0, 0},   // 107   24,  3         //  5.696ms
    { 48,  2,  4,  3, 2, 0, 0, 0, 0},   // 108   24,  6, 2      //  6.886ms
    { 48,  2,  4,  2, 3, 0, 0, 0, 0},   // 109   24,  6, 3      //  6.003ms
    { 48,  2,  3,  8, 0, 0, 0, 0, 0},   // 110   24,  8         //  4.762ms
    { 48,  2,  2,  4, 3, 0, 0, 0, 0},   // 111   24, 12, 3      //  5.926ms
    { 48,  2,  2,  3, 4, 0, 0, 0, 0},   // 112   24, 12, 4      //  5.786ms
    { 48,  2,  2,  2, 2, 3, 0, 0, 0},   // 113   24, 12, 6, 3   //  6.285ms
    //----------------------------------------
    { 40,  8,  5,  0, 0, 0, 0, 0, 0},   // 114   5              //  4.045ms*
    { 40,  5,  8,  0, 0, 0, 0, 0, 0},   // 115   8              //  4.787ms
    //----------------------------------------
    { 36,  4,  3,  3, 0, 0, 0, 0, 0},   // 116    9, 3          //  4.378ms*
    { 36,  3,  4,  3, 0, 0, 0, 0, 0},   // 117   12, 3          //  4.518ms
    { 36,  3,  3,  4, 0, 0, 0, 0, 0},   // 118   12, 4          //  4.442ms
    //----------------------------------------
    { 32,  8,  4,  0, 0, 0, 0, 0, 0},   // 119    4             //  2.701ms
    { 32,  4,  8,  0, 0, 0, 0, 0, 0},   // 120    8             //  2.214ms
    { 32,  2, 16,  0, 0, 0, 0, 0, 0},   // 121   16             //  1.818ms*
    //----------------------------------------
    { 28,  7,  4,  0, 0, 0, 0, 0, 0},   // 122   4              //  4.339ms
    //----------------------------------------
    { 24,  8,  3,  0, 0, 0, 0, 0, 0},   // 123   3              //  2.240ms
    { 24,  3,  8,  0, 0, 0, 0, 0, 0},   // 124   8              //  1.792ms*
    //----------------------------------------
    { 20,  5,  4,  0, 0, 0, 0, 0, 0},   // 125   4              //  2.253ms
    { 20,  4,  5,  0, 0, 0, 0, 0, 0},   // 126   5              //  1.613ms*
    //----------------------------------------
    { 16, 16,  0,  0, 0, 0, 0, 0, 0},   // 127                  //  0.512ms
    //--------------------------------------------------------------------------------
    {  0,  0,  0,  0, 0, 0, 0, 0, 0}    // END
};

static void fft_real_radix_tbl_clear(void)
{
    memset(radix_tbl, 0, sizeof(radix_tbl));
}

//
// Assumed: npoint%4=0
//
//  0: NA(Not Assigned), else: assgined
//
static int  fft_real_radix_top(int npoint, int f_first)
{
int     i, m;
int     fft_n_thresh = FFT_N_MAX*17/20;   // =0.85=:=108/128
int     f_small_prohibit;

    // Top radix is to be decided by having decided the teriminal.
    // Terminal radix should be one of 
    static int  trm[5] = {16, 5, 8, 3, 4};
    // for better performance in the order(more left, more better).
    // (There is no terminal 7. 2 is the worst.)

    f_small_prohibit = (f_first==1 && npoint > fft_n_thresh) ? 1 : 0;
    for (i=0;i<5;i++) {
        if (npoint%trm[i]==0) {
            m=npoint/trm[i];
            if (m==1) return trm[i];    // The terminal.
            if (m%8==0) return 8;
            if (m%4==0) return 4;
            if (m%7==0) return 7;
            if (m%5==0) return 5;
            if (f_small_prohibit==1) {  // Select 4 for smaller fft output buf size.
                if (m%2==0) return 4;
                if (m%3==0) return 4;
            } else {
                if (m%2==0) return 2;
                if (m%3==0) return 3;
            }
            return 0;   // NA.
        }
    }
    return 0;   // NA.
}

//
//  0/1: fail/success
//
// radix_tbl[k], k<=16, assumed to have been set.
static int  fft_real_radix_tbl_set_auto(int npoint)
{
int i, m, n;
int r;
int r_tmp[32];
int f_first     = 1;
int f_success   = 0;

    if (npoint%4!=0) return f_success;  // NG.  (from cos/sin table nature)

    // Decide radices
    memset(r_tmp, 0, sizeof(r_tmp));
    m=npoint;
    i=0;
    while(m>1 && i<32) {
        r = fft_real_radix_top(m, f_first);
        if (r==0) break;    // Fail
        r_tmp[i] = r;
        f_first = 0;
        m/=r_tmp[i];
        i++;
        if (m==1) {
            f_success = 1;  // Success
            break;
        }
    }
    n=i;
    if (f_success==0) return f_success;

    // Set the radices in the table
    m=npoint;
    for (i=0;i<n;i++) {
        radix_tbl[m]=r_tmp[i];
        m/=radix_tbl[m];
    }
    return f_success;
}

// Manually chosen data set
// only for npoint<128.
//
// The table will be changed by running fft_real_radix_tbl_set_auto().
//
static void fft_real_radix_tbl_set_default(void)
{
    fft_real_radix_tbl_clear();

    // Data: Optimum radix combination
    // {128,  8, 16,  0, 0, 0, 0, 0, 0},   // 0    16              // 12.42 ms*
    // {120,  8,  3,  5, 0, 0, 0, 0, 0},   // 24   15,  5          // 18.125ms*
    // {112,  7, 16,  0, 0, 0, 0, 0, 0},   // 38   16              // 19.52ms*
    // {108,  4,  3,  3, 3, 0, 0, 0, 0},   // 41   27,  3, 3       // 18.65 ms*
    // {100,  4,  5,  5, 0, 0, 0, 0, 0},   // 47   25, 5           // 17.54 ms*
    // { 96,  2,  3, 16, 0, 0, 0, 0, 0},   // 55   48, 16          // 10.227ms*
    //           ~~~~~~~-->48
    // { 84,  4,  7,  3, 0, 0, 0, 0, 0},   // 58   21, 3           // 17.50ms*
    // { 80,  2,  8,  5, 0, 0, 0, 0, 0},   // 66   40, 5           // 10.125ms*
    //           ~~~~~~~-->40
    // { 72,  3,  3,  8, 0, 0, 0, 0, 0},   // 79   24, 8           //  9.075ms*
    //           ~~~~~~~-->24
    // { 64,  4, 16,  0, 0, 0, 0, 0, 0},   // 83   16              //  4.915ms*
    // { 60,  4,  3,  5, 0, 0, 0, 0, 0},   // 88   15, 5           //  7.808ms*
    // { 56,  7,  8,  0, 0, 0, 0, 0, 0},   // 91   8               //  9.267ms
    // { 48,  3, 16,  0, 0, 0, 0, 0, 0},   // 97   16              //  3.942ms*
    // { 40,  8,  5,  0, 0, 0, 0, 0, 0},   // 106   5              //  4.045ms*
    // { 36,  4,  3,  3, 0, 0, 0, 0, 0},   // 108    9, 3          //  4.378ms*
    // { 32,  2, 16,  0, 0, 0, 0, 0, 0},   // 113   16             //  1.818ms*
    // { 28,  7,  4,  0, 0, 0, 0, 0, 0},   // 114   4              //  4.339ms
    // { 24,  3,  8,  0, 0, 0, 0, 0, 0},   // 116   8              //  1.792ms*
    // { 20,  4,  5,  0, 0, 0, 0, 0, 0},   // 118   5              //  1.613ms*

    // Set default radix table
    radix_tbl[128]  = 8;    // 16
    radix_tbl[120]  = 8;    // 15
    radix_tbl[112]  = 7;    // 16
    radix_tbl[108]  = 4;    // 27
    radix_tbl[100]  = 4;    // 25
    radix_tbl[ 96]  = 2;    // 48
    radix_tbl[ 84]  = 4;    // 21
    radix_tbl[ 80]  = 2;    // 40
    radix_tbl[ 72]  = 3;    // 36
    radix_tbl[ 64]  = 4;    // 16
    radix_tbl[ 60]  = 4;    // 15
    radix_tbl[ 56]  = 7;    //  8
    radix_tbl[ 48]  = 3;    // 16
    radix_tbl[ 40]  = 8;    //  5
    radix_tbl[ 36]  = 4;    //  9
    radix_tbl[ 32]  = 2;    // 16
    radix_tbl[ 28]  = 7;    //  4
    radix_tbl[ 24]  = 3;    //  8
    radix_tbl[ 20]  = 4;    //  5
    //
    radix_tbl[27]  = 3;
    radix_tbl[25]  = 5;
    radix_tbl[21]  = 7;
    radix_tbl[15]  = 3;
    radix_tbl[ 9]  = 3;
    //
    radix_tbl[12]  = 3;
    radix_tbl[10]  = 2;
    radix_tbl[ 6]  = 2;
    //
    // Below are just for monitoring purpose.
    radix_tbl[16]  = 16;
    radix_tbl[ 8]  = 8;
    radix_tbl[ 5]  = 5;
    radix_tbl[ 4]  = 4;
    radix_tbl[ 3]  = 3;
    radix_tbl[ 2]  = 2;
}

static void fft_real_init_primitives(void)
{
    w_n4_re = cos(-M_PI/4);
    fft_data_heap_ptr = 0;
    fft_data_heap_ptr_max = 0;

    fft_n_samll_init();

        // DEBUG
        // SCIprintf("FFT_N_MAX                    = %d\n", FFT_N_MAX);
        // SCIprintf("FFT_N_STG1               = %d\n", FFT_N_STG1);
        // SCIprintf("FFT_N_STG2               = %d\n", FFT_N_STG2);
        // SCIprintf("FFT_OUT_LEN_STG1         = %d\n", FFT_OUT_LEN_STG1);
        // SCIprintf("FFT_OUT_LEN_STG2         = %d\n", FFT_OUT_LEN_STG2);
        // SCIprintf("FFT_DATA_HEAP_LEN_STRICT = %d\n", FFT_DATA_HEAP_LEN_STRICT);
        // SCIprintf("FFT_DATA_HEAP_LEN            = %d\n", FFT_DATA_HEAP_LEN);
}

//--------------------------------------------------------------------------------
// dft_init() must be done in advance for sin_cos table init.
//--------------------------------------------------------------------------------
void    fft_real_n_init(int npoint)     // npoint must be
                                        //  =k*4                ; from cos/sin table constraint
                                        //  <=128               ; from memory size
                                        //  prime factor <=7    ; from memory size at radix_odd function
{
int status;

    g_fft_init_status = 0;
    g_fft_status      = 0;
    fft_n = npoint;

    fft_real_init_primitives();

    if (fft_n<=1) {
        g_fft_init_status = -100;   // error.
        return;
    }
    if (fft_n>FFT_N_MAX) {
        g_fft_init_status = -101;   // error.
        return;
    }
    if ((fft_n%4)!=0) {
        g_fft_init_status = -300;   // sin/cos table constraint
        return;
    }
    switch(1) { 
        case 0:
            fft_real_radix_tbl_set_default();   // <=128, manually optimized.
            break;
        default:
        case 1:
            status=fft_real_radix_tbl_set_auto(npoint); // auto set
            if (status!=1) {
                g_fft_init_status = -301;       // no valid radices combination
            }
            break;
    }
}

void    fft_real_init_by_preset(int no)
{
int i;
int m;
int rdx;

    g_fft_init_status = 0;
    g_fft_status      = 0;
    fft_real_radix_tbl_clear();
    fft_real_init_primitives();

    fft_radix_preset_no = no;
    fft_n = radix_preset[no][0];
    if (fft_n<=1) {
        g_fft_init_status = -100;   // error.
        return;
    }
    if (fft_n>FFT_N_MAX) {
        g_fft_init_status = -101;   // error.
        return;
    }
    if ((fft_n%4)!=0) {
        g_fft_init_status = -300;   // error. sin/cos table constraint
        return;
    }

    // Set radix table according to the preset number, radix_preset[no]
    m = fft_n;
    for (i=1;m>1;i++) {
        rdx = radix_preset[no][i];
        if (!(rdx>1)) {
            g_fft_init_status = -101;   // error, unexpected preset parameter
            return;
        }
        if ((m % rdx)!=0) {
            g_fft_init_status = -102;   // error, not a multiple of a radix
            return;
        }
        radix_tbl[m] = (m==rdx) ? 0: rdx;     // 0: end
        m /= rdx;
    }
}

int     fft_real_get_npoint(void)
{
    return fft_n;
}

// fft_real_radix_tbl_set_auto(int npoint) will replace.
//
// OBSOLETE:
//
// assumed: radix_tbl[] has valid values.
int     fft_real_default_n_available(int npoint)
{
    if (npoint%4!=0 || npoint>FFT_N_MAX) {
        return 0;   // NG
    }
    if (radix_tbl[npoint]>1) {
        return 1;   // OK
    }
    return 0;       // NG
}

void    fft_real_radix_tbl_mon(void)
{
int i;

    SCIprintf("radix_tbl[] --------------------\n");
    for (i=FFT_N_MAX;i>=2;i--) {
        if (radix_tbl[i]>0) {
            SCIprintf(" radix_tbl[%4d] = %3d\n", i, radix_tbl[i]);
        }
    }
    SCIprintf("--------------------------------\n");
}

// Monitor radices just about to be used
void    fft_real_radix_mon(int cr_sw)   // 0: CR off, 1: CR on
{
int i, m;
int rdx;

    m=fft_n;
    rdx = radix_tbl[m];
    SCIprintf("npoint=%3d, Radix=", m);
    for (i=0;rdx>1;i++) {
        if (i==0) {
            SCIprintf("%d", rdx);
        } else {
            SCIprintf("x%d", rdx);
        }
        m /= rdx;
        rdx = radix_tbl[m];
    }
    if (cr_sw==1) {
        SCIprintf("\n");
    }
}

// Monitor radix preset table
void    fft_real_radix_preset_tbl_mon(
            int no,     // >=0: show the number, -1: show all
            int cr_sw)  // 0: CR off, 1: CR on
{
int i, j;
unsigned char c;

    for (j=0;radix_preset[j][0]>0;j++) {
        if (no<0 || j==no) {
            c = (j==fft_radix_preset_no) ? '*' : ' ';
            SCIprintf("%c[%3d] npoint=%3d, Radix=%d", c, j, radix_preset[j][0], radix_preset[j][1]);
            for (i=2; radix_preset[j][i]>0;i++) {
                SCIprintf("x%d", radix_preset[j][i]);
            }
            if (cr_sw==1) {
                SCIprintf("\n");
            }
        }
    }
}

static double *data_allocate(int n)
{
double *dp;
int i;

    dp = &(fft_data_heap[fft_data_heap_ptr]);
    if (g_fft_status!=0) return dp;     // Error
    fft_data_heap_ptr += n;

    if (fft_data_heap_ptr>fft_data_heap_ptr_max) {
        fft_data_heap_ptr_max = fft_data_heap_ptr;
    }
    if (fft_data_heap_ptr>FFT_DATA_HEAP_LEN) {      // heap overflow
            SCIprintf("Fatal error: data_allocate: fft_data_heap overflow.\n");
            SCIprintf("  requested n = %d\n", n);
            SCIprintf("  fft_data_heap_ptr = %d\n", fft_data_heap_ptr);
            SCIprintf("  ptr > FFT_DATA_HEAP_LEN(=%d)\n", FFT_DATA_HEAP_LEN);
        fft_data_heap_ptr -= n;
        for (i=0;i<FFT_DATA_HEAP_LEN;i++) {
            fft_data_heap[i] = 0.0;                 // Destroy the heap area with zero to alarm.
        }
        g_fft_status = -9;
    }
    return dp;
}

//--------------------------------------------------------------------------------
// Complex Multiply and Accumulate for FFT of Radix-n
//
//  input data: one stage buffers
//--------------------------------------------------------------------------------
static void c_mac_stage(    const   int n,          // count for MAC
                            const   int j,
                            int     *w_idx,         // twiddle factor index before % fft_n.
                                                    // content: 0: blank, {1, 2, ..., (n-1)} used.
                            double  *w_re,          // work for twiddle factors, real
                            double  *w_im,          //                           imaginary
                            double  *pa_out_re[],   // sub-FFT result real for phase=0, 1, 2, ..., (n-1)
                            double  *pa_out_im[],   // sub-FFT result real for phase=0, 1, 2, ..., (n-1)
                            double  *out_rep,
                            double  *out_imp)
{
int     i, ph;
int     wi;
double  sum_re, sum_im;

    for (i=1;i<n;i++) {
        wi = w_idx[i] % fft_n;          // to table index
        w_re[i] =       g_cos_tbl[wi];  // Twiddle factors: real
        w_im[i] = g_minus_sin_tbl[wi];  //                  imaginary
    }
    // Multiply twiddle factor and Accumulate
    //  (a+ib)(x+iy)=(ax-by)+i(ay+bx)
    sum_re = pa_out_re[0][j];
    sum_im = pa_out_im[0][j];
    for (ph=1;ph<n;ph++) {
        sum_re+=w_re[ph]*pa_out_re[ph][j] - w_im[ph]*pa_out_im[ph][j];
        sum_im+=w_re[ph]*pa_out_im[ph][j] + w_im[ph]*pa_out_re[ph][j];
    }

    *out_rep = sum_re;
    *out_imp = sum_im;
}

//--------------------------------------------------------------------------------
// Complex Multiply and Accumulate for FFT of Radix-n: Conjugate version
//  Use input data after conjugate.
//
//  input data: one stage buffers
//--------------------------------------------------------------------------------
static void c_mac_stage_conj(
                            const   int n,         // count for MAC
                            const   int j,
                            int     *w_idx,         // twiddle factor index before % fft_n.
                                                    // content: 0: blank, {1, 2, ..., (n-1)} used.
                            double  *w_re,          // work for twiddle factors, real
                            double  *w_im,          //                           imaginary
                            double  *pa_out_re[],   // sub-FFT result real for phase=0, 1, 2, ..., (n-1)
                            double  *pa_out_im[],   // sub-FFT result real for phase=0, 1, 2, ..., (n-1)
                            double  *out_rep,
                            double  *out_imp)
{
int     i, ph;
int     wi;
double  sum_re, sum_im;

    for (i=1;i<n;i++) {
        wi = w_idx[i] % fft_n;          // to table index
        w_re[i] =       g_cos_tbl[wi];  // Twiddle factors: real
        w_im[i] = g_minus_sin_tbl[wi];  //                  imaginary
    }
    // Multiply twiddle factor and Accumulate, with pa_out_im[ph][j] conjugated
    //  (a+ib)(x+iy)=(ax-by)+i(ay+bx)
    sum_re =  pa_out_re[0][j];
    sum_im = -pa_out_im[0][j];
    for (ph=1;ph<n;ph++) {
        sum_re+=  w_re[ph]*pa_out_re[ph][j] + w_im[ph]*pa_out_im[ph][j];
        sum_im+= -w_re[ph]*pa_out_im[ph][j] + w_im[ph]*pa_out_re[ph][j];
    }

    *out_rep = sum_re;
    *out_imp = sum_im;
}

//--------------------------------------------------------------------------------
// Complex Multiply and Accumulate for FFT of Radix-2
//
//  input data: one stage buffers
//--------------------------------------------------------------------------------
static void c_mac_radix2(   const   int j,
                            int     w1_idx,         // twiddle factor index before % fft_n.
                            double  *pa_out_re[],   // sub-FFT result real for phase=0, 1
                            double  *pa_out_im[],   // sub-FFT result real for phase=0, 1
                            double  *out_rep,
                            double  *out_imp)
{
int     wi;
double  w1_re,  w1_im;
double  sum_re, sum_im;

    wi   = w1_idx % fft_n;
    w1_re=       g_cos_tbl[wi];  // Twiddle factors: real
    w1_im= g_minus_sin_tbl[wi];  //                  imaginary
    // Multiply twiddle factor and Accumulate
    //  (a+ib)(x+iy)=(ax-by)+i(ay+bx)
    sum_re = pa_out_re[0][j];
    sum_im = pa_out_im[0][j];
    sum_re+=w1_re*pa_out_re[1][j] - w1_im*pa_out_im[1][j];
    sum_im+=w1_re*pa_out_im[1][j] + w1_im*pa_out_re[1][j];

    *out_rep = sum_re;
    *out_imp = sum_im;
}

//--------------------------------------------------------------------------------
// Radix-4 butterfly
//
//  input data: 4 sub-fft result
//  output:     upper 2 ports   (omit lower 2ports, since FFT input is assumed real.)
//--------------------------------------------------------------------------------
static void radix4_butterfly(
                            int     j,
                            int     flg_conj,       // Flag for taking input as conjugate
                            int     *w_idx,         // twiddle factor index before % fft_n.
                                                    // content: 0: blank, {1, 2, 3} used.
                            double  *pa_out_re[],   // sub-FFT result real      for phase=0, 1, 2, 3
                            double  *pa_out_im[],   // sub-FFT result imaginary for phase=0, 1, 2, 3
                            double  *out0_rep,
                            double  *out0_imp,
                            double  *out1_rep,
                            double  *out1_imp)
{
int     i, ph;
int     wi;
                    // w_re[0]=1.0; w_im[0]=0.0.
double  w_re[4];    // 0: dummy, {1, 2, 3}: used.
double  w_im[4];    // 0: dummy, {1, 2, 3}: used.

double  tmp_re[4];
double  tmp_im[4];

    for (i=1;i<4;i++) {
        wi = w_idx[i] % fft_n;          // to table index
        w_re[i] =       g_cos_tbl[wi];  // Twiddle factors: real
        w_im[i] = g_minus_sin_tbl[wi];  // Twiddle factors: imaginary
    }
    // Multiply twiddle factor
    tmp_re[0] = pa_out_re[0][j];
    if (flg_conj==1) {
        tmp_im[0] = -pa_out_im[0][j];
        for (ph=1;ph<4;ph++) {
            tmp_re[ph] =  w_re[ph]*pa_out_re[ph][j] + w_im[ph]*pa_out_im[ph][j];
            tmp_im[ph] = -w_re[ph]*pa_out_im[ph][j] + w_im[ph]*pa_out_re[ph][j];
        }
    } else {
        tmp_im[0] = pa_out_im[0][j];
        for (ph=1;ph<4;ph++) {
            tmp_re[ph] = w_re[ph]*pa_out_re[ph][j] - w_im[ph]*pa_out_im[ph][j];
            tmp_im[ph] = w_re[ph]*pa_out_im[ph][j] + w_im[ph]*pa_out_re[ph][j];
        }
    }
    //---
    // port 0: tmp0 +      tmp1 +      tmp2 +      tmp3
    // port 1: tmp0 + (-i)*tmp1 + (-1)*tmp2 + (+i)*tmp3     // i: imaginary unit
    // port 2: omit
    // port 3: omit
    //---
    // Output
    *out0_rep =  tmp_re[0] + tmp_re[1] + tmp_re[2] + tmp_re[3];
    *out1_rep =  tmp_re[0] + tmp_im[1] - tmp_re[2] - tmp_im[3];
    *out0_imp =  tmp_im[0] + tmp_im[1] + tmp_im[2] + tmp_im[3];
    *out1_imp =  tmp_im[0] - tmp_re[1] - tmp_im[2] + tmp_re[3];
}

//--------------------------------------------------------------------------------
// Radix-2 butterfly
//
//  input data: 2 sub-fft result
//  output:     upper 1 port    (omit the lower port, since FFT input is assumed real.)
//--------------------------------------------------------------------------------
static void radix2_butterfly(
                            int     j,
                            int     flg_conj,       // Flag for taking input as conjugate
                            int     w1_idx,         // twiddle factor index before % fft_n.
                            double  *pa_out_re[],   // sub-FFT result real      for phase=0, 1, 2, 3
                            double  *pa_out_im[],   // sub-FFT result imaginary for phase=0, 1, 2, 3
                            double  *out0_rep,
                            double  *out0_imp)
{
int     wi;
double  w1_re;
double  w1_im;

double  tmp_re[4];
double  tmp_im[4];

    wi = w1_idx % fft_n;        // to table index
    w1_re=       g_cos_tbl[wi]; // Twiddle factors: real
    w1_im= g_minus_sin_tbl[wi]; // Twiddle factors: imaginary
    // Multiply twiddle factor
    tmp_re[0] = pa_out_re[0][j];
    if (flg_conj==1) {
        tmp_im[0] = -pa_out_im[0][j];
        tmp_re[1] =  w1_re*pa_out_re[1][j] + w1_im*pa_out_im[1][j];
        tmp_im[1] = -w1_re*pa_out_im[1][j] + w1_im*pa_out_re[1][j];
    } else {
        tmp_im[0] = pa_out_im[0][j];
        tmp_re[1] = w1_re*pa_out_re[1][j] - w1_im*pa_out_im[1][j];
        tmp_im[1] = w1_re*pa_out_im[1][j] + w1_im*pa_out_re[1][j];
    }
    //---
    // port 0: tmp0 + tmp1
    // port 1: tmp0 - tmp1      // Omit
    //---
    // Output
    *out0_rep =  tmp_re[0] + tmp_re[1];
    *out0_imp =  tmp_im[0] + tmp_im[1];
}

//--------------------------------------------------------------------------------
// Radix-8 butterfly
//
//  input data: 8 sub-fft result
//  output:     upper 4 ports   (omit lower 4ports, since FFT input is assumed real.)
//--------------------------------------------------------------------------------
static void radix8_butterfly(
                            int     j,
                            int     flg_conj,       // Flag for taking input as conjugate
                            int     *w_idx,         // twiddle factor index before % fft_n.
                                                    // content: 0: blank, {1, 2, ..., 7} used.
                            double  *pa_out_re[],   // sub-FFT result real      for phase=0, 1, 2, ..., 7
                            double  *pa_out_im[],   // sub-FFT result imaginary for phase=0, 1, 2, ..., 7
                            double  *out0_rep,
                            double  *out0_imp,
                            double  *out1_rep,
                            double  *out1_imp,
                            double  *out2_rep,
                            double  *out2_imp,
                            double  *out3_rep,
                            double  *out3_imp)
{
int     i, ph;
int     wi;
                    // w_re[0]=1.0; w_im[0]=0.0.
double  w_re[8];    // 0: dummy, {1, 2, ..., 7}: used.
double  w_im[8];    // 0: dummy, {1, 2, ..., 7}: used.

double  tmp_re[8];
double  tmp_im[8];

double s_04_re, s_04_im;
double s_26_re, s_26_im;
double s_15_re, s_15_im;
double s_37_re, s_37_im;

double d_04_re, d_04_im;
double d_26_re, d_26_im;
double d_15_re, d_15_im;
double d_37_re, d_37_im;

double d_15_re_p_d_37_re;
double d_15_re_m_d_37_re;
double d_15_im_p_d_37_im;
double d_15_im_m_d_37_im;

double p1_w_re;
double p1_w_im;

double p3_w_re;
double p3_w_im;

    for (i=1;i<8;i++) {
        wi = w_idx[i] % fft_n;          // to table index
        w_re[i] =       g_cos_tbl[wi];  // Twiddle factors: real
        w_im[i] = g_minus_sin_tbl[wi];  // Twiddle factors: imaginary
    }
    // Multiply twiddle factor
    tmp_re[0] = pa_out_re[0][j];
    if (flg_conj==1) {
        tmp_im[0] = -pa_out_im[0][j];
        for (ph=1;ph<8;ph++) {
            tmp_re[ph] =  w_re[ph]*pa_out_re[ph][j] + w_im[ph]*pa_out_im[ph][j];
            tmp_im[ph] = -w_re[ph]*pa_out_im[ph][j] + w_im[ph]*pa_out_re[ph][j];
        }
    } else {
        tmp_im[0] = pa_out_im[0][j];
        for (ph=1;ph<8;ph++) {
            tmp_re[ph] = w_re[ph]*pa_out_re[ph][j] - w_im[ph]*pa_out_im[ph][j];
            tmp_im[ph] = w_re[ph]*pa_out_im[ph][j] + w_im[ph]*pa_out_re[ph][j];
        }
    }
    //---
    // port 0: (tmp0 + tmp4) +   (tmp2 + tmp6) +          (tmp1 + tmp5) +       (tmp3 + tmp7)
    // port 1: (tmp0 - tmp4) - i*(tmp2 - tmp6) + a*[(1-i)*(tmp1 - tmp5) - (1+i)*(tmp3 - tmp7)]
    // port 3: (tmp0 - tmp4) + i*(tmp2 - tmp6) - a*[(1+i)*(tmp1 - tmp5) - (1-i)*(tmp3 - tmp7)]
    // port 2: (tmp0 + tmp4) -   (tmp2 + tmp6) - i*[      (tmp1 + tmp5) -       (tmp3 + tmp7)]
    // port 4: omit
    // port 5: omit
    // port 6: omit
    // port 7: omit
    //---
    s_04_re = tmp_re[0] + tmp_re[4];
    s_04_im = tmp_im[0] + tmp_im[4];
    s_26_re = tmp_re[2] + tmp_re[6];
    s_26_im = tmp_im[2] + tmp_im[6];
    s_15_re = tmp_re[1] + tmp_re[5];
    s_15_im = tmp_im[1] + tmp_im[5];
    s_37_re = tmp_re[3] + tmp_re[7];
    s_37_im = tmp_im[3] + tmp_im[7];

    d_04_re = tmp_re[0] - tmp_re[4];
    d_04_im = tmp_im[0] - tmp_im[4];
    d_26_re = tmp_re[2] - tmp_re[6];
    d_26_im = tmp_im[2] - tmp_im[6];
    d_15_re = tmp_re[1] - tmp_re[5];
    d_15_im = tmp_im[1] - tmp_im[5];
    d_37_re = tmp_re[3] - tmp_re[7];
    d_37_im = tmp_im[3] - tmp_im[7];

    d_15_re_p_d_37_re = d_15_re + d_37_re;
    d_15_re_m_d_37_re = d_15_re - d_37_re;
    d_15_im_p_d_37_im = d_15_im + d_37_im;
    d_15_im_m_d_37_im = d_15_im - d_37_im;

    p1_w_re = w_n4_re*(d_15_re_m_d_37_re + d_15_im_p_d_37_im);
    p1_w_im = w_n4_re*(d_15_im_m_d_37_im - d_15_re_p_d_37_re);

    p3_w_re = w_n4_re*(d_15_re_m_d_37_re - d_15_im_p_d_37_im);
    p3_w_im = w_n4_re*(d_15_im_m_d_37_im + d_15_re_p_d_37_re);

    // port 0:
    *out0_rep = s_04_re + s_26_re + s_15_re + s_37_re;
    *out0_imp = s_04_im + s_26_im + s_15_im + s_37_im;
    // port 2:
    *out2_rep = s_04_re - s_26_re + s_15_im - s_37_im;
    *out2_imp = s_04_im - s_26_im - s_15_re + s_37_re;
    // port 1:
    *out1_rep = d_04_re + d_26_im + p1_w_re;
    *out1_imp = d_04_im - d_26_re + p1_w_im;
    // port 3:
    *out3_rep = d_04_re - d_26_im - p3_w_re;
    *out3_imp = d_04_im + d_26_re - p3_w_im;
}

//--------------------------------------------------------------------------------
// FFT, Radix-2
//  input:  real
//          N=even,  4<= N <=128
//  output: positive frequencies only: {0, 1, ..., fNYQ}.   // fNYQ = N/2
//          negative frequencies {(fNYQ+1), (fNYQ+2), ..., (N-1)} are omitted.
//--------------------------------------------------------------------------------
static void fft_real_radix2(    double      din_re[],   // size=npoint*interval
                                const int   npoint,
                                const int   interval,   // interval for valid data
                                double      out_re[],   // 0..fNYQ  // fNYQ=npoint/2
                                double      out_im[])   // 0..fNYQ
{
int     fft_data_heap_ptr_save;
int     ph;
int     j, j_src, j_nyq;
int     flg_conj;
int     fi_0;
int     w1_idx;
double  *pa_out_re[2];  // Pointer Array for FFT result of sub sampled data
double  *pa_out_im[2];

const   int fi_max      = npoint/2;
const   int seg_npoint  = npoint/2;
const   int interval_x2 = interval*2;

const   int sub_fft_out_len = seg_npoint/2 + 1;     // e.g. (128/2)/2+1=64/2+1=33;  (32/2)/2+1=9

    // Allocate sub fft output memory
    fft_data_heap_ptr_save = fft_data_heap_ptr;
    for (ph=0;ph<2;ph++) {
        pa_out_re[ph] = data_allocate(sub_fft_out_len);
        pa_out_im[ph] = data_allocate(sub_fft_out_len);
    }

    // FFT each (npoint/2)-point sequence
    for (ph=0;ph<2;ph++) {
        fft_real_recursive(&(din_re[ph*interval]), seg_npoint, interval_x2, &(pa_out_re[ph][0]), &(pa_out_im[ph][0]));
    }
    // Twiddle-factor-Multiply and Accumulate
    //  e.g.: npoint=128
    //      for fi = {0, 1, ..., 64}        // 64=Nyquist frequency
    //      (1) for fi = {0, 1, ..., 31}, {32, 33, ..., 63}
    fi_0  = 0;
    j_nyq = seg_npoint/2;       // e.g. 1=2/2
    for (j=0;j<seg_npoint;j++) {    // 2 = {0, 1}
        w1_idx = (fi_0*interval);   // twiddle factor index
        // Decide normal or mirror access
        if (j>j_nyq) {
            flg_conj = 1;
            j_src    = seg_npoint - j;  // Mirror access. e.g. 32-31=1.
        } else {
            flg_conj = 0;
            j_src    = j;               // Normal access.
        }
        // Radix-2 butterfly
        radix2_butterfly(j_src, flg_conj, w1_idx, pa_out_re, pa_out_im, &(out_re[fi_0]), &(out_im[fi_0]));
        fi_0++;
    }
    //---
    //  e.g.: npoint=128
    //      (2) for fi = 64, Nyquist frequency
    fi_0 = fi_max;
    j=0;
    w1_idx = (fi_0*interval);   // twiddle factor index
    // Complex Multiply and Accumulate
    c_mac_radix2(j, w1_idx, pa_out_re, pa_out_im, &(out_re[fi_0]), &(out_im[fi_0]));

    // Restore data pointer
    fft_data_heap_ptr = fft_data_heap_ptr_save;
}

//--------------------------------------------------------------------------------
// FFT, Radix-4
//  input:  real
//          N=128, 64, 32, 16, 8.
//  output: positive frequencies only: {0, 1, ..., fNYQ}.   // fNYQ = N/2
//          negative frequencies {(fNYQ+1), (fNYQ+2), ..., (N-1)} are omitted.
//--------------------------------------------------------------------------------
static void fft_real_radix4(    double      din_re[],   // size=npoint*interval
                                const int   npoint,
                                const int   interval,   // interval for valid data
                                double      out_re[],   // 0..fNYQ  // fNYQ=npoint/2
                                double      out_im[])   // 0..fNYQ
{
int     fft_data_heap_ptr_save;
int     ph;
int     i;
int     j, j_src, j_nyq;
int     flg_conj;
int     fi_0, fi_1;
int     w_idx[4];       // 0: dummy, {1, 2, 3}: used.
double  w_re[4];        // work for twiddle factors, real
double  w_im[4];        //                           imaginary
double  *pa_out_re[4];  // Pointer Array for FFT result of sub sampled data
double  *pa_out_im[4];

const   int fi_max      = npoint/2;
const   int seg_npoint  = npoint/4;
const   int interval_x4 = interval*4;

const   int sub_fft_out_len = seg_npoint/2 + 1;     // e.g. (128/4)/2+1=32/2+1=17;  (32/4)/2+1=5

    // Allocate sub fft output memory
    fft_data_heap_ptr_save = fft_data_heap_ptr;
    for (ph=0;ph<4;ph++) {
        pa_out_re[ph] = data_allocate(sub_fft_out_len);
        pa_out_im[ph] = data_allocate(sub_fft_out_len);
    }

    // FFT each (npoint/4)-point sequence
    for (ph=0;ph<4;ph++) {
        fft_real_recursive(&(din_re[ph*interval]), seg_npoint, interval_x4, &(pa_out_re[ph][0]), &(pa_out_im[ph][0]));
    }
    // Twiddle-factor-Multiply and Accumulate
    //  e.g.: npoint=128
    //      for fi = {0, 1, ..., 64}        // 64=Nyquist frequency
    //      (1) for fi = {0, 1, ..., 31}, {32, 33, ..., 63}
    fi_0  = 0;
    fi_1  = seg_npoint;
    j_nyq = seg_npoint/2;
    for (j=0;j<seg_npoint;j++) {
        for (i=1;i<4;i++) {
            w_idx[i] = (i*fi_0*interval);   // twiddle factor index
        }
        // Decide normal or mirror access
        if (j>j_nyq) {
            flg_conj = 1;
            j_src    = seg_npoint - j;  // Mirror access. e.g. 32-31=1.
        } else {
            flg_conj = 0;
            j_src    = j;               // Normal access.
        }
        // Radix-4 butterfly
        radix4_butterfly(j_src, flg_conj, w_idx, pa_out_re, pa_out_im, &(out_re[fi_0]), &(out_im[fi_0]),
                                                                       &(out_re[fi_1]), &(out_im[fi_1]));
        fi_0++;
        fi_1++;
    }
    //---
    //  e.g.: npoint=128
    //      (2) for fi = 64, Nyquist frequency
    fi_0 = fi_max;
    j=0;
    for (i=1;i<4;i++) {
        w_idx[i] = (i*fi_0*interval);   // twiddle factor index
    }
    // Complex Multiply and Accumulate
    c_mac_stage(4, j, w_idx, w_re, w_im, pa_out_re, pa_out_im, &(out_re[fi_0]), &(out_im[fi_0]));

    // Restore data pointer
    fft_data_heap_ptr = fft_data_heap_ptr_save;
}

//--------------------------------------------------------------------------------
// FFT, Radix-8
//  input:  real
//          N=8*k
//  output: positive frequencies only: {0, 1, ..., fNYQ}.   // fNYQ = N/2
//          negative frequencies {(fNYQ+1), (fNYQ+2), ..., (N-1)} are omitted.
//--------------------------------------------------------------------------------
static void fft_real_radix8(    double      din_re[],   // size=npoint*interval
                                const int   npoint,
                                const int   interval,   // interval for valid data
                                double      out_re[],   // 0..fNYQ  // fNYQ=npoint/2
                                double      out_im[])   // 0..fNYQ
{
int     fft_data_heap_ptr_save;
int     ph;
int     i;
int     j, j_src, j_nyq;
int     flg_conj;
int     fi_0, fi_1, fi_2, fi_3;
int     w_idx[8];       // 0: dummy, {1, 2, ...,7}: used.
double  w_re[8];        // work for twiddle factors, real
double  w_im[8];        //                           imaginary
double  *pa_out_re[8];  // Pointer Array for FFT result of sub sampled data
double  *pa_out_im[8];

const   int fi_max      = npoint/2;
const   int seg_npoint  = npoint/8;
const   int interval_x8 = interval*8;

const   int sub_fft_out_len = seg_npoint/2 + 1;     // e.g. (128/8)/2+1=16/2+1=9;  (64/8)/2+1=5

    // Allocate sub fft output memory
    fft_data_heap_ptr_save = fft_data_heap_ptr;
    for (ph=0;ph<8;ph++) {
        pa_out_re[ph] = data_allocate(sub_fft_out_len);
        pa_out_im[ph] = data_allocate(sub_fft_out_len);
    }

    // FFT each (npoint/8)-point sequence
    for (ph=0;ph<8;ph++) {
        fft_real_recursive(&(din_re[ph*interval]), seg_npoint, interval_x8, &(pa_out_re[ph][0]), &(pa_out_im[ph][0]));
    }
    // Twiddle-factor-Multiply and Accumulate
    //  e.g.: npoint=128
    //      for fi = {0, 1, ..., 64}        // 64=Nyquist frequency
    //      (1) for fi = {0, 1, ..., 15}, {16, 17, ..., 31}, {32, 33, ..., 47}, {48, 49, ..., 63}
    fi_0  = 0;
    fi_1  = 1*seg_npoint;
    fi_2  = 2*seg_npoint;
    fi_3  = 3*seg_npoint;
    j_nyq = seg_npoint/2;
    for (j=0;j<seg_npoint;j++) {
        for (i=1;i<8;i++) {
            w_idx[i] = (i*fi_0*interval);   // twiddle factor index
        }
        // Decide normal or mirror access
        if (j>j_nyq) {
            flg_conj = 1;
            j_src    = seg_npoint - j;  // Mirror access.
        } else {
            flg_conj = 0;
            j_src    = j;               // Normal access.
        }
        // Radix-8 butterfly
        radix8_butterfly(j_src, flg_conj, w_idx, pa_out_re, pa_out_im, &(out_re[fi_0]), &(out_im[fi_0]),
                                                                       &(out_re[fi_1]), &(out_im[fi_1]),
                                                                       &(out_re[fi_2]), &(out_im[fi_2]),
                                                                       &(out_re[fi_3]), &(out_im[fi_3]));
        fi_0++; fi_1++; fi_2++; fi_3++;
    }
    //---
    //  e.g.: npoint=128
    //      (2) for fi = 64, Nyquist frequency
    fi_0 = fi_max;
    j=0;
    for (i=1;i<8;i++) {
        w_idx[i] = (i*fi_0*interval);   // twiddle factor index
    }
    // Complex Multiply and Accumulate
    c_mac_stage(8, j, w_idx, w_re, w_im, pa_out_re, pa_out_im, &(out_re[fi_0]), &(out_im[fi_0]));

    // Restore data pointer
    fft_data_heap_ptr = fft_data_heap_ptr_save;
}

//--------------------------------------------------------------------------------
// FFT, Radix-Odd   // 3, 5 or 7
//  input:  real
//          N=Odd*k; e.g. 12(=3*4), 28=(7*4), 120(=5*24), etc.
//  output: positive frequencies only: {0, 1, ..., fNYQ}.   // fNYQ = floor(N/2)
//          negative frequencies {(fNYQ+1), (fNYQ+2), ..., (N-1)} are omitted.
//--------------------------------------------------------------------------------
static void fft_real_radix_odd( int         radix_no,   // radix number(odd): 3, 5 or 7
                                double      din_re[],   // size=npoint*interval
                                const int   npoint,
                                const int   interval,   // interval for valid data
                                double      out_re[],   // 0..fNYQ  // fNYQ=npoint/2
                                double      out_im[])   // 0..fNYQ
{
int     fft_data_heap_ptr_save;
int     ph;
int     i;
int     j, j_src, j_nyq;
int     fi;
int     seg_k;
int     w_idx[7];       // 0: dummy, {1, 2, .., (radix_no-1)}: used.
double  w_re[7];        // work for twiddle factors, real
double  w_im[7];        //                           imaginary
double  *pa_out_re[7];  // Pointer Array for FFT result of sub sampled data
double  *pa_out_im[7];

const   int fi_max      = npoint/2;
const   int seg_npoint  = npoint/radix_no;
const   int interval_xn = interval*radix_no;
const   int seg_mid     = radix_no/2;       // e.g. 3/2=1

const   int sub_fft_out_len = seg_npoint/2 + 1;     // e.g. (60/3)/2+1=20/2+1=11;

    // Allocate sub fft output memory
    fft_data_heap_ptr_save = fft_data_heap_ptr;
    for (ph=0;ph<radix_no;ph++) {
        pa_out_re[ph] = data_allocate(sub_fft_out_len);
        pa_out_im[ph] = data_allocate(sub_fft_out_len);
    }
    // FFT each segment
    for (ph=0;ph<radix_no;ph++) {
        fft_real_recursive(&(din_re[ph*interval]), seg_npoint, interval_xn, &(pa_out_re[ph][0]), &(pa_out_im[ph][0]));
    }
    // Twiddle-factor-Multiply and Accumulate
    //  e.g.: npoint=60
    //      for fi = {0, 1, ..., 30}    // 30=Nyquist frequency
    fi  = 0;
    j_nyq = seg_npoint/2;
    for (seg_k=0;seg_k<=seg_mid;seg_k++) {
        for (j=0;j<seg_npoint;j++) {
            for (i=1;i<radix_no;i++) {
                w_idx[i] = (i*fi*interval); // twiddle factor index
            }
            // Complex Multiply and Accumulate, deciding normal or mirror access
            if (j>j_nyq) {
                j_src    = seg_npoint - j;  // Mirror access.
                c_mac_stage_conj(radix_no, j_src, w_idx, w_re, w_im, pa_out_re, pa_out_im, &(out_re[fi]), &(out_im[fi]));
            } else {
                j_src    = j;               // Normal access.
                c_mac_stage     (radix_no, j_src, w_idx, w_re, w_im, pa_out_re, pa_out_im, &(out_re[fi]), &(out_im[fi]));
            }
            if (fi>=fi_max) goto FFT_REAL_RADIX_ODD_RET;    // Terminate at Nyquist frequency.
            fi++;
        }
    }
FFT_REAL_RADIX_ODD_RET:
    // Restore data pointer
    fft_data_heap_ptr = fft_data_heap_ptr_save;
}

void fft_real_mon(          double      din_re[],   // size=npoint*interval
                            const int   npoint,     // e.g. 128, 120, 32, 8, 2, etc.
                            const int   interval,   // interval for valid data
                            double      out_re[],   // range={0, 1, 2, ..., fNYQ},  fNYQ=npoint/2
                            double      out_im[])   // same above
{
int i;
int fi, fi_max;

    if (npoint!=8) return;

    SCIprintf("----------------------------------------\n");
    SCIprintf("DEBUG MON: fft_real_mon: npoint=%d, interval=%d\n", npoint, interval);
    SCIprintf("---input---\n");
    for (i=0;i<npoint;i++) {
        SCIprintf(" i=%3d: re = %12.3le\n", i, din_re[i]);
    }
    SCIprintf("---output---\n");
    fi_max = npoint/2;
    for (fi=0;fi<=fi_max;fi++) {
        SCIprintf(" fi=%2d: re = %12.3le, im = %12.3le\n", fi, out_re[fi], out_im[fi]);
    }
    SCIprintf("----------------------------------------\n");
}

void fft_real_recursive(    double      din_re[],   // size=npoint*interval
                            const int   npoint,     // e.g. 128, 120, 32, 8, 2, etc.
                            const int   interval,   // interval for valid data
                            double      out_re[],   // range={0, 1, 2, ..., fNYQ},  fNYQ=npoint/2
                            double      out_im[])   // same above
{
int     radixn;

    // Terminal stage
    switch(npoint) {
        case 16:
             fft_real_n16(&(din_re[0]), interval, &(out_re[0]), &(out_im[0]));
             return;
        case 8:
            fft_real_n8(&(din_re[0]), interval, &(out_re[0]), &(out_im[0]));
            return;
        case 5:
            fft_real_n5(&(din_re[0]), interval, &(out_re[0]), &(out_im[0]));
            return;
        case 4:
            fft_real_n4(&(din_re[0]), interval, &(out_re[0]), &(out_im[0]));
            return;
        case 3:
            fft_real_n3(&(din_re[0]), interval, &(out_re[0]), &(out_im[0]));
            return;
        case 2:
            fft_real_n2(&(din_re[0]), interval, &(out_re[0]), &(out_im[0]));
            return;
        default:
            // continue below.
            break;
    }
    if (g_fft_status!=0) {      // Error has occured.
        memset(out_re, 0, sizeof(double)*(npoint/2+1));
        memset(out_im, 0, sizeof(double)*(npoint/2+1));
        return;
    }
    // Intermedate stages
    radixn = radix_tbl[npoint];
    switch(radixn) {
        case 8:
            fft_real_radix8(&(din_re[0]), npoint, interval, &(out_re[0]), &(out_im[0]));
            return;
        case 4:
            fft_real_radix4(&(din_re[0]), npoint, interval, &(out_re[0]), &(out_im[0]));
            return;
        case 2:
            fft_real_radix2(&(din_re[0]), npoint, interval, &(out_re[0]), &(out_im[0]));
            return;
        case 3:
        case 5:
        case 7:
            fft_real_radix_odd(radixn, &(din_re[0]), npoint, interval, &(out_re[0]), &(out_im[0]));
            return;
        default:
            // radix-n undefined.
            SCIprintf("ERROR: fft_real_recursive(): radixn=%d, not supported.\n", radixn);
            SCIprintf(" npoint=%d, interval=%d\n", npoint, interval);

            g_fft_status = -2;
            return;
    }
}

//--------------------------------------------------------------------------------
// Entry: FFT(N=128, 120, etc.)
//  input:  real
//  output: positive frequencies only.          eg. N=128, {0 ..., 64}
//          negative frequencies are omitted.   e.g.       {65, ..., 127}
//--------------------------------------------------------------------------------
// dft_init(), fft_init() must be done in advance.
//--------------------------------------------------------------------------------
void    fft_real_entry( double din_re[],        // npoint
                        double out_re[],        // size = npoint/2+1
                        double out_im[])        // size = npoint/2+1
{
int     intrv;
int     fi_max;     // frequency index

    if (g_fft_init_status!=0) return;

    g_fft_status = 0;
    intrv   = 1;
    fi_max  = fft_n/2;
    fft_real_recursive(&(din_re[0]), fft_n, intrv, &(out_re[0]), &(out_im[0]));
    if (g_fft_status==0) {
        fft_mul(fi_max, (1.0)/(double)fft_n, &(out_re[0]), &(out_im[0]));
    } else {    // Error.
        memset(out_re, 0, sizeof(double)*(fi_max+1));
        memset(out_im, 0, sizeof(double)*(fi_max+1));
    }
}

void fft_address_show(void)
{
    SCIprintf("fft_address_show():--------------------\n");
    SCIprintf("FFT_N_MAX            = %d\n",    FFT_N_MAX);
    SCIprintf("FFT_DATA_HEAP_LEN    = %d\n",    FFT_DATA_HEAP_LEN);
    SCIprintf("&fft_data_heap[]     = %08Xh - %08Xh\n", &fft_data_heap[0],  &fft_data_heap[FFT_DATA_HEAP_LEN-1]);
    SCIprintf("&radix_tbl[]         = %08Xh - %08Xh\n", &(radix_tbl[0]),    &radix_tbl[RADIX_TBL_SIZE-1]);
    SCIprintf("&radix_preset[]      = %08Xh - %08Xh\n", &(radix_preset[0][0]),
                                                            ((char*)&(radix_preset[0][0]))+sizeof(radix_preset)-1);
}
