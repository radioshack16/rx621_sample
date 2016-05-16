/***********************************************************************/
/* Utility
/***********************************************************************/
#include "util.h"
#include "sci.h"            // SCIprintf()
#include "mtu.h"
#include "globals.h"        // g_ms_abs_count

//------------------------------------------------------------
// bit 4 to hex char
//------------------------------------------------------------
unsigned char   b4toc(unsigned int d)
{
    if (d<10) {
        return '0'+d;
    } else if (d<16) {
        return 'A'+(d-10);
    } else {
        return 'F';
    }
}

//------------------------------------------------------------
// time measurement unit=1ms
//------------------------------------------------------------
static int  ms_abs_count_marked = 0;
void    ms_abs_count_mark(void)
{
    ms_abs_count_marked = g_ms_abs_count;
}

int     ms_abs_count_get_diff(void)
{
int ms_capt, ms_diff;

    ms_capt  = g_ms_abs_count;
    ms_diff = ms_capt - ms_abs_count_marked;
    return ms_diff;
}

void    ms_abs_count_show(char *s)
{
    int ms_capt, ms_diff;

    ms_capt  = g_ms_abs_count;
    ms_diff = ms_capt - ms_abs_count_marked;
    SCIprintf("%s: g_ms_abs_count = %d[ms], diff = %d[ms]\n", s, ms_capt, ms_diff);     // Time consuming
}

//------------------------------------------------------------
// time measurement by MTU ch3
//  unit=5.33us, 12.8us, etc.. See MTU2_ch3_unit_in_sec().
//------------------------------------------------------------
void    tcnt_restart(void)
{
    MTU2_ch3_restart();
}

double  tcnt_get_in_ms(void)
{
int     n;
double  t;

    n=MTU2_ch3_count();
    t = MTU2_ch3_unit_in_sec()*1000.0*(double)n;    // [ms]
    return t;
}

void    tcnt_show(  char *s,
                    int cr_sw)  // 0: CR off, 1: CR on
{
double  t;

    t = tcnt_get_in_ms();
    SCIprintf("%s: elapsed = %8.3lf[ms]", s, t);
    if (cr_sw==1) {
        SCIprintf("\n");
    }
}
