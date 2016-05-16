//
// MTU functions.
//
#ifndef __MTU__
#define __MTU__

void    MTU2_init(void);

void    MTU2_ch0_stop(void);
void    MTU2_ch0_clear(void);
void    MTU2_ch0_start(void);
void    MTU2_ch0_restart(void);
int     MTU2_ch0_if_half_passed(void);

void    MTU2_ch3_stop(void);
void    MTU2_ch3_start(void);
void    MTU2_ch3_restart(void);
int     MTU2_ch3_count(void);
double  MTU2_ch3_unit_in_sec(void);

#endif
