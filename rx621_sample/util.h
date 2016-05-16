#ifndef __UTL__
#define __UTL__

unsigned char   b4toc(unsigned int d);

void    ms_abs_count_mark(void);
int     ms_abs_count_get_diff(void);
void    ms_abs_count_show(char *s);

void    tcnt_restart(void);
double  tcnt_get_in_ms(void);
void    tcnt_show(  char *s,
                    int cr_sw); // 0: CR off, 1: CR on
#endif
