#include    <stdio.h>
#include    "sci.h"
#include    "stack_var_adrs.h"

static int* adrs_min  = (int*)0xFFFFFFFF;

void stack_var_adrs_min_update(void)
{
int dummy;

    adrs_min = (&dummy<adrs_min) ? &dummy : adrs_min;
}

int* stack_var_adrs_min_get(void)
{
    return adrs_min;
}

void stack_var_adrs_min_show(void)
{
    SCIprintf("stack_var_adrs_min_show(): adrs_min = %08Xh\n", adrs_min);
}
