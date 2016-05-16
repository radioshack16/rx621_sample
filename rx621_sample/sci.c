//
// SCI
// RX621
//
// 20160508
//
#include    "iodefine.h"

#include    <stdio.h>
#include    <stdarg.h>

#include    "sci.h"

static  char    vsbuff[SCI_BUF_LEN];    /* 文字列展開用バッファ(必要なら増やす) */

//============================================================
// SCI
//============================================================
// +-----+---------+-----------+
// |SCI5 | CPU-pin | Board pin |
// +-----+---------+-----------+
// |TXD  | P49     | CN2-49    |
// |RXD  | P50     | CN2-48    |
// +-----+---------+-----------+
//------------------------------------------------------------

void SCI_init(void) {
    //------------------------------------------------------------
    // Release module
    //-----------------------------------------------------------
    MSTP_SCI5 = 0;  // 0: release, 1: stop
    //------------------------------------------------------------

    //------------------------------------------------------------
    // SCI5: Serial Control Register
    //------------------------------------------------------------
    SCI5.SCR.BIT.CKE    = 0;    // Clock Enable:                        use INTERNAL(*1) baud rate generator
    SCI5.SCR.BIT.TEIE   = 0;    // Transmit End Interrupt Enable:       Disable
    SCI5.SCR.BIT.MPIE   = 0;    // Multi-Processor Interrupt Enable:    Normal
    SCI5.SCR.BIT.RE     = 0;    // Receive Enable:                      Disable
    SCI5.SCR.BIT.TE     = 0;    // Transmit Enale:                      Disable     // Enable later.
    SCI5.SCR.BIT.RIE    = 0;    // Receive Interrupt Enable:            Disable
    SCI5.SCR.BIT.TIE    = 0;    // Transmit Interrupt Enable:           Disable
    //------------------------------------------------------------

    //------------------------------------------------------------
    // SCI5: Seaial Moder Register
    //------------------------------------------------------------
    SCI5.SMR.BIT.CKS    = 0;    // Clock select:        PCLK/1 clock (n=0; BRR)
    SCI5.SMR.BIT.MP     = 0;    // Multi Processor Communication:   Prohibit
    SCI5.SMR.BIT.STOP   = 0;    // Stop bit length:     1 bit
    SCI5.SMR.BIT.PM     = 0;    // Parity mode:         Even parity
    SCI5.SMR.BIT.PE     = 0;    // Parity Enable;       Paritny none
    SCI5.SMR.BIT.CHR    = 0;    // Character length:    8 bit
    SCI5.SMR.BIT.CM     = 0;    // Communication mode:  Asynchronous
    //------------------------------------------------------------

    //------------------------------------------------------------
    // SCI5: Smart Card Mode Register
    //    SCI5.TDR = 'A';
    //------------------------------------------------------------
    SCI5.SCMR.BIT.SMIF  = 0;    // Smart Card Interface mode select:    Serial Communication
    SCI5.SCMR.BIT.SINV  = 0;    // Data Invert bit:     not invert, the same as TDR/RDR
    SCI5.SCMR.BIT.SDIR  = 0;    // Direction:           LSB first
    // SCI5.SCMR.BIT.BCP2 none
    //------------------------------------------------------------
    //------------------------------------------------------------
    // SCI5: Serial Extended Mode Register (Don't care, since (*1) precedes.)
    //------------------------------------------------------------
    SCI5.SEMR.BIT.ACS0  = 0;    // Asynchronous clock source select: External clock
    SCI5.SEMR.BIT.ABCS  = 0;    // Clock pulse count select for (asynchronouns mode) 1bit:  16
    //------------------------------------------------------------

    //------------------------------------------------------------
    // SCI5: Bit Rate Register
    //------------------------------------------------------------
    SCI5.BRR = SCI_BRR_9600baud_48MHz;
    //------------------------------------------------------------

    //------------------------------------------------------------
    // SCI5: Serial Control Register, again.
    //------------------------------------------------------------
    SCI5.SCR.BIT.TE     = 1;    // Transmit Enale:  Enable
}

//--------------------------------------------------------------------------------
void SCIput(char c)
//--------------------------------------------------------------------------------
{
    while(SCI5.SSR.BIT.TEND!=1);    // Wait until transmit ends.
    SCI5.TDR = c;
    while(SCI5.SSR.BIT.TEND!=0);    // Wait until transmit starts.
}

static  int     len = 0;
//--------------------------------------------------------------------------------
// SCIprintf
//--------------------------------------------------------------------------------
void SCIprintf(char *fmt, ...)
{
int     i;
va_list arg;

    va_start(arg, fmt);
    *vsbuff = '\0';
    vsprintf(vsbuff, fmt, arg);
    va_end(arg);

    for(i=0;i<SCI_BUF_LEN;i++) {
        if( vsbuff[i] == 0) break;
        if( vsbuff[i] == '\n') {    // if LF,
            SCIput('\r');           // insert CR before LF for MS-DOS.
            len++;
        }
        SCIput(vsbuff[i]);
        len++;
    }
}

int SCIprintf_len(void)
{
    SCIprintf("SCIprintf_len=%d\n", len);
    return  len;
}


