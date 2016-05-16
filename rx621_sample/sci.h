//
// SCI functions.
//

#ifndef __SCI__
#define __SCI__

#define SCI_BUF_LEN                 80
                                                        // PCLK:
#define SCI_BRR_9600baud_16MHz      51                  // 16MHz.
#define SCI_BRR_9600baud_20MHz      64  // =65.1-1      // 20MHz.
#define SCI_BRR_9600baud_25MHz      80                  // 25MHz.
#define SCI_BRR_9600baud_48MHz     155                  // 48MHz.

#define SCI_BRR_19200baud_16MHz     25                  // 16MHz.

void    SCI_init(void);
void    SCIput(char c);                     // Tx 1 char.
void    SCIprintf(char *fmt, ...);

// debug.
int SCIprintf_len(void);

#endif
