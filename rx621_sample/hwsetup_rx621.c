//----------------------------------------
// Akizuki RX621 CPU board
// CPU chip:
// ROM 512KB
// RAM 96KB
//----------------------------------------
//
// Kazuyuki Hirooka
//
// 20160507
//
#include "iodefine.h"
#include "hwsetup_rx621.h"
#include "mtu.h"

#define     MTU0_N_PERIOD       48000               // Period = N*(1/48MHz);    48000 -->1ms
// #define     MTU0_N_PERIOD       20000               // Period = N*(1/20MHz);    20000 -->1ms
#define     MTU0_N_PERIOD_HALF  (MTU0_N_PERIOD>>1)  // *(1/2)

#define     MTU2_CH3_UNIT_SEC   (256.0/48e6)        // =5.333[us]=(1/(48MHz/256))=PRESCALE/PLCK

//======================================================================
// Port
//======================================================================
void PORT_init(void)
{
    PORT2.DR.BYTE=0xff;     // Port2/8bit: output
    PORT2.DDR.BYTE=0x00;    // Port2/8bit: all 0

                            // PORT2/bit0, CN2/25: 0/1: LED OFF/ON
}

//------------------------------------------------------------
// ??? Set port direction as input not to affect ADC performance
//------------------------------------------------------------
//  port ???
// void PORT_set_input_for_ADC(void)
// {
//     // PORT???.PDR.BYTE  = 0x00;
// }

//============================================================
// MTU2
//============================================================

//////////////////////////////////////////////////////////////
// MTU2/ch0: stop, clear, start, restart
//////////////////////////////////////////////////////////////
void MTU2_ch0_stop(void)
{
    MTUA.TSTR.BIT.CST0 = 0;     // ch0 stop
}

void MTU2_ch0_clear(void)
{
    MTU0.TCNT = 0x0000;         // ch0 clear
}

void MTU2_ch0_start(void)
{
    MTUA.TSTR.BIT.CST0 = 1;     // ch0 start
}

void MTU2_ch0_restart(void)
{
    MTU2_ch0_stop();
    MTU2_ch0_clear();
    MTU2_ch0_start();
}

//--------------------------------
// Check MTU2/ ch0 counter half passed?
//--------------------------------
int MTU2_ch0_if_half_passed(void)
{
    return (MTU0.TCNT >= MTU0_N_PERIOD_HALF);
}

//////////////////////////////////////////////////////////////
// MTU2/ch3: stop, clear, start, restart
//////////////////////////////////////////////////////////////
void MTU2_ch3_stop(void)
{
    MTUA.TSTR.BIT.CST3 = 0;     // ch3 stop
}

void MTU2_ch3_clear(void)
{
    MTU3.TCNT = 0x0000;         // ch3 clear
}

void MTU2_ch3_start(void)
{
    MTUA.TSTR.BIT.CST3 = 1;     // ch3 start
}

void MTU2_ch3_restart(void)
{
    MTU2_ch3_stop();
    MTU2_ch3_clear();
    MTU2_ch3_start();
}

int MTU2_ch3_count(void)
{
    return MTU3.TCNT;
}

double  MTU2_ch3_unit_in_sec(void)
{
    return MTU2_CH3_UNIT_SEC;
}

//============================================================
// MTU2 init
//  PCLK: 48MHz at AKI-RX621
//    (differ from 20MHz at AKI-RX220)
//============================================================
// Memo: MTU module name:
//  MTU2:   RX621, SH2-Tiny(7125)   // with no 'a' appended.
//  MTU2a:  RX220
//---
// ch0:     for 1ms-interrupt
//  input:  PCLK=48MHz xtal clock
//  output: periodic 1000Hz, T=1ms
//          TIOCOA, 1 on compare match
//---
// ch3:     for elapsed time measurement
//  input:  PCLK=48MHz xtal clock
//          Prescale=1/256 ==> T=5.33us         // NOT 12.8us at the RX220.
//  output: none
//---
void MTU2_init(void)
{
    //------------------------------------------------------------
    // Release module to run
    //------------------------------------------------------------
    MSTP_MTUA = 0;  // Release MTUA module(MTU/0-5)     // SYSTEM.MSTPCRA.BIT.MSTPA9
    MSTP_MTUB = 1;  // Stop    MTUB module(MTU/6-11)    // SYSTEM.MSTPCRA.BIT.MSTPA8

    //------------------------------------------------------------
    // Stop all channels' counter   // count start  0: stop, 1: start
    //------------------------------------------------------------
    MTUA.TSTR.BYTE  = 0x00;
    MTUB.TSTR.BYTE  = 0x00;
    //------------------------------------------------------------
    // No synchronization among channels
    //------------------------------------------------------------
    MTUA.TSYR.BYTE  = 0x00;
    MTUB.TSYR.BYTE  = 0x00;

    ////////////////////////////////////////////////////////////////////////////////
    // ch0: for 1ms-system-tick generation
    ////////////////////////////////////////////////////////////////////////////////
    //----------------------------------------
    // Timer Control Register
    //----------------------------------------
    MTU0.TCR.BIT.TPSC = 0;          // Timer count clock / Prescaler select:
                                    //   Internal clock 0:*PCLK/1, 1: PCLK/4, 2: PCLK/16, 3: PCLK/64,
                                    //   External clock 4: TCLKA,  5: TCLKB,  6: TCLKC,   7: TCLKD
    MTU0.TCR.BIT.CKEG = 0;          // Clock edge by:    0:*rising edge, 1: falling edge, 2, 3: both edges
    MTU0.TCR.BIT.CCLR = 1;          // Counter clear by: 1: TGRA compare match
    //----------------------------------------
    // Timer Mode Register
    //----------------------------------------
    MTU0.TMDR.BIT.MD  = 0;          // Mode of operation: 0: normal, 2: PWM, 4: phase measure mode1, etc.
    MTU0.TMDR.BIT.BFA = 0;          // Mode for TGRA: 0:*normal, 1: TGRC is a buffer of TGRA
    MTU0.TMDR.BIT.BFB = 0;          // Mode for TGRB: 0:*normal, 1: TGRD is a buffer of TGRB
    MTU0.TMDR.BIT.BFE = 0;          // Mode for TGRE: 0:*normal, 1: TRGF is a buffer of TGRE
    //----------------------------------------
    // Timer I/O Control Register
    //----------------------------------------
    // High
    //    MTU0.TIORH.BIT.IOA = 2;         // 2: Output from TIOC0A: 0 at initial, 1 on compare match
    // No output until it is required.
    MTU0.TIORH.BIT.IOA = 0;         // 2: Output from TIOC0A: 0 at initial, 1 on compare match
    MTU0.TIORH.BIT.IOB = 0;         // 0: Prohibit ouput from TIOC0B
    // Low
    MTU0.TIORL.BIT.IOC = 0;         // 0: Prohibit ouput from TIOC0C
    MTU0.TIORL.BIT.IOD = 0;         // 0: Prohibit ouput from TIOC0D
    //----------------------------------------
    // Timer Buffer Transfer Mode Register
    //----------------------------------------
    MTU0.TBTM.BYTE = 0;             // 0: Transfer at compare match A/B/E
    //----------------------------------------
    // Timer Counter Register
    //----------------------------------------
    MTU0.TCNT = 0x0000;
    //----------------------------------------
    // Timer General Register
    //----------------------------------------
    // $$$ SET PERIOD
    MTU0.TGRA = (MTU0_N_PERIOD-1);  // TGRA for compare match
                                    // -1 NECESSARY, since count {0, 1, 2, ..., (N-1)}
    MTU0.TGRE = MTU0.TGRA;          // the same timing as TGRA  // used for ADC start trigger
    MTU0.TGRB = 0xFFFF;             // the same as init value
    MTU0.TGRC = 0xFFFF;             // the same as init value
    MTU0.TGRD = 0xFFFF;             // the same as init value
    MTU0.TGRE = 0xFFFF;             // the same as init value
    MTU0.TGRF = 0xFFFF;             // the same as init value
    //----------------------------------------
    // Timer Interrupt Enable Register
    //----------------------------------------
    MTU0.TIER.BYTE  = 0;            // disable all MTU0-related interrupt
    MTU0.TIER2.BYTE = 0;            //
                                    // except:
    // $$$ SET INTERRUPT
    MTU0.TIER.BIT.TGIEA = 1;        // Interrupt Enable by compare-match-A: 0: disable, 1: enable
    //----------------------------------------
    // MTU0 remains stopped here.
    //----------------------------------------

    ////////////////////////////////////////////////////////////////////////////////
    // ch1-2
    ////////////////////////////////////////////////////////////////////////////////
    // Blank

    ////////////////////////////////////////////////////////////////////////////////
    // ch3: for processing time measurement
    ////////////////////////////////////////////////////////////////////////////////
    //----------------------------------------
    // Timer Control Register
    //----------------------------------------
    MTU3.TCR.BIT.TPSC = 4;          // Timer count clock / Prescaler select:
                                    //   Internal clock 0: PCLK/1,   1: PCLK/4, 2: PCLK/16, 3: PCLK/64,
                                    //                 *4: PCLK/256, 5: PCLK/1024,
                                    //   External clock 6: TCLKA,    7: TCLKB
    MTU3.TCR.BIT.CKEG = 0;          // Clock edge by:    0:*rising edge, 1: falling edge, 2, 3: both edges
    MTU3.TCR.BIT.CCLR = 0;          // Counter clear by: 0: free run
    //----------------------------------------
    // Timer Mode Register
    //----------------------------------------
    MTU3.TMDR.BIT.MD  = 0;          // Mode of operation: 0: normal, 2: PWM, 4: phase measure mode1, etc.
    MTU3.TMDR.BIT.BFA = 0;          // Mode for TGRA: 0:*normal, 1: TGRC is a buffer of TGRA
    MTU3.TMDR.BIT.BFB = 0;          // Mode for TGRB: 0:*normal, 1: TGRD is a buffer of TGRB
    //            BFE none.
    //----------------------------------------
    // Timer I/O Control Register
    //----------------------------------------
    // High
    MTU3.TIORH.BIT.IOA = 0;         // 0: Prohibit ouput from TIOC0A
    MTU3.TIORH.BIT.IOB = 0;         // 0: Prohibit ouput from TIOC0B
    // Low
    MTU3.TIORL.BIT.IOC = 0;         // 0: Prohibit ouput from TIOC0C
    MTU3.TIORL.BIT.IOD = 0;         // 0: Prohibit ouput from TIOC0D
    //----------------------------------------
    // Timer Buffer Transfer Mode Register
    //----------------------------------------
    MTU3.TBTM.BYTE = 0;             // 0: Transfer at compare match A/B/E
    //----------------------------------------
    // Timer Counter Register
    //----------------------------------------
    MTU3.TCNT = 0x0000;
    //----------------------------------------
    // Timer General Register
    //----------------------------------------
    MTU3.TGRA = 0xFFFF;             // the same as init value
    MTU3.TGRB = 0xFFFF;             // the same as init value
    MTU3.TGRC = 0xFFFF;             // the same as init value
    MTU3.TGRD = 0xFFFF;             // the same as init value
    //----------------------------------------
    // Timer Interrupt Enable Register
    //----------------------------------------
    MTU3.TIER.BYTE  = 0;            // disable all MTU3-related interrupt
    //   TIER2 none
    //----------------------------------------
    // MTU3 remains stopped here.
    //----------------------------------------

    ////////////////////////////////////////////////////////////////////////////////
    // ch4
    ////////////////////////////////////////////////////////////////////////////////
    //----------------------------------------
    // Timer AD Conversion Request Controll Register
    //----------------------------------------
    MTU4.TADCR.WORD = 0;        // Disable all.

    ////////////////////////////////////////////////////////////////////////////////
    // ch5-11
    ////////////////////////////////////////////////////////////////////////////////
    // Blank
}

/*
// ??? TO BE MODIFIED TO USE IN RX621
//
// RX220

//============================================================
// ADC
//============================================================
// Data:                12bit
// Ch count:            12
//                      4 of 16ch are not assigned for 64-pin chip.
//---
// AD Data Register:    ADDR0-ADDR15
//---
// Data Alignment:      right aligned in 16bit
// ADC start trigger:   TRG0EN(MTU0.TGRE compare match),
//                      1ms period
// Mode:                single scan
// Interrupt:           generate interrupt request S12ADI0
//                      on the end of scan
//---
// registers not used:
//  ADRD:   AD Self Diagnose Data register
//  ADOCDR: AD Internal Reference Data Register
//---
// LQFP64:
//  AN00:   pin60:  (PortXX)
//  AN01:   pin58:  (PortXX)
//  AN02:   pin57:  (PortXX)
//  AN03:   pin56:  (PortXX)
//  AN04:   pin55:  (PortXX)
//  AN05:   NA
//  AN06:   pin53:  (PortXX)
//  AN07:   NA
//  AN08:   pin51:  (PortXX)
//  AN09:   pin50:  (PortXX)
//  AN10:   pin49:  (PortXX)
//  AN11:   pin48:  (PortXX)
//  AN12:   pin47:  (PortXX)
//  AN13:   pin46:  (PortXX)
//  AN14:   NA
//  AN15:   NA
//------------------------------------------------------------
void ADC_init(void) {
    int n;

    MSTP_S12AD = 0;     // Release module stop

    //------------------------------------------------------------
    // AD Control (Scan) Register
    //------------------------------------------------------------
    S12AD.ADCSR.BIT.ADST    = 0;    // AD Start bit: 0: ADC stop, 1: ADC start
                                    // Stop while register setting.
    S12AD.ADCSR.BIT.DBLANS  = 0;    // Don't care in single scan.
    S12AD.ADCSR.BIT.GBADIE  = 0;    // Disable interrupt on the end of group B scan
    S12AD.ADCSR.BIT.DBLE    = 0;    // Not double trigger mode
    S12AD.ADCSR.BIT.EXTRG   = 0;    // Select synchronous trigger by MTU, ELC for ADC start
    S12AD.ADCSR.BIT.TRGE    = 0;    // Disable trigger for ADC start. Enable later.
    S12AD.ADCSR.BIT.ADIE    = 1;    // Enable Interrupt S12ADI0 on the end of scan
    S12AD.ADCSR.BIT.ADCS    = 0;    // Scan mode: 0:*single scan, 1: group scan, 2: continuous scan
    // S12AD.ADCSR.BIT.ADST         // AD Start bit: auto set/reset by the ADC circuit.
    //------------------------------------------------------------

    //------------------------------------------------------------
    // ADANSA: AD Analog channel select A
    //------------------------------------------------------------
    S12AD.ADANSA.WORD  = 0xffff;    //  Enable all 16ch, except:
    S12AD.ADANSA.BIT.ANSA5  = 0;    //  Disable
    S12AD.ADANSA.BIT.ANSA7  = 0;    //  Disable
    S12AD.ADANSA.BIT.ANSA14 = 0;    //  Disable
    S12AD.ADANSA.BIT.ANSA15 = 0;    //  Disable
    //------------------------------------------------------------

    //------------------------------------------------------------
    // ADANSB: AD Analog channel select B
    //------------------------------------------------------------
    S12AD.ADANSB.WORD = 0x0000;     // Not used in single scan.
    //------------------------------------------------------------

    //------------------------------------------------------------
    // ADADS: AD Add-mode channel select
    //------------------------------------------------------------
    S12AD.ADADS.WORD = 0x0000;      // Disable add-mode for all ch.
    //------------------------------------------------------------

    //------------------------------------------------------------
    // ADADC: AD Add count
    //------------------------------------------------------------
    S12AD.ADADC.BIT.ADC = 0;        // No addition.
    //------------------------------------------------------------

    //------------------------------------------------------------
    // ADCER: AD Control Extention Register
    //------------------------------------------------------------
    S12AD.ADCER.BIT.ADRFMT  = 0;    // AD Data Register Format:     right aligned
    S12AD.ADCER.BIT.ACE     = 1;    // Auto Clear Enable:           Enable
    S12AD.ADCER.BIT.DIAGVAL = 2;    // Self diagnose voltage value: (Vref x 1/2) (Not used)
    S12AD.ADCER.BIT.DIAGLD  = 1;    // Self diagnose voltage:       fixed (Not used)
    S12AD.ADCER.BIT.DIAGM   = 0;    // Self diagnose:               Disable
    //------------------------------------------------------------

    //------------------------------------------------------------
    // ADSTRGR: AD Start Trigger select Register
    //------------------------------------------------------------
    S12AD.ADSTRGR.BIT.TRSA  = 4;    // TRG0EN (MTU0/TRGE compare match)
    S12AD.ADSTRGR.BIT.TRSB  = 0;    // Group B trigger:             ADTRG0# (Not used)
    //------------------------------------------------------------

    //------------------------------------------------------------
    // ADEXICR: AD Extended Input Control Register
    //------------------------------------------------------------
    S12AD.ADEXICR.BIT.OCSAD = 0;    // Addtion mode of internal Vref:   No
    S12AD.ADEXICR.BIT.TSS   = 0;    // (Although RX220 does not have TSS bit.)
    S12AD.ADEXICR.BIT.OCS   = 0;    // AD convert internal Vref:        No  // Must be 0 for single scan mode.
    //------------------------------------------------------------

    //------------------------------------------------------------
    // ADSSTR: AD Sampling State Register
    //------------------------------------------------------------
    // Necessary:
    //  12 <= n <= 255
    //  T>0.4us  ;Sampling time T=n*(1/fclk)=n*50ns     // fclk=20MHz   ???
    //------------------------------------------------------------
    // n=255, T=12.75us     // Tentative
    //------------------------------------------------------------
    n=255;
    S12AD.ADSSTR0   = n;    // AN000
    //
    S12AD.ADSSTRL   = n;    // AN008--AN015
    S12AD.ADSSTRT   = n;    // (Although RX220 does not have Temperature sensor.)
    S12AD.ADSSTRO   = n;    // internal vref
    //
    S12AD.ADSSTR1   = n;    // AN001
    S12AD.ADSSTR2   = n;    // AN002
    S12AD.ADSSTR3   = n;    // AN003
    S12AD.ADSSTR4   = n;    // AN004
    S12AD.ADSSTR5   = n;    // AN005
    S12AD.ADSSTR6   = n;    // AN006
    S12AD.ADSSTR7   = n;    // AN007
    //------------------------------------------------------------

    //------------------------------------------------------------
    // ADDISCR: AD Disconnect Detect Control Register
    //------------------------------------------------------------
    S12AD.ADDISCR.BIT.ADNDIS = 0;   // Assist Disconnect Detect: No
    //------------------------------------------------------------

}

void ADC_trigger_enable(int trig_enable) {
    if (trig_enable==1) {
        S12AD.ADCSR.BIT.TRGE  = 1;  // Enable  trigger for ADC start
    } else {
        S12AD.ADCSR.BIT.TRGE  = 0;  // Disable trigger for ADC start
    }
}
*/
