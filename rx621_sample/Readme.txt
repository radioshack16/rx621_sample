-------- PROJECT GENERATOR --------
PROJECT NAME :	rx621_sample
PROJECT DIRECTORY :	C:\hewwork\rx621_sample\rx621_sample
CPU SERIES :	RX600
CPU TYPE :	RX62N
TOOLCHAIN NAME :	Renesas RX Standard Toolchain
TOOLCHAIN VERSION :	1.2.1.0
GENERATION FILES :
    C:\hewwork\rx621_sample\rx621_sample\dbsct.c
        Setting of B,R Section
    C:\hewwork\rx621_sample\rx621_sample\typedefine.h
        Aliases of Integer Type
    C:\hewwork\rx621_sample\rx621_sample\lowlvl.src
        Program of Low level
    C:\hewwork\rx621_sample\rx621_sample\lowsrc.c
        Program of I/O Stream
    C:\hewwork\rx621_sample\rx621_sample\sbrk.c
        Program of sbrk
    C:\hewwork\rx621_sample\rx621_sample\iodefine.h
        Definition of I/O Register
    C:\hewwork\rx621_sample\rx621_sample\intprg.c
        Interrupt Program
    C:\hewwork\rx621_sample\rx621_sample\vecttbl.c
        Initialize of Vector Table
    C:\hewwork\rx621_sample\rx621_sample\vect.h
        Definition of Vector
    C:\hewwork\rx621_sample\rx621_sample\resetprg.c
        Reset Program
    C:\hewwork\rx621_sample\rx621_sample\hwsetup.c
        Hardware Setup file
    C:\hewwork\rx621_sample\rx621_sample\rx621_sample.c
        Main Program
    C:\hewwork\rx621_sample\rx621_sample\lowsrc.h
        Header file of I/O Stream file
    C:\hewwork\rx621_sample\rx621_sample\sbrk.h
        Header file of sbrk file
    C:\hewwork\rx621_sample\rx621_sample\stacksct.h
        Setting of Stack area
START ADDRESS OF SECTION :
 H'1000	B_1,R_1,B_2,R_2,B,R,SU,SI
 H'FFFF8000	PResetPRG
 H'FFFF8100	C_1,C_2,C,C$*,D_1,D_2,D,P,PIntPRG,W*,L
 H'FFFFFFD0	FIXEDVECT

* When the user program is executed,
* the interrupt mask has been masked.
* 
* Program start 0xFFFF8000.
* RAM start 0x1000.

DATE & TIME : 2016/05/07 0:22:49
