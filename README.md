# rx621_sample
Small sample C program as a HEW project for Akizuki RX621 development board
------
# Functions
## DFT/FFT and self test
On power up, DFT and FFT functions are tested and processing time is shown via serial port.  
 - FFT function spec:  
  - Input:          real only  
  - Output:         0 to Nuquist frequency  
  - Precision:      double(8 or 4byte)  
    (Renesas tool HEW can handle double type as either 8 or 4 byte.)  
  - Number of input data supported: **N=4(2^p)(3^q)(5^r)(7^s)**  
    （4<=N<=2048 and p, q, r and s is an integer >= 0）  
    (N can be extended to 4096, if double is handled as 4byte)
  - Desimation type: decimation-in-time  
  - Radix:          mixed radix  
      - intermediate stage: 8, 7, 5, 4, 3, 2  
      - terminal stage:     16, 8, 5, 4, 3, 2  
  - Implementation: recursively implemented:  
        a core function calls one of radix-specific functions depending on its input data count parameter M, and the called radix-specific function calls the core function
        to perform sub-M FFT.  

  - performance summary
    ![FFT(real) performance summary](./FFT_real_performance.png)
    For further details, see the [20160514_FFT_performance_on_RX621_RX220.pdf](20160514_FFT_performance_on_RX621_RX220.pdf) file.  

------
# License: MIT License  
- See the [LICENSE.txt](LICENSE.txt) file for details.

-----
Blog : <http://solar-club.jp/member/radioshack16/>
-----
