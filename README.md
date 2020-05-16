# cache-contoller
The goal of this project is to determine the best architecture for a cache controller that interfaces to a 32- bit microprocessor with a 32-bit data bus.

## Overview:
While the microprocessor is general purpose in design and performs number of functions, it is desired to speed up certain signal processing functions, such as a fast Fourier transform (FFT) routine.

### Detailed Requirements:
1. We are given the task of determining the architecture of a cache memory to speed up a microprocessor system. After the program being executed was profiled, the largest hit in performance was seen in a function called Radix2FFT which is shown in the program.
2. The goal is to determine the best architecture for the 256 kB cache including the associative set size (N), number of cache lines (L), and burst length (BL), write strategy (write-back or write-through), and replacement strategy (round robin vs LRU).
3. The project seeks to allow 4GB of DRAM to be interfaced using a 32-bit data bus (arranged as 1G x 32-bits) and the cache should be limited to 256 kB in size. The size of the FFT will be 32768 points (512kB) in normal operation.
4. By adding extra code to the FFT program, it is possible to determine exactly how memory is accessed. Direct observations include: -

* How often and in what order the loop variables and data elements accessed.
* How can thrashing of cache lines be limited?
* How does the organization of cache, write strategy, and replacement strategy affect performance?
