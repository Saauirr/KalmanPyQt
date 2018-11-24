# KalmanOpenMP

Simulation of Kalman filtering of mass-spring-damper systems using OpenMP with C++ and Python using SWIG.

### Table of Contents
**[Getting Started](#getting-started)**<br>
**[Prerequisites](#prerequisites)**<br>

## Getting Started

Development in processor manufacturing technology in the last half-decade allowed one of the few modern chip manufacturers such as Intel to produce very compact in size but can deliver more computing power. Founding fathers of early computers would not even dreamed of in terms of size.

Modern day computers, even laptops which are limited in space for hardware, can be purchased with at least four CPU cores inside. Of course there are higher-end CPU's that have more than four cores: this [CPU](https://www.intel.com/content/www/us/en/products/processors/core/x-series/i9-9980xe.html) from Intel has whopping 16 cores!  
By distributing the work to two or more cores, loads can be computed in parallel to speed up the process and as a result, the time it takes to compute can be minimized.

OpenMP is an API(Application Programming Interface) for homogeneous, shared-memory system such as a multi-core processor. For the use on laptops or a single motherboard computers at home, a RAM(Random Access Memory) is an example for such shared-memory. For this project, I will be exploring into a

## Prerequisites

To be able to benefit from parallel computing using OpenMP, first you need a CPU with more than a single core. Assuming the computers nowadays are at least working with dual cores, you can check the number of CPU cores by typing `lscpu` into the terminal.
Below is the output from the machine I will be using for this project:
```console
alexk@AW15R2:~$ lscpu
Architecture:        x86_64
CPU op-mode(s):      32-bit, 64-bit
Byte Order:          Little Endian
CPU(s):              8
On-line CPU(s) list: 0-7
Thread(s) per core:  2
Core(s) per socket:  4
Socket(s):           1
NUMA node(s):        1
Vendor ID:           GenuineIntel
CPU family:          6
Model:               94
Model name:          Intel(R) Core(TM) i7-6700HQ CPU @ 2.60GHz
Stepping:            3
CPU MHz:             800.093
CPU max MHz:         3500.0000
CPU min MHz:         800.0000
BogoMIPS:            5184.00
Virtualization:      VT-x
L1d cache:           32K
L1i cache:           32K
L2 cache:            256K
L3 cache:            6144K
NUMA node0 CPU(s):   0-7
```

