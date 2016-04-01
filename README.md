About
----

JADE++ is a high performance implementation of adaptive differential
evolution optimization algorithm from Jingqiao Zhang and Arthur
C. Sanderson book 'Adaptive Differential Evolution. A Robust Approach
to Multimodal Problem Optimization' Springer, 2009.  JADE++ is
designed to run efficiently in parallel on multicore processors,
multiprocessor systems, clusters and supercomputers with help of MPI.

JADE++ needs MPI and Cmake installed to compile and run. It also needs
C++11 compatible complier.

Usage
-----

For Debian/Ubuntu systems single line install with

    # apt-get install openmpi-bin openmpi-doc libopenmpi-dev cmake

and to use LLVM Clang as a compiler

    # apt-get install clang libc++-dev

Use jade.cc and jade.h as a C++ library.

Download
-------

Checkout with the [released version](https://github.com/kostyfisik/jade/releases/tag/1.0), used in papers below!

Papers
------

The optimaizer was used to obtain results in the following papers:

1. "Reduction of scattering using thin all-dielectric shells designed by stochastic optimizer"
   Konstantin Ladutenko, Ovidio Peña-Rodríguez, Irina Melchakova, Ilya
   Yagupov, and Pavel Belov  J. Appl. Phys., vol. 116, pp. 184508,
   2014 http://dx.doi.org/10.1063/1.4900529

2. "Superabsorption of light by nanoparticles" Konstantin Ladutenko,
   Pavel Belov, Ovidio Peña-Rodríguez, Ali Mirzaei, Andrey
   E. Miroshnichenko and Ilya V. Shadrivov  Nanoscale, 2015,7,
   18897-18901 http://dx.doi.org/10.1039/C5NR05468K

Self-tests
----------

Edit go.sh to run JADE++ on your number of processes.
 
    ./go.sh single

normaly should compile JADE++ and run a single test with Rosenbrock
function (f5). On success it will finish with (almost) zero mean value of
global minima positioned at (1.0, 1.0, ..., 1.0) coordinate.
https://en.wikipedia.org/wiki/Rosenbrock_function
All individuals (candidate solutions) are shown as
evaluated.

The souce code of this test can be used as a `Hello world` example
with JADE++, you can find it in file [test-jade-single-function.cc](https://github.com/kostyfisik/jade/blob/master/src/test-jade-single-function.cc)

     ./go.sh test

to run optimization of all standard test functions (in 30D and 100D cases), will last much longer.
Example for

    ./go.sh test |grep runs

value of final best fitness function found - mean (stddev)

```
30D
f1	gen1500	5.8e-53 (4.4e-52) runs(60) at (-100,100)
f2	gen2000	6.6e-23 (4.4e-22) runs(60) at (-10,10)
f3	gen5000	2.1e-92 (7.6e-92) runs(60) at (-100,100)
f4	gen5000	7.8e-07 (2.9e-07) runs(60) at (-100,100)
f5	gen3000	6.6e-02 (5.1e-01) runs(60) at (-30,30)
f6	gen100	6.9e+00 (1.8e+00) runs(60) at (-100,100)
f7	gen3000	1.9e-02 (7.6e-03) runs(60) at (-1.28,1.28)
f8	gen1000	7.6e+01 (5.7e+02) runs(60) at (-500,500)
f9	gen1000	2.3e-04 (9.2e-05) runs(60) at (-5.12,5.12)
f10	gen500	4.3e-09 (3.2e-09) runs(60) at (-32,32)
f11	gen500	1.5e-07 (8.1e-07) runs(60) at (-600,600)
f12	gen500	1.1e-15 (5.4e-15) runs(60) at (-50,50)
f13	gen500	2.5e-15 (7.2e-15) runs(60) at (-50,50)
100D
f1	gen2000	6.0e-67 (1.8e-66) runs(60) at (-100,100)
f2	gen3000	1.5e-50 (8.3e-50) runs(60) at (-10,10)
f3	gen8000	4.4e-38 (6.3e-38) runs(60) at (-100,100)
f4	gen15000	2.4e-02 (4.7e-03) runs(60) at (-100,100)
f5	gen6000	4.0e-01 (1.2e+00) runs(60) at (-30,30)
f6	gen100	1.4e+02 (1.7e+01) runs(60) at (-100,100)
f7	gen6000	9.8e-03 (3.4e-03) runs(60) at (-1.28,1.28)
f8	gen1000	8.7e+03 (3.4e+02) runs(60) at (-500,500)
f9	gen3000	3.0e-01 (6.4e-02) runs(60) at (-5.12,5.12)
f10	gen500	4.8e-07 (1.3e-07) runs(60) at (-32,32)
f11	gen500	1.2e-04 (9.5e-04) runs(60) at (-600,600)
f12	gen500	1.8e-13 (1.6e-13) runs(60) at (-50,50)
f13	gen500	1.2e-11 (2.9e-11) runs(60) at (-50,50)

with PMCRADE feature on
30D
f1	gen1500	9.4e-79 (3.2e-78) runs(60) at (-100,100)
f2	gen2000	3.5e-52 (4.5e-52) runs(60) at (-10,10)
f3	gen5000	1.1e-90 (8.4e-90) runs(60) at (-100,100)
f4	gen5000	1.2e-38 (6.9e-38) runs(60) at (-100,100)
f5	gen3000	2.7e-01 (9.9e-01) runs(60) at (-30,30)
f6	gen100	4.0e+00 (1.5e+00) runs(60) at (-100,100)
f7	gen3000	5.7e-04 (2.0e-04) runs(60) at (-1.28,1.28)
f8	gen1000	2.0e+00 (1.5e+01) runs(60) at (-500,500)
f9	gen1000	7.8e-06 (4.0e-05) runs(60) at (-5.12,5.12)
f10	gen500	7.5e-12 (3.4e-12) runs(60) at (-32,32)
f11	gen500	0.0e+00 (0.0e+00) runs(60) at (-600,600)
f12	gen500	1.2e-22 (2.5e-22) runs(60) at (-50,50)
f13	gen500	1.2e-21 (1.7e-21) runs(60) at (-50,50)
100D
f1	gen2000	9.9e-72 (9.7e-72) runs(60) at (-100,100)
f2	gen3000	1.2e-45 (3.6e-45) runs(60) at (-10,10)
f3	gen8000	2.0e-38 (2.8e-38) runs(60) at (-100,100)
f4	gen15000	6.9e-62 (4.7e-61) runs(60) at (-100,100)
f5	gen6000	3.3e-01 (1.1e+00) runs(60) at (-30,30)
f6	gen100	9.2e+01 (1.4e+01) runs(60) at (-100,100)
f7	gen6000	8.6e-04 (1.8e-04) runs(60) at (-1.28,1.28)
f8	gen1000	8.8e+03 (3.4e+02) runs(60) at (-500,500)
f9	gen3000	9.4e-02 (6.0e-02) runs(60) at (-5.12,5.12)
f10	gen500	1.9e-08 (6.6e-09) runs(60) at (-32,32)
f11	gen500	6.6e-04 (2.2e-03) runs(60) at (-600,600)
f12	gen500	1.7e-16 (1.6e-16) runs(60) at (-50,50)
f13	gen500	9.5e-15 (1.7e-14) runs(60) at (-50,50)

```


