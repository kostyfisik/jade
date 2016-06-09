About
----

JADE++ is a high performance C++ implementation of adaptive differential
evolution optimization algorithm from Jingqiao Zhang and Arthur
C. Sanderson book 'Adaptive Differential Evolution. A Robust Approach
to Multimodal Problem Optimization' Springer, 2009.  JADE++ is
designed to run efficiently in parallel on multicore processors,
multiprocessor systems, clusters and supercomputers with help of MPI.

JADE++ needs MPI and Cmake installed to compile and run. It also needs
C++11 compatible complier.

Feel free to contact me with questions about JADE++ via e-mail
k.ladutenko@metalab.ifmo.ru!

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

value of final best fitness function found - mean value (and
stddev). Ideal value is to be zero and JADE is usually very
close to it. However, some functions (like f6 and f8) are really hard
to opimize.

``` C++
/// %brief Discontinuous step function
double f6(std::vector<double> x) {
  double sum = 0;
  for (auto x_i : x) sum += pow2(floor(x_i + 0.5));
  return sum;
}

double f8(std::vector<double> x) {
  double sum = 0;
  for (auto x_i : x) sum += -x_i * sin(sqrt(std::abs(x_i)));
  double D = static_cast<double>(x.size()); 
  return sum + D*418.98288727243369;
}
```

Test results
------------
Results from ./go.sh at [revision](
https://github.com/kostyfisik/jade/commit/27ebf553682405e8ee18bcaf66a5a835da21b112 )
See
[test-jade.cc](https://github.com/kostyfisik/jade/blob/master/src/test-jade.cc)
for more details.

```
dim 30, repeats 50
func, gen, mean, (sigma)
      f1        f2        f3        f4        f5        f6        f7        f8        f9       f10       f11       f12       f13 
    1500      2000      5000      5000      3000       100      3000      1000      1000       500       500       500       500 
 5.7e-79   5.7e-52   4.1e-93   3.8e-34   1.6e-01   4.3e+00   5.4e-04   -8.0e-13   3.3e-06   7.4e-12   3.5e-04   1.1e-22   1.0e-21 
(1.6e-78) (9.8e-52) (1.9e-92) (2.6e-33) (7.8e-01) (1.6e+00) (1.8e-04) (7.8e-12) (4.0e-06) (3.9e-12) (1.7e-03) (2.6e-22) (1.1e-21)

dim 100, repeats 50
func, gen, mean, (sigma)
      f1        f2        f3        f4        f5        f6        f7        f8        f9       f10       f11       f12       f13 
    2000      3000      8000     15000      6000       100      6000      1000      3000       500       500       500       500 
 1.2e-71   6.5e-46   3.3e-38   2.6e-61   6.4e-01   9.2e+01   8.6e-04   8.8e+03   8.0e-02   1.8e-08   2.2e-14   1.9e-03   6.0e-15 
(1.6e-71) (2.5e-45) (4.5e-38) (1.6e-60) (1.5e+00) (1.3e+01) (2.0e-04) (3.7e+02) (5.4e-02) (4.5e-09) (1.4e-14) (7.4e-03) (5.2e-15)
```


