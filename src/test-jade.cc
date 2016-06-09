// Using Doxygen 1.8.0 (with Markdown)
///
/// @file   test-jade.cc
/// @author Ladutenko Konstantin <kostyfisik at gmail (.) com>
/// @copyright 2013 Ladutenko Konstantin
/// @section LICENSE
/// This file is part of JADE++.
///
/// JADE++ is free software: you can redistribute it and/or modify
/// it under the terms of the GNU General Public License as published by
/// the Free Software Foundation, either version 3 of the License, or
/// (at your option) any later version.
///
/// JADE++ is distributed in the hope that it will be useful,
/// but WITHOUT ANY WARRANTY; without even the implied warranty of
/// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
/// GNU General Public License for more details.
///
/// You should have received a copy of the GNU General Public License
/// along with JADE++.  If not, see <http://www.gnu.org/licenses/>.
///
/// @date   Wed Aug 14 13:40:37 2013
/// @brief  Test of JADE++ lib, Doxygen mainpage description.
#include <mpi.h>
#include <vector>
#include <string>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include "./jade.h"
#include <stdarg.h>  // For va_start, etc.
#include <memory>    // For std::unique_ptr

std::string string_format(const std::string fmt_str, ...) {
  // Erik Aronesty answer at
  // http://stackoverflow.com/questions/2342162/stdstring-formatting-like-sprintf
  int final_n, n = ((int)fmt_str.size()) * 2; /* Reserve two times as much as the length of the fmt_str */
  std::string str;
  std::unique_ptr<char[]> formatted;
  va_list ap;
  while(1) {
    formatted.reset(new char[n]); /* Wrap the plain char array into the unique_ptr */
    strcpy(&formatted[0], fmt_str.c_str());
    va_start(ap, fmt_str);
    final_n = vsnprintf(&formatted[0], n, fmt_str.c_str(), ap);
    va_end(ap);
    if (final_n < 0 || final_n >= n)
      n += abs(final_n - n + 1);
    else
      break;
  }
  return std::string(formatted.get());
}

//const int x100 = 100;
//debug
const int x100 = 100;
const int total_repeats = 50;
template<class T> inline T pow2(const T value) {return value*value;}
/// @brief Fitness test functions f1-f13 for benchmarks
///
/// @param x Position vector to check
///
/// @return Value to be minimized by changing x
double f1(std::vector<double> x) {
  double accumed = 0;
  for (auto x_i : x) accumed += pow2(x_i);
  return accumed;
}
double f2(std::vector<double> x) {
  double sum = 0, product = 1.0;
  for (auto x_i : x) {
    product *= std::abs(x_i);
    sum += std::abs(x_i);
  }
  return sum+product;
}
double f3(std::vector<double> x) {
  long D = x.size();
  double sum1 = 0;
  for (int i = 0; i < D; ++i) {
    double sum2 = 0;
    for (int j = 0; j < i; ++j) sum2 += x[j];
    sum1 += pow2(sum2);
  }    
  return sum1;
}
double f4(std::vector<double> x) {
  double y = std::abs(x.front());
  for (auto x_i : x) if (std::abs(x_i) > y) y = std::abs(x_i);
  return y;
}
/// %brief Rosenbrock function
double f5(std::vector<double> x) {
  long D = x.size();
  double sum = 0.0;
  for (int i = 0; i < D-1; ++i)
    sum += 100.0*pow2(x[i+1]-pow2(x[i])) + pow2(x[i] - 1.0);
  return sum;
}
/// %brief Discontinuous step function
double f6(std::vector<double> x) {
  double sum = 0;
  for (auto x_i : x) sum += pow2(floor(x_i + 0.5));
  return sum;
}
std::random_device rd;
std::mt19937_64 generator(rd());
std::uniform_real_distribution<double> distribution(0.0, 1.0);
/// %brief Noisy quartic function
double f7(std::vector<double> x) {
  long D = x.size();
  double sum = 0;
  for (int i = 0; i < D; ++i)
    sum += (i+1)*pow2(pow2(x[i]));
  return sum + distribution(generator);
}
double f8(std::vector<double> x) {
  double sum = 0;
  for (auto x_i : x) sum += -x_i * sin(sqrt(std::abs(x_i)));
  double D = static_cast<double>(x.size()); 
  return sum + D*418.98288727243369;
}
const double PI=3.14159265358979323846;
double f9(std::vector<double> x) {
  double sum = 0;
  for (auto x_i : x) sum += pow2(x_i) - 10*cos(2.0*PI*x_i) + 10.0;  
  return sum;
}
double f10(std::vector<double> x) {
  double exp1_sum = 0, exp2_sum = 0;
  for (auto x_i : x) {
    exp1_sum += pow2(x_i);
    exp2_sum += cos(2.0*PI*x_i);
  }
  double D = static_cast<double>(x.size()); 
  return -20.0 * exp(-0.2 * sqrt(exp1_sum/D))
    - exp(exp2_sum/D) + 20.0 + 2.718281828459045235;
}
double f11(std::vector<double> x) {
  double sum = 0, product = 1.0;
  long D = x.size();
  for (int i = 0; i < D; ++i) {
    sum += pow2(x[i]);
    product *= cos(x[i] / sqrt(i + 1));
  }
  return sum/4000.0 - product + 1.0;
}
double y(double x_i) {return 1.0 + (x_i + 1.0)/4.0;}
double u(double x_i, double a, double k) {
  if (x_i > a) return k*pow2(pow2(x_i - a));
  if (x_i < -a) return k*pow2(pow2(-x_i - a));
  return 0;
}
double f12(std::vector<double> x) {
  long D = x.size();
  double sum_y = 0;
  for (int i = 0; i < D-1; ++i)
    sum_y += pow2(y(x[i]) - 1.0) * (1.0 + 10.0*pow2(sin(PI * y(x[i + 1]))));
  double sum_u = 0;
  for (auto x_i : x) sum_u += u(x_i, 10.0, 100.0);
  return PI/static_cast<double>(D)
    * (10.0*pow2(sin(PI*y(x[0]))) + sum_y + pow2(y(x[D-1]) - 1.0 ))
    + sum_u;                                      ;
}
double f13(std::vector<double> x) {
  long D = x.size();
  double sum_y = 0;
  for (int i = 0; i < D-1; ++i)
    sum_y += pow2(x[i] - 1.0) * (1.0 + pow2(sin(3.0 * PI * x[i + 1])));
  double sum_u = 0;
  for (auto x_i : x) sum_u += u(x_i, 5.0, 100.0);
  return 0.1 * (pow2(sin(3.0 * PI* x[0])) + sum_y
                + pow2(x[D-1] - 1.0 )*(1.0 + pow2(sin(2.0 * PI * x[D-1]))) )
    + sum_u;
}
/// @brief Vector of pointers to test function
std::vector<double (*)(std::vector<double>)> f =
  {&f1, &f2, &f3, &f4, &f5, &f6, &f7, &f8, &f9, &f10, &f11, &f12, &f13};
std::vector<std::string> comment =
  {"f1", "f2", "f3", "f4", "f5", "f6", "f7", "f8", "f9", "f10", "f11", "f12", "f13"};
std::vector<double> lbound =
  {-100, -10, -100, -100, -30, -100, -1.28, -500, -5.12, -32, -600, -50, -50};
std::vector<double> ubound =
  { 100,  10,  100,  100,  30,  100,  1.28,  500,  5.12,  32,  600,  50,  50};
/// @brief Generations for each test function
//orig
std::vector<long> gen30D =
  {15, 20, 50, 50, 30, 1, 30, 10, 10, 5, 5, 5, 5};
std::vector<long> gen100D =
  {20, 30, 80, 150, 60, 1, 60, 10, 30, 5, 5, 5, 5};
// //mod
// std::vector<long> gen30D =
//   {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
// std::vector<long> gen100D =
//   {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};

std::vector<std::vector <double> > fitness30D, fitness100D;
std::vector<double> mean30D, mean100D;
std::vector<double> sigma30D, sigma100D;

/// @brief Run tests of JADE++.
///
/// @param argc
/// @param argv well known input parameters.
///
/// @return Zero by default.
int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  int done_status = jade::kDone;
  const int number_of_test_functions = 13;
  fitness30D.resize(number_of_test_functions);
  fitness100D.resize(number_of_test_functions);
  mean30D.resize(number_of_test_functions);
  mean100D.resize(number_of_test_functions);
  sigma30D.resize(number_of_test_functions);
  sigma100D.resize(number_of_test_functions);
   /// Settings for optimization algorithm;
  //int total_population = 3 * dimenstion;  /// Total number of individuals in population.
  for (int counts = 0; counts < total_repeats; ++counts) {
    if (fitness30D[0].size() >= total_repeats) break;
    if (rank == 0) printf("30D\n");
    int dimension = 30;  /// Number of parameters to optimize.
    int total_population = 100;  /// Total number of individuals in population.
    // int dimension = 3;  /// Number of parameters to optimize.
    // int total_population = 10;  /// Total number of individuals in population.
    for (int i = 0; i < number_of_test_functions; ++i) {
      // if (i != 1) continue;
      jade::SubPopulation sub_population;
      if (sub_population.Init(total_population, dimension) == jade::kDone) {

	// std::vector<double> feed(dimension, 0.0);
	// std::vector< std::vector<double> > feed_vector(2, feed);
	// sub_population.SetFeed(feed_vector);

        sub_population.FitnessFunction = f[i];
        /// Low and upper bound for all dimenstions;
        sub_population.SetAllBounds(lbound[i], ubound[i]);
        sub_population.SetTargetToMinimum();
        // sub_population.SetTargetToMaximum();
        sub_population.SetTotalGenerationsMax(gen30D[i]*x100);
        sub_population.SetBestShareP(0.05);
        sub_population.SetAdapitonFrequencyC(0.1);
        sub_population.SetDistributionLevel(0);
	sub_population.SwitchOffPMCRADE();
        //sub_population.PrintParameters("f1");
        sub_population.RunOptimization();
        auto current = sub_population.GetFinalFitness();
        if (sub_population.ErrorStatus()) continue;
        for (auto val : current) fitness30D[i].push_back(val);
        double sum = 0;
        for (auto x : fitness30D[i]) sum += x;
        double size = static_cast<double>(fitness30D[i].size());
        double mean = sum/size;
        double sigma = 0;
        for (auto x : fitness30D[i]) sigma += pow2(x - mean);
        sigma = sqrt(sigma/size);
        if (rank == 0)
          printf("%s\tgen%li\t%4.1e (%4.1e) runs(%g) at (%g,%g)\n",
                 comment[i].c_str(), gen30D[i]*x100, mean,sigma,size, lbound[i], ubound[i]);
	mean30D[i] = mean;
	sigma30D[i] = sigma;

      } else {
        printf("Some error!\n");
      }
    }  // end of for all test functions
    if (rank == 0) printf("100D\n");
    dimension = 100;  /// Number of parameters to optimize.
    total_population = 400;  /// Total number of individuals in population.
    // dimension = 5;  /// Number of parameters to optimize.
    // total_population = 40;  /// Total number of individuals in population.

    for (int i = 0; i < number_of_test_functions; ++i) {
      //if (i != 1) continue;
      jade::SubPopulation sub_population;
      if (sub_population.Init(total_population, dimension) == jade::kDone) {

	// std::vector<double> feed(dimension, 0.0);
	// std::vector< std::vector<double> > feed_vector(2, feed);
	// sub_population.SetFeed(feed_vector);

        sub_population.FitnessFunction = f[i];
        /// Low and upper bound for all dimenstions;
        sub_population.SetAllBounds(lbound[i], ubound[i]);
        sub_population.SetTargetToMinimum();
        // sub_population.SetTargetToMaximum();
        sub_population.SetTotalGenerationsMax(gen100D[i]*x100);
        sub_population.SetBestShareP(0.05);
        sub_population.SetAdapitonFrequencyC(0.1);
        sub_population.SetDistributionLevel(0);
	sub_population.SwitchOffPMCRADE();
        //sub_population.PrintParameters("f1");
        sub_population.RunOptimization();
        auto current = sub_population.GetFinalFitness();
        if (sub_population.ErrorStatus()) continue;
        for (auto val : current) fitness100D[i].push_back(val);
        double sum = 0;
        for (auto x : fitness100D[i]) sum += x;
        double size = static_cast<double>(fitness100D[i].size());
        double mean = sum/size;
        double sigma = 0;
        for (auto x : fitness100D[i]) sigma += pow2(x - mean);
        sigma = sqrt(sigma/size);
        if (rank == 0)
          printf("%s\tgen%li\t%4.1e (%4.1e) runs(%g) at (%g,%g)\n",
                 comment[i].c_str(), gen100D[i]*x100, mean,sigma,size, lbound[i], ubound[i]);
	mean100D[i] = mean;
	sigma100D[i] = sigma;
      } else {
        printf("Some error!\n");
      }
    }  // end of for all test functions
  }  // end of collecting runs
  // Output
  if (rank == 0) {
    FILE *fp;
    std::string fname = "test-jade.txt", out;
    fp = fopen(fname.c_str(), "w");
    double size = static_cast<double>(fitness30D[0].size());
    out = "dim 30, repeats "+string_format("%g",size)+"\nfunc, gen, mean, (sigma)\n";
    for (auto com : comment) out += string_format(" %7s  ",com.c_str());
    out.pop_back();    out += "\n";
    for (auto gen:gen30D) out += string_format(" %7li  ",gen*x100);
    out.pop_back();    out += "\n";
    for (auto mean:mean30D) out += string_format(" %4.1e  ",mean);
    out.pop_back();    out += "\n";
    for (auto sigma:sigma30D) out += string_format("(%4.1e) ",sigma);
    out.pop_back();    out += "\n";
    
    size = static_cast<double>(fitness100D[0].size());
    out += "\ndim 100, repeats "+string_format("%g",size)+"\nfunc, gen, mean, (sigma)\n";
    for (auto com : comment) out += string_format(" %7s  ",com.c_str());
    out.pop_back();    out += "\n";
    for (auto gen:gen100D) out += string_format(" %7li  ",gen*x100);
    out.pop_back();    out += "\n";
    for (auto mean:mean100D) out += string_format(" %4.1e  ",mean);
    out.pop_back();    out += "\n";
    for (auto sigma:sigma100D) out += string_format("(%4.1e) ",sigma);
    out.pop_back();    out += "\n";
    printf("%s", out.c_str());
    fprintf(fp, "%s", out.c_str());
    fclose(fp);      
  }
  MPI_Finalize();
  return done_status;
}  // end of main
// ************************************************************************* //
/// @page ChangeLog
/// ## ChangeLog
/// ### Version 0.0.1
/// - CMake files are configured for ease of use.
// ************************************************************************* //
// ************************************************************************* //
/// @mainpage JADE++ Documentation
///
/// [OpenMPI]: http://www.open-mpi.org
/// [Metamaterials Laboratory]: http://phoi.ifmo.ru/metamaterials
/// [NRU ITMO]: http://en.ifmo.ru "National Research University ITMO"
/// [CMake]: http://www.cmake.org "CMake"
/// [Ioffe Institute]: http://www.ioffe.ru/index_en.html
///                    "Ioffe Physical Technical Instute"
///
/// ### Parallel adaptive differential evolution software
///
/// JADE++ is a free (GPLv3+) high performance implementation of
/// adaptive differential evolution optimization algorithm from
/// Jingqiao Zhang and Arthur C. Sanderson book 'Adaptive Differential
/// Evolution. A Robust Approach to Multimodal Problem Optimization'
/// Springer, 2009.
///
/// JADE++ was designed to run efficiently with:
/// - multi-core processors
/// - multiprocessor systems
/// - clusters and supercomputers
///
/// JADE++ software is developed by [Metamaterials Laboratory] of
/// Photonics and Optical Informatics Department of [NRU ITMO]. It is
/// also supported by [Ioffe Institute].
///
/// Contact Ladutenko Konstantin <kostyfisik at gmail (.)  com>
/// with any questions about JADE++.
///
/// ### Prerequisites
/// - **MPI** (to use JADE++, to compile it) <br><br>
///   JADE++ needs some MPI realization to be compiled and to be
///   run. It was developed using [OpenMPI]. For Debian/Ubuntu systems
///   it can be installed with:
///
///         # apt-get install openmpi-bin openmpi-doc libopenmpi-dev
///
///   <br>
/// - **And [CMake] build system** (to compile JADE++) <br><br>
///   Use it on Linux/Unix and other operation systems.
