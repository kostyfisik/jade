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
#include <cmath>
#include <cstdio>
#include "./jade.h"
//const int x100 = 100;
//debug
const int x100 = 100;
const int total_repeats = 5;
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
  double sum = 0, product = 0;
  for (auto x_i : x) {
    product *= abs(x_i);
    sum += abs(x_i);
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
  double y = abs(x.front());
  for (auto x_i : x) if (abs(x_i) > y) y = abs(x_i);
  return y;
}
/// %brief Rosenbrock function
double f5(std::vector<double> x) {
  long D = x.size();
  double sum = 0;
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
  for (auto x_i : x) sum += -x_i * sin(sqrt(abs(x_i)));
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
  double sum = 0, product = 0;
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
std::vector<long> gen30 =
  {15, 20, 50, 50, 30, 1, 30, 10, 10, 5, 5, 5, 5};
std::vector<std::vector <double> > fitness30, fitness100;
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
  fitness30.resize(number_of_test_functions);
  fitness100.resize(number_of_test_functions);
   /// Settings for optimization algorithm;
  int dimenstion = 30;  /// Number of parameters to optimize.
  int total_population = 100;  /// Total number of individuals in population.
  // int dimenstion = 100;  /// Number of parameters to optimize.
  // int total_population = 400;  /// Total number of individuals in population.
  //int total_population = 3 * dimenstion;  /// Total number of individuals in population.
  for (int i = 0; i < number_of_test_functions; ++i) {
    // if (i == 8) continue;
    jade::SubPopulation sub_population;
    if (sub_population.Init(total_population, dimenstion) == jade::kDone) {
      sub_population.FitnessFunction = f[i];
      /// Low and upper bound for all dimenstions;
      sub_population.SetAllBounds(lbound[i], ubound[i]);
      sub_population.SetTargetToMinimum();
      // sub_population.SetTargetToMaximum();
      sub_population.SetTotalGenerationsMax(gen30[i]*x100);
      sub_population.SetBestShareP(0.05);
      sub_population.SetAdapitonFrequencyC(0.1);
      sub_population.SetDistributionLevel(0);
      //sub_population.PrintParameters("f1");
      sub_population.RunOptimization();
      auto current = sub_population.GetFinalFitness();
      for (auto val : current) fitness30[i].push_back(val);
      double sum = 0;
      for (auto x : fitness30[i]) sum += x;
      double size = static_cast<double>(fitness30[i].size());
      double mean = sum/size;
      double sigma = 0;
      for (auto x : fitness30[i]) sigma += pow2(x - mean);
      sigma = sqrt(sigma/size);
      if (rank == 0)
        printf("%s\tgen%li\t%4.1e (%4.1e) runs(%g) at (%g,%g)\n",
               comment[i].c_str(), gen30[i]*x100, mean,sigma,size, lbound[i], ubound[i]);
    } else {
      printf("Some error!\n");
    }
  }
  MPI_Finalize();
  return done_status;
}  // end of main
// ************************************************************************* //
/// @page ChangeLog
/// #ChangeLog
/// ## Version 0.0.1
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
/// #Parallel adaptive differential evolution software
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
/// #Prerequisites
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
