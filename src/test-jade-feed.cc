// Using Doxygen 1.8.0 (with Markdown)
///
///
/// @file   test-jade-single-function.cc
/// @author Ladutenko Konstantin <kostyfisik at gmail (.) com>
/// @date   Mon Oct  7 17:23:35 2013
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
template<class T> inline T pow2(const T value) {return value*value;}
/// @brief Fitness test function
///
/// @param x Position vector to check
///
/// @return Value to be minimized by changing x
// double func_single(std::vector<double> x) {
//   // double w=x.front();
//   // double ww = w*w;
//   // double f = 1.0/2.0 *
//   //   ( 8.0 - 6.0*ww - 2.0*(3.0 - 3.0*ww)*ww/(-1.0 + ww)
//   //     )/
//   //   (-12.0 + 9.0*ww + 3.0*(3.0 - 3.0*ww)*ww/(-1.0 + ww)
//   //    -4.0*(4.0 - 4.0*ww)*ww/(pow2(-1.0 + ww))     
//   //    );
//   double w2 = x[0];
//   double wD = x[1];
//   double w22= w2*w2, wD2 = wD*wD, wD3 = wD2*wD, wD4 = wD3*wD;
//   double
//   f = 1.0/2.0* ( 2.0*(4.0*wD2 - 3.0*w22)/wD3
//                 - 2.0*(3.0*wD2 - 3.0*w22)*w22/
//                   ((-wD2 + w22)*wD3)
//                )/(
//                   -3.0*(4.0*wD2 - 3.0*w22)/wD4
//                   +3.0*(3.0*wD2 - 3.0*w22)*w22/
//                    ((-wD2 + w22)*wD4)
//                   -4.0*(4.0*wD2 - 4.0*w22)*w22/
//                    ((-wD2 + w22)*wD2)
//                  );
//   return f;
// }

/// %brief Rosenbrock function (f5)
double func_single(std::vector<double> x) {
  long D = x.size();
  double sum = 0.0;
  for (int i = 0; i < D-1; ++i)
    sum += 100.0*pow2(x[i+1]-pow2(x[i])) + pow2(x[i] - 1.0);
  return sum;
}

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  try {
    
    /// Settings for optimization algorithm;
    int dimension = 30;  /// Number of parameters to optimize.
    int total_population = 5 * dimension;  /// Total number of individuals in population.
    //int total_population = 100;  /// Total number of individuals in population.
    jade::SubPopulation sub_population;
    sub_population.Init(total_population, dimension);
    sub_population.SetFeed({{0.0,0.0}});
    sub_population.FitnessFunction = &func_single;
    /// Low and upper bound for all dimensions;
    // double lbound = 0.0;
    // double ubound = 1.1;
    //  sub_population.SetAllBoundsVectors({0, 0.1}, {1.5, 1.1});
    sub_population.SetAllBounds(-30.0, 30.0);
    
    sub_population.SetAllBoundsVectors({0, 0.1}, {1.5, 1.1});
    sub_population.SetTargetToMinimum();
    //sub_population.SetTargetToMaximum();
    sub_population.SetTotalGenerationsMax(10);
    sub_population.SetBestShareP(0.05);
    sub_population.SetAdapitonFrequencyC(0.1);
    sub_population.SetDistributionLevel(0);
    //sub_population.PrintParameters("f1");
    sub_population.RunOptimization();
    //  sub_population.PrintResult("Alena");
    sub_population.PrintResult("Rosenbrock function");
    // if (rank == 0) {
    //   //for (double x = lbound; x < ubound; x+= (ubound - lbound)/100.0) {
    //   double step = 0.0000001;
    //   for (double w2 = 0.866025; w2 < 0.86604+step; w2+= step) {
    //     for (double wD = 1.0999; wD < 1.1+step; wD+= step) {
    //       std::vector<double> input = {w2, wD};
    //       if ( func_single(input) > 1650162
    //             && std::abs(w2 - wD) > 0)
    //         printf("f(%20.7f,%20.7f)=%+5.3f\n",w2,wD, func_single(input));
    //     }
    //   }
    //   printf("\n");
    // }
  } catch( const std::invalid_argument& ia ) {
    // Will catch if  multi_layer_mie fails or other errors.
    std::cerr << "Invalid argument: " << ia.what() << std::endl;
    MPI_Abort(MPI_COMM_WORLD, 1);  
  }  
  MPI_Finalize();
  return 0;
}  // end of main
