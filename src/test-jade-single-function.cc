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
double func_single(std::vector<double> x) {
  double w=x.front();
  double ww = w*w;
  double f = 1.0/2.0 *
    ( 8.0 - 6.0*ww - 2.0*(3.0 - 3.0*ww)*ww/(-1.0 + ww)
      )/
    (-12.0 + 9.0*ww + 3.0*(3.0 - 3.0*ww)*ww/(-1.0 + ww)
     -4.0*(4.0 - 4.0*ww)*ww/(pow2(-1.0 + ww))     
     );
  return f;
}

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   /// Settings for optimization algorithm;
  int dimenstion = 1;  /// Number of parameters to optimize.
  int total_population = 3 * dimenstion;  /// Total number of individuals in population.
  //int total_population = 100;  /// Total number of individuals in population.
  jade::SubPopulation sub_population;
  sub_population.Init(total_population, dimenstion);
  sub_population.FitnessFunction = &func_single;
  /// Low and upper bound for all dimenstions;
  double lbound = 0.0;
  double ubound = 2.0/std::sqrt(3.0);
  sub_population.SetAllBounds(lbound, ubound);
  //sub_population.SetTargetToMinimum();
  sub_population.SetTargetToMaximum();
  sub_population.SetTotalGenerationsMax(30);
  sub_population.SetBestShareP(0.05);
  sub_population.SetAdapitonFrequencyC(0.1);
  sub_population.SetDistributionLevel(0);
  //sub_population.PrintParameters("f1");
  sub_population.RunOptimization();
  sub_population.PrintResult("Alena");
  if (rank == 0) {
    //for (double x = lbound; x < ubound; x+= (ubound - lbound)/100.0) {
    for (double x = lbound; x < ubound; x+= 0.01) {
      std::vector<double> input = {x};
      printf("f(%g)=%g\n",x, func_single(input));
    }
    std::vector<double> input = {ubound};
    printf("final f(ubound)=%g\n", func_single(input));
  }
  MPI_Finalize();
  return 0;
}  // end of main
