///
/// @file   optimize-cloak.cc
/// @author Ladutenko Konstantin <kostyfisik at gmail (.) com>
/// @date   Tue Sep  3 00:37:05 2013
/// @copyright 2013 Ladutenko Konstantin
///
/// optimize-cloak is free software: you can redistribute it and/or modify
/// it under the terms of the GNU General Public License as published by
/// the Free Software Foundation, either version 3 of the License, or
/// (at your option) any later version.
///
/// optimize-cloak is distributed in the hope that it will be useful,
/// but WITHOUT ANY WARRANTY; without even the implied warranty of
/// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
/// GNU General Public License for more details.
///
/// You should have received a copy of the GNU General Public License
/// along with optimize-cloak.  If not, see <http://www.gnu.org/licenses/>.
///
/// optimize-cloak uses nmie.c from scattnlay by Ovidio Pena
/// <ovidio@bytesfall.com> as a linked library. He has an additional condition to 
/// his library:
//    The only additional condition is that we expect that all publications         //
//    describing  work using this software , or all commercial products             //
//    using it, cite the following reference:                                       //
//    [1] O. Pena and U. Pal, "Scattering of electromagnetic radiation by           //
//        a multilayered sphere," Computer Physics Communications,                  //
//        vol. 180, Nov. 2009, pp. 2348-2354.                                       //
///
/// @brief Simulate scattering from dielectric sphere covered with gold/dielectric
/// double shell using scattnlay lib (or gold shell inside dielectric ball)
/// 
#include <mpi.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <cstdio>
#include "./jade.h"
#include <cstdlib>
#include <cstdio>
#include <ctime>
#include <stdexcept>
#include <string>
#include "./nmie/ucomplex.h"
#include "./nmie/nmie-wrapper.h"
#include "./nmie/Au-dispersion.h"
const double pi=3.14159265358979323846;
template<class T> inline T pow2(const T value) {return value*value;}
nmie::MultiLayerMie multi_layer_mie;
bool isUsingPEC = true;
//bool isUsingPEC = false;
// Semouchkina APPLIED PHYSICS LETTERS 102, 113506 (2013)
double lambda_work = 3.75; // cm
//    double f_work = 30/lambda_work; // 8 GHz
double a = 0.75*lambda_work;  // 2.8125 cm
double b = pi*pow2(a);
//size param = 2 pi r/wl = 2pi0.75 = 4.71
double layer_thickness = 0.015*a;
int number_of_layers = 5;
void SetTarget();
void SetThickness();

double EvaluateScatter(std::vector<double> index) {
  double Qext, Qsca, Qabs, Qbk;
  std::vector<complex> cindex;
  cindex.clear();
  double k = 0.0;
  for (auto n : index) cindex.push_back({n, k});
  multi_layer_mie.SetCoatingIndex(cindex);
  multi_layer_mie.RunMie(&Qext, &Qsca, &Qabs, &Qbk);
  double total_r = multi_layer_mie.GetTotalRadius();
  return Qsca*pi*pow2(total_r);
}
int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  try {
    // Set common parameters for all wavelengths.
    SetTarget();
    multi_layer_mie.SetQfaild(1000.0);  // Searching for minima
    multi_layer_mie.SetWavelength(lambda_work);
    double Qext, Qsca, Qabs, Qbk;
    multi_layer_mie.RunMie(&Qext, &Qsca, &Qabs, &Qbk);
    double total_r = multi_layer_mie.GetTotalRadius();
    if (rank == 0) printf("Initial RCS: %g\n", Qsca*pi*pow2(total_r));

    SetThickness();
    // Set optimizer
    jade::SubPopulation sub_population;
    double total_population = number_of_layers * 5;
    double dimension = number_of_layers;
    sub_population.Init(total_population, dimension);
    sub_population.FitnessFunction = &EvaluateScatter;
    /// Low and upper bound for all dimenstions;
    double from_n = 1.0, to_n = 5.0;
    sub_population.SetAllBounds(from_n, to_n);
    sub_population.SetTargetToMinimum();
    // sub_population.SetTargetToMaximum();
    sub_population.SetTotalGenerationsMax(100);
        //sub_population.PrintParameters("f1");
    sub_population.RunOptimization();
    auto current = sub_population.GetFinalFitness();
    if (rank == 0) {
      for (auto c : current) printf("%g ",c);
      printf("\n");
    }  // end of if first process
  } catch( const std::invalid_argument& ia ) {
    // Will catch if  multi_layer_mie fails or other errors.
    std::cerr << "Invalid argument: " << ia.what() << std::endl;
  }  
  MPI_Finalize();
  return 0;
}
void SetTarget() {
  if (isUsingPEC) {
    double target_shell_share = 0.01;
    double eps_re = 1.0;
    double eps_im = 1250.0;
    double n = sqrt(0.5*(sqrt(pow2(eps_re) + pow2(eps_im)) + eps_re ));
    double k = sqrt(0.5*(sqrt(pow2(eps_re) + pow2(eps_im)) - eps_re ));
    multi_layer_mie.AddTargetLayer((1.0-target_shell_share)*a, {1.0, 0.0000000});
    multi_layer_mie.AddTargetLayer(target_shell_share*a, {n, k});
  } else {      
    multi_layer_mie.AddTargetLayer(a, {2.0, 0.0001});
  }
}
void SetThickness() {
  std::vector<double> thickness;
  thickness.clear();
  if (number_of_layers < 0)
    throw std::invalid_argument("Number of coating layers should be >= 0!");
  for (int i = 0; i < number_of_layers; ++i) thickness.push_back(layer_thickness);
  multi_layer_mie.SetCoatingThickness(thickness);
}
