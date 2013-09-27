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
#include "./gnuplot-wrapper/gnuplot-wrapper.h"
#include "./nmie/ucomplex.h"
#include "./nmie/nmie-wrapper.h"
#include "./nmie/Au-dispersion.h"
const double pi=3.14159265358979323846;
template<class T> inline T pow2(const T value) {return value*value;}
nmie::MultiLayerMie multi_layer_mie;  // Mie model.
jade::SubPopulation sub_population;  // Optimizer of parameters for Mie model.
// ********************************************************************** //
// ********************************************************************** //
// ********************************************************************** //
bool isUsingPEC = true;
//bool isUsingPEC = false;
bool isOnlyIndexOptimization = true;
//bool isOnlyIndexOptimization = false;
// Semouchkina APPLIED PHYSICS LETTERS 102, 113506 (2013)
double lambda_work = 3.75; // cm
//    double f_work = 30/lambda_work; // 8 GHz
double a = 0.75*lambda_work;  // 2.8125 cm
//double a = lambda_work;  // 
//double b = pi*pow2(a);
//size param = 2 pi r/wl = 2pi0.75 = 4.71
int mul = 1;
double layer_thickness = 0.015*a/static_cast<double>(mul);
int number_of_layers = 8 * mul;
int total_generations = 20;
void SetTarget();
void SetThickness();
double SetInitialModel();
void SetOptimizer();
double EvaluateScatterOnlyIndex(std::vector<double> input);
double EvaluateScatter(std::vector<double> input);
void PrintCoating(std::vector<double> current, double initial_RCS,
                    jade::SubPopulation sub_population);
// ********************************************************************** //
// ********************************************************************** //
// ********************************************************************** //
int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  try {
    double initial_RCS = SetInitialModel();
    SetOptimizer();
    sub_population.RunOptimization();
    auto current = sub_population.GetFinalFitness();
    //Output results
    int output_rank = 0;
    for (unsigned int i = 0; i < current.size(); ++i)
      if (current[output_rank] > current[i]) output_rank = i;
    if (rank == output_rank) {
      for (auto c : current) printf("All %g\n",c);
      PrintCoating(current, initial_RCS, sub_population);
    }  // end of if first process
    sub_population.PrintResult("-- ");
  } catch( const std::invalid_argument& ia ) {
    // Will catch if  multi_layer_mie fails or other errors.
    std::cerr << "Invalid argument: " << ia.what() << std::endl;
  }  
  MPI_Finalize();
  return 0;
}
// ********************************************************************** //
// ********************************************************************** //
// ********************************************************************** //
double SetInitialModel() {
  // Set common parameters for all wavelengths.
  SetTarget();
  multi_layer_mie.SetQfaild(1000.0);  // Searching for minima
  multi_layer_mie.SetWavelength(lambda_work);
  double Qext, Qsca, Qabs, Qbk;
  multi_layer_mie.RunMie(&Qext, &Qsca, &Qabs, &Qbk);
  double total_r = multi_layer_mie.GetTotalRadius();
  double initial_RCS = Qsca*pi*pow2(total_r);
  return initial_RCS;
}
// ********************************************************************** //
// ********************************************************************** //
// ********************************************************************** //
void SetOptimizer() {
  long dimension = 0;
  if (isOnlyIndexOptimization) {
    SetThickness();
    dimension = number_of_layers;
    sub_population.FitnessFunction = &EvaluateScatterOnlyIndex;
  } else {
    dimension = number_of_layers * 2;
    sub_population.FitnessFunction = &EvaluateScatter;
  }
  long total_population = dimension * 3;
  sub_population.Init(total_population, dimension);
  /// Low and upper bound for all dimenstions;
  double from_n = 1.0, to_n = 8.0;
  sub_population.SetAllBounds(from_n, to_n);
  sub_population.SetTargetToMinimum();
  sub_population.SetTotalGenerationsMax(total_generations);
}
// ********************************************************************** //
// ********************************************************************** //
// ********************************************************************** //
void SetTarget() {
  if (isUsingPEC) {
    double target_shell_share = 0.01;
    double eps_re = 1.0;
    //double eps_im = 1250.0; // a = 0.75 lambda
    double eps_im = 900.0;
    double n = sqrt(0.5*(sqrt(pow2(eps_re) + pow2(eps_im)) + eps_re ));
    double k = sqrt(0.5*(sqrt(pow2(eps_re) + pow2(eps_im)) - eps_re ));
    multi_layer_mie.AddTargetLayer((1.0-target_shell_share)*a, {1.0, 0.0000000});
    multi_layer_mie.AddTargetLayer(target_shell_share*a, {n, k});
  } else {      
    multi_layer_mie.AddTargetLayer(a, {2.0, 0.0001});
  }
}
// ********************************************************************** //
// ********************************************************************** //
// ********************************************************************** //
void SetThickness() {
  std::vector<double> thickness;
  thickness.clear();
  if (number_of_layers < 0)
    throw std::invalid_argument("Number of coating layers should be >= 0!");
  for (int i = 0; i < number_of_layers; ++i) thickness.push_back(layer_thickness);
  multi_layer_mie.SetCoatingThickness(thickness);
}
// ********************************************************************** //
// ********************************************************************** //
// ********************************************************************** //
double EvaluateScatterOnlyIndex(std::vector<double> input) {
  double Qext, Qsca, Qabs, Qbk;
  std::vector<complex> cindex;
  cindex.clear();
  double k = 0.0;
  for (auto n : input) cindex.push_back({n, k});
  multi_layer_mie.SetCoatingIndex(cindex);
  try {
    multi_layer_mie.RunMie(&Qext, &Qsca, &Qabs, &Qbk);
  } catch( const std::invalid_argument& ia ) {
    printf(".");
    // Will catch if  multi_layer_mie fails or other errors.
    //std::cerr << "Invalid argument: " << ia.what() << std::endl;
  }  
  double total_r = multi_layer_mie.GetTotalRadius();
  return Qsca*pi*pow2(total_r);
}
// ********************************************************************** //
// ********************************************************************** //
// ********************************************************************** //
double EvaluateScatter(std::vector<double> input) {
  std::vector<double> thickness;
  thickness.clear();
  std::vector<complex> cindex;
  cindex.clear();
  double k = 0.0;
  for (int i = 0; i < number_of_layers; ++i) {
    cindex.push_back({input[i], k});
    if (input[i+number_of_layers] < 1.0) input[i+number_of_layers] = 1.0; 
    thickness.push_back(input[i+number_of_layers]*layer_thickness);
  }
  multi_layer_mie.SetCoatingIndex(cindex);
  multi_layer_mie.SetCoatingThickness(thickness);
  double Qext, Qsca, Qabs, Qbk;
  try {
    multi_layer_mie.RunMie(&Qext, &Qsca, &Qabs, &Qbk);
  } catch( const std::invalid_argument& ia ) {
    printf(".");
    // Will catch if  multi_layer_mie fails or other errors.
    //std::cerr << "Invalid argument: " << ia.what() << std::endl;
  }  
  double total_r = multi_layer_mie.GetTotalRadius();
  return Qsca*pi*pow2(total_r);
}
// ********************************************************************** //
// ********************************************************************** //
// ********************************************************************** //
void PrintCoating(std::vector<double> current, double initial_RCS,
                    jade::SubPopulation sub_population) {
  double best_RCS = 0.0;
  auto best_x = sub_population.GetBest(&best_RCS);
  printf("Target R=%g, WL=%g\n", a, lambda_work);
  printf("Initial RCS: %g\n", initial_RCS);
  printf("Final RCS: %g (%4.1f%%)\n", best_RCS, (best_RCS/initial_RCS-1.0)*100.0);
  printf ("Layer:\t");
  for (int i = 0; i < number_of_layers; ++i)
    printf("% 5i\t",i+1);
  printf ("\n");
  printf ("Index:\t");
  for (int i = 0; i < number_of_layers; ++i)
    printf("%5.4g\t",best_x[i]);
  printf("\n");
  double total_coating_width = 0.0;
  if (!isOnlyIndexOptimization) {
    printf ("Width:\t");
    for (int i = 0; i < number_of_layers; ++i) {
      double width = best_x[i+number_of_layers]*layer_thickness;
      printf("%5.4g\t",width);
      total_coating_width += width;
    }
    printf("\n");
  } else {
    total_coating_width = layer_thickness*number_of_layers;
  }
  printf("Total coating width: %g\n", total_coating_width);
}
// ********************************************************************** //
// ********************************************************************** //
// ********************************************************************** //
