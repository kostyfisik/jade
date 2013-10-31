// * Посчитать от нулевой толщины до d крит.
// * Нужна зависимость d критического от максимального eps
// * сделать для d крит разбиение на 128 слоев
// *Только ветка 1
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
// Mie model. Used in fitness function for optimization of
// sub_population.
nmie::MultiLayerMie multi_layer_mie;  
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
//double a = 1; // Krasnok PEC
double a = 0.75*lambda_work;  // 2.8125 cm
//double a = lambda_work;  // 
//double b = pi*pow2(a);
//size param = 2 pi r/wl = 2pi0.75 = 4.71
//double layer_thickness = 0.015*a;
double layer_thickness = 0.0;
double n = 4;
double k = 0;
int number_of_layers = 8;
// Production parameters
//int total_generations = 1200;
//double thickness_step = 0.02;
// Test parameters
int total_generations = 120;
double thickness_step = 0.2;
void SetTarget(double n, double k);
void SetThickness();
double SetInitialModel(double n, double k);
void SetOptimizer();
double EvaluateScatterOnlyIndex(std::vector<double> input);
double EvaluateScatter(std::vector<double> input);
std::vector< std::vector<double> > EvaluateSpectraForBestDesign();
void PrintCoating(std::vector<double> current, double initial_RCS,
                    jade::SubPopulation sub_population);
void PrintGnuPlotIndex(double initial_RCS,
                  jade::SubPopulation sub_population);
void PrintGnuPlotSpectra(std::vector< std::vector<double> > spectra,
                         double initial_RCS);
// ********************************************************************** //
// ********************************************************************** //
// ********************************************************************** //
double loss_index = 1e-11;
int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  try {
    if (isUsingPEC) {n = -1.0; k = -1.0;}
    // for (k = 0.0; k < 1.0; k = (k+0.0001)*2.0)  {
      double initial_RCS = SetInitialModel(n, k);
      for (double total_thickness = 0.4; total_thickness < 0.5;
           total_thickness += thickness_step) {
        //double total_thickness = 0.45;
        for (number_of_layers = 32; number_of_layers < 40; number_of_layers *=2) {
          layer_thickness = total_thickness / number_of_layers;
          SetOptimizer();
          sub_population.RunOptimization();
          auto current = sub_population.GetFinalFitness();
          // Output results
          int output_rank = 0;
          for (unsigned int i = 0; i < current.size(); ++i)
            if (current[output_rank] > current[i]) output_rank = i;
          if (rank == output_rank) {
            for (auto c : current) printf("All %g\n",c);
            PrintCoating(current, initial_RCS, sub_population);
          }  // end of output for process with best final fitness
          PrintGnuPlotIndex(initial_RCS, sub_population);
          PrintGnuPlotSpectra(EvaluateSpectraForBestDesign(), initial_RCS);
          sub_population.PrintResult("-- ");
        }  // end of changing number of layers
      }  // end of total coating thickness sweep
    // }  // end of k sweep
  } catch( const std::invalid_argument& ia ) {
    // Will catch if  multi_layer_mie fails or other errors.
    std::cerr << "Invalid argument: " << ia.what() << std::endl;
    MPI_Abort(MPI_COMM_WORLD, 1);  
  }  
  MPI_Finalize();
  return 0;
}
// ********************************************************************** //
// ********************************************************************** //
// ********************************************************************** //
void PrintGnuPlotIndex(double initial_RCS,
                  jade::SubPopulation sub_population) {
  gnuplot::GnuplotWrapper wrapper;
  double best_RCS = 0.0;
  auto best_x = sub_population.GetBest(&best_RCS);
  double total_coating_width = layer_thickness*number_of_layers;
  double index_sum = 0.0;
  for (auto i : best_x) index_sum+=i;
  char plot_name [300];
  snprintf(plot_name, 300,
           "TargetR%g-index-n%gk%g-CoatingW%06.3f-FinalRCS%7.4fdiff%+4.1f%%-n%lu-s%015.12f-index",
           a, n, k, total_coating_width,
           best_RCS, (best_RCS/initial_RCS-1.0)*100.0, best_x.size(), index_sum);
  wrapper.SetPlotName(plot_name);
  wrapper.SetXLabelName("Layer #");
  wrapper.SetYLabelName("Index");
  wrapper.SetDrawStyle("with histeps lw 2");
  wrapper.SetXRange({0.51, number_of_layers+0.49});
  for (int i = 0; i < number_of_layers; ++i) 
    wrapper.AddMultiPoint({i+1.0, best_x[i]});
  wrapper.AddColumnName("Layer N");
  wrapper.AddColumnName("Index");
  wrapper.MakeOutput();
}
// ********************************************************************** //
// ********************************************************************** //
// ********************************************************************** //
double SetInitialModel(double n, double k) {
  // Set common parameters for all wavelengths.
  SetTarget(n, k);
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
  sub_population.SwitchOffPMCRADE();

  sub_population.SetBestShareP(0.02);
  sub_population.SetAdapitonFrequencyC(1.0/20.0);

}
// ********************************************************************** //
// ********************************************************************** //
// ********************************************************************** //
void SetTarget(double n, double k) {
  multi_layer_mie.AddTargetLayer(a, {n, k});
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
  double k = loss_index;
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
  std::vector<complex> cindex;
  double k = loss_index;
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
std::vector< std::vector<double> > EvaluateSpectraForBestDesign() {
  double best_RCS;
  auto best_x = sub_population.GetBest(&best_RCS);
  // Setting Mie model to the best state.
  sub_population.FitnessFunction(best_x);
  return multi_layer_mie.GetSpectra(lambda_work*0.5, lambda_work*1.5, 1000);
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
void PrintGnuPlotSpectra(std::vector< std::vector<double> > spectra,
                         double initial_RCS) {
  gnuplot::GnuplotWrapper wrapper;
  double best_RCS = 0.0;
  auto best_x = sub_population.GetBest(&best_RCS);
  double total_coating_width = layer_thickness*number_of_layers;
  double index_sum = 0.0;
  for (auto i : best_x) index_sum+=i;
  char plot_name [300];
  snprintf(plot_name, 300,
           "TargetR%g-index-n%gk%g-CoatingW%06.3f-FinalRCS%07.4fdiff%+4.1f%%-n%lu-s%015.12f-spectra",
           a, n, k, total_coating_width,
           best_RCS, (best_RCS/initial_RCS-1.0)*100.0, best_x.size(), index_sum);
  wrapper.SetPlotName(plot_name);
  wrapper.SetXLabelName("WL");
  wrapper.SetYLabelName("Mie efficiency");
  wrapper.SetDrawStyle("w l lw 2");
  wrapper.SetXRange({spectra.front()[0], spectra.back()[0]});
  for (auto multi_point : spectra) wrapper.AddMultiPoint(multi_point);
  wrapper.AddColumnName("WL");
  wrapper.AddColumnName("Qext");
  wrapper.AddColumnName("Qsca");
  wrapper.AddColumnName("Qabs");
  wrapper.AddColumnName("Qbk");
  wrapper.MakeOutput();
}
// ********************************************************************** //
// ********************************************************************** //
// ********************************************************************** //
