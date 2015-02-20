/// @file   optimize-feed-cloak.cc
/// @author Ladutenko Konstantin <kostyfisik at gmail (.) com>
/// @date   Tue Sep  3 00:37:05 2013
/// @copyright 2013 Ladutenko Konstantin
///
/// optimize-feed-cloak is free software: you can redistribute it and/or modify
/// it under the terms of the GNU General Public License as published by
/// the Free Software Foundation, either version 3 of the License, or
/// (at your option) any later version.
///
/// optimize-feed-cloak is distributed in the hope that it will be useful,
/// but WITHOUT ANY WARRANTY; without even the implied warranty of
/// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
/// GNU General Public License for more details.
///
/// You should have received a copy of the GNU General Public License
/// along with optimize-feed-cloak.  If not, see <http://www.gnu.org/licenses/>.
///
/// optimize-feed-cloak uses nmie.c from scattnlay by Ovidio Pena
/// <ovidio@bytesfall.com> as a linked library. He has an additional condition to 
/// his library:
//    The only additional condition is that we expect that all publications         //
//    describing  work using this software , or all commercial products             //
//    using it, cite the following reference:                                       //
//    [1] O. Pena and U. Pal, "Scattering of electromagnetic radiation by           //
//        a multilayered sphere," Computer Physics Communications,                  //
//        vol. 180, Nov. 2009, pp. 2348-2354.                                       //
///
/// @brief Simulate scattering from dielectric sphere covered with
///  shell using scattnlay lib 
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
//#include "./nmie/ucomplex.h"
#include "./nmie/nmie-wrapper.h"
//#include "./nmie/Au-dispersion.h"
const double pi=3.14159265358979323846;
template<class T> inline T pow2(const T value) {return value*value;}

void SetTarget(double n, double k);
void SetThickness();
double SetInitialModel();
void SetOptimizer();
double EvaluateScatterOnlyIndex(std::vector<double> input);
std::vector< std::vector<double> > EvaluateSpectraForBestDesign();
void PrintCoating(std::vector<double> current, double initial_RCS,
                    jade::SubPopulation sub_population);
void PrintGnuPlotIndex(double initial_RCS,
                  jade::SubPopulation sub_population);
void PrintGnuPlotSpectra(std::vector< std::vector<double> > spectra,
                         double initial_RCS);
void Output();
// ********************************************************************** //
// ********************************************************************** //
// ********************************************************************** //
// Mie model. Used in fitness function for optimization of
// sub_population_.
nmie::MultiLayerMie multi_layer_mie_;  
jade::SubPopulation sub_population_;  // Optimizer of parameters for Mie model.

double lambda_work = 3.75; // cm
double a = 0.75*lambda_work;  // 2.8125 cm - size of PEC core
int total_generations = 20;
double layer_thickness = 0.2;
int max_number_of_layers_ = 8;

double min_index = 1e-11;
int number_of_layers_ = 0;
double initial_RCS_ = 0.0;
// ********************************************************************** //
// ********************************************************************** //
// ********************************************************************** //
int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);
  try {
    initial_RCS_ = SetInitialModel();
    for (number_of_layers_ = 1; number_of_layers_ < max_number_of_layers_; ++number_of_layers_) {
      SetOptimizer();
      sub_population_.RunOptimization();
      Output();         // Output results
    }  // end of changing number of layers
  } catch( const std::invalid_argument& ia ) {
    // Will catch if  multi_layer_mie_ fails or other errors.
    std::cerr << "Invalid argument: " << ia.what() << std::endl;
    MPI_Abort(MPI_COMM_WORLD, 1);  
  }  
  MPI_Finalize();
  return 0;
}
// ********************************************************************** //
// ********************************************************************** //
// ********************************************************************** //
void Output() {
        auto current = sub_population_.GetFinalFitness();
        int output_rank = 0;
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        for (unsigned int i = 0; i < current.size(); ++i)
          if (current[output_rank] > current[i]) output_rank = i;
        if (rank == output_rank) {
          for (auto c : current) printf("All %g\n",c);
          PrintCoating(current, initial_RCS_, sub_population_);
        }  // end of output for process with best final fitness
        PrintGnuPlotIndex(initial_RCS_, sub_population_);
        PrintGnuPlotSpectra(EvaluateSpectraForBestDesign(), initial_RCS_);
        sub_population_.PrintResult("-- ");
} 
// ********************************************************************** //
// ********************************************************************** //
// ********************************************************************** //
void PrintGnuPlotIndex(double initial_RCS,
                  jade::SubPopulation sub_population) {
  gnuplot::GnuplotWrapper wrapper;
  double best_RCS = 0.0;
  auto best_x = sub_population.GetBest(&best_RCS);
  double total_coating_width = layer_thickness*number_of_layers_;
  double index_sum = 0.0;
  for (auto i : best_x) index_sum+=i;
  char plot_name [300];
  snprintf(plot_name, 300,
           "TargetR%6.4f-CoatingW%06.3f-FinalRCS%7.4fdiff%+4.1f%%-n%lu-s%015.12f-index",
           a, total_coating_width,
           best_RCS, (best_RCS/initial_RCS-1.0)*100.0, best_x.size(), index_sum);
  wrapper.SetPlotName(plot_name);
  wrapper.SetXLabelName("Layer #");
  wrapper.SetYLabelName("Index");
  wrapper.SetDrawStyle("with histeps lw 2");
  wrapper.SetXRange({0.51, number_of_layers_+0.49});
  for (int i = 0; i < number_of_layers_; ++i) 
    wrapper.AddMultiPoint({i+1.0, best_x[i]});
  wrapper.AddColumnName("Layer N");
  wrapper.AddColumnName("Index");
  wrapper.MakeOutput();
}
// ********************************************************************** //
// ********************************************************************** //
// ********************************************************************** //
double SetInitialModel() {
  // Set common parameters for all wavelengths.
  multi_layer_mie_.AddTargetLayer(a, {-1.0, -1.0});  // Set PEC core
  multi_layer_mie_.SetQfaild(1000.0);  // Searching for minima
  multi_layer_mie_.SetWavelength(lambda_work);
  double Qext, Qsca, Qabs, Qbk;
  multi_layer_mie_.RunMie(&Qext, &Qsca, &Qabs, &Qbk);
  double total_r = multi_layer_mie_.GetTotalRadius();
  double initial_RCS = Qsca*pi*pow2(total_r);
  return initial_RCS;
}
// ********************************************************************** //
// ********************************************************************** //
// ********************************************************************** //
void SetOptimizer() {
  long dimension = 0;
  SetThickness();
  dimension = number_of_layers_;
  sub_population_.FitnessFunction = &EvaluateScatterOnlyIndex;
  long total_population = dimension * 3;
  sub_population_.Init(total_population, dimension);
  /// Low and upper bound for all dimensions;
  double from_n = 1.0, to_n = 20.0;
  sub_population_.SetAllBounds(from_n, to_n);
  sub_population_.SetTargetToMinimum();
  sub_population_.SetTotalGenerationsMax(total_generations);
  sub_population_.SwitchOffPMCRADE();

  sub_population_.SetBestShareP(0.1);
  sub_population_.SetAdapitonFrequencyC(1.0/20.0);

}
// ********************************************************************** //
// ********************************************************************** //
// ********************************************************************** //
void SetTarget(double n, double k) {
  multi_layer_mie_.AddTargetLayer(a, {n, k});
}
// ********************************************************************** //
// ********************************************************************** //
// ********************************************************************** //
void SetThickness() {
  std::vector<double> thickness;
  thickness.clear();
  if (number_of_layers_ < 0)
    throw std::invalid_argument("Number of coating layers should be >= 0!");
  for (int i = 0; i < number_of_layers_; ++i) thickness.push_back(layer_thickness);
  multi_layer_mie_.SetCoatingThickness(thickness);
}
// ********************************************************************** //
// ********************************************************************** //
// ********************************************************************** //
double EvaluateScatterOnlyIndex(std::vector<double> input) {
  double Qext, Qsca, Qabs, Qbk;
  std::vector<complex> cindex;
  cindex.clear();
  double k = min_index;
  for (auto n : input) cindex.push_back({n, k});
  multi_layer_mie_.SetCoatingIndex(cindex);
  try {
    multi_layer_mie_.RunMie(&Qext, &Qsca, &Qabs, &Qbk);
  } catch( const std::invalid_argument& ia ) {
    printf(".");
    // Will catch if  multi_layer_mie_ fails or other errors.
    //std::cerr << "Invalid argument: " << ia.what() << std::endl;
  }  
  double total_r = multi_layer_mie_.GetTotalRadius();
  return Qsca*pi*pow2(total_r);
}
// ********************************************************************** //
// ********************************************************************** //
// ********************************************************************** //
std::vector< std::vector<double> > EvaluateSpectraForBestDesign() {
  double best_RCS;
  auto best_x = sub_population_.GetBest(&best_RCS);
  // Setting Mie model to the best state.
  sub_population_.FitnessFunction(best_x);
  return multi_layer_mie_.GetSpectra(lambda_work*0.5, lambda_work*1.5, 1000);
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
  for (int i = 0; i < number_of_layers_; ++i)
    printf("% 5i\t",i+1);
  printf ("\n");
  printf ("Index:\t");
  for (int i = 0; i < number_of_layers_; ++i)
    printf("%5.4g\t",best_x[i]);
  printf("\n");
  double total_coating_width = 0.0;
  total_coating_width = layer_thickness*number_of_layers_;
  printf("Total coating width: %g\n", total_coating_width);
}
// ********************************************************************** //
// ********************************************************************** //
// ********************************************************************** //
void PrintGnuPlotSpectra(std::vector< std::vector<double> > spectra,
                         double initial_RCS) {
  gnuplot::GnuplotWrapper wrapper;
  double best_RCS = 0.0;
  auto best_x = sub_population_.GetBest(&best_RCS);
  double total_coating_width = layer_thickness*number_of_layers_;
  double index_sum = 0.0;
  for (auto i : best_x) index_sum+=i;
  char plot_name [300];
  snprintf(plot_name, 300,
           "TargetR%6.4f-CoatingW%06.3f-FinalRCS%07.4fdiff%+4.1f%%-n%lu-s%015.12f-spectra",
           a, total_coating_width,
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
