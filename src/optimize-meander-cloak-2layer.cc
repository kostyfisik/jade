///
/// @file   optimize-meander-cloak.cc
/// @author Ladutenko Konstantin <kostyfisik at gmail (.) com>
/// @date   Thu Jan 30 21:14:37 2014
/// @copyright 2014 Ladutenko Konstantin
///
/// optimize-meander-cloak is free software: you can redistribute it and/or modify
/// it under the terms of the GNU General Public License as published by
/// the Free Software Foundation, either version 3 of the License, or
/// (at your option) any later version.
///
/// optimize-meander-cloak is distributed in the hope that it will be useful,
/// but WITHOUT ANY WARRANTY; without even the implied warranty of
/// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
/// GNU General Public License for more details.
///
/// You should have received a copy of the GNU General Public License
/// along with optimize-meander-cloak.  If not, see <http://www.gnu.org/licenses/>.
///
/// optimize-meander-cloak uses nmie.c from scattnlay by Ovidio Pena
/// <ovidio@bytesfall.com> as a linked library. He has an additional condition to 
/// his library:
// The only additional condition is that we expect that all publications         //
// describing  work using this software , or all commercial products             //
// using it, cite the following reference:                                       //
// [1] O. Pena and U. Pal, "Scattering of electromagnetic radiation by           //
//     a multilayered sphere," Computer Physics Communications,                  //
//     vol. 180, Nov. 2009, pp. 2348-2354.                                       //
///
/// @brief Simulate scattering from PEC sphere covered with dielectric
/// multilayered shell
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
// Mie model. Used in fitness function for optimization of
// sub_population.
nmie::MultiLayerMie multi_layer_mie;  
jade::SubPopulation sub_population;  // Optimizer of parameters for Mie model.
// ********************************************************************** //
// ********************************************************************** //
// ********************************************************************** //
bool isUsingPEC = true;
//bool isUsingPEC = false;
bool isStratFromLowIndex = true;
//bool isStratFromLowIndex = false;
double low_index = 0.67;
double hi_index = 8;
// Semouchkina APPLIED PHYSICS LETTERS 102, 113506 (2013)
double lambda_work = 0.532; // cm
//    double f_work = 30/lambda_work; // 8 GHz
double a = 0.75*lambda_work;  // 2.8125 cm
int number_of_layers = 2;
double total_thickness = 0.001;
double layer_thickness = total_thickness /
  static_cast<double>(number_of_layers);
double min_layer_thickness = 0.0;
double max_layer_thickness =0.2;
double n = 4;
double k = 0;
// // Production parameters
// int total_generations = 100;
// double thickness_step = 0.02;
// Test parameters
int total_generations = 1500;
int population_multiplicator = 60;
double thickness_step = 0.003;

void SetTarget(double n, double k);
void SetMeanderIndex();
double SetInitialModel(double n, double k);
void SetOptimizer();

double EvaluateScatterOnlyThickness(std::vector<double> input);
std::vector< std::vector<double> > EvaluateSpectraForBestDesign();
void PrintCoating(std::vector<double> current, double initial_RCS,
                    jade::SubPopulation sub_population);
void PrintGnuPlotIndex(double initial_RCS,
                  jade::SubPopulation sub_population);
void PrintGnuPlotSpectra(std::vector< std::vector<double> > spectra,
                         double initial_RCS);
void PrintGnuPlotThickness(std::vector< std::vector<double> > spectra,
			   double meta_n, std::string ylabel);
// ********************************************************************** //
// ********************************************************************** //
// ********************************************************************** //
double loss_index = 1e-7;
int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  try {
    if (isUsingPEC) {n = -1.0; k = -1.0;}
    // for (k = 0.0; k < 1.0; k = (k+0.0001)*2.0)  {
    double initial_RCS = SetInitialModel(n, k);
    std::vector< std::vector<double> > spectra1;
    std::vector< std::vector<double> > spectra2;
    int output_rank = 0;
        
    for (total_thickness = 0.003; total_thickness < 0.12;
	 total_thickness += thickness_step) {
      //for (number_of_layers = 8; number_of_layers < 15; number_of_layers *=2) {
      //          number_of_layers = 4;
      // layer_thickness = total_thickness /
      //   static_cast<double>(number_of_layers);
      //minimal_thickness = layer_thickness * 0.001;
      SetOptimizer();
      sub_population.RunOptimization();
      auto current = sub_population.GetFinalFitness();
      double best_RCS = 0.0;
      auto best_x = sub_population.GetBest(&best_RCS);
      spectra1.push_back({total_thickness, best_x[0]});
      spectra2.push_back({total_thickness, best_RCS});
      // Output results
      for (unsigned int i = 0; i < current.size(); ++i)
	if (current[output_rank] > current[i]) output_rank = i;
      if (rank == output_rank) {
	for (auto c : current) printf("All %g\n",c);
	printf("####### total_thickness = %g\n",total_thickness);
	PrintCoating(current, initial_RCS, sub_population);
      }  // end of output for process with best final fitness
      //PrintGnuPlotIndex(initial_RCS, sub_population);
      // PrintGnuPlotSpectra(EvaluateSpectraForBestDesign(), initial_RCS);
      //sub_population.PrintResult("-- ");
      //}  // end of changing number of layers
    }  // end of total coating thickness sweep
    PrintGnuPlotThickness(spectra1, 0.1, "1st layer share");
    PrintGnuPlotThickness(spectra2, 0.2, "Total RCS");
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
void PrintGnuPlotThickness(std::vector< std::vector<double> > spectra,
			   double meta_n, std::string ylabel) {
  gnuplot::GnuplotWrapper wrapper;
  //  auto best_x = sub_population.GetBest(&best_RCS);
  double total_coating_width = layer_thickness*number_of_layers;
  double index_sum = 0.0;
  //for (auto i : best_x) index_sum+=i;
  char plot_name [300];
  snprintf(plot_name, 300,
           "TargetR%6.4f-CoatingW%06.3f-index%g-spectra",
           a, total_coating_width,meta_n);
  wrapper.SetPlotName(plot_name);
  wrapper.SetXLabelName("Thickness");
  wrapper.SetYLabelName(ylabel);
  wrapper.SetDrawStyle("w l lw 2");
  wrapper.SetXRange({spectra.front()[0], spectra.back()[0]});
  for (auto multi_point : spectra) wrapper.AddMultiPoint(multi_point);
  wrapper.AddColumnName("Thickness");
  wrapper.AddColumnName("Qsca");
  wrapper.MakeOutput();
}
// ********************************************************************** //
// ********************************************************************** //
// ********************************************************************** //
void PrintGnuPlotIndex(double initial_RCS,
                  jade::SubPopulation sub_population) {
  gnuplot::GnuplotWrapper wrapper;
  double best_RCS = 0.0;
  auto best_x = sub_population.GetBest(&best_RCS);  
  double total_coating_width = 0.0;
  for (auto w : best_x) total_coating_width += w;
  double index_sum = 0.0;
  for (auto i : best_x) index_sum+=i;
  char plot_name [300];
  snprintf(plot_name, 300,
           "TargetR%g-index-n%gk%g-FinalRCS%7.4fdiff%+4.1f%%-CoatingW%06.3f-n%lu-s%015.12f-index",
           a, n, k,
           best_RCS, (best_RCS/initial_RCS-1.0)*100.0, total_coating_width, best_x.size(), index_sum);
  
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
  SetMeanderIndex();
  dimension = number_of_layers;
  //Two layers only!!
  dimension = 1;
  sub_population.FitnessFunction = &EvaluateScatterOnlyThickness;
  long total_population = dimension * population_multiplicator;
  sub_population.Init(total_population, dimension);
  /// Low and upper bound for all dimenstions;
  sub_population.SetAllBounds(0.001,
                              1.0);
  // sub_population.SetAllBounds(min_layer_thickness,
  //                             max_layer_thickness);
  sub_population.SetTargetToMinimum();
  sub_population.SetTotalGenerationsMax(total_generations);
  sub_population.SwitchOffPMCRADE();

  sub_population.SetBestShareP(0.1);
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
void SetMeanderIndex() {
  std::vector<complex> cindex;
  double k = loss_index;
  int shift = isStratFromLowIndex ? 1 : 0;
  if (number_of_layers < 0)
    throw std::invalid_argument("Number of coating layers should be >= 0!");
  if (low_index > hi_index)
    throw std::invalid_argument("low_index should be less the hi_index!");
  for (int i = 0; i < number_of_layers; ++i) {
    if ( (i+shift) % 2 ) cindex.push_back({low_index,k});
    else cindex.push_back({hi_index, k});
  }  // end of for each layer
  multi_layer_mie.SetCoatingIndex(cindex);
  for (auto i : cindex) printf("%g ", i.r);
  printf("\n");
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
double EvaluateScatterOnlyThickness(std::vector<double> input) {
  std::vector<double> thickness;
  
  //Two layers only!!
  if (number_of_layers != 2) 
        throw std::invalid_argument("Number of coating layers should be = 2!");
  if (input[0]<0.001) input[0]=0.001;
  if (input[0]>0.999) input[0]=0.999;
  
  thickness.push_back(input[0]*total_thickness);
  thickness.push_back((1.0-input[0])*total_thickness);
  
// for (int i = 0; i < number_of_layers; ++i) {
  //   thickness.push_back(input[i]);
  // }
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
  // printf("total_r = %g, Qsca = %g\n", total_r, Qsca);
  // throw std::invalid_argument("Break point!");
  return Qsca*pi*pow2(total_r);
}
// ********************************************************************** //
// ********************************************************************** //
// ********************************************************************** //
// double EvaluateScatter(std::vector<double> input) {
//   std::vector<double> thickness;
//   std::vector<complex> cindex;
//   double k = loss_index;
//   for (int i = 0; i < number_of_layers; ++i) {
//     cindex.push_back({input[i], k});
//     if (input[i+number_of_layers] < 1.0) input[i+number_of_layers] = 1.0; 
//     thickness.push_back(input[i+number_of_layers]*layer_thickness);
//   }
//   multi_layer_mie.SetCoatingIndex(cindex);
//   multi_layer_mie.SetCoatingThickness(thickness);
//   double Qext, Qsca, Qabs, Qbk;
//   try {
//     multi_layer_mie.RunMie(&Qext, &Qsca, &Qabs, &Qbk);
//   } catch( const std::invalid_argument& ia ) {
//     printf(".");
//     // Will catch if  multi_layer_mie fails or other errors.
//     //std::cerr << "Invalid argument: " << ia.what() << std::endl;
//   }  
//   double total_r = multi_layer_mie.GetTotalRadius();
//   return Qsca*pi*pow2(total_r);
// }
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
  printf ("Width:\t");
  for (int i = 0; i < number_of_layers; ++i)
    printf("%5.4g\t",best_x[i]);
  printf("\n");
  int shift = isStratFromLowIndex ? 1 : 0;
  printf ("Index:\t");
  for (int i = 0; i < number_of_layers; ++i) {
    if ( (i+shift) % 2 ) printf("%5.4g\t",low_index);
    else printf("%5.4g\t", hi_index);
  }  // end of for each layer
  printf("\n");  
  double total_coating_width = 0.0;
  for (auto w : best_x) total_coating_width += w;
  printf("Total coating width: %g\n", total_coating_width);
  printf("Layer width limits min/max: %g/%g\n", min_layer_thickness, max_layer_thickness);
  printf("Layer index limits min/max: %g/%g\n", low_index, hi_index);
  
}
// ********************************************************************** //
// ********************************************************************** //
// ********************************************************************** //
void PrintGnuPlotSpectra(std::vector< std::vector<double> > spectra,
                         double initial_RCS) {
  gnuplot::GnuplotWrapper wrapper;
  double best_RCS = 0.0;
  auto best_x = sub_population.GetBest(&best_RCS);
  double total_coating_width = 0.0;
  for (auto w : best_x) total_coating_width += w;
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
