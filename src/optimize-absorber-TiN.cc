/**
 * @file   optimize-absorber-TiN.cc
 * @author Konstantin Ladutenko <kostyfisik at gmail (.) com>
 * @date   Wed Mar 11 10:59:03 2015
**/
/// @copyright 2015  Konstantin Ladutenko
///
/// optimize-absorber-TiN is free software: you can redistribute it and/or modify
/// it under the terms of the GNU General Public License as published by
/// the Free Software Foundation, either version 3 of the License, or
/// (at your option) any later version.
///
/// optimize-absorber-TiN is distributed in the hope that it will be useful,
/// but WITHOUT ANY WARRANTY; without even the implied warranty of
/// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
/// GNU General Public License for more details.
///
/// You should have received a copy of the GNU General Public License
/// along with optimize-absorber-TiN.  If not, see <http://www.gnu.org/licenses/>.
///
/// optimize-absorber-TiN uses nmie.cc from scattnlay by Ovidio Pena
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
#include "./nmie/nmie.h"
#include "./read-spectra/read-spectra.h"
const double pi=3.14159265358979323846;
template<class T> inline T pow2(const T value) {return value*value;}
// Mie model. Used in fitness function for optimization of
// sub_population.
void SetTarget(double n, double k);
void SetOptimizer();

double EvaluateFitness(std::vector<double> input);
std::vector< std::vector<double> > EvaluateSpectraForBestDesign();
void PrintCoating(std::vector<double> current, double initial_RCS,
                    jade::SubPopulation sub_population);
void PrintGnuPlotSpectra(std::vector< std::vector<double> > spectra);
// ********************************************************************** //
// ********************************************************************** //
// ********************************************************************** //
nmie::MultiLayerMie multi_layer_mie;  
jade::SubPopulation sub_population;  // Optimizer of parameters for Mie model.
read_spectra::ReadSpectra core_index_, TiN_;
// ********************************************************************** //
// ********************************************************************** //
// ********************************************************************** //
double lambda_best_ = 0.0; // nm
double core_share_ = 0.0, TiN_share_ = 0.0;
double Qabs_=0.0, initial_Qabs_=0.0;
double total_r_ = 0.0; 
// ********************************************************************** //
// Set model: core->TiN->shell
double max_r_ = 150; // nm
double max_TiN_width_ = 10; // nm
// Set dispersion
double from_wl_ = 300.0, to_wl_ = 900.0;
int samples_ = 3;
bool isGaAs = false; // Select Si of GaAs as a material for core and shell
// Set optimizer
int total_generations_ = 150;
int population_multiplicator_ = 6;
double step_r_ = max_r_ / 3.0;
// ********************************************************************** //
// ********************************************************************** //
// ********************************************************************** //
int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  try {
    if (isGaAs)
      core_index_.ReadFromFile("GaAs.txt")
    else
      core_index_.ReadFromFile("Si.txt")
    core_index_.ResizeToComplex(from_wl, to_wl, samples_).ToIndex();
    TiN_.ReadFromFile("TiN.txt").ResizeToComplex(from_wl, to_wl, samples_).ToIndex();
    if (core_index_.GetIndex().size()
	!= TiN_.GetIndex().size()) throw invalid_argument("Unexpected sampling of dispersion!/n");
    // if (rank == 0) {  printf("\nSi:\n");  Si_.PrintData(); printf("\nGaAs:\n");
    //   GaAs_.PrintData();  printf("\nTiN:\n");   TiN_.PrintData(); }
    initial_Qabs_ = FitnessFunction({1.0,0.0});  // Only core
    int output_rank = 0;
    for (total_r_ = step_r_; total_r+=step_r_; total_r < max_r_*1.00001) {
    //   SetOptimizer();
    //   sub_population.RunOptimization();
    //   auto current = sub_population.GetFinalFitness();
    //   double best_RCS = 0.0;
    //   auto best_x = sub_population.GetBest(&best_RCS);
    //   spectra1.push_back({total_thickness, best_x[0]});
    //   spectra2.push_back({total_thickness, best_RCS});
    //   // Output results
    //   for (unsigned int i = 0; i < current.size(); ++i)
    // 	if (current[output_rank] > current[i]) output_rank = i;
    //   if (rank == output_rank) {
    // 	for (auto c : current) printf("All %g\n",c);
    // 	printf("####### total_thickness = %g\n",total_thickness);
    // 	PrintCoating(current, initial_RCS, sub_population);
    //   }  // end of output for process with best final fitness
    //   //PrintGnuPlotIndex(initial_RCS, sub_population);
    //   // PrintGnuPlotSpectra(EvaluateSpectraForBestDesign(), initial_RCS);
    //   //sub_population.PrintResult("-- ");
    }  // end of total coating thickness sweep
    // PrintGnuPlotThickness(spectra1, 0.1, "1st layer share");
    // PrintGnuPlotThickness(spectra2, 0.2, "Total RCS");
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
void SetOptimizer() {
  //Width is optimized for two layers only!!
  //The third one fills to total_r_
  long dimension = 2;
  sub_population.FitnessFunction = &EvaluateFitness;
  long total_population = dimension * population_multiplicator_;
  sub_population.Init(total_population, dimension);
  /// Low and upper bound for all dimenstions;
  sub_population.SetAllBounds(0.0, 1.0);
  sub_population.SetTargetToMaximum();
  sub_population.SetTotalGenerationsMax(total_generations_);
  //sub_population.SwitchOffPMCRADE();

  sub_population.SetBestShareP(0.1);
  sub_population.SetAdapitonFrequencyC(1.0/20.0);

}
// ********************************************************************** //
// ********************************************************************** //
// ********************************************************************** //
double EvaluateFitness(std::vector<double> input) {
  double Qsca;//, Qabs, Qbk;
  auto core = core_index_.GetIndex();
  auto TiN = TiN_.GetIndex();
  std::vector<std::complex<double> > cindex;
  cindex.clear();
  double k = loss_index;
  for (auto n : input) cindex.push_back({n, k});
  multi_layer_mie.SetCoatingIndex(cindex);
  multi_layer_mie.SetWavelength(lambda_work);
  try {
    multi_layer_mie.RunMieCalculations();
    Qsca = multi_layer_mie.GetQsca();
  } catch( const std::invalid_argument& ia ) {
    auto best_x = sub_population_.GetWorst(&Qabs_);
    printf("#");
    // Will catch if  multi_layer_mie_ fails or other errors.
    //std::cerr << "Invalid argument: " << ia.what() << std::endl;
  }  
  return Qabs_;
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
  // auto best_x = sub_population.GetBest(&best_RCS);
  // printf("Target R=%g, WL=%g\n", a, lambda_work);
  // printf("Initial RCS: %g\n", initial_RCS);
  //  printf("Final RCS: %g (%4.1f%%)\n", best_RCS, (best_RCS/initial_RCS-1.0)*100.0);
  // printf ("Layer:\t");
  // for (int i = 0; i < number_of_layers; ++i)
  //   printf("% 5i\t",i+1);
  // printf ("\n");
  // printf ("Width:\t");
  // for (int i = 0; i < number_of_layers; ++i)
  //   printf("%5.4g\t",best_x[i]);
  // printf("\n");
  // int shift = isStratFromLowIndex ? 1 : 0;
  // printf ("Index:\t");
  // for (int i = 0; i < number_of_layers; ++i) {
  //   if ( (i+shift) % 2 ) printf("%5.4g\t",low_index);
  //   else printf("%5.4g\t", hi_index);
  // }  // end of for each layer
  // printf("\n");  
  // double total_coating_width = 0.0;
  // for (auto w : best_x) total_coating_width += w;
  // printf("Total coating width: %g\n", total_coating_width);
  // printf("Layer width limits min/max: %g/%g\n", min_layer_thickness, max_layer_thickness);
  // printf("Layer index limits min/max: %g/%g\n", low_index, hi_index);
  
}
// ********************************************************************** //
// ********************************************************************** //
// ********************************************************************** //
void PrintGnuPlotSpectra(std::vector< std::vector<double> > spectra) {
  gnuplot::GnuplotWrapper wrapper;
  char plot_name [300];
  snprintf(plot_name, 300,
           "TotalR%06.2f-Qabs%06.3f:%+4.1f%%-core%07.2f:%+4.1f-TiN%07.2f:%+4.1fshell%07.2f:%+4.1f-spectra",
           total_r_, Qabs_, (Qabs_/initial_Qabs_-1.0)*100.0,
	   (max_r_ - max_TiN_width_*TiN_share_)*core_share_, core_share_,
	   max_TiN_width_*TiN_share_, TiN_share_
	   (max_r_ - max_TiN_width_*TiN_share_)*(1.0-core_share_), (1.0-core_share_));
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
