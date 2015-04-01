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
#include "./nmie/nmie-wrapper.h"
#include "./read-spectra/read-spectra.h"
const double pi=3.14159265358979323846;
template<class T> inline T pow2(const T value) {return value*value;}
// Mie model. Used in fitness function for optimization of
// sub_population.
void SetOptimizer();

double EvaluateFitness(std::vector<double> input);
double EvaluateFitnessChannel(std::vector<double> input);
std::vector< std::vector<double> > EvaluateSpectraForBestDesign();
std::vector< std::vector<double> > EvaluateSpectraForChannels();
void Print();
void PrintCoating(std::vector<double> current, double initial_RCS,
                    jade::SubPopulation sub_population);
void PrintGnuPlotSpectra(std::vector< std::vector<double> > spectra);
void PrintGnuPlotChannels(std::vector< std::vector<double> > spectra);
void PrintGnuPlotChannelSweep(std::vector< std::vector<double> > spectra);
void PrintGnuPlotEpsilon(std::vector< std::vector<double> > spectra);
// ********************************************************************** //
// ********************************************************************** //
// ********************************************************************** //
nmie::MultiLayerMie multi_layer_mie_;  
jade::SubPopulation sub_population_;  // Optimizer of parameters for Mie model.
std::string sign_, full_sign_;
double eps_re_ = 0.0, eps_im_ = 0.0;
// ********************************************************************** //
// ********************************************************************** //
// ********************************************************************** //
double lambda_best_ = 0.0; // nm
int fails_ = 0;
double Qabs_=0.0, initial_Qabs_=0.0;
double total_r_ = 0.0; 
// ********************************************************************** //
// Set model: core->TiN->shell
const double max_r_ = 159.0; // nm
// Set dispersion
double at_wl_ = 500.0;
double from_wl_ = at_wl_, to_wl_ = at_wl_;
int samples_ = 1;
// double from_wl_ = 300.0, to_wl_ = 900.0;
// int samples_ = 151;
double plot_from_wl_ = 300.0, plot_to_wl_ = 900.0;
int plot_samples_ = 351;
double plot_step_wl_ = (plot_to_wl_-plot_from_wl_)/static_cast<double>(plot_samples_);
// Set optimizer
int total_generations_ = 250;
int population_multiplicator_ = 160;
double bound = 10.0;
double step_r_ = 1.0; //max_r_ / 159.0;
int channel_ = 3;
// ********************************************************************** //
// ********************************************************************** //
// ********************************************************************** //
const double eps_num_=1e-11;
int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  try {
    sub_population_.FitnessFunction = &EvaluateFitness;
    sign_ = "bulk-bn"+std::to_string(channel_);
    std::cout << "Sign: " << sign_ << std::endl;
    SetOptimizer();
    if (step_r_ <=0.0) throw std::invalid_argument("Radius step should be positive!/n");
    if (channel_ < 1) throw std::invalid_argument("Channel number starts from 1!/n");
    auto best_x = sub_population_.GetBest(&Qabs_);
    std::vector< std::vector<double> > channel_sweep;
    std::vector< std::vector<double> > epsilon_sweep;
    int least_size_sweep = 15;
    //double best_Qabs = 0.0, best_total_r = 0.0;
    // ***************************************************
    // **************  Main loop   ***********************
    // ***************************************************  
    for (total_r_ = step_r_; total_r_ < max_r_*1.00001; total_r_+=step_r_) {
      //for (total_r_ = 145.0; total_r_ < 147.0; total_r_+=0.05) {
      if (rank == 0) printf("\nTotal R = %g\n", total_r_);    
      sub_population_.RunOptimization();
      // Plot spectra from each process
      PrintGnuPlotSpectra(EvaluateSpectraForBestDesign());
      PrintGnuPlotChannels(EvaluateSpectraForChannels());
      if (rank == 0) {
	Print();
	std::vector<std::complex<double> > an = multi_layer_mie_.GetAn();
	std::vector<std::complex<double> > bn = multi_layer_mie_.GetBn();
	for (int i = 0; i < 3; ++i)
	  printf("an[%d]=%g, bn[%d]=%g\n",i,std::abs(an[i]),i, std::abs(bn[i]));
      }
      std::vector<double> tmp({total_r_});
      //std::vector<double> channels(multi_layer_mie_.GetQabs_channel());
      std::vector<std::complex<double> > an = multi_layer_mie_.GetAn();
      std::vector<std::complex<double> > bn = multi_layer_mie_.GetBn();
      for (int i = 0; i < an.size(); ++i) {
	tmp.push_back(an[i].real() - pow2(std::abs(an[i])));
	tmp.push_back(bn[i].real() - pow2(std::abs(bn[i])));
      }
      channel_sweep.push_back(tmp);
      if (least_size_sweep > tmp.size()) least_size_sweep = tmp.size();
      epsilon_sweep.push_back({total_r_, eps_re_, eps_im_});
    }  // end of total coating thickness sweep
    for (auto& row : channel_sweep) row.resize(least_size_sweep);
    PrintGnuPlotChannelSweep(channel_sweep);
    PrintGnuPlotEpsilon(epsilon_sweep);
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
void Print() {
  printf("Spectra: %s\n", full_sign_.c_str());
  // TiN_share_ and core_share_ are set in EvaluateSpectraForBestDesign()
  printf("eps_re:%.19g\neps_im:%.19g\n\n",eps_re_, eps_im_);  
}
// ********************************************************************** //
// ********************************************************************** //
// ********************************************************************** //
void SetOptimizer() {
  //Width is optimized for two layers only!!
  //The third one fills to total_r_
  long dimension = 2;
  //sub_population_.FitnessFunction = &EvaluateFitness;
  //sub_population_.FitnessFunction = &EvaluateFitnessChannel;
  long total_population = dimension * population_multiplicator_;
  sub_population_.Init(total_population, dimension);
  /// Low and upper bound for all dimenstions;
  
  sub_population_.SetAllBounds(-bound, bound);
  sub_population_.SetTargetToMaximum();
  sub_population_.SetTotalGenerationsMax(total_generations_);
  //sub_population.SwitchOffPMCRADE();

  sub_population_.SetBestShareP(0.1);
  sub_population_.SetAdapitonFrequencyC(1.0/20.0);

}
// ********************************************************************** //
// ********************************************************************** //
// ********************************************************************** //
double EvaluateFitness(std::vector<double> input) {
  if (input.size() != 2) throw std::invalid_argument("Wrong input dimension!/n");
  eps_re_ = input[0];
  eps_im_ = input[1];
  std::complex<double> epsilon(eps_re_, eps_im_);
  multi_layer_mie_.ClearTarget();
  multi_layer_mie_.AddTargetLayer(total_r_, std::sqrt(epsilon));
  multi_layer_mie_.SetWavelength(at_wl_);
  try {
    multi_layer_mie_.RunMieCalculations();
    //std::vector<double> channels(multi_layer_mie_.GetQabs_channel_normalized());
    //Qabs_ = multi_layer_mie_.GetQabs();
    std::complex<double> an = multi_layer_mie_.GetBn()[channel_ -1];
    Qabs_ = an.real() - pow2(std::abs(an));
  } catch( const std::invalid_argument& ia ) {
    printf(".");
    sub_population_.GetWorst(&Qabs_);
  }
  return Qabs_;
}
// ********************************************************************** //
// ********************************************************************** //
// ********************************************************************** //
std::vector< std::vector<double> > EvaluateSpectraForBestDesign() {
  fails_ = 0;
  auto best_x = sub_population_.GetBest(&Qabs_);
  // Setting Mie model to the best state.
  if (best_x.size() != 2) throw std::invalid_argument("Wrong input dimension!/n");
  eps_re_ = best_x[0];
  eps_im_ = best_x[1];
  std::complex<double> epsilon(eps_re_, eps_im_);
  multi_layer_mie_.ClearTarget();
  multi_layer_mie_.AddTargetLayer(total_r_, std::sqrt(epsilon));
  double Qabs = 0.0, Qext=0.0, Qsca=0.0, Qbk =0.0;
  std::vector< std::vector<double> > spectra;
  for (double wl = plot_from_wl_; wl < plot_to_wl_; wl+=plot_step_wl_) {
    multi_layer_mie_.SetWavelength(wl);
    try {
      multi_layer_mie_.RunMieCalculations();
      Qabs = multi_layer_mie_.GetQabs();
      Qext = multi_layer_mie_.GetQext();
      Qsca = multi_layer_mie_.GetQsca();
      Qabs = multi_layer_mie_.GetQabs();
      spectra.push_back({wl,Qext,Qsca,Qabs,Qbk});
    } catch( const std::invalid_argument& ia ) {
      printf(".");
      ++fails_;
    }
  }  // end of for all points of the spectrum
  return spectra;
}
// ********************************************************************** //
// ********************************************************************** //
// ********************************************************************** //
std::vector< std::vector<double> > EvaluateSpectraForChannels() {
  fails_ = 0;
  auto best_x = sub_population_.GetBest(&Qabs_);
  // Setting Mie model to the best state.
  if (best_x.size() != 2) throw std::invalid_argument("Wrong input dimension!/n");
  eps_re_ = best_x[0];
  eps_im_ = best_x[1];
  std::complex<double> epsilon(eps_re_, eps_im_);

  multi_layer_mie_.ClearTarget();
  multi_layer_mie_.AddTargetLayer(total_r_, std::sqrt(epsilon));
  double max_Qabs = 0.0, Qabs = 0.0, Qext=0.0, Qsca=0.0, Qbk =0.0;
  std::vector< std::vector<double> > spectra;
  int least_size = 10000;
  // Scan all wavelengths
  for (double wl = plot_from_wl_; wl < plot_to_wl_; wl+=plot_step_wl_) {
    multi_layer_mie_.SetWavelength(wl);
    try {
      multi_layer_mie_.RunMieCalculations();
       std::vector<double> tmp({wl});
      //std::vector<double> channels(multi_layer_mie_.GetQabs_channel());
      std::vector<double> channels(multi_layer_mie_.GetQabs_channel_normalized());
      tmp.insert(tmp.end(), channels.begin(), channels.end());
      spectra.push_back(tmp);
      if (least_size > tmp.size()) least_size = tmp.size();
   } catch( const std::invalid_argument& ia ) {
      printf(".");
      ++fails_;
    }
  }  // end of for all points of the spectrum
  for (auto& row : spectra) row.resize(least_size);
  return spectra;
}
// ********************************************************************** //
// ********************************************************************** //
// ********************************************************************** //
void PrintGnuPlotSpectra(std::vector< std::vector<double> > spectra) {
  gnuplot::GnuplotWrapper wrapper;
  char plot_name [300];
  auto best_x = sub_population_.GetBest(&Qabs_);
  EvaluateFitness(best_x);
  Qabs_ = multi_layer_mie_.GetQabs();
  snprintf(plot_name, 300,
           "%s-TotalR%06.2fnm-Qabs%016.13f--epsRe%07.2f--epsIm%07.2f-fails%d-spectra",
	   sign_.c_str(), total_r_, Qabs_,  eps_re_, eps_im_, fails_);
  full_sign_ = std::string(plot_name);
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
void PrintGnuPlotChannels(std::vector< std::vector<double> > spectra) {
  gnuplot::GnuplotWrapper wrapper;
  char plot_name [300];
  auto best_x = sub_population_.GetBest(&Qabs_);
  EvaluateFitness(best_x);
  Qabs_ = multi_layer_mie_.GetQabs();
  std::vector<std::complex<double> > an = multi_layer_mie_.GetAn();
  std::vector<std::complex<double> > bn = multi_layer_mie_.GetBn();
  snprintf(plot_name, 300,
           "%s-ch-TotalR%06.2fnm-Qabs%016.13f--epsRe%07.4f--epsIm%07.4f--an%g--bn%g--fails%d-spectra",  sign_.c_str(), total_r_, Qabs_,  eps_re_, eps_im_, std::abs(an[0]), std::abs(bn[0]), fails_);
  full_sign_ = std::string(plot_name);
  wrapper.SetPlotName(plot_name);
  wrapper.SetXLabelName("WL");
  wrapper.SetYLabelName("NACS");
  wrapper.SetDrawStyle("w l lw 2");
  wrapper.SetXRange({spectra.front()[0], spectra.back()[0]});
  for (auto multi_point : spectra) wrapper.AddMultiPoint(multi_point);
  wrapper.AddColumnName("WL");
  for (int i = 1; i < spectra.front().size(); ++i) {
    char column_name[10];
    snprintf(column_name, 10, "[%d]",i);
    wrapper.AddColumnName(column_name);
  }
  wrapper.MakeOutput();
}
// ********************************************************************** //
// ********************************************************************** //
// ********************************************************************** //
void PrintGnuPlotChannelSweep(std::vector< std::vector<double> > spectra) {
  gnuplot::GnuplotWrapper wrapper;
  char plot_name [300];
  auto best_x = sub_population_.GetBest(&Qabs_);
  EvaluateFitness(best_x);
  Qabs_ = multi_layer_mie_.GetQabs();
  std::vector<std::complex<double> > an = multi_layer_mie_.GetAn();
  std::vector<std::complex<double> > bn = multi_layer_mie_.GetBn();
  snprintf(plot_name, 300,
           "%s-sweep-ch",  sign_.c_str());
  full_sign_ = std::string(plot_name);
  wrapper.SetPlotName(plot_name);
  wrapper.SetXLabelName("Total_R");
  wrapper.SetYLabelName("NACS");
  wrapper.SetDrawStyle("w l lw 2");
  wrapper.SetXRange({spectra.front()[0], spectra.back()[0]});
  wrapper.SetYRange({-0.02, 0.27});
  for (auto multi_point : spectra) wrapper.AddMultiPoint(multi_point);
  wrapper.AddColumnName("WL");
  for (int i = 1; i < spectra.front().size(); ++i) {
    char column_name[10];
    std::string type;
    if ( (i-1)%2 == 0) type ="an";
    else type="bn";
    snprintf(column_name, 10, "%s[%d]",type.c_str(), (i-1)/2+1);
    wrapper.AddColumnName(column_name);
  }
  wrapper.MakeOutput();
}
// ********************************************************************** //
// ********************************************************************** //
// ********************************************************************** //

void PrintGnuPlotEpsilon(std::vector< std::vector<double> > spectra) {
  gnuplot::GnuplotWrapper wrapper;
  char plot_name [300];
  auto best_x = sub_population_.GetBest(&Qabs_);
  EvaluateFitness(best_x);
  Qabs_ = multi_layer_mie_.GetQabs();
  std::vector<std::complex<double> > an = multi_layer_mie_.GetAn();
  std::vector<std::complex<double> > bn = multi_layer_mie_.GetBn();
  snprintf(plot_name, 300,
           "%s-epsilon",  sign_.c_str());
  full_sign_ = std::string(plot_name);
  wrapper.SetPlotName(plot_name);
  wrapper.SetXLabelName("Total_R");
  wrapper.SetYLabelName("Epsilon");
  wrapper.SetDrawStyle("w l lw 2");
  wrapper.SetXRange({spectra.front()[0], spectra.back()[0]});
  for (auto multi_point : spectra) wrapper.AddMultiPoint(multi_point);
  wrapper.AddColumnName("WL");
  wrapper.AddColumnName("Re");
  wrapper.AddColumnName("Im");
  wrapper.MakeOutput();
}
// ********************************************************************** //
// ********************************************************************** //
// ********************************************************************** //
