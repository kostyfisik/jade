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
std::vector< std::vector<double> > EvaluateSpectraForChannels(std::vector<double>& best_x,
							      double total_r);
void Print();
void PrintCoating(std::vector<double> current, double initial_RCS,
                    jade::SubPopulation sub_population);
void PrintGnuPlotSpectra(std::vector< std::vector<double> > spectra);
void PrintGnuPlotChannels(std::vector< std::vector<double> > spectra);
// ********************************************************************** //
// ********************************************************************** //
// ********************************************************************** //
nmie::MultiLayerMie multi_layer_mie_;  
jade::SubPopulation sub_population_;  // Optimizer of parameters for Mie model.
read_spectra::ReadSpectra core_index_, TiN_;
read_spectra::ReadSpectra plot_core_index_, plot_TiN_;
std::string sign_, full_sign_;
// ********************************************************************** //
// ********************************************************************** //
// ********************************************************************** //
double lambda_best_ = 0.0; // nm
int fails_ = 0;
double core_share_ = 0.0, TiN_share_ = 0.0;
double Qabs_=0.0, initial_Qabs_=0.0;
double total_r_ = 0.0; 
// ********************************************************************** //
// Set model: core->TiN->shell
const double max_r_ = 159.0; // nm
const double max_TiN_width_ = max_r_; // nm
//const double max_TiN_width_ = 10; // nm
// Set dispersion
double at_wl_ = 500.0;
double from_wl_ = at_wl_, to_wl_ = at_wl_;
int samples_ = 1;
// double from_wl_ = 300.0, to_wl_ = 900.0;
// int samples_ = 151;
double plot_from_wl_ = 300.0, plot_to_wl_ = 900.0;
int plot_samples_ = 1501;
//bool isGaAs = false; // Select Si of GaAs as a material for core and shell
bool isGaAs = true;
// Set optimizer
int total_generations_ = 150;
int population_multiplicator_ = 160;
double step_r_ = 1.0; //max_r_ / 159.0;
// ********************************************************************** //
// ********************************************************************** //
// ********************************************************************** //
const double eps_=1e-11;
int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  try {
    sub_population_.FitnessFunction = &EvaluateFitness;
    //sub_population_.FitnessFunction = &EvaluateFitnessChannel;
    //std::string core_filename("GaAs.txt");
    std::string core_filename("Si.txt");
    //std::string core_filename("Ag.txt");
    //std::string TiN_filename("TiN.txt");
    std::string TiN_filename("Ag.txt");
    //std::string TiN_filename("Si.txt");
    std::string shell_filename(core_filename);
    sign_ = core_filename.substr(0, core_filename.find("."))+"-"+
      TiN_filename.substr(0, core_filename.find("."))+"-"+
      shell_filename.substr(0, core_filename.find("."));
    std::cout << "Sign: " << sign_ << std::endl;
    core_index_.ReadFromFile(core_filename);
    plot_core_index_.ReadFromFile(core_filename);
    core_index_.ResizeToComplex(from_wl_, to_wl_, samples_).ToIndex();
    plot_core_index_.ResizeToComplex(plot_from_wl_, plot_to_wl_, plot_samples_)
      .ToIndex();
    TiN_.ReadFromFile(TiN_filename).ResizeToComplex(from_wl_, to_wl_, samples_)
      .ToIndex();
    plot_TiN_.ReadFromFile(TiN_filename)
      .ResizeToComplex(plot_from_wl_, plot_to_wl_, plot_samples_).ToIndex();
    if (core_index_.GetIndex().size()
	!= TiN_.GetIndex().size()) throw std::invalid_argument("Unexpected sampling of dispersion!/n");
    // if (rank == 0) {  printf("\ncore:\n");  core_index_.PrintData(); printf("\nTiN:\n");   TiN_.PrintData(); }
    SetOptimizer();
    if (step_r_ <=0.0) throw std::invalid_argument("Radius step should be positive!/n");
    auto best_x = sub_population_.GetBest(&Qabs_);
    double best_Qabs = 0.0, best_total_r = 0.0;
    // ***************************************************
    // **************  Main loop   ***********************
    // ***************************************************  
    for (total_r_ = step_r_; total_r_ < max_r_*1.00001; total_r_+=step_r_) {
      //for (total_r_ = 145.0; total_r_ < 147.0; total_r_+=0.05) {
      if (rank == 0) printf("\nTotal R = %g\n", total_r_);    
      sub_population_.RunOptimization();
      // Plot spectra from each process
      PrintGnuPlotSpectra(EvaluateSpectraForBestDesign());
      if (rank == 0) Print();
      auto best_local_x = sub_population_.GetBest(&Qabs_);
      if (rank == 0) {
	std::vector<std::complex<double> > an = multi_layer_mie_.GetAn();
	std::vector<std::complex<double> > bn = multi_layer_mie_.GetBn();
	for (int i = 0; i < 3; ++i)
	  printf("an[%d]=%g, bn[%d]=%g\n",i,std::abs(an[i]),i, std::abs(bn[i]));
      }
      if (Qabs_ > best_Qabs) {
	best_Qabs = Qabs_;
	best_total_r = total_r_;
	best_x = best_local_x;	
      }
      PrintGnuPlotChannels(EvaluateSpectraForChannels(best_local_x, total_r_));
    }  // end of total coating thickness sweep
    PrintGnuPlotChannels(EvaluateSpectraForChannels(best_x, best_total_r));
    if (rank == 0) {
      printf("==========best=========\n");
      for (auto x : best_x) printf("besttt %g,  ", x);
      printf("Q=%g, bbest Q = %g\n", EvaluateFitness(best_x), best_Qabs);
    }
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
  const double TiN_width = (total_r_>max_TiN_width_ ? max_TiN_width_ : total_r_) * TiN_share_;
  const double core_width = (total_r_ - TiN_width) * core_share_;
  const double shell_width = total_r_ - core_width - TiN_width;
  printf("core_width:%.19g\nfirst_shell_width:%.19g\nouter_shell_width:%.19g\n",
	 core_width, TiN_width, shell_width);
  
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
  sub_population_.SetAllBounds(eps_, 1.0-eps_);
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
  core_share_ = input[0];
  TiN_share_ = input[1];
  const double TiN_width = (total_r_>max_TiN_width_ ? max_TiN_width_ : total_r_) * TiN_share_;
  const double core_width = (total_r_ - TiN_width) * core_share_;
  const double shell_width = total_r_ - core_width - TiN_width;
  if (TiN_width <= 0.0 || core_width <= 0.0 || shell_width <= 0.0) {
    printf("core_share = %g\n",core_share_);
    printf("TiN_share = %g\n",TiN_share_);
    printf("total_r_ = %g\n",total_r_);
    
    printf("TiN <=0:   %g\n", TiN_width);
    printf("core <=0:   %g\n",core_width);
    printf("shell <=0:   %g\n",shell_width);
  }
  
  auto core_data = core_index_.GetIndex();
  auto TiN_data = TiN_.GetIndex();
  double max_Qabs = 0.0, Qabs = 0.0;
  for (int i=0; i < core_data.size(); ++i) {
    const double& wl = core_data[i].first;
    const std::complex<double>& core = core_data[i].second;
    const std::complex<double>& TiN = TiN_data[i].second;
    const std::complex<double>& shell = core;
    multi_layer_mie_.ClearTarget();
    //    double min_share = 0.00001;
    //    if (core_share_ > min_share) 
      multi_layer_mie_.AddTargetLayer(core_width, core);
    //    if (TiN_share_ > min_share) 
      multi_layer_mie_.AddTargetLayer(TiN_width, TiN);
    //    if (shell_width/total_r_ > min_share)
      multi_layer_mie_.AddTargetLayer(shell_width, shell);
    multi_layer_mie_.SetWavelength(wl);
    try {
      multi_layer_mie_.RunMieCalculations();
      Qabs = multi_layer_mie_.GetQabs();
    } catch( const std::invalid_argument& ia ) {
      printf(".");
      sub_population_.GetWorst(&Qabs);
    }
    if (Qabs > max_Qabs) max_Qabs = Qabs;
  }  // end of for all points of the spectrum
  Qabs_ = Qabs;
  return Qabs_;
}
// ********************************************************************** //
// ********************************************************************** //
// ********************************************************************** //
double EvaluateFitnessChannel(std::vector<double> input) {
  if (input.size() != 2) throw std::invalid_argument("Wrong input dimension!/n");
  core_share_ = input[0];
  TiN_share_ = input[1];
  const double TiN_width = (total_r_>max_TiN_width_ ? max_TiN_width_ : total_r_) * TiN_share_;
  const double core_width = (total_r_ - TiN_width) * core_share_;
  const double shell_width = total_r_ - core_width - TiN_width;
  if (TiN_width <= 0.0 || core_width <= 0.0 || shell_width <= 0.0) {
    printf("core_share = %g\n",core_share_);
    printf("TiN_share = %g\n",TiN_share_);
    printf("total_r_ = %g\n",total_r_);
    
    printf("TiN <=0:   %g\n", TiN_width);
    printf("core <=0:   %g\n",core_width);
    printf("shell <=0:   %g\n",shell_width);
  }
  
  auto core_data = core_index_.GetIndex();
  auto TiN_data = TiN_.GetIndex();
  double max_Qabs = 0.0, Qabs = 0.0;
  for (int i=0; i < core_data.size(); ++i) {
    const double& wl = core_data[i].first;
    const std::complex<double>& core = core_data[i].second;
    const std::complex<double>& TiN = TiN_data[i].second;
    const std::complex<double>& shell = core;
    multi_layer_mie_.ClearTarget();
    multi_layer_mie_.AddTargetLayer(core_width, core);
    multi_layer_mie_.AddTargetLayer(TiN_width, TiN);
    multi_layer_mie_.AddTargetLayer(shell_width, shell);
    multi_layer_mie_.SetWavelength(wl);
    try {
      multi_layer_mie_.RunMieCalculations();
      //Qabs = multi_layer_mie_.GetQsca();
      std::vector<double> channels(multi_layer_mie_.GetQabs_channel_normalized());
      //if (channels.size() > 2) Qabs = channels[0];
      //if (channels.size() > 2) Qabs = channels[0]+channels[1];
      //if (channels.size() > 4) Qabs = channels[0]+channels[1]+channels[2];
      if (channels.size() > 4) Qabs = channels[0]*channels[1]*channels[2];
    } catch( const std::invalid_argument& ia ) {
      printf(".");
    }
    if (Qabs > max_Qabs) max_Qabs = Qabs;
  }  // end of for all points of the spectrum
  Qabs_ = Qabs;
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
  core_share_ = best_x[0];
  TiN_share_ = best_x[1];
  const double TiN_width = (total_r_>max_TiN_width_ ? max_TiN_width_ : total_r_) * TiN_share_;
  const double core_width = (total_r_ - TiN_width) * core_share_;
  const double shell_width = total_r_ - core_width - TiN_width;
  if (TiN_width <= 0.0 || core_width <= 0.0 || shell_width <= 0.0) {
    printf("core_share = %g\n",core_share_);
    printf("TiN_share = %g\n",TiN_share_);
    printf("total_r_ = %g\n",total_r_);
    
    printf("TiN <=0:   %g\n", TiN_width);
    printf("core <=0:   %g\n",core_width);
    printf("shell <=0:   %g\n",shell_width);
  }
  
  auto core_data = plot_core_index_.GetIndex();
  auto TiN_data = plot_TiN_.GetIndex();
  double max_Qabs = 0.0, Qabs = 0.0, Qext=0.0, Qsca=0.0, Qbk =0.0;
  std::vector< std::vector<double> > spectra;
  for (int i=0; i < core_data.size(); ++i) {
    const double& wl = core_data[i].first;
    const std::complex<double>& core = core_data[i].second;
    const std::complex<double>& TiN = TiN_data[i].second;
    const std::complex<double>& shell = core;
    multi_layer_mie_.ClearTarget();
    multi_layer_mie_.AddTargetLayer(core_width, core);
    multi_layer_mie_.AddTargetLayer(TiN_width, TiN);
    multi_layer_mie_.AddTargetLayer(shell_width, shell);
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
    if (Qabs > max_Qabs) max_Qabs = Qabs;
  }  // end of for all points of the spectrum
  return spectra;
}
// ********************************************************************** //
// ********************************************************************** //
// ********************************************************************** //
std::vector< std::vector<double> > EvaluateSpectraForChannels 
(std::vector<double>& best_x, double best_total_r) {
  fails_ = 0;
  // Setting Mie model to the best state.
  total_r_ = best_total_r;
  Qabs_ = EvaluateFitness(best_x);
  if (best_x.size() != 2) throw std::invalid_argument("Wrong input dimension!/n");
  core_share_ = best_x[0];
  TiN_share_ = best_x[1];
  const double TiN_width = (total_r_>max_TiN_width_ ? max_TiN_width_ : total_r_) * TiN_share_;
  const double core_width = (total_r_ - TiN_width) * core_share_;
  const double shell_width = total_r_ - core_width - TiN_width;
  if (TiN_width <= 0.0 || core_width <= 0.0 || shell_width <= 0.0) {
    printf("core_share = %g\n",core_share_);
    printf("TiN_share = %g\n",TiN_share_);
    printf("total_r_ = %g\n",total_r_);
    
    printf("TiN <=0:   %g\n", TiN_width);
    printf("core <=0:   %g\n",core_width);
    printf("shell <=0:   %g\n",shell_width);
  }
  
  auto core_data = plot_core_index_.GetIndex();
  auto TiN_data = plot_TiN_.GetIndex();
  //double max_Qabs = 0.0, Qabs = 0.0, Qext=0.0, Qsca=0.0, Qbk =0.0;
  std::vector< std::vector<double> > spectra;
  int least_size = 10000;
  // Scan all wavelengths
  for (int i=0; i < core_data.size(); ++i) {
    const double& wl = core_data[i].first;
    const std::complex<double>& core = core_data[i].second;
    const std::complex<double>& TiN = TiN_data[i].second;
    const std::complex<double>& shell = core;
    multi_layer_mie_.ClearTarget();
    multi_layer_mie_.AddTargetLayer(core_width, core);
    multi_layer_mie_.AddTargetLayer(TiN_width, TiN);
    multi_layer_mie_.AddTargetLayer(shell_width, shell);
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
  // for (auto row : spectra) {
  //   printf(" WL=");
  //   for (double value:row) printf(" %g,\t",value);
  //   printf("\n\n");
  // }
  return spectra;
}
// ********************************************************************** //
// ********************************************************************** //
// ********************************************************************** //
void PrintCoating(std::vector<double> current, double initial_RCS,
                    jade::SubPopulation sub_population) {
  //double best_RCS = 0.0;
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
  const double TiN_width = (total_r_>max_TiN_width_ ? max_TiN_width_ : total_r_) * TiN_share_;
  const double core_width = (total_r_ - TiN_width) * core_share_;
  const double shell_width = total_r_ - core_width - TiN_width;  
  snprintf(plot_name, 300,
           "%s-TotalR%06.2fnm-Qabs%016.13f--core%07.2fnm--inshell%07.2fnm--outshell%07.2fnm-fails%d-spectra", sign_.c_str(),
           total_r_, Qabs_,  core_width, TiN_width, shell_width, fails_);
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
  const double TiN_width = (total_r_>max_TiN_width_ ? max_TiN_width_ : total_r_) * TiN_share_;
  const double core_width = (total_r_ - TiN_width) * core_share_;
  const double shell_width = total_r_ - core_width - TiN_width;  
  snprintf(plot_name, 300,
           "o-spectra-%s-channels-TotalR%06.2fnm-Qabs%016.13f--core%07.2fnm--inshell%07.2fnm--outshell%07.2fnm-fails%d", sign_.c_str(),
           total_r_, Qabs_,  core_width, TiN_width, shell_width, fails_);
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

