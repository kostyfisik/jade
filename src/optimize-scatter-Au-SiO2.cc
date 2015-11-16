/**
 * @file   optimize-scatter-Au-SiO2.cc
 * @author Konstantin Ladutenko <kostyfisik at gmail (.) com>
 * @date   Sun Nov 15 17:10:31 2015
 * 
 * @brief  
**/
/// @copyright 2015  Konstantin Ladutenko
///
/// optimize-scatter-Au-SiO2 is free software: you can redistribute it and/or modify
/// it under the terms of the GNU General Public License as published by
/// the Free Software Foundation, either version 3 of the License, or
/// (at your option) any later version.
///
/// optimize-scatter-Au-SiO2 is distributed in the hope that it will be useful,
/// but WITHOUT ANY WARRANTY; without even the implied warranty of
/// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
/// GNU General Public License for more details.
///
/// You should have received a copy of the GNU General Public License
/// along with optimize-scatter-Au-SiO2.  If not, see <http://www.gnu.org/licenses/>.
///
/// optimize-scatter-Au-SiO2 uses nmie.cc from scattnlay by Ovidio Pena
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
#include "./nmie/nmie-applied.h"
#include "./read-spectra/read-spectra.h"
const double pi=3.14159265358979323846;
template<class T> inline T pow2(const T value) {return value*value;}
// Mie model. Used in fitness function for optimization of
// sub_population.
void SetOptimizer();

double EvaluateFitness(std::vector<double> input);
std::vector< std::vector<double> > EvaluateSpectraForBestDesign();
std::vector< std::vector<double> > EvaluateSpectraForChannels(std::vector<double>& best_x,
							      double total_r);
void Print();
void PrintCoating(std::vector<double> current, double initial_RCS,
                    jade::SubPopulation sub_population);
void PrintGnuPlotSpectra(std::vector< std::vector<double> > spectra);
void PrintGnuPlotChannels(std::vector< std::vector<double> > spectra);
void PrintGnuPlotChannelSweep(std::vector< std::vector<double> > spectra);
void ReadDesignSpectra();
std::vector<double> ShareToWidth(double total_r, std::vector<double> share);

// ********************************************************************** //
// ********************************************************************** //
// ********************************************************************** //
nmie::MultiLayerMieApplied multi_layer_mie_;  
jade::SubPopulation sub_population_;  // Optimizer of parameters for Mie model.
// You should provide a Material.txt file for each "Material" in index desing;
// File format:
// comment line starts with #,
// data line is tab separated:  WL,nm  re(epsilon)  im(epsilon)
std::vector< std::string> index_design_ = {"Au", "SiO2", "Au"};
std::vector< read_spectra::ReadSpectra > index_spectra_, index_plot_spectra_;
std::string sign_, full_sign_;
// ********************************************************************** //
// ********************************************************************** //
// ********************************************************************** //
int fails_ = 0;
double core_share_ = 0.0, mid_share_ = 0.0;
double Q_=0.0, initial_Q_=0.0;
double total_r_ = 0.0; 
// ********************************************************************** //
// Set model: core->mid->shell
const double max_r_ = 100.0; // nm
const double max_mid_width_ = max_r_; // nm
//const double max_mid_width_ = 10; // nm
// Set dispersion
double at_wl_ = 800.0;
double from_wl_ = at_wl_, to_wl_ = at_wl_;
int samples_ = 1;
// double from_wl_ = 300.0, to_wl_ = 900.0;
// int samples_ = 151;
double plot_from_wl_ = at_wl_-100.0, plot_to_wl_ = at_wl_+100.0;
int plot_samples_ = 501;
bool isQsca = true;
//bool isQsca = false;
// Set optimizer
bool isFindMax = true;
//bool isFindMax = false;
int total_generations_ = 150;
int population_multiplicator_ = 160;
int dim_=3;
double step_r_ = 15.0; //max_r_ / 159.0;
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
    ReadDesignSpectra();
    std::cout << "Sign: " << sign_ << std::endl;

    SetOptimizer();
    if (step_r_ <=0.0) throw std::invalid_argument("Radius step should be positive!/n");
    auto best_x = sub_population_.GetBest(&Q_);
    double best_Q = 0.0, best_total_r = 0.0;
    std::vector< std::vector<double> > channel_sweep;
    int least_size_sweep = 15;
    // ***************************************************
    // **************  Main loop   ***********************
    // ***************************************************  
    //for (total_r_ = 62.0; total_r_ < 65.0; total_r_+=step_r_) {
    for (total_r_ = step_r_; total_r_ < max_r_*1.00001; total_r_+=step_r_) {
      // if (
      // 	  (total_r_ > 80-0.01 && total_r_ < 83-0.01)
      // 	  ||
      // 	  (total_r_ > 64-0.01 && total_r_ < 66-0.01)
      // 	  ||
      // 	  (total_r_ > 55-0.01 && total_r_ < 57-0.01)
      // 	  )
      // 	  step_r_ = 0.1;
      // else step_r_ = 1.0;
      //for (total_r_ = 145.0; total_r_ < 147.0; total_r_+=0.05) {
      if (rank == 0) printf("\nTotal R = %g\n", total_r_);    
      sub_population_.RunOptimization();
      // Plot spectra from each process
      PrintGnuPlotSpectra(EvaluateSpectraForBestDesign());
      if (rank == 0) Print();
      auto best_local_x = sub_population_.GetBest(&Q_);
      if (rank == 0) {
	std::vector<std::complex<double> > an = multi_layer_mie_.GetAn();
	std::vector<std::complex<double> > bn = multi_layer_mie_.GetBn();
	for (int i = 0; i < 3; ++i)
	  printf("an[%d]=%g, bn[%d]=%g\n",i,std::abs(an[i]),i, std::abs(bn[i]));
      }
      if (Q_ > best_Q) {
	best_Q = Q_;
	best_total_r = total_r_;
	best_x = best_local_x;	
      }
      PrintGnuPlotChannels(EvaluateSpectraForChannels(best_local_x, total_r_));
      sub_population_.FitnessFunction(best_local_x);
      std::vector<double> tmp({total_r_});
      //std::vector<double> channels(multi_layer_mie_.GetQabs_channel());
      std::vector<std::complex<double> > an = multi_layer_mie_.GetAn();
      std::vector<std::complex<double> > bn = multi_layer_mie_.GetBn();
      for (int i = 0; i < an.size(); ++i) {
	if (isQsca) {
	  tmp.push_back(pow2(std::abs(an[i])));
	  tmp.push_back(pow2(std::abs(bn[i])));	  
	}else{
	  tmp.push_back(an[i].real() - pow2(std::abs(an[i])));
	  tmp.push_back(bn[i].real() - pow2(std::abs(bn[i])));
	}
      }
      channel_sweep.push_back(tmp);
      if (least_size_sweep > tmp.size()) least_size_sweep = tmp.size();
    }  // end of total coating thickness sweep
    PrintGnuPlotChannels(EvaluateSpectraForChannels(best_x, best_total_r));
    for (auto& row : channel_sweep) row.resize(least_size_sweep);
    PrintGnuPlotChannelSweep(channel_sweep);
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
  // mid_share_ and core_share_ are set in EvaluateSpectraForBestDesign()
  const double mid_width = (total_r_>max_mid_width_ ? max_mid_width_ : total_r_) * mid_share_;
  const double core_width = (total_r_ - mid_width) * core_share_;
  const double shell_width = total_r_ - core_width - mid_width;
  printf("core_width:%.19g\nfirst_shell_width:%.19g\nouter_shell_width:%.19g\n",
	 core_width, mid_width, shell_width);
  
}
// ********************************************************************** //
// ********************************************************************** //
// ********************************************************************** //
void SetOptimizer() {
  //The last layer fills to total_r_
  long dimension = index_design_.size() - 1;
  //sub_population_.FitnessFunction = &EvaluateFitness;
  //sub_population_.FitnessFunction = &EvaluateFitnessChannel;
  long total_population = dimension * population_multiplicator_;
  sub_population_.Init(total_population, dimension);
  /// Low and upper bound for all dimenstions;
  sub_population_.SetAllBounds(0.0, 1.0);
  if (isFindMax) sub_population_.SetTargetToMaximum();
  else  sub_population_.SetTargetToMinimum();
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
  std::vector<double> width = ShareToWidth(total_r_, input);

  double Q = 0.0;
  for (int i=0; i < index_spectra_[0].GetIndex().size(); ++i) {
    multi_layer_mie_.ClearTarget();
    const auto core_data = index_spectra_[0].GetIndex();
    const double& wl = core_data[i].first;
    multi_layer_mie_.SetWavelength(wl);
    for (int l = 0; l < width.size(); ++l) {
      if (width[l] > 0.0) {
    	const auto index_data = index_spectra_[l].GetIndex();
    	const std::complex<double>& index = index_data[i].second;	
    	multi_layer_mie_.AddTargetLayer(width[l], index);
      }
    }  // end of for each layer
    try {
      multi_layer_mie_.RunMieCalculation();
      if (isQsca)  Q = multi_layer_mie_.GetQsca();
      else Q = multi_layer_mie_.GetQabs();
    } catch( const std::invalid_argument& ia ) {
      printf(".");
      sub_population_.GetWorst(&Q);
    }
  }  // end of for all points of the spectrum
  Q_ = Q;
  return Q_;
}
// ********************************************************************** //
// ********************************************************************** //
// ********************************************************************** //
std::vector< std::vector<double> > EvaluateSpectraForBestDesign() {
  fails_ = 0;
  auto best_x = sub_population_.GetBest(&Q_);
  // Setting Mie model to the best state.
  if (best_x.size() != 2) throw std::invalid_argument("Wrong input dimension!/n");
  core_share_ = best_x[0];
  mid_share_ = best_x[1];
  const double mid_width = (total_r_>max_mid_width_ ? max_mid_width_ : total_r_) * mid_share_;
  const double core_width = (total_r_ - mid_width) * core_share_;
  const double shell_width = total_r_ - core_width - mid_width + eps_;
  if (mid_width <= 0.0 || core_width <= 0.0 || shell_width <= 0.0) {
    printf("core_share = %g\n",core_share_);
    printf("mid_share = %g\n",mid_share_);
    printf("total_r_ = %g\n",total_r_);
    
    printf("mid <=0:   %g\n", mid_width);
    printf("core <=0:   %g\n",core_width);
    printf("shell <=0:   %g\n",shell_width);
  }

  auto core_data = index_plot_spectra_[0].GetIndex();
  auto mid_data =  index_plot_spectra_[1].GetIndex();

  double Qabs = 0.0, Qext=0.0, Qsca=0.0, Qbk =0.0;
  std::vector< std::vector<double> > spectra;
  for (int i=0; i < core_data.size(); ++i) {
    const double& wl = core_data[i].first;
    const std::complex<double>& core = core_data[i].second;
    const std::complex<double>& mid = mid_data[i].second;
    const std::complex<double>& shell = core;
    multi_layer_mie_.ClearTarget();
    multi_layer_mie_.AddTargetLayer(core_width, core);
    multi_layer_mie_.AddTargetLayer(mid_width, mid);
    multi_layer_mie_.AddTargetLayer(shell_width, shell);
    multi_layer_mie_.SetWavelength(wl);
    try {
      multi_layer_mie_.RunMieCalculation();
      Qabs = multi_layer_mie_.GetQabs();
      Qext = multi_layer_mie_.GetQext();
      Qsca = multi_layer_mie_.GetQsca();
      Qbk = multi_layer_mie_.GetQbk();
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
std::vector< std::vector<double> > EvaluateSpectraForChannels 
(std::vector<double>& best_x, double best_total_r) {
  fails_ = 0;
  // Setting Mie model to the best state.
  total_r_ = best_total_r;
  Q_ = sub_population_.FitnessFunction(best_x);
  if (best_x.size() != 2) throw std::invalid_argument("Wrong input dimension!/n");
  core_share_ = best_x[0];
  mid_share_ = best_x[1];
  const double mid_width = (total_r_>max_mid_width_ ? max_mid_width_ : total_r_) * mid_share_;
  const double core_width = (total_r_ - mid_width) * core_share_;
  const double shell_width = total_r_ - core_width - mid_width + eps_;
  if (mid_width <= 0.0 || core_width <= 0.0 || shell_width <= 0.0) {
    printf("core_share = %g\n",core_share_);
    printf("mid_share = %g\n",mid_share_);
    printf("total_r_ = %g\n",total_r_);
    
    printf("mid <=0:   %g\n", mid_width);
    printf("core <=0:   %g\n",core_width);
    printf("shell <=0:   %g\n",shell_width);
  }

  auto core_data = index_plot_spectra_[0].GetIndex();
  auto mid_data =  index_plot_spectra_[1].GetIndex();
  //double max_Qabs = 0.0, Qabs = 0.0, Qext=0.0, Qsca=0.0, Qbk =0.0;
  std::vector< std::vector<double> > spectra;
  int least_size = 15;
  // Scan all wavelengths
  for (int i=0; i < core_data.size(); ++i) {
    const double& wl = core_data[i].first;
    const std::complex<double>& core = core_data[i].second;
    const std::complex<double>& mid = mid_data[i].second;
    const std::complex<double>& shell = core;
    multi_layer_mie_.ClearTarget();
    multi_layer_mie_.AddTargetLayer(core_width, core);
    multi_layer_mie_.AddTargetLayer(mid_width, mid);
    multi_layer_mie_.AddTargetLayer(shell_width, shell);
    multi_layer_mie_.SetWavelength(wl);
    try {
      multi_layer_mie_.RunMieCalculation();
      std::vector<double> tmp({wl});
      //std::vector<double> channels(multi_layer_mie_.GetQabs_channel());
      std::vector<std::complex<double> > an = multi_layer_mie_.GetAn();
      std::vector<std::complex<double> > bn = multi_layer_mie_.GetBn();
      for (int i = 0; i < an.size(); ++i) {
	if (isQsca) {
	  tmp.push_back(pow2(std::abs(an[i])));
	  tmp.push_back(pow2(std::abs(bn[i])));
	}else{
	  tmp.push_back(an[i].real() - pow2(std::abs(an[i])));
	  tmp.push_back(bn[i].real() - pow2(std::abs(bn[i])));
	}
      }
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
void PrintGnuPlotSpectra(std::vector< std::vector<double> > spectra) {
  auto best_x = sub_population_.GetBest(&Q_);
  sub_population_.FitnessFunction(best_x);
  if (isQsca) Q_=multi_layer_mie_.GetQsca();
  else Q_=multi_layer_mie_.GetQabs();
  gnuplot::GnuplotWrapper wrapper;
  char plot_name [300];
  const double mid_width = (total_r_>max_mid_width_ ? max_mid_width_ : total_r_) * mid_share_;
  const double core_width = (total_r_ - mid_width) * core_share_;
  const double shell_width = total_r_ - core_width - mid_width;  
  if (isQsca) snprintf(plot_name, 300,
           "%s-TotalR%06.2fnm-Qsca%016.13f--core%07.2fnm--inshell%07.2fnm--outshell%07.2fnm-fails%d-spectra", sign_.c_str(), total_r_, Q_,  core_width, mid_width, shell_width, fails_);
  else snprintf(plot_name, 300,
           "%s-TotalR%06.2fnm-Qabs%016.13f--core%07.2fnm--inshell%07.2fnm--outshell%07.2fnm-fails%d-spectra", sign_.c_str(), total_r_, Q_,  core_width, mid_width, shell_width, fails_);
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
  const double mid_width = (total_r_>max_mid_width_ ? max_mid_width_ : total_r_) * mid_share_;
  const double core_width = (total_r_ - mid_width) * core_share_;
  const double shell_width = total_r_ - core_width - mid_width;  
  if (isQsca) snprintf(plot_name, 300,
           "o-spectra-%s-channels-TotalR%06.2fnm-Qsca%016.13f--core%07.2fnm--inshell%07.2fnm--outshell%07.2fnm-fails%d", sign_.c_str(),
           total_r_, Q_,  core_width, mid_width, shell_width, fails_);
  else snprintf(plot_name, 300,
           "o-spectra-%s-channels-TotalR%06.2fnm-Qabs%016.13f--core%07.2fnm--inshell%07.2fnm--outshell%07.2fnm-fails%d", sign_.c_str(),
           total_r_, Q_,  core_width, mid_width, shell_width, fails_);
  full_sign_ = std::string(plot_name);
  wrapper.SetPlotName(plot_name);
  wrapper.SetXLabelName("WL");
  if (isQsca)   wrapper.SetYLabelName("NRCS");
  else  wrapper.SetYLabelName("NACS");
  wrapper.SetDrawStyle("w l lw 2");
  wrapper.SetXRange({spectra.front()[0], spectra.back()[0]});
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
void PrintGnuPlotChannelSweep(std::vector< std::vector<double> > spectra) {
  gnuplot::GnuplotWrapper wrapper;
  char plot_name [300];
  snprintf(plot_name, 300,
           "sweep-ch");
  full_sign_ = std::string(plot_name);
  wrapper.SetPlotName(plot_name);
  wrapper.SetXLabelName("Total_R");
  if (isQsca) wrapper.SetYLabelName("NRCS");
  else wrapper.SetYLabelName("NACS");
  wrapper.SetDrawStyle("w l lw 2");
  wrapper.SetXRange({spectra.front()[0], spectra.back()[0]});
  //wrapper.SetYRange({-0.02, 0.27});
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
void ReadDesignSpectra() {
  sign_ = "";
  for (auto name : index_design_) {
    sign_ += name+"-";
    read_spectra::ReadSpectra tmp;
    tmp.ReadFromFile(name+".txt");
    tmp.ResizeToComplex(from_wl_, to_wl_, samples_).ToIndex();
    index_spectra_.push_back(tmp);
    tmp.ResizeToComplex(plot_from_wl_, plot_to_wl_, plot_samples_).ToIndex();
    index_plot_spectra_.push_back(tmp);
  }
  sign_.pop_back();  
}
// ********************************************************************** //
// ********************************************************************** //
// ********************************************************************** //
std::vector<double> ShareToWidth(double total_r, std::vector<double> share) {
  std::vector<double> width;
  double remainder = total_r;
  for (auto layer : share) {
    if (layer < eps_) layer = 0.0;
    if (layer > 1.0 - eps_) layer = 1.0;
    double current_width = remainder*layer;
    width.push_back(current_width);
    remainder -= current_width;
    // To avoid roundoff errors with outer layers width is set to be zero.
    if (remainder < 0.0) remainder = 0.0;
  }
  width.push_back(remainder);
  return width;
}
// ********************************************************************** //
// ********************************************************************** //
// ********************************************************************** //

