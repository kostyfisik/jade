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
#include <complex>
#include "./jade.h"
#include <cstdlib>
#include <cstdio>
#include <ctime>
#include <stdexcept>
#include <string>
#include "./gnuplot-wrapper/gnuplot-wrapper.h"
#include "./nmie/nmie-wrapper.h"
//#include "./read-spectra/read-spectra.h"
const double pi=3.14159265358979323846;
const double speed_of_light = 299792458;
template<class T> inline T pow2(const T value) {return value*value;}
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
jade::SubPopulation sub_population_;  // Optimizer of parameters for Mie model.
// ********************************************************************** //
// ********************************************************************** //
// ********************************************************************** //
nmie::MultiLayerMie multi_layer_mie_;  
double w2l( double l) {return 2.0 * pi * speed_of_light/l;};
double l2w( double w) {return 2.0 * pi * speed_of_light/w;};
// ********************************************************************** //
// ********************************************************************** //
// ********************************************************************** //
const double omega_p_ = 1.0e9;
const double lambda_p_ = w2l(omega_p_);
double r1_ = 0.4749*lambda_p_;
double r2_ = 0.6404*lambda_p_;
double r3_ = 0.8249*lambda_p_;
double core_width_ = r1_;
double inshell_width_ = r2_ - r1_;
double outshell_width_ = r3_ - r2_;
double from_omega_ = 0.288*omega_p_;
double to_omega_ = 0.3*omega_p_;
double epsilon_d_ = 12.96;
double inshell_index_ = std::sqrt(epsilon_d_);
double gamma_d_ = 0.0;
double gamma_bulk_ = 0.002*omega_p_;
bool isLossy_ = false;

std::complex<double> epsilon_m(double omega) {
  // std::cout << "omega:" << omega << "  omega_p:"<<omega_p_ <<std::endl;
  return 1.0 - pow2(omega_p_)/(pow2(omega)+std::complex<double>(0,1)*gamma_d_*omega);
}
// Set dispersion
double from_wl_ = w2l(from_omega_);
double to_wl_ = w2l(to_omega_);
int samples_ = 300;
double plot_from_wl_ = from_wl_, plot_to_wl_ = to_wl_;
int plot_samples_ = samples_;
// Set optimizer
int total_generations_ = 50;
int population_multiplicator_ = 16;
// ********************************************************************** //
// ********************************************************************** //
// ********************************************************************** //
const double eps_=1e-11;
int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  try {
    const double lossless = 0.0;
    int least_size = 10000;
    std::vector< std::vector<double> > spectra;
    if (from_omega_ > to_omega_) throw std::invalid_argument("Wrong omega range!");
    SetOptimizer();
    double omega_step = (to_omega_ - from_omega_)/samples_;
    for (double omega = from_omega_; omega < to_omega_; omega += omega_step) {
      //double metal_index_re = std::sqrt(epsilon_m(omega));
      multi_layer_mie_.ClearTarget();
      double A = 1.0;
      double V_f = 7.37e-4*lambda_p_*omega_p_;
      double l_r = r1_;
      if (isLossy_) gamma_d_ = gamma_bulk_ + A*V_f/l_r;
      multi_layer_mie_.AddTargetLayer(core_width_, std::sqrt(epsilon_m(omega)));
      multi_layer_mie_.AddTargetLayer(inshell_width_, {inshell_index_, lossless});
      l_r = r3_-r2_;
      if (isLossy_) gamma_d_ = gamma_bulk_ + A*V_f/l_r;
      multi_layer_mie_.AddTargetLayer(outshell_width_, std::sqrt(epsilon_m(omega)));
      multi_layer_mie_.SetWavelength(w2l(omega));
      try {
	multi_layer_mie_.RunMieCalculations();
	double Qsca = multi_layer_mie_.GetQsca();
	double norm = (Qsca * pi*pow2(r3_)) / ( pow2(w2l(omega)) / (2.0*pi) );
	std::vector<double> tmp({omega/omega_p_, norm});
	std::vector<double> channels(multi_layer_mie_.GetQsca_channel());
	//std::vector<double> channels(multi_layer_mie_.GetQsca_channel_normalized());
	//std::vector<double> channels(multi_layer_mie_.GetQabs_channel_normalized());
	//std::vector<double> channels(multi_layer_mie_.GetQabs_channel());
	tmp.insert(tmp.end(), channels.begin(), channels.end());
	spectra.push_back(tmp);
	if (least_size > tmp.size()) least_size = tmp.size();
      } catch( const std::invalid_argument& ia ) {
	std::cerr << "Invalid argument: " << ia.what() << std::endl;
	printf(".");
      }
    }  // end of omega sweep
    //PrintGnuPlotSpectra(spectra);
    for (auto& row : spectra) {
      //if (rank==0) std::cout<< least_size<< " -- "<<row.size()<<std::endl;
      row.resize(least_size);
    }
   
    PrintGnuPlotChannels(spectra);
  } catch( const std::invalid_argument& ia ) {
    // Will catch if  multi_layer_mie_ fails or other errors.
    std::cerr << "Invalid argument: " << ia.what() << std::endl;
    MPI_Abort(MPI_COMM_WORLD, 1);  
  }  
  MPI_Finalize();
  return 0;
}
    // try {
    //   multi_layer_mie_.RunMieCalculations();
    //   std::vector<double> tmp({wl});
    //   //std::vector<double> channels(multi_layer_mie_.GetQabs_channel());
    //   std::vector<double> channels(multi_layer_mie_.GetQabs_channel_normalized());
    //   tmp.insert(tmp.end(), channels.begin(), channels.end());
    //   spectra.push_back(tmp);
    //   if (least_size > tmp.size()) least_size = tmp.size();
    // } catch( const std::invalid_argument& ia ) {
    //   printf(".");
    //   ++fails_;
    // }
// ********************************************************************** //
// ********************************************************************** //
// ********************************************************************** //
void SetOptimizer() {
  //Width is optimized for two layers only!!
  //The third one fills to total_r_
  long dimension = 3;
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
void PrintGnuPlotSpectra(std::vector< std::vector<double> > spectra) {
  gnuplot::GnuplotWrapper wrapper;
  char plot_name [300];
  snprintf(plot_name, 300,
           "TotalR%06.f-spectra", r3_);
  wrapper.SetPlotName(plot_name);
  wrapper.SetXLabelName("Omega");
  wrapper.SetYLabelName("Norm RCS");
  wrapper.SetDrawStyle("w l lw 2");
  wrapper.SetXRange({spectra.front()[0], spectra.back()[0]});
  for (auto multi_point : spectra) wrapper.AddMultiPoint(multi_point);
  wrapper.AddColumnName("Omega/Omega_p");
  wrapper.AddColumnName("Qsca");
  wrapper.MakeOutput();
}
// ********************************************************************** //
// ********************************************************************** //
// ********************************************************************** //
void PrintGnuPlotChannels(std::vector< std::vector<double> > spectra) {
  gnuplot::GnuplotWrapper wrapper;
  char plot_name [300];
  snprintf(plot_name, 300,
           "TotalR%06.f-spectra", r3_);
  wrapper.SetPlotName(plot_name);
  wrapper.SetXLabelName("\\omega/\\omega_p");
  wrapper.SetYLabelName("Norm ACS");
  wrapper.SetDrawStyle("w l lw 2");
  wrapper.SetXRange({spectra.front()[0], spectra.back()[0]});
  for (auto multi_point : spectra) wrapper.AddMultiPoint(multi_point);
  wrapper.AddColumnName("Omega/Omega_p");
  wrapper.AddColumnName("total");
  for (int i = 2; i < spectra.front().size(); ++i) {
    char column_name[10];
    snprintf(column_name, 10, "[%d]",i-1);
    wrapper.AddColumnName(column_name);
  }
  wrapper.MakeOutput();
}
// ********************************************************************** //
// ********************************************************************** //
// ********************************************************************** //

