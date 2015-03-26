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
//#include "./jade.h"
#include <cstdlib>
#include <cstdio>
#include <ctime>
#include <stdexcept>
#include <string>
#include "./gnuplot-wrapper/gnuplot-wrapper.h"
#include "./nmie/nmie-wrapper.h"
#include "./read-spectra/read-spectra.h"
const double pi=3.14159265358979323846;
const double speed_of_light = 299792458;
template<class T> inline T pow2(const T value) {return value*value;}
// Mie model. Used in fitness function for optimization of
// sub_population.
//void SetOptimizer();

double EvaluateFitness(std::vector<double> input);
void PrintGnuPlotSpectra(std::vector< std::vector<double> > spectra);
void PrintGnuPlotChannels(std::vector< std::vector<double> > spectra);
// ********************************************************************** //
// ********************************************************************** //
// ********************************************************************** //
nmie::MultiLayerMie multi_layer_mie_;  
double eV2nm( double l) {return 1239.84/l;}
double nm2eV( double w) {return 1239.84/w;}
read_spectra::ReadSpectra core_index_;
// ********************************************************************** //
// ********************************************************************** //
// ********************************************************************** //
double core_width_ = 10;
double from_omega_ = 2.0;
double to_omega_ = 5.0;
double glass_index_ = 1.5;
// Set dispersion
double from_wl_ = eV2nm(to_omega_);
double to_wl_ = eV2nm(from_omega_);
int samples_ = 300;
double plot_from_wl_ = from_wl_, plot_to_wl_ = to_wl_;
int plot_samples_ = samples_;
// ********************************************************************** //
// ********************************************************************** //
// ********************************************************************** //
int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  try {
    std::string core_filename("Ag.txt");
    core_index_.ReadFromFile(core_filename);
    core_index_.ResizeToComplex(from_wl_, to_wl_, samples_).ToIndex();
    auto core_data = core_index_.GetIndex();
    const double lossless = 0.0;
    int least_size = 10000;
    std::vector< std::vector<double> > spectra;
    if (from_omega_ > to_omega_) throw std::invalid_argument("Wrong omega range!");
    for (int i=core_data.size()-1; i>-1;--i) {
      const double& wl = core_data[i].first;
      const std::complex<double>& core_index_value = core_data[i].second;
      const std::complex<double> index_norm={core_index_value.real()/glass_index_,
	    core_index_value.imag()};
      //double metal_index_re = std::sqrt(epsilon_m(omega));
      multi_layer_mie_.ClearTarget();
      multi_layer_mie_.AddTargetLayer(core_width_, index_norm);
      multi_layer_mie_.SetWavelength(wl/glass_index_);
      try {
	multi_layer_mie_.RunMieCalculations();
	double Qsca = multi_layer_mie_.GetQabs();
	// double norm = (Qsca * pi*pow2(r3_)) / ( pow2(w2l(omega)) / (2.0*pi) );
	// std::vector<double> tmp({omega/omega_p_, norm});
	// //std::vector<double> channels(multi_layer_mie_.GetQsca_channel_normalized());
	// //std::vector<double> channels(multi_layer_mie_.GetQabs_channel_normalized());
	// std::vector<double> channels(multi_layer_mie_.GetQabs_channel());
	// tmp.insert(tmp.end(), channels.begin(), channels.end());
	// spectra.push_back(tmp);
	// if (least_size > tmp.size()) least_size = tmp.size();
	spectra.push_back({nm2eV(wl), Qsca*2.0*pi*pow2(core_width_)});
      } catch( const std::invalid_argument& ia ) {
	std::cerr << "Invalid argument: " << ia.what() << std::endl;
	printf(".");
      }
    }  // end of omega sweep
    PrintGnuPlotSpectra(spectra);
    // for (auto& row : spectra) {
    //   //if (rank==0) std::cout<< least_size<< " -- "<<row.size()<<std::endl;
    //   row.resize(least_size);
    // }
   
    // PrintGnuPlotChannels(spectra);
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
void PrintGnuPlotSpectra(std::vector< std::vector<double> > spectra) {
  gnuplot::GnuplotWrapper wrapper;
  char plot_name [300];
  snprintf(plot_name, 300,
           "TotalR%06.f-spectra", core_width_);
  wrapper.SetPlotName(plot_name);
  wrapper.SetXLabelName("Omega, eV");
  wrapper.SetYLabelName("ACS");
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
           "TotalR%06.f-spectra", core_width_);
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

