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
void SetMie();
void SetGeometry();

double EvaluateFitness(std::vector<double> input);
double EvaluateFitnessMod(std::vector<double> input);
double EvaluateFitnessChannel(std::vector<double> input);
std::vector< std::vector<double> > EvaluateSpectraForChannels();
void Print();
void PrintCoating(std::vector<double> current, double initial_RCS,
                    jade::SubPopulation sub_population);
//void PrintGnuPlotSpectra(std::vector< std::vector<double> > spectra);
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
//bool isPreset = true;
bool isPreset = false;
bool isQabs = true;
//bool isQabs = false;
bool isOuterR = true;
//bool isOuterR = false;
double outer_r_ = 0.5;
int dim_=5;
//S.Fan params
std::vector<double> input_ = {0.4749, 0.6404, 0.8249, 0.2932, 1.0};
//Custom 
//std::vector<double> input_ = {0.25, 0.30, 0.00, 0.2932}; //12.92
//std::vector<double> input_ = {+0.12,    +0.08,    +0.02,    +0.54}; //21.35
//std::vector<double> input_ = {0.6235, 0.7489, 0.8249, 0.2932}; //17.61
//std::vector<double> input_ = {0.45464, 0.520026, 0.0, 0.2932}; //29.69
double r1_ = input_[0]*lambda_p_;
double r2_ = input_[1]*lambda_p_;
double r3_ = input_[2]*lambda_p_;
double omega_resonance_ = input_[3]*omega_p_;
double core_width_ = r1_;
double inshell_width_ = r2_ - r1_;
double outshell_width_ = r3_ - r2_;
double from_omega_ = 0.288*omega_p_;
double to_omega_ = 0.3*omega_p_;
double omega_ = 0.0;
double epsilon_d_ = 12.96;
double inshell_index_ = std::sqrt(epsilon_d_);
double gamma_d_ = 0.0;
double gamma_bulk_ = 0.002*omega_p_;
//bool isLossy_ = false;
bool isLossy_ = true;

std::complex<double> epsilon_m(double omega) {
  // std::cout << "omega:" << omega << "  omega_p:"<<omega_p_ <<std::endl;
  return 1.0 - pow2(omega_p_)/(pow2(omega)
			       +std::complex<double>(0,1)*gamma_d_*omega
			       *input_[4]);
}
// Set dispersion
double from_wl_ = w2l(from_omega_);
double to_wl_ = w2l(to_omega_);
int samples_ = 1300;
double plot_from_wl_ = from_wl_, plot_to_wl_ = to_wl_;
int plot_samples_ = samples_;
double plot_xshare_ = 0.1;
// Set optimizer
int total_generations_ = 1500;
int population_multiplicator_ = 160;
double Qsca_best_ = 0.0;
// ********************************************************************** //
// ********************************************************************** //
// ********************************************************************** //
const double eps_=1e-11;
int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  try {
    std::vector< std::vector<double> > spectra;
    if (from_omega_ > to_omega_) throw std::invalid_argument("Wrong omega range!");
    
    if (!isPreset) {
      SetOptimizer();
      // std::vector<double> feed = {input_[0],input_[1]};
      // sub_population_.SetFeed({feed});
      sub_population_.RunOptimization();
      auto best_x  = sub_population_.GetBest(&Qsca_best_);
      for (int i = 0; i< best_x.size(); ++i)
	input_[i] = best_x[i];
    }
    
    SetGeometry();    SetMie();
    if (rank ==0) {printf("Input_:"); for (auto value : input_) printf(" %g,", value);  }
    multi_layer_mie_.RunMieCalculations();
    std::vector<double> channels;
    if (isQabs){
      Qsca_best_ = multi_layer_mie_.GetQabs();
      channels = multi_layer_mie_.GetQabs_channel();
    } else {
      Qsca_best_ = multi_layer_mie_.GetQsca();
      channels = multi_layer_mie_.GetQsca_channel();
    }
    if (rank ==0) {
	printf("\nQsca_best: %g\n",Qsca_best_);
	for (auto value : channels) printf(" %g,", value);    
	PrintGnuPlotChannels(EvaluateSpectraForChannels());
    }
  } catch( const std::invalid_argument& ia ) {
    // Will catch if  multi_layer_mie_ fails or other errors.
    std::cerr << "Invalid argument: " << ia.what() << std::endl;
    //MPI_Abort(MPI_COMM_WORLD, 1);  
  }  
  MPI_Finalize();
  return 0;
}
// ********************************************************************** //
// ********************************************************************** //
// ********************************************************************** //
double EvaluateFitness(std::vector<double> input) {
  if (input.size() > 5) throw std::invalid_argument("Wrong input dimension!/n");
  for (int i = 0; i< input.size(); ++i)
    input_[i] = input[i];
  SetGeometry();
  SetMie();
  double Qsca = 0.0;
  try {
    multi_layer_mie_.RunMieCalculations();
    Qsca = multi_layer_mie_.GetQsca();
  } catch( const std::invalid_argument& ia ) {
    printf(".");
    sub_population_.GetWorst(&Qsca);
  }
  return Qsca;
}
// ********************************************************************** //
// ********************************************************************** //
// ********************************************************************** //
double EvaluateFitnessQabs(std::vector<double> input) {
  if (input.size() > 5) throw std::invalid_argument("Wrong input dimension!/n");
  for (int i = 0; i< input.size(); ++i)
    input_[i] = input[i];
  SetGeometry();
  SetMie();
  double Qabs = 0.0;
  try {
    multi_layer_mie_.RunMieCalculations();
    Qabs = multi_layer_mie_.GetQabs();
  } catch( const std::invalid_argument& ia ) {
    printf(".");
    sub_population_.GetWorst(&Qabs);
  }
  return Qabs;
}
// ********************************************************************** //
// ********************************************************************** //
// ********************************************************************** //
double EvaluateFitnessMod(std::vector<double> input) {
  if (input.size() > 5) throw std::invalid_argument("Wrong input dimension!/n");
  for (int i = 0; i< input.size(); ++i)
    input_[i] = input[i];
  SetGeometry();
  SetMie();
  double Qsca = 0.0;
  try {
    multi_layer_mie_.RunMieCalculations();
    Qsca = multi_layer_mie_.GetQsca();
    std::vector<double> channels(multi_layer_mie_.GetQsca_channel());
    // no loss
    //Qsca = channels[0]*channels[1]*channels[2];
    // with losses
    Qsca = channels[0]*channels[1];
  } catch( const std::invalid_argument& ia ) {
    printf(".");
    Qsca = 0.0;
    sub_population_.GetWorst(&Qsca);
  }
  return Qsca;
}
// ********************************************************************** //
// ********************************************************************** //
// ********************************************************************** //
std::vector< std::vector<double> > EvaluateSpectraForChannels() {
    int least_size = 10000;
    std::vector< std::vector<double> > spectra;
    double best_3 = input_[3];
    omega_resonance_ = best_3;
    // double omega_step = (to_omega_ - from_omega_)/samples_;
    // for (input_[3] = from_omega_/omega_p_; input_[3] < to_omega_/omega_p_;
    // 	 input_[3] += omega_step/omega_p_) {
    for (input_[3] = best_3*(1.0 - plot_xshare_);
    	 input_[3] < best_3*(1.0 + plot_xshare_);
    	 input_[3] += best_3*2.0*plot_xshare_/samples_) {
      SetGeometry();
      SetMie();
      try {
	multi_layer_mie_.RunMieCalculations();
	double Qsca;
	std::vector<double> channels;
	if (isQabs) {
	  Qsca = multi_layer_mie_.GetQabs();
	  channels = multi_layer_mie_.GetQabs_channel();
	} else {
	  Qsca = multi_layer_mie_.GetQsca();
	  channels = multi_layer_mie_.GetQsca_channel();
	}
	double norm = (Qsca * pi*pow2(r3_)) / ( pow2(w2l(omega_)) / (2.0*pi) );
	std::vector<double> tmp({omega_/omega_p_, norm});
	for (auto& value : channels) 
	  value *= pi*pow2(r3_) / ( pow2(w2l(omega_)) / (2.0*pi) );
	tmp.insert(tmp.end(), channels.begin(), channels.end());
	spectra.push_back(tmp);
	if (least_size > tmp.size()) least_size = tmp.size();
	if (std::abs((input_[3]-0.2932)/input_[3])<eps_ || 
	    std::abs((input_[3]-0.5426))<0.00001) {
	  printf("\n\nFrom spectra:\n");
	  printf("Input_:");
	  for (auto value : input_) printf(" %g,", value);      
	  multi_layer_mie_.RunMieCalculations();
	  double Qsca;
	  std::vector<double> channels;
	  if (isQabs) {
	    Qsca = multi_layer_mie_.GetQabs();
	    channels = multi_layer_mie_.GetQabs_channel();
	  } else {
	    Qsca = multi_layer_mie_.GetQsca();
	    channels = multi_layer_mie_.GetQsca_channel();
	  }
	  printf("\nQsca_best: %g\n",Qsca);
	  for (auto value : channels) printf(" %g,", value);
     	}
      } catch( const std::invalid_argument& ia ) {
	std::cerr << "Invalid argument: " << ia.what() << std::endl;
	printf(".");
      }
    }  // end of omega sweep
    for (auto& row : spectra) {
      //if (rank==0) std::cout<< least_size<< " -- "<<row.size()<<std::endl;
      row.resize(least_size);
    }
    return spectra;
}
// ********************************************************************** //
// ********************************************************************** //
// ********************************************************************** //
void SetGeometry() {
  if (input_[4] < 0.2) input_[4] = 0.2;
  if (isOuterR) {
    input_[2] = outer_r_;
    if (input_[0] > input_[2]) input_[0] = input_[2]-2.0*eps_;
    if (input_[1] > input_[2]) input_[1] = input_[2]-eps_;
  }
  if (input_[1] < input_[0]) input_[1] = input_[0]+eps_;
  if (input_[2] < input_[1]) input_[2] = input_[1]+eps_;
  r1_ = input_[0]*lambda_p_;
  r2_ = input_[1]*lambda_p_;
  r3_ = input_[2]*lambda_p_;
  omega_ = input_[3]*omega_p_;
  core_width_ = r1_;
  inshell_width_ = r2_ - r1_;
  outshell_width_ = r3_ - r2_;
}
// ********************************************************************** //
// ********************************************************************** //
// ********************************************************************** //
void SetMie() {
      multi_layer_mie_.ClearTarget();
      double A = 1.0;
      double V_f = 7.37e-4*lambda_p_*omega_p_;
      double l_r = r1_;
      if (isLossy_) gamma_d_ = gamma_bulk_ + A*V_f/l_r;
      multi_layer_mie_.AddTargetLayer(core_width_, std::sqrt(epsilon_m(omega_)));
      const double lossless = 0.0;    
      multi_layer_mie_.AddTargetLayer(inshell_width_, {inshell_index_, lossless});
      l_r = r3_-r2_;
      if (isLossy_) gamma_d_ = gamma_bulk_ + A*V_f/l_r;
      multi_layer_mie_.AddTargetLayer(outshell_width_, std::sqrt(epsilon_m(omega_)));
      multi_layer_mie_.SetWavelength(w2l(omega_));
}
// ********************************************************************** //
// ********************************************************************** //
// ********************************************************************** //
void SetOptimizer() {
  //Width is optimized for two layers only!!
  //The third one fills to total_r_
  long dimension = dim_;
  if (isQabs)   sub_population_.FitnessFunction = &EvaluateFitnessQabs;
  else {
    sub_population_.FitnessFunction = &EvaluateFitness;
    //sub_population_.FitnessFunction = &EvaluateFitnessMod;
  }
  long total_population = dimension * population_multiplicator_;
  sub_population_.Init(total_population, dimension);
  /// Low and upper bound for all dimenstions;
sub_population_.SetAllBounds(eps_, 2.0-eps_);
  if (dimension == 2) sub_population_.SetAllBounds(eps_, 0.8249-eps_);

  sub_population_.SetTargetToMaximum();
  sub_population_.SetTotalGenerationsMax(total_generations_);
  //sub_population.SwitchOffPMCRADE();

  sub_population_.SetBestShareP(0.1);
  sub_population_.SetAdapitonFrequencyC(1.0/20.0);

}
// ********************************************************************** //
// ********************************************************************** //
// ********************************************************************** //
// void PrintGnuPlotSpectra(std::vector< std::vector<double> > spectra) {
//   gnuplot::GnuplotWrapper wrapper;
//   char plot_name [300];
//   snprintf(plot_name, 300,
//            "TotalR%06.f-spectra", r3_);
//   wrapper.SetPlotName(plot_name);
//   wrapper.SetXLabelName("Omega");
//   wrapper.SetYLabelName("Norm RCS");
//   wrapper.SetDrawStyle("w l lw 2");
//   wrapper.SetXRange({spectra.front()[0], spectra.back()[0]});
//   for (auto multi_point : spectra) wrapper.AddMultiPoint(multi_point);
//   wrapper.AddColumnName("Omega/Omega_p");
//   wrapper.AddColumnName("Qsca");
//   wrapper.MakeOutput();
// }
// ********************************************************************** //
// ********************************************************************** //
// ********************************************************************** //
void PrintGnuPlotChannels(std::vector< std::vector<double> > spectra) {
  gnuplot::GnuplotWrapper wrapper;
  char plot_name [300];
  std::string type;
  if (isQabs) type="Qabs";
  else type = "Qsca";
  snprintf(plot_name, 300,
           "%s%07.4f-Rcore%06.4f-Rin%06.4f-Rout%06.4f-omega%06.4f-natural%06.4f-spectra-dim%d",type.c_str(),  Qsca_best_,
	   input_[0],input_[1],input_[2],omega_resonance_,input_[4],dim_);
  wrapper.SetPlotName(plot_name);
  wrapper.SetXLabelName("\\omega/\\omega_p");
  if (isQabs)   wrapper.SetYLabelName("Norm ACS");
  else wrapper.SetYLabelName("Norm RCS");
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

