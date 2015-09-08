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
#include "./nmie/nmie-applied.h"
//#include "./read-spectra/read-spectra.h"
const double pi=3.14159265358979323846;
const double speed_of_light = 299792458;
template<class T> inline T pow2(const T value) {return value*value;}
void SetOptimizer();
void SetMie();
void SetGeometry();

double EvaluateFitness(std::vector<double> input);
jade::SubPopulation sub_population_;  // Optimizer of parameters for Mie model.
// ********************************************************************** //
// ********************************************************************** //
// ********************************************************************** //
nmie::MultiLayerMieApplied multi_layer_mie_;  
double l2w( double l) {return 2.0 * pi * speed_of_light/l;};  // lambda to omega
double w2l( double w) {return 2.0 * pi * speed_of_light/w;};
// ********************************************************************** //
// ********************************************************************** //
// ********************************************************************** //
const double lambda_0_ = 500.0e-9;  // 500 nm
const double omega_0_ = l2w(lambda_0_);
bool isPreset = true;
//bool isPreset = false;
//bool isOuterR = true;
bool isOuterR = false;
double outer_r_ = 0.5;
int dim_=3;
//Alu params
//std::vector<double> input_ = {0.13, 0.16, 0.194,1.0};
//std::vector<double> input_ = {0.142, 0.166, 0.5194, 1.0};
//std::vector<double> input_ = {0.142, 0.166, 0.194, 1.0};
std::vector<double> input_ = {0.00635, 0.00747, 0.00747000001, 1.0};
//0.142135
//std::vector<double> input_ = {6.3535e-3, 7.4765e-3, 7.4766e-3, 1};
double r1_ = input_[0]*lambda_0_;
double r2_ = input_[1]*lambda_0_;
double r3_ = input_[2]*lambda_0_;
double omega_ = input_[3]*omega_0_;
double core_width_ = r1_;
double inshell_width_ = r2_ - r1_;
double outshell_width_ = r3_ - r2_;
double from_omega_ = 0.1*omega_0_;
double to_omega_ = 2.0*omega_0_;
std::complex<double> inshell_index_(0,0);
std::complex<double> core_index_ = std::sqrt(std::complex<double>(1.29,0.01));
std::complex<double> outshell_index_ = std::sqrt(std::complex<double>(8.4,2.33));
// std::complex<double> core_index_ = std::sqrt(std::complex<double>(1.25,0.03));
// std::complex<double> outshell_index_ = std::sqrt(std::complex<double>(8.6,0.96)) ;



const double gamma_d_ = 2.0*pi*17.64*1.0e12;
const double omega_p_ = 2.0*pi*2069.0*1.0e12;

std::complex<double> epsilon_m(double omega) {
  //return 1.53 - pow2(omega_p_)/(omega*(omega +std::complex<double>(0,1)*gamma_d_));
  return -10.37+ 0.35*std::complex<double>(0,1);
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
int population_multiplicator_ = 250;
double Qsca_best_ = 0.0;
double Qabs_best_ = 0.0;
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
    if (rank ==0) {printf("Input_:"); for (auto value : input_) printf(" %24.22f,", value);  }
    multi_layer_mie_.RunMieCalculation();
    Qsca_best_ = multi_layer_mie_.GetQsca();
    Qabs_best_ = multi_layer_mie_.GetQabs();
    if (rank ==0) {
      printf("\nQabs: %24.22f\nQsca: %24.22f\nZeta=%24.22f\n",Qabs_best_,Qsca_best_, Qabs_best_/Qsca_best_);
      double Cabs = Qabs_best_*pi*pow2(r3_);
      double A = 3.0*pow2(lambda_0_)/(8.0*pi);
      printf("Cabs = %g\nA = %g\n",Cabs,A);
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
  // Copy releveant part of optimizer input to global var input_ for Mie calculations
  for (int i = 0; i< input.size(); ++i)
    input_[i] = input[i];
  SetGeometry();
  SetMie();
  double Zeta = 0.0;
  try {
    multi_layer_mie_.RunMieCalculation();
    Zeta = multi_layer_mie_.GetQabs()/multi_layer_mie_.GetQsca();
  } catch( const std::invalid_argument& ia ) {
    printf(".");
    sub_population_.GetWorst(&Zeta);
  }
  double r_outer = input_[2]*lambda_0_;
  double Qabs = multi_layer_mie_.GetQabs();
  double Cabs = Qabs*pi*pow2(r3_);
  double A = 3.0*pow2(lambda_0_)/(8.0*pi);
  double Q0 = 5.0, Z0=1000.0;

  //return Qabs > Q0 ? Zeta*pow2(pow2(Q0/Qabs)) : Zeta*pow2(pow2(Qabs/Q0));
  //return Cabs > A ? Zeta*pow2(pow2(A/Cabs)) : Zeta*pow2(pow2(Cabs/A));
  return Zeta > Z0 ? Z0/Zeta*pow2(Qabs) : Zeta/Z0*pow2(Qabs);
  //return Cabs > A ? Zeta : 0.0;
}
// ********************************************************************** //
// ********************************************************************** //
// ********************************************************************** //
void SetGeometry() {
  if (isOuterR) {
    input_[2] = outer_r_;
    if (input_[0] > input_[2]) input_[0] = input_[2]-2.0*eps_;
    if (input_[1] > input_[2]) input_[1] = input_[2]-eps_;
  }
  if (input_[1] < input_[0]) input_[1] = input_[0]+eps_;
  if (input_[2] < input_[1]) input_[2] = input_[1]+eps_;
  r1_ = input_[0]*lambda_0_;
  r2_ = input_[1]*lambda_0_;
  r3_ = input_[2]*lambda_0_;
  omega_ = input_[3]*omega_0_;
  core_width_ = r1_;
  inshell_width_ = r2_ - r1_;
  outshell_width_ = r3_ - r2_;
}
// ********************************************************************** //
// ********************************************************************** //
// ********************************************************************** //
void SetMie() {
      multi_layer_mie_.ClearTarget();
      multi_layer_mie_.AddTargetLayer(core_width_, core_index_);
      multi_layer_mie_.AddTargetLayer(inshell_width_, std::sqrt(epsilon_m(omega_)));
      // if (omega_ > omega_0_*0.999 && omega_ < omega_0_*1.001)
      // 	printf("eps = %g, %gj",epsilon_m(omega_).real(), epsilon_m(omega_).imag());
      multi_layer_mie_.AddTargetLayer(outshell_width_, outshell_index_);
      multi_layer_mie_.SetWavelength(w2l(omega_));
}
// ********************************************************************** //
// ********************************************************************** //
// ********************************************************************** //
void SetOptimizer() {
  //Width is optimized for two layers only!!
  //The third one fills to total_r_
  long dimension = dim_;
  sub_population_.FitnessFunction = &EvaluateFitness;
  long total_population = dimension * population_multiplicator_;
  sub_population_.Init(total_population, dimension);
  /// Low and upper bound for all dimenstions;
  sub_population_.SetAllBounds(eps_, 2.0-eps_);
  //sub_population_.SetAllBounds(eps_, input_[2]-eps_);
  sub_population_.SetTargetToMaximum();
  sub_population_.SetTotalGenerationsMax(total_generations_);
  //sub_population.SwitchOffPMCRADE();

  sub_population_.SetBestShareP(0.1);
  sub_population_.SetAdapitonFrequencyC(1.0/20.0);

}
// ********************************************************************** //
// ********************************************************************** //
// ********************************************************************** //

