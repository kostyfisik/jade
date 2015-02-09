/**
 * @file   coating-w-sweep.cc
 * @author Konstantin Ladutenko <kostyfisik at gmail (.) com>
 * @date   Sun Feb  8 13:48:43 2015
 */
/// @copyright 2013 Konstantin Ladutenko 
///
/// coating-w-sweep is free software: you can redistribute it and/or modify
/// it under the terms of the GNU General Public License as published by
/// the Free Software Foundation, either version 3 of the License, or
/// (at your option) any later version.
///
/// coating-w-sweep is distributed in the hope that it will be useful,
/// but WITHOUT ANY WARRANTY; without even the implied warranty of
/// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
/// GNU General Public License for more details.
///
/// You should have received a copy of the GNU General Public License
/// along with coating-w-sweep.  If not, see <http://www.gnu.org/licenses/>.
///
/// coating-w-sweep uses nmie.c from scattnlay by Ovidio Pena
/// <ovidio@bytesfall.com> as a linked library. He has an additional condition to 
/// his library:
//    The only additional remark is that we expect that all publications            //
//    describing  work using this software , or all commercial products             //
//    using it, cite the following reference:                                       //
//    [1] O. Pena and U. Pal, "Scattering of electromagnetic radiation by           //
//        a multilayered sphere," Computer Physics Communications,                  //
//        vol. 180, Nov. 2009, pp. 2348-2354.                                       //
///
/// @brief Simulate scattering from  a core-shell sphere 
/// 
#include <mpi.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstdio>
#include <ctime>
#include <stdexcept>
#include <string>
#include "./gnuplot-wrapper/gnuplot-wrapper.h"
#include "./nmie/ucomplex.h"
#include "./nmie/nmie-wrapper.h"
#include "./nmie/Au-dispersion.h"
const double pi=3.14159265358979323846;
template<class T> inline T pow2(const T value) {return value*value;}
// Mie model. Used in fitness function for optimization of
// sub_population.
nmie::MultiLayerMie multi_layer_mie;  
// ********************************************************************** //
// ********************************************************************** //
// ********************************************************************** //
bool isUsingPEC = true;
//bool isUsingPEC = false;
bool isOnlyIndexOptimization = true;
//bool isOnlyIndexOptimization = false;
double lambda_work = 0.532; // nm
//    double f_work = 30/lambda_work; // 8 GHz
//double a = 1; // Krasnok PEC
double a = 0.75*lambda_work;  // 2.8125 cm
//double a = lambda_work;  // 
//double b = pi*pow2(a);
//size param = 2 pi r/wl = 2pi0.75 = 4.71
//double layer_thickness = 0.015*a;
double layer_thickness = 0.0;
double n = 4;
double k = 0;
int number_of_layers = 1;
// Production parameters
//int total_generations = 1200;
//double thickness_step = 0.02;
// Test parameters
int total_generations = 120;
double thickness_step = 0.2;
void SetTarget(double n, double k);
void SetThickness();
double SetInitialModel(double n, double k);
void SetOptimizer();
double EvaluateScatterOnlyIndex(std::vector<double> input);
double EvaluateScatter(std::vector<double> input);
std::vector< std::vector<double> > EvaluateSpectraForBestDesign();
void PrintGnuPlotSpectra(std::vector< std::vector<double> > spectra,
                         double initial_RCS);
void PrintGnuPlotThickness(std::vector< std::vector<double> > spectra,
			   double meta_n);
// ********************************************************************** //
// ********************************************************************** //
// ********************************************************************** //
double loss_index = 1e-11;
int main(int argc, char *argv[]) {
  int rank;
  try {
    if (isUsingPEC) {n = -1.0; k = -1.0;}
    double r_begin = a;
    lambda_work = 0.532; // nm
    double meta_n = 0.8, meta_k = 1e-11;
    for (meta_n = 1.0; meta_n > 0.6; meta_n-=0.01) {
      double step_thickness = 0.003;
      std::vector< std::vector<double> > spectra;
      for (layer_thickness = step_thickness; layer_thickness < 0.13; layer_thickness+=step_thickness) {
	multi_layer_mie.ClearTarget();
	double initial_RCS = SetInitialModel(n, k);
	multi_layer_mie.AddTargetLayer(layer_thickness, {meta_n,meta_k});
	double Qext, Qsca, Qabs, Qbk;
	try {
	  multi_layer_mie.RunMie(&Qext, &Qsca, &Qabs, &Qbk);
      } catch( const std::invalid_argument& ia ) {
	  printf(".");
	  // Will catch if  multi_layer_mie fails or other errors.
	  //std::cerr << "Invalid argument: " << ia.what() << std::endl;
      }  
	double total_r = multi_layer_mie.GetTotalRadius();
	double total_RCS = Qsca*pi*pow2(total_r);
	//if (total_RCS < 100) 
	spectra.push_back({layer_thickness, total_RCS});      
      } // end of thickness sweep
      PrintGnuPlotThickness(spectra, meta_n);
    } // end of meta_n sweep
    // auto spectra = multi_layer_mie.GetSpectra(lambda_work*0.1, lambda_work*15.0, 1000);
    // auto totalR =  multi_layer_mie.GetTotalRadius();
    // printf("spectra size %li, initRCS=%g, totalR=%g\n", spectra.size(), initial_RCS, totalR);
    // for (auto &point:spectra) {
    //   //printf("%5g\t\%g\n", point[0], point[2]);      
    //   point[0]= 1/point[0];
    //   // for (int i=1; i<point.size(); ++i)
    //   // 	point[i] *= pi*pow2(totalR);
    // }
    // PrintGnuPlotSpectra(spectra, initial_RCS);
      
    //}  // end of total coating thickness sweep
    // }  // end of k sweep
  } catch( const std::invalid_argument& ia ) {
    // Will catch if  multi_layer_mie fails or other errors.
    std::cerr << "Invalid argument: " << ia.what() << std::endl;
  }  
  return 0;
}

// ********************************************************************** //
// ********************************************************************** //
// ********************************************************************** //
void PrintGnuPlotThickness(std::vector< std::vector<double> > spectra,
			   double meta_n) {
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
  wrapper.SetYLabelName("Total RCS");
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
void PrintGnuPlotSpectra(std::vector< std::vector<double> > spectra,
                         double initial_RCS) {
  gnuplot::GnuplotWrapper wrapper;
  double best_RCS = 0.0;
  //  auto best_x = sub_population.GetBest(&best_RCS);
  double total_coating_width = layer_thickness*number_of_layers;
  double index_sum = 0.0;
  //for (auto i : best_x) index_sum+=i;
  char plot_name [300];
  snprintf(plot_name, 300,
           "TargetR%6.4f-CoatingW%06.3f-FinalRCS%07.4fdiff%+4.1f%%-spectra",
           a, total_coating_width,
           best_RCS, (best_RCS/initial_RCS-1.0)*100.0);
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
double SetInitialModel(double n, double k) {
  // Set common parameters for all wavelengths.
  SetTarget(n, k);
  multi_layer_mie.SetQfaild(1000.0);  // Searching for minima
  multi_layer_mie.SetWavelength(lambda_work);
  double Qext, Qsca, Qabs, Qbk;
  try {
    multi_layer_mie.RunMie(&Qext, &Qsca, &Qabs, &Qbk);
  } catch( const std::invalid_argument& ia ) {
    // Will catch if  multi_layer_mie fails or other errors.
    // std::cerr << "Invalid argument: " << ia.what() << std::endl;
  }  
  double total_r = multi_layer_mie.GetTotalRadius();
  double initial_RCS = Qsca*pi*pow2(total_r);
  return initial_RCS;
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
void SetThickness() {
  std::vector<double> thickness;
  thickness.clear();
  if (number_of_layers < 0)
    throw std::invalid_argument("Number of coating layers should be >= 0!");
  for (int i = 0; i < number_of_layers; ++i) thickness.push_back(layer_thickness);
  multi_layer_mie.SetCoatingThickness(thickness);
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
double EvaluateScatter(std::vector<double> input) {
  std::vector<double> thickness;
  std::vector<complex> cindex;
  double k = loss_index;
  for (int i = 0; i < number_of_layers; ++i) {
    cindex.push_back({input[i], k});
    if (input[i+number_of_layers] < 1.0) input[i+number_of_layers] = 1.0; 
    thickness.push_back(input[i+number_of_layers]*layer_thickness);
  }
  multi_layer_mie.SetCoatingIndex(cindex);
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
  return Qsca*pi*pow2(total_r);
}
// ********************************************************************** //
// ********************************************************************** //
// ********************************************************************** //
