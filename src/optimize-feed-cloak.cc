/// @file   optimize-feed-cloak.cc
/// @author Ladutenko Konstantin <kostyfisik at gmail (.) com>
/// @date   Tue Sep  3 00:37:05 2013
/// @copyright 2013 Ladutenko Konstantin
///
/// optimize-feed-cloak is free software: you can redistribute it and/or modify
/// it under the terms of the GNU General Public License as published by
/// the Free Software Foundation, either version 3 of the License, or
/// (at your option) any later version.
///
/// optimize-feed-cloak is distributed in the hope that it will be useful,
/// but WITHOUT ANY WARRANTY; without even the implied warranty of
/// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
/// GNU General Public License for more details.
///
/// You should have received a copy of the GNU General Public License
/// along with optimize-feed-cloak.  If not, see <http://www.gnu.org/licenses/>.
///
/// optimize-feed-cloak uses nmie.c from scattnlay by Ovidio Pena
/// <ovidio@bytesfall.com> as a linked library. He has an additional condition to 
/// his library:
//    The only additional condition is that we expect that all publications         //
//    describing  work using this software , or all commercial products             //
//    using it, cite the following reference:                                       //
//    [1] O. Pena and U. Pal, "Scattering of electromagnetic radiation by           //
//        a multilayered sphere," Computer Physics Communications,                  //
//        vol. 180, Nov. 2009, pp. 2348-2354.                                       //
///
/// @brief Simulate scattering from dielectric sphere covered with
///  shell using scattnlay lib 
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
//#include "./nmie/Au-dispersion.h"
const double pi=3.14159265358979323846;
template<class T> inline T pow2(const T value) {return value*value;}

void SetWidth();
double SetInitialModel();
void SetOptimizer();
double EvaluateScatterOnlyIndex(std::vector<double> input);
std::vector< std::vector<double> > EvaluateSpectraForBestDesign();
void PrintCoating(std::vector<double> current, double initial_RCS, jade::SubPopulation sub_population);
void PrintGnuPlotIndex(double initial_RCS, jade::SubPopulation sub_population);
void PrintGnuPlotSpectra(std::vector< std::vector<double> > spectra, double initial_RCS);
void Output();
// ********************************************************************** //
// ********************************************************************** //
// ********************************************************************** //
// Mie model. Used in fitness function for optimization of
// sub_population_.
nmie::MultiLayerMie multi_layer_mie_;  
jade::SubPopulation sub_population_;  // Optimizer of parameters for Mie model.
// ********************************************************************** //
// ********************************************************************** //

double lambda_work_ = 3.75; // cm
double a_ = 0.75*lambda_work_;  // 2.8125 cm - size of PEC core
double min_index_ = 1e-11;
//double min_index_ = 1.0;
double default_index = 1.0;
int number_of_layers_ = 0;
double initial_RCS_ = 0.0;

double Qfailed_ = 1000;
int total_generations_ = 1200;
int population_multiplicator_ = 60;
double layer_width_ = 0.2;
int max_number_of_layers_ = 24;
double from_epsilon_ = -100.0, to_epsilon_ = 100.0;
// ********************************************************************** //
int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);
  try {
    initial_RCS_ = SetInitialModel();
    //multi_layer_mie_.SetMaxTermsNumber(15);
    multi_layer_mie_.SetCoatingWidth({0.1,0.1});
    printf("With %g  coating= (26.24).\n",
	   EvaluateScatterOnlyIndex({-0.29, 24.6}));
    multi_layer_mie_.SetCoatingWidth({0.1,0.1,0.1});
    printf("With %g  coating> (26.24).\n",
	   EvaluateScatterOnlyIndex({-0.29, 24.6, 1.0}));
    //multi_layer_mie_.SetMaxTermsNumber(-1);

      // 26.24: 25||   -0.29   +24.60 
      // 28.48: 38||   -0.29   +24.60    +1.00

    std::vector< std::vector<double> > feed_vector(1, std::vector<double>(1,1.0));

    for (number_of_layers_ = 1; number_of_layers_ < max_number_of_layers_; ++number_of_layers_) {
      SetOptimizer();
      for (auto &feed : feed_vector) {
	while (feed.size()<number_of_layers_) feed.push_back(default_index);
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        if (rank == 0) {
	  for (auto index:feed) printf(" %+7.2f", index);
	  printf("\n");
	}
      }
      sub_population_.SetFeed(feed_vector);
      sub_population_.RunOptimization();
      double best_RCS = 0.0;
      feed_vector.clear();
      feed_vector.push_back(sub_population_.GetBest(&best_RCS));
      Output();         // Output results

      
    }  // end of changing number of layers
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
void Output() {
        auto current = sub_population_.GetFinalFitness();
        int output_rank = 0;
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        for (unsigned int i = 0; i < current.size(); ++i)
          if (current[output_rank] > current[i]) output_rank = i;
        if (rank == output_rank) {
          for (auto c : current) printf("All %g\n",c);
          PrintCoating(current, initial_RCS_, sub_population_);
        }  // end of output for process with best final fitness
        PrintGnuPlotIndex(initial_RCS_, sub_population_);
        PrintGnuPlotSpectra(EvaluateSpectraForBestDesign(), initial_RCS_);
        sub_population_.PrintResult("-- ");
} 
// ********************************************************************** //
// ********************************************************************** //
// ********************************************************************** //
void PrintGnuPlotIndex(double initial_RCS, jade::SubPopulation sub_population) {
  gnuplot::GnuplotWrapper wrapper;
  double best_RCS = 0.0;
  auto best_x = sub_population.GetBest(&best_RCS);
  double total_coating_width = layer_width_ * number_of_layers_;
  double index_sum = 0.0;
  for (auto i : best_x) index_sum+=i;
  char plot_name [300];
  snprintf(plot_name, 300,
           "TargetR%6.4f-CoatingW%06.3f-FinalRCS%07.4fdiff%+05.1f%%-n%02lu-s%015.12f-index",
           a_, total_coating_width,
           best_RCS, (best_RCS/initial_RCS-1.0)*100.0, best_x.size(), index_sum);
  wrapper.SetPlotName(plot_name);
  wrapper.SetXLabelName("Layer #");
  wrapper.SetYLabelName("Epsilon");
  wrapper.SetDrawStyle("with histeps lw 2");
  wrapper.SetXRange({0.51, number_of_layers_+0.49});
  for (int i = 0; i < number_of_layers_; ++i) 
    wrapper.AddMultiPoint({i+1.0, best_x[i]});
  wrapper.AddColumnName("Layer N");
  wrapper.AddColumnName("Epsilon");
  wrapper.MakeOutput();
}
// ********************************************************************** //
// ********************************************************************** //
// ********************************************************************** //
double SetInitialModel() {
  // Set PEC core
  multi_layer_mie_.SetTargetPEC(a_);
  multi_layer_mie_.SetWavelength(lambda_work_);
  multi_layer_mie_.RunMieCalculations();
  double Qsca = multi_layer_mie_.GetQsca();
  double total_r = multi_layer_mie_.GetTotalRadius();
  double initial_RCS = Qsca*pi*pow2(total_r);
  multi_layer_mie_.SetCoatingWidth({0.1});
  multi_layer_mie_.SetCoatingIndex({{1.0,0.0}});
  multi_layer_mie_.RunMieCalculations();
  double Qsca1 = multi_layer_mie_.GetQsca();
  double total_r1 = multi_layer_mie_.GetTotalRadius();
  double initial_RCS1 = Qsca1*pi*pow2(total_r1);
  printf("With %g (r=%g) and without %g (r=%g)air coating.\n",
	 initial_RCS1, total_r1, 
	 initial_RCS, total_r);
  return initial_RCS;
}
// ********************************************************************** //
// ********************************************************************** //
// ********************************************************************** //
void SetOptimizer() {
  long dimension = 0;
  SetWidth();
  dimension = number_of_layers_;
  sub_population_.FitnessFunction = &EvaluateScatterOnlyIndex;
  long total_population = dimension * population_multiplicator_;
  sub_population_.Init(total_population, dimension);
  /// Low and upper bound for all dimensions;
  sub_population_.SetAllBounds(from_epsilon_, to_epsilon_);
  sub_population_.SetTargetToMinimum(); //Check SetQFailed!!!
  sub_population_.SetTotalGenerationsMax(total_generations_);
  sub_population_.SwitchOffPMCRADE();

  sub_population_.SetBestShareP(0.1);
  sub_population_.SetAdapitonFrequencyC(1.0/20.0);

}
// ********************************************************************** //
// ********************************************************************** //
// ********************************************************************** //
void SetWidth() {
  std::vector<double> width;
  width.clear();
  if (number_of_layers_ < 0)
    throw std::invalid_argument("Number of coating layers should be >= 0!");
  for (int i = 0; i < number_of_layers_; ++i) width.push_back(layer_width_);
  multi_layer_mie_.SetCoatingWidth(width);
}
// ********************************************************************** //
// ********************************************************************** //
// ********************************************************************** //
double EvaluateScatterOnlyIndex(std::vector<double> input) {
  double Qsca;
  std::vector<std::complex<double>> cindex;
  cindex.clear();
  double k = min_index_, n=min_index_;
  for (double epsilon : input) {
    k = min_index_, n=min_index_;
    // sqrt(epsilon) = n + i*k
    if (epsilon > 0.0) n=std::sqrt(epsilon);
    else k = std::sqrt(-epsilon);
    if (n < min_index_) n = min_index_;
    if (k < min_index_) k = min_index_;
    //printf("eps= %g, n=%g, k=%g\n", epsilon, n, k);
    cindex.push_back(std::complex<double>(n, k));
  }
  multi_layer_mie_.SetCoatingIndex(cindex);
  try {
    multi_layer_mie_.RunMieCalculations();
    Qsca = multi_layer_mie_.GetQsca();
  } catch( const std::invalid_argument& ia ) {
    auto best_x = sub_population_.GetWorst(&Qfailed_);
    Qsca = Qfailed_;
    printf(".");
    // Will catch if  multi_layer_mie_ fails or other errors.
    //std::cerr << "Invalid argument: " << ia.what() << std::endl;
  }  
  double total_r = multi_layer_mie_.GetTotalRadius();
  return Qsca*pi*pow2(total_r);
}
// ********************************************************************** //
// ********************************************************************** //
// ********************************************************************** //
std::vector< std::vector<double> > EvaluateSpectraForBestDesign() {
  double best_RCS;
  auto best_x = sub_population_.GetBest(&best_RCS);
  // Setting Mie model to the best state.
  sub_population_.FitnessFunction(best_x);
  return multi_layer_mie_.GetSpectra(lambda_work_*0.5, lambda_work_*1.5, 1000);
}
// ********************************************************************** //
// ********************************************************************** //
// ********************************************************************** //
void PrintCoating(std::vector<double> current, double initial_RCS, jade::SubPopulation sub_population) {
  double best_RCS = 0.0;
  auto best_x = sub_population.GetBest(&best_RCS);
  printf("Target R=%g, WL=%g\n", a_, lambda_work_);
  printf("Initial RCS: %g\n", initial_RCS);
  printf("Final RCS: %g (%4.1f%%)\n", best_RCS, (best_RCS/initial_RCS-1.0)*100.0);
  printf ("Layer:\t");
  for (int i = 0; i < number_of_layers_; ++i)
    printf("% 7i\t",i+1);
  printf ("\n");
  printf ("Epsil:\t");
  for (int i = 0; i < number_of_layers_; ++i)
    printf("%+7.2f\t",best_x[i]);
  printf("\n");
  double total_coating_width = 0.0;
  total_coating_width = layer_width_*number_of_layers_;
  printf("Total coating width: %g\n", total_coating_width);
}
// ********************************************************************** //
// ********************************************************************** //
// ********************************************************************** //
void PrintGnuPlotSpectra(std::vector< std::vector<double> > spectra, double initial_RCS) {
  gnuplot::GnuplotWrapper wrapper;
  double best_RCS = 0.0;
  auto best_x = sub_population_.GetBest(&best_RCS);
  double total_coating_width = layer_width_*number_of_layers_;
  double index_sum = 0.0;
  for (auto i : best_x) index_sum+=i;
  char plot_name [300];
  snprintf(plot_name, 300,
           "TargetR%6.4f-CoatingW%06.3f-FinalRCS%07.4fdiff%+05.1f%%-n%02lu-s%015.12f-spectra",
           a_, total_coating_width,
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
