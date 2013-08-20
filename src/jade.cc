///
/// @file   jade.cc
/// @author Ladutenko Konstantin <kostyfisik at gmail (.) com>
/// @date   Thu Aug 15 19:26:53 2013
/// @copyright 2013 Ladutenko Konstantin
/// @section LICENSE
/// This file is part of JADE++.
///
/// JADE++ is free software: you can redistribute it and/or modify
/// it under the terms of the GNU General Public License as published by
/// the Free Software Foundation, either version 3 of the License, or
/// (at your option) any later version.
///
/// JADE++ is distributed in the hope that it will be useful,
/// but WITHOUT ANY WARRANTY; without even the implied warranty of
/// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
/// GNU General Public License for more details.
///
/// You should have received a copy of the GNU General Public License
/// along with JADE++.  If not, see <http://www.gnu.org/licenses/>.
///
/// @brief JADE++ is a free (GPLv3+) high performance implementation of
/// adaptive differential evolution optimization algorithm from
/// Jingqiao Zhang and Arthur C. Sanderson book 'Adaptive Differential
/// Evolution. A Robust Approach to Multimodal Problem Optimization'
/// Springer, 2009.
#include "./jade.h"
#include <mpi.h>
#include <random>
#include <cstdio>
#include <cmath>
#include <algorithm>
#include <map>
#include <string>
namespace jade {
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  int SubPopulation::RunOptimization() {
    CreateInitialPopulation();
    EvaluateCurrentVectors();
    printf("Optimize\n");

    return kDone;
  }  // end of int SubPopulation::RunOptimization()
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  /// %todo Change returned kError to meangfull error code.
  int SubPopulation::Init(long long total_population, long long dimension) {// NOLINT
    total_population_  = total_population;
    if (total_population_ < 1) return kError;
    dimension_ = dimension;
    if (dimension_ < 1) return kError;
    MPI_Comm_rank(MPI_COMM_WORLD, &process_rank_);
    MPI_Comm_size(MPI_COMM_WORLD, &number_of_processes_);
    if (process_rank_ < 0) return kError;
    if (number_of_processes_ < 1) return kError;
    std::random_device rd;
    generator_.seed(rd());
    // //debug
    // CheckRandom();
    if (number_of_processes_ > total_population_) return kError;
    if (number_of_processes_ == total_population_) subpopulation_ = 1;
    // //debug section
    // if (process_rank_ == 0) {
    //   for (long long popul = 1; popul < 10000; popul++) {
    //     total_population_ = popul;
    //     for (int procs = 1; procs < 300; procs++) {
    //       if (procs >= popul) break;
    //       number_of_processes_ = procs;
    //       long long popul_eval = 0;
    //       for (int rank = 0; rank < procs; rank++) {
    //         process_rank_ = rank;
    // //end of debug section
    double subpopulation_size = static_cast <double> (total_population_)
      / static_cast<double>(number_of_processes_);
    double subpopulation_start = static_cast<double>(process_rank_)
      * subpopulation_size;
    double subpopulation_finish = static_cast<double>(process_rank_ + 1)
      * subpopulation_size;
    // Evaluate index!
    index_first_ = ceil(subpopulation_start);
    index_last_ = ceil(subpopulation_finish) - 1;
    // HACK! try to deal with double rounding unstability.
    if (process_rank_ + 1 == number_of_processes_)
      index_last_ = total_population_ - 1;
    subpopulation_ = index_last_ - index_first_ + 1;
    if (subpopulation_ == 0) return kError;
    // //debug section
    //         printf("%lli-%lli ", index_first_, index_last_);
    //         popul_eval += index_last_ - index_first_ + 1;
    //       }
    //       if (popul != popul_eval)
    //         printf("procs %i for popul %lli (%lli)\n",
    //                procs, popul, popul_eval);
    //     }
    //   }
    // }
    // //end of debug section
    current_generation_ = 0;
    x_current_vectors_.resize(subpopulation_);
    evaluated_fitness_for_current_vectors_.resize(subpopulation_);
    for (auto &x : x_current_vectors_) x.resize(dimension_);
    // //debug
    // printf("%i, x1 size = %li \n", process_rank_, x_current_vectors_.size());
    x_lbound_.resize(dimension_);
    x_ubound_.resize(dimension_);
    mutation_F_.resize(subpopulation_);
    crossover_CR_.resize(subpopulation_);
    return kDone;
  }  // end of void SubPopulation::Test()
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  int SubPopulation::CreateInitialPopulation() {
    for (auto &x : x_current_vectors_)
      for (int i = 0; i < dimension_; ++i) {
        if (x_lbound_[i] > x_ubound_[i]) return kError;
        x[i] = rand(x_lbound_[i], x_ubound_[i]);                            // NOLINT
      }
    // //debug
    // for (auto x : x_current_vectors_[0]) printf("%g ",x);
    return kDone;
  }  // end of int SubPopulation::CreateInitialPopulation()
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  int SubPopulation::EvaluateCurrentVectors() {
    evaluated_fitness_for_current_vectors_.clear();
    for (long long i = 0; i < subpopulation_; ++i) {                        // NOLINT
      auto tmp = std::make_pair(FitnessFunction(x_current_vectors_[i]), i);
      evaluated_fitness_for_current_vectors_.push_back(tmp);
    }
    if (find_minimum_)
      evaluated_fitness_for_current_vectors_
        .sort([](const std::pair<double, long long>& a,                     // NOLINT
                 const std::pair<double, long long>& b) {                   // NOLINT
                return a.first < b.first;
              });
    else
      evaluated_fitness_for_current_vectors_
        .sort([](const std::pair<double, long long>& a,                     // NOLINT
                 const std::pair<double, long long>& b) {                   // NOLINT
                return a.first > b.first;
              });
    // //debug
    // printf("\n After ");
    // for (auto val : evaluated_fitness_for_current_vectors_)
    //   printf("%g ", val.first);
    return kDone;
  }  // end of int SubPopulation::EvaluateCurrentVectors()
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  int SubPopulation::SetAllBounds(double lbound, double ubound) {
    if (lbound >= ubound) return kError;
    for (auto &x : x_lbound_) x = lbound;
    for (auto &x : x_ubound_) x = ubound;
    return kDone;
  }  // end of int SubPopulation::SetAllBounds(double lbound, double ubound)
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  double SubPopulation::randn(double mean, double stddev) {
    std::normal_distribution<double> distribution(mean, stddev);
    return distribution(generator_);
  }  // end of double SubPopulation::randn(double mean, double stddev)
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  double SubPopulation::randc(double location, double scale) {
    std::cauchy_distribution<double> distribution(location, scale);
    return distribution(generator_);
  }  // end of double SubPopulation::randc(double location, double scale)
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  long long SubPopulation::randint(long long lbound, long long ubound) {    // NOLINT
    std::uniform_int_distribution<long long> distribution(lbound, ubound);  // NOLINT
    return distribution(generator_);
  }
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  double SubPopulation::rand(double lbound, double ubound) {                // NOLINT
    std::uniform_real_distribution<double> distribution(lbound, ubound);
    return distribution(generator_);
  }  // end of double rand(double lbound, double ubound)
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  void SubPopulation::CheckRandom() {
    std::map<int, int> hist_n, hist_c, hist_int, hist;
    int factor = 1000;
    for (int n = 0; n < 100*factor; ++n) {
      ++hist_n[std::round(randn(0, 2))];
      ++hist_c[std::round(randc(0, 2))];
      ++hist_int[std::round(randint(0, 10))];
      ++hist[std::round(rand(0, 15))];                                      // NOLINT 
    }
    if (process_rank_ == 0) {
      printf("Normal (0,2)\n");
      for (auto p : hist_n) {
        if (p.second > factor)
          printf("%i: % 4i %s\n", process_rank_, p.first,
                 std::string(p.second/factor, '*').c_str());
      }  // end of for p in hist_n
      printf("Cauchy (0,2)\n");
      for (auto p : hist_c) {
        if (p.second > factor)
          printf("%i: % 4i %s\n", process_rank_, p.first,
                 std::string(p.second/factor, '*').c_str());
      }  // end of for p in hist_c
      printf("Uniform int (0,10)\n");
      for (auto p : hist_int) {
        if (p.second > factor)
          printf("%i: % 4i %s\n", process_rank_, p.first,
                 std::string(p.second/factor, '*').c_str());
      }  // end of for p in hist_int
      printf("Uniform double (0,15)\n");
      for (auto p : hist) {
        if (p.second > factor)
          printf("%i: % 4i %s\n", process_rank_, p.first,
                 std::string(p.second/factor, '*').c_str());
      }  // end of for p in hist
    }  //  end of if current process_rank_ == 0
  }
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //

}  // end of namespace jade
