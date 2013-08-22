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
  /// @todo Replace all simple kError returns with something meangfull.
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  int SubPopulation::Selection(std::vector<double> crossover_u,
                               long long i)  {
    bool is_evaluated = false;
    double f_current;
    for (auto f : evaluated_fitness_for_current_vectors_) {
      if (f.second == i) {
        f_current = f.first;
        is_evaluated = true;
      }
    }  // end of searching of pre-evaluated fitness for current individual
    if (!is_evaluated) error_status_ = kError;
    double f_best = evaluated_fitness_for_current_vectors_.front().first;
    double f_crossover_u = FitnessFunction(crossover_u);
    bool is_success = f_crossover_u > f_current
        || f_crossover_u == f_best;  //Selected for maxima search
    if (is_find_minimum_) is_success = !is_success;
    if (!is_success) {
      // Case of current x and f were new for current generation.
      x_vectors_next_generation_[i] = x_vectors_current_[i];
      for (auto &f : evaluated_fitness_for_next_generation_) {
        if (f.second == i) f.first = f_current;
      }  // end of saving old fitness value in new generation.
    } else {  // if is_success == true
      x_vectors_next_generation_[i] = crossover_u;
      for (auto &f : evaluated_fitness_for_next_generation_) 
        if (f.second == i) f.first = f_crossover_u;
      to_be_archived_best_A_.push_back(x_vectors_current_[i]);
      successful_mutation_parameters_S_F_.push_back(mutation_F_[i]);
      successful_crossover_parameters_S_CR_.push_back(crossover_CR_[i]);
      if (process_rank_ == kOutput)
        printf("n%lli f_new=%4.2f\n",i,f_crossover_u);
      PrintSingleVector(crossover_u);
    }  // end of dealing with success crossover
    return kDone;
  } // end of int SubPopulation::Selection(std::vector<double> crossover_u);
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  int SubPopulation::ArchiveCleanUp() {
    return kDone;
  } // end of int SubPopulation:: ArchiveCleanUp();
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  int SubPopulation::Adaption()  {
    return kDone;
  } // end of int SubPopulation::Adaption();
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  int SubPopulation::RunOptimization() {
    if (process_rank_ == kOutput) printf("Start optimization..\n");
    if (error_status_) return error_status_;
    adaptor_mutation_mu_F_ = 0.5;
    adaptor_crossover_mu_CR_ = 0.5;
    archived_best_A_.clear();
    CreateInitialPopulation();
    x_vectors_next_generation_ = x_vectors_current_;
    EvaluateCurrentVectors();
    evaluated_fitness_for_next_generation_ =
      evaluated_fitness_for_current_vectors_;
    for (long long g = 0; g < total_generations_max_; ++g) {
      to_be_archived_best_A_.clear();
      successful_mutation_parameters_S_F_.clear();
      successful_crossover_parameters_S_CR_.clear();        
      //debug section
      if (process_rank_ == kOutput)
        printf("==============  Generation %lli =============\n", g);
      PrintPopulation();      
      PrintEvaluated();
      //end of debug section
      for (long long i = 0; i < subpopulation_; ++i) {
        SetCRiFi(i);
        std::vector<double> mutated_v, crossover_u;
        mutated_v = Mutation(i);
        crossover_u = Crossover(mutated_v, i);
        Selection(crossover_u, i);
      }  // end of for all individuals in subpopulation
      ArchiveCleanUp();
      Adaption();
      x_vectors_current_.swap(x_vectors_next_generation_);
      evaluated_fitness_for_current_vectors_
        .swap(evaluated_fitness_for_next_generation_);
      SortEvaluatedCurrent();
      if (error_status_) return error_status_;
    }  // end of stepping generations
    return kDone;
  }  // end of int SubPopulation::RunOptimization()
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  std::vector<double> SubPopulation::Crossover(std::vector<double> mutation_v,
                                                long long i) {
    const double CR_i = crossover_CR_[i];
    std::vector<double> crossover_u, x_current;
    crossover_u.resize(dimension_);
    x_current = x_vectors_current_.at(i);
    long long j_rand = randint(0, dimension_ - 1);
    for (long long c = 0; c < dimension_; ++c) {
      if (c == j_rand || rand(0,1) < CR_i)
        crossover_u[c] = mutation_v[c];
      else
        crossover_u[c] = x_current[c];
    }
    // //debug section
    // if (process_rank_ == kOutput)
    //   printf("x -> v -> u with CR_i=%4.2f j_rand=%lli\n", CR_i, j_rand);
    // PrintSingleVector(mutation_v);
    // PrintSingleVector(x_current);
    // PrintSingleVector(crossover_u);
    // //end of debug section
    return crossover_u;
  } // end of  std::vector<double> SubPopulation::Crossover();
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  std::vector<double> SubPopulation::Mutation(long long i) {
    std::vector<double> mutation_v, x_best_current, x_random_current,
      x_random_archive_and_current, x_current;    
    x_current = x_vectors_current_.at(i);
    x_best_current = GetXpBestCurrent();
    // //debug
    // if (process_rank_ == kOutput) printf("x_best: ");
    // PrintSingleVector(x_best_current);
    long long index_of_random_current = -1;
    x_random_current = GetXRandomCurrent(&index_of_random_current, i);
    // //debug
    // if (process_rank_ == kOutput) printf("x_random: ");
    // PrintSingleVector(x_random_current);    
    x_random_archive_and_current = GetXRandomArchiveAndCurrent(index_of_random_current, i);
    // //debug
    // if (process_rank_ == kOutput) printf("x_random with archive: ");
    // PrintSingleVector(x_random_archive_and_current);
    mutation_v.resize(dimension_);
    double F_i = mutation_F_[i];
    for (long long c = 0; c < dimension_; ++c) {
      // Mutation
      mutation_v[c] = x_current[c]
        + F_i * (x_best_current[c] - x_current[c])
        + F_i * (x_random_current[c]
                 - x_random_archive_and_current[c]);
      // Bounds control
      if (mutation_v[c] > x_ubound_[c])
        mutation_v[c] = (x_ubound_[c] + x_current[c])/2;
      if (mutation_v[c] < x_lbound_[c])
        mutation_v[c] = (x_lbound_[c] + x_current[c])/2;
    }
    // //debug section
    // int isSame = 888, isSame2 = 7777, isSame3 =11111;
    // for (long long c = 0; c < dimension_; ++c) {
    //   double norm = std::abs(mutation_v[c] - x_current[c]);
    //   double norm2 = std::abs(x_random_current[c] - x_current[c]);
    //   double norm3 = std::abs(x_random_current[c]
    //                           - x_random_archive_and_current[c]);
    //   if ( norm  > 0.0001) isSame = 0;
    //   if ( norm2  > 0.0001) isSame2 = 0;
    //   if ( norm3  > 0.0001) isSame3 = 0;
    // }
    // if (process_rank_ == kOutput) printf("mutation_v%i%i%i:  ",isSame, isSame2, isSame3);
    // PrintSingleVector(mutation_v);
    // if (process_rank_ == kOutput) printf("current _v%i%i%i:  ",isSame, isSame2, isSame3);
    // PrintSingleVector(x_current);    
    // if (process_rank_ == kOutput)
    //   printf("  -> f = %4.2f                                    F_i=%4.2f\n",
    //          FitnessFunction(mutation_v), F_i);
    // //end of debug section
    return mutation_v;
  } // end of std::vector<double> SubPopulation::Mutation();
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  std::vector<double> SubPopulation::GetXpBestCurrent() {
    const long long n_best_total = static_cast<long long>
      (floor(subpopulation_ * best_share_p_ ));
    if (n_best_total == subpopulation_) error_status_ = kError;
    long long best_n = randint(0, n_best_total);
    long long best_n_index = -1, i = 0;
    for (auto x : evaluated_fitness_for_current_vectors_) {
      if (i == best_n) {
        best_n_index = x.second;
        break;
      }
      ++i;
    }
    if (best_n_index >= x_vectors_current_.size()) error_status_ = kError; 
    return x_vectors_current_.at(best_n_index);
  }  // end of std::vector<double> SubPopulation::GetXpBestCurrent();
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  std::vector<double> SubPopulation::GetXRandomCurrent(long long *index,
                                              long long forbidden_index) {
    long long random_n = randint(0, subpopulation_-1);
    while (random_n == forbidden_index) random_n = randint(0, subpopulation_-1);
    (*index) = random_n;
    if (random_n >= x_vectors_current_.size()) error_status_ = kError; 
    return x_vectors_current_.at(random_n);
  }  // end of std::vector<double> SubPopulation::GetXRandomCurrent()
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  std::vector<double> SubPopulation::GetXRandomArchiveAndCurrent(
                   long long forbidden_index1, long long forbidden_index2) {
    long long archive_size = archived_best_A_.size();
    long long random_n = randint(0, subpopulation_ + archive_size - 1);
    while (random_n == forbidden_index1 || random_n == forbidden_index2)
      random_n = randint(0, subpopulation_ + archive_size - 1);
    if (random_n < subpopulation_) return x_vectors_current_.at(random_n);
    random_n -= subpopulation_;
    long long i = 0;
    for (auto x : archived_best_A_) {
      if (i == random_n) {
        //debug
        if (process_rank_ == kOutput) printf("Using Archive!!\n");
        return x;
      }
      ++i;
    }
    error_status_ = kError;
    std::vector<double> x;
    return x;    
  }  // end of std::vector<double> SubPopulation::GetXRandomArchiveAndCurrent()
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  int SubPopulation::SetCRiFi(long long i) {
    while (1) {
      mutation_F_[i] = randc(adaptor_mutation_mu_F_, 0.1);
      if (mutation_F_[i] > 1) {
        mutation_F_[i] = 1;
        break;
      }
      if (mutation_F_[i] > 0) break;
    }
    crossover_CR_[i] = randn(adaptor_crossover_mu_CR_,0.1);
    if (crossover_CR_[i] > 1) crossover_CR_[i] = 1;
    if (crossover_CR_[i] < 0) crossover_CR_[i] = 0;    
    return kDone;
  }  // end of int SubPopulation::SetCRiFi(long long i)
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  /// %todo Change returned kError to meangfull error code.
  int SubPopulation::Init(long long total_population, long long dimension) {// NOLINT
    total_population_  = total_population;
    if (total_population_ < 1) {
      error_status_ = kError;
      return kError;
    }
    dimension_ = dimension;
    if (dimension_ < 1) {
      error_status_ = kError;
      return kError;
    };
    MPI_Comm_rank(MPI_COMM_WORLD, &process_rank_);
    MPI_Comm_size(MPI_COMM_WORLD, &number_of_processes_);
    if (process_rank_ < 0) {
      error_status_ = kError;
      return kError;
    };
    if (number_of_processes_ < 1) {
      error_status_ = kError;
      return kError;
    };
    std::random_device rd;
    generator_.seed(rd());
    // //debug
    // CheckRandom();
    if (number_of_processes_ > total_population_) {
      error_status_ = kError;
      return kError;
    };
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
    if (subpopulation_ == 0) {
      error_status_ = kError;
      return kError;
    }
    // //debug section
    //         if (process_rank_ == kOutput) printf("%lli-%lli ", index_first_, index_last_);
    //         popul_eval += index_last_ - index_first_ + 1;
    //       }
    //       if (popul != popul_eval)
    //         if (process_rank_ == kOutput) printf("procs %i for popul %lli (%lli)\n",
    //                procs, popul, popul_eval);
    //     }
    //   }
    // }
    // //end of debug section
    current_generation_ = 0;
    x_vectors_current_.resize(subpopulation_);
    for (auto &x : x_vectors_current_) x.resize(dimension_);
    x_vectors_next_generation_.resize(subpopulation_);
    for (auto &x : x_vectors_next_generation_) x.resize(dimension_);
    evaluated_fitness_for_current_vectors_.resize(subpopulation_);
    // //debug
    // if (process_rank_ == kOutput) printf("%i, x1 size = %li \n", process_rank_, x_vectors_current_.size());
    x_lbound_.resize(dimension_);
    x_ubound_.resize(dimension_);
    mutation_F_.resize(subpopulation_);
    crossover_CR_.resize(subpopulation_);
    return kDone;
  }  // end of void SubPopulation::Test()
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  int SubPopulation::PrintPopulation() {
    if (process_rank_ == kOutput) {
      for (int i = 0; i < subpopulation_; ++i) {
        printf("n%i:", i);
        for (int c = 0; c < dimension_; ++c) {
          printf(" %5.2f ", x_vectors_current_[i][c]);
        }
        printf("\n");
      }  // end of for each individual
    }  // end of if output
    return kDone;
  }  // end of int SubPupulation::PrintPopulation()
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  int SubPopulation::PrintEvaluated() {
    if (process_rank_ == kOutput) {
      for (auto x : evaluated_fitness_for_current_vectors_)
        printf("%lli:%4.2f  ", x.second, x.first);
      printf("\n");
    }  // end of if output
    return kDone;
  }  // end of int SubPopulation::PrintEvaluated()
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  int SubPopulation::PrintSingleVector(std::vector<double> x) {
    if (process_rank_ == kOutput) {
      for (auto c : x) printf("%5.2f ", c);
      printf("\n");
    }  // end of output
    return kDone;
  }  // end of int SubPopulation::PrintSingleVector(std::vector<double> x)
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  int SubPopulation::CreateInitialPopulation() {
    for (auto &x : x_vectors_current_)
      for (int i = 0; i < dimension_; ++i) {
        if (x_lbound_[i] > x_ubound_[i]) {
          error_status_ = kError;
          return kError;
        }
        x[i] = rand(x_lbound_[i], x_ubound_[i]);                            // NOLINT
      }  // end of for each dimension
    // //debug
    // for (auto x : x_vectors_current_[0]) if (process_rank_ == kOutput) printf("%g ",x);
    return kDone;
  }  // end of int SubPopulation::CreateInitialPopulation()
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  int SubPopulation::SortEvaluatedCurrent() {
    evaluated_fitness_for_current_vectors_
      .sort([=](const std::pair<double, long long>& a,                     // NOLINT
               const std::pair<double, long long>& b) {                   // NOLINT
              bool cmp = a.first < b.first;
              if (is_find_minimum_) return cmp;
              return !cmp;
            });
    return kDone;
  }  // end of int SubPopulation::SortEvaluatedCurrent()
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  int SubPopulation::EvaluateCurrentVectors() {
    evaluated_fitness_for_current_vectors_.clear();
    for (long long i = 0; i < subpopulation_; ++i) {                        // NOLINT
      auto tmp = std::make_pair(FitnessFunction(x_vectors_current_[i]), i);
      evaluated_fitness_for_current_vectors_.push_back(tmp);
    }
    SortEvaluatedCurrent();
    // //debug
    // if (process_rank_ == kOutput) printf("\n After ");
    // for (auto val : evaluated_fitness_for_current_vectors_)
    //   if (process_rank_ == kOutput) printf("%g ", val.first);
    return kDone;
  }  // end of int SubPopulation::EvaluateCurrentVectors()
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  int SubPopulation::SetAllBounds(double lbound, double ubound) {
    if (lbound >= ubound) {
      error_status_ = kError;
      return kError;
    }
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
      if (process_rank_ == kOutput) printf("Normal (0,2)\n");
      for (auto p : hist_n) {
        if (p.second > factor)
          if (process_rank_ == kOutput) printf("%i: % 4i %s\n", process_rank_, p.first,
                 std::string(p.second/factor, '*').c_str());
      }  // end of for p in hist_n
      if (process_rank_ == kOutput) printf("Cauchy (0,2)\n");
      for (auto p : hist_c) {
        if (p.second > factor)
          if (process_rank_ == kOutput) printf("%i: % 4i %s\n", process_rank_, p.first,
                 std::string(p.second/factor, '*').c_str());
      }  // end of for p in hist_c
      if (process_rank_ == kOutput) printf("Uniform int (0,10)\n");
      for (auto p : hist_int) {
        if (p.second > factor)
          if (process_rank_ == kOutput) printf("%i: % 4i %s\n", process_rank_, p.first,
                 std::string(p.second/factor, '*').c_str());
      }  // end of for p in hist_int
      if (process_rank_ == kOutput) printf("Uniform double (0,15)\n");
      for (auto p : hist) {
        if (p.second > factor)
          if (process_rank_ == kOutput) printf("%i: % 4i %s\n", process_rank_, p.first,
                 std::string(p.second/factor, '*').c_str());
      }  // end of for p in hist
    }  //  end of if current process_rank_ == 0
  }
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //

}  // end of namespace jade
