#ifndef SRC_JADE_H_
#define SRC_JADE_H_
///
/// @file   jade.h
/// @author Ladutenko Konstantin <kostyfisik at gmail (.) com>
/// @date   Thu Aug 15 19:21:57 2013
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

/// @brief JADE++ is a free (GPLv3+) high performance implementation of
/// adaptive differential evolution optimization algorithm from
/// Jingqiao Zhang and Arthur C. Sanderson book 'Adaptive Differential
/// Evolution. A Robust Approach to Multimodal Problem Optimization'
/// Springer, 2009.
#include <list>
#include <random>
#include <vector>
namespace jade {
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  /// @brief Population controlled by single MPI process.
  class SubPopulation {
   public:
    int CreateInitialPopulation();
    /// @brief Vizualize used random distributions (to do manual check).
    double (*fitness_function)(std::vector<double> x) = nullptr;
    void CheckRandom();
    int Init(long long total_population, long long dimension);
    int SetAllBounds(double lbound, double ubound);
   private:
    /// @name Population, individuals and algorithm .
    // @{
    /// @brief Total number of individuals in all subpopulations.
    long long total_population_ = 0;
    /// @brief Number of individuals in subpopulation
    long long subpopulation_ = 0;
    /// @brief All individuals are indexed. First and last index of
    /// individuals in subpopulations.
    long long index_first_ = -1, index_last_ = -1;
    /// @brief Dimension of the optimization task (number of variables to optimize).
    long long dimension_ = -1;
    /// @brief Current generation of evalution process;
    long long current_generation_ = -1;
    /// @brief Current state vectors of all individuals in subpopulation.
    std::vector<std::vector<double> > x_current_vectors_;
    /// @brief Archived best solutions (state vactors)
    std::list<std::vector<double> > archived_best_;
    /// @brief Low and upper bounds for x vectors.
    std::vector<double> x_lbound_;
    std::vector<double> x_ubound_;
    /// @brief JADE+ adaption parameter for mutation factor
    double adaptor_mutation_mu_F_ = 0.5;
    /// @brief JADE+ adaption parameter for crossover probability
    double adaptor_crossover_mu_CR_ = 0.5;
    /// @brief Individual mutation and crossover parameters for each individual.
    std::vector<double> mutation_F_, crossover_CR_;
    std::list<double> successful_mutation_parameters_S_F_;
    std::list<double> successful_crossover_parameters_S_CR_;
    // @}
    /// @name Random generation
    /// Names are in notation from Jingqiao Zhang and Arthur C. Sanderson book.
    // @{
    /// @todo Select random generator enginge for best results in DE!
    std::mt19937_64 generator_;
    // std::ranlux48 generator_;
    /// @brief randn(\mu, \sigma^2 ) denotes a random value from a normal
    /// distribution of mean \mu and variance \sigma^2
    double randn(double mean, double stddev);
    /// @brief randc(\mu, \delta ) a random value from a Cauchy distribution
    /// with location and scale parameters \mu and \delta
    double randc(double location, double scale);
    /// @brief randint(1, D) is an integer randomly chosen from 1 to D
    long long randint(long long lbound, long long ubound);
    /// @brief rand(a, b) is an uniform random number chosen from a to b
    double rand(double lbound, double ubound);
    // @}
    int process_rank_;
    int number_of_processes_;
  };  // end of class SubPopulation
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  /// @brief Error codes
  ///
  /// Error codes used with jade
  enum Errors {
    /// no error
    kDone = 0,
    /// Unspecified (pending to be described).
    kError
  };  // end of enum Errors
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
}  // end of namespace jade
#endif  // SRC_JADE_H_
