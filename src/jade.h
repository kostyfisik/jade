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
#include<random>
namespace jade {
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  /// @brief Population controlled by single MPI process.
  class SubPopulation {
   public:
    void CheckRandom();
    void Init();
    void SetTotalPopulation(int population) {
      total_population_ = population;
    };

   private:
    /// @brief Fast access (without MPI call) to process rank of
    /// containing object. Should be set in init()
    int process_rank_;
    /// @brief Total number of individuals in population (sum of all
    /// sub-pupulations)
    int total_population_;
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
  };  // End of class SubPopulation
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
  };  // End of enum Errors
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
}  // End of namespace jade
#endif  // SRC_JADE_H_
