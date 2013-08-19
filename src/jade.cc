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
#include <map>
#include <cstdio>
#include <string>
#include <cmath>
namespace jade {
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
      ++hist[std::round(rand(0, 15))];
    }
    if (process_rank_ == 0) {
      printf("Normal (0,2)\n");
      for (auto p : hist_n) {
        if (p.second > factor) printf("%i: % 4i %s\n", process_rank_, p.first,
                                      std::string(p.second/factor, '*').c_str());
      }  // end of for p in hist_n
      printf("Cauchy (0,2)\n");
      for (auto p : hist_c) {
        if (p.second > factor) printf("%i: % 4i %s\n", process_rank_, p.first,
                                      std::string(p.second/factor, '*').c_str());
      }  // end of for p in hist_c
      printf("Uniform int (0,10)\n");
      for (auto p : hist_int) {
        if (p.second > factor) printf("%i: % 4i %s\n", process_rank_, p.first,
                                      std::string(p.second/factor, '*').c_str());
      }  // end of for p in hist_int
      printf("Uniform double (0,15)\n");
      for (auto p : hist) {
        if (p.second > factor) printf("%i: % 4i %s\n", process_rank_, p.first,
                                      std::string(p.second/factor, '*').c_str());
      }  // end of for p in hist
    }  //  end of if current process_rank_ == 0
  }
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  void SubPopulation::Init() {
    MPI_Comm_rank(MPI_COMM_WORLD, &process_rank_);
    std::random_device rd;
    generator_.seed(rd());
    CheckRandom();
  }  // End of void SubPopulation::Test()
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  double SubPopulation::randn(double mean, double stddev) {
    std::normal_distribution<double> distribution(mean, stddev);
    return distribution(generator_);
  }  // End of double SubPopulation::randn(double mean, double stddev)
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  double SubPopulation::randc(double location, double scale) {
    std::cauchy_distribution<double> distribution(location, scale);
    return distribution(generator_);
  }  // End of double SubPopulation::randc(double location, double scale)
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  long long SubPopulation::randint(long long lbound, long long ubound) {
    std::uniform_int_distribution<long long> distribution(lbound, ubound);
    return distribution(generator_);
  }
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  double SubPopulation::rand(double lbound, double ubound) {
    std::uniform_real_distribution<double> distribution(lbound, ubound);
    return distribution(generator_);
  }  // End of double rand(double lbound, double ubound)
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //

}  // End of namespace jade
