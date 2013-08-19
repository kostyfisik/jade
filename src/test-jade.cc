// Using Doxygen 1.8.0 (with Markdown)
///
/// @file   test-jade.cc
/// @author Ladutenko Konstantin <kostyfisik at gmail (.) com>
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
/// @date   Wed Aug 14 13:40:37 2013
/// @brief  Test of JADE++ lib, Doxygen mainpage description.
#include <mpi.h>
#include <cstdio>
#include "./jade.h"
/// @brief Run tests of JADE++.
///
/// @param argc
/// @param argv well known input parameters.
///
/// @return Zero by default.
int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);
  int done_status = jade::kDone;
  jade::SubPopulation sub_population;
  // Settings for optimization algorithm;
  int total_population = 10;  // Total number of individuals in population.
  if (sub_population.Init(total_population) == jade::kDone) {
    printf("Optimize\n");
  }
  // while (1) {  // use break to report error with done_status.
  //   done_status = halo_exchange_process.Init(argc, argv);
  //   if (done_status != jade::kDone) break;
  //   // Split total simualation domain into subdomains for parallel
  //   // execution. Each subdomain is simulated on its own MPI processes
  //   // (physicaly it would be a single core, processor, node etc.).
  //   done_status = halo_exchange_process.RunDecomposition();
  //   if (done_status != jade::kDone) break;
  //   done_status = halo_exchange_process.InitSimulation();
  //   if (done_status != jade::kDone) break;
  //   done_status = halo_exchange_process.RunSimulation();
  //   break;
  // }  // end of while breaked with errors
  // if (done_status != jade::kDone
  //     && done_status != jade::kErrorProcessNotInGrid)
  //   MPI_Abort(MPI_COMM_WORLD, done_status);
  MPI_Finalize();
  return done_status;
}  // end of main
// ************************************************************************* //
/// @page ChangeLog
/// #ChangeLog
/// ## Version 0.0.1
/// - CMake files are configured for ease of use.
// ************************************************************************* //
// ************************************************************************* //
/// @mainpage JADE++ Documentation
///
/// [OpenMPI]: http://www.open-mpi.org
/// [Metamaterials Laboratory]: http://phoi.ifmo.ru/metamaterials
/// [NRU ITMO]: http://en.ifmo.ru "National Research University ITMO"
/// [CMake]: http://www.cmake.org "CMake"
/// [Ioffe Institute]: http://www.ioffe.ru/index_en.html
///                    "Ioffe Physical Technical Instute"
///
/// #Parallel adaptive differential evolution software
///
/// JADE++ is a free (GPLv3+) high performance implementation of
/// adaptive differential evolution optimization algorithm from
/// Jingqiao Zhang and Arthur C. Sanderson book 'Adaptive Differential
/// Evolution. A Robust Approach to Multimodal Problem Optimization'
/// Springer, 2009.
///
/// JADE++ was designed to run efficiently with:
/// - multi-core processors
/// - multiprocessor systems
/// - clusters and supercomputers
///
/// JADE++ software is developed by [Metamaterials Laboratory] of
/// Photonics and Optical Informatics Department of [NRU ITMO]. It is
/// also supported by [Ioffe Institute].
///
/// Contact Ladutenko Konstantin <kostyfisik at gmail (.)  com>
/// with any questions about JADE++.
///
/// #Prerequisites
/// - **MPI** (to use JADE++, to compile it) <br><br>
///   JADE++ needs some MPI realization to be compiled and to be
///   run. It was developed using [OpenMPI]. For Debian/Ubuntu systems
///   it can be installed with:
///
///         # apt-get install openmpi-bin openmpi-doc libopenmpi-dev
///
///   <br>
/// - **And [CMake] build system** (to compile JADE++) <br><br>
///   Use it on Linux/Unix and other operation systems.
