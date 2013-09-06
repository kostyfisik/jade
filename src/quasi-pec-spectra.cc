///
/// @file   optimize-cloak.cc
/// @author Ladutenko Konstantin <kostyfisik at gmail (.) com>
/// @date   Tue Sep  3 00:37:05 2013
/// @copyright 2013 Ladutenko Konstantin
///
/// optimize-cloak is free software: you can redistribute it and/or modify
/// it under the terms of the GNU General Public License as published by
/// the Free Software Foundation, either version 3 of the License, or
/// (at your option) any later version.
///
/// optimize-cloak is distributed in the hope that it will be useful,
/// but WITHOUT ANY WARRANTY; without even the implied warranty of
/// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
/// GNU General Public License for more details.
///
/// You should have received a copy of the GNU General Public License
/// along with optimize-cloak.  If not, see <http://www.gnu.org/licenses/>.
///
/// optimize-cloak uses nmie.c from scattnlay by Ovidio Pena
/// <ovidio@bytesfall.com> as a linked library. He has an additional condition to 
/// his library:
//    The only additional condition is that we expect that all publications         //
//    describing  work using this software , or all commercial products             //
//    using it, cite the following reference:                                       //
//    [1] O. Pena and U. Pal, "Scattering of electromagnetic radiation by           //
//        a multilayered sphere," Computer Physics Communications,                  //
//        vol. 180, Nov. 2009, pp. 2348-2354.                                       //
///
/// @brief Simulate scattering from dielectric sphere covered with gold/dielectric
/// double shell using scattnlay lib (or gold shell inside dielectric ball)
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
#include <string>
#include "./nmie/ucomplex.h"
#include "./nmie/nmie-wrapper.h"
#include "./nmie/Au-dispersion.h"
template<class T> inline T pow2(const T value) {return value*value;}
int main(int argc, char *argv[]) {
  try {
    // Semouchkina APPLIED PHYSICS LETTERS 102, 113506 (2013)
    double lambda_work = 3.75; // cm
    //    double f_work = 30/lambda_work; // 8 GHz
    double a = 0.75*lambda_work;  // 2.8125 cm
    //size param = 2 pi r/wl = 2pi0.75 = 4.71
    // Qsca
    // pec lum = 2.4
    // nmie/lum
    // n5 = 1.917/2.02
    // n2 = 2.43/2.52
    // n1.1 =  0.42/0.43
    // n1.2 = 1.59/1.62
    // Set common parameters for all wavelengths.
    nmie::MultiLayerMie multi_layer_mie;
    double target_shell_share = 0.01;
    double eps_re = 1.0;
    double eps_im = 1250.0;
    double n = sqrt(0.5*(sqrt(pow2(eps_re)+
                              pow2(eps_im)
                              ) + eps_re ));
    double k = sqrt(0.5*(sqrt(pow2(eps_re)+
                              pow2(eps_im)
                              ) - eps_re ));
    // multi_layer_mie.AddTargetLayer(a, {2.0, 0.0});
    multi_layer_mie.AddTargetLayer((1.0-target_shell_share)*a, {1.0, 0.0000000});
    multi_layer_mie.SetCoatingThickness({target_shell_share*a});
    multi_layer_mie.SetCoatingIndex({{n, k }});
    //multi_layer_mie.SetCoatingIndex({{index, 0.0}});
    // multi_layer_mie.SetWavelength(lambda_work);
    double Qext, Qsca, Qabs, Qbk;
    for (int i = 0; i<100; ++i) {
      lambda_work = 3.0 + (5.0-3.0)/100.0*i;
      multi_layer_mie.SetWavelength(lambda_work);
      multi_layer_mie.RunMie(&Qext, &Qsca, &Qabs, &Qbk);
      // try {
      //   multi_layer_mie.RunMie(&Qext, &Qsca, &Qabs, &Qbk);
      // } catch(const std::invalid_argument& ia) {continue;}
      printf("%g\t%g\t%g\t%g\t%g\n", lambda_work, Qext, Qsca,Qabs,Qbk);        
    }
  } catch( const std::invalid_argument& ia ) {
    // Will catch if  multi_layer_mie fails or other errors.
    std::cerr << "Invalid argument: " << ia.what() << std::endl;
  }  
  return 0;
}
