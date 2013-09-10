#ifndef SRC_NMIE_NMIE_WRAPPER_H_
#define SRC_NMIE_NMIE_WRAPPER_H_
///
/// @file   nmie-wrapper.h
/// @author Ladutenko Konstantin <kostyfisik at gmail (.) com>
/// @date   Tue Sep  3 00:40:47 2013
/// @copyright 2013 Ladutenko Konstantin
///
/// nmie-wrapper is free software: you can redistribute it and/or modify
/// it under the terms of the GNU General Public License as published by
/// the Free Software Foundation, either version 3 of the License, or
/// (at your option) any later version.
///
/// nmie-wrapper is distributed in the hope that it will be useful,
/// but WITHOUT ANY WARRANTY; without even the implied warranty of
/// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
/// GNU General Public License for more details.
///
/// You should have received a copy of the GNU General Public License
/// along with nmie-wrapper.  If not, see <http://www.gnu.org/licenses/>.
///
/// nmie-wrapper uses nmie.c from scattnlay by Ovidio Pena
/// <ovidio@bytesfall.com> as a linked library. He has an additional condition to 
/// his library:
//    The only additional condition is that we expect that all publications         //
//    describing  work using this software , or all commercial products             //
//    using it, cite the following reference:                                       //
//    [1] O. Pena and U. Pal, "Scattering of electromagnetic radiation by           //
//        a multilayered sphere," Computer Physics Communications,                  //
//        vol. 180, Nov. 2009, pp. 2348-2354.                                       //
///
/// @brief  Wrapper class around nMie function for ease of use
/// 
///
#include "ucomplex.h"
#include <vector>
namespace nmie {
  const double PI=3.14159265358979323846;
  class MultiLayerMie {
   public:
    void SetWavelength(double wavelength) {wavelength_ = wavelength;};
    void AddTargetLayer(double thickness, complex layer_index);
    void SetCoatingThickness(std::vector<double> thickness);
    void SetCoatingIndex(std::vector<complex> index);
    void SetQfaild (double Q) {Qfaild_ = Q;};
    void RunMie(double *Qext, double *Qsca, double *Qabs, double *Qbk);
    void RunMieDebug(int L, double x[], complex m[], int nTheta, double Theta[], double *Qext, double *Qsca, double *Qabs, double *Qbk, double *Qpr, double *g, double *Albedo, complex S1[], complex S2[]);

  private:
    void GenerateSizeParameter();
    void GenerateIndex();
    double wavelength_ = 1.0;
    double Qfaild_ = 100.0;
    std::vector<double> target_thickness_, coating_thickness_;
    std::vector<complex> target_index_, coating_index_;
    std::vector<double> size_parameter_;
    std::vector<complex> index_;
  };  // end of class MultiLayerMie
}  // end of namespace nmie
#endif  // SRC_NMIE_NMIE_WRAPPER_H_
