///
/// @file   nmie-wrapper.cc
/// @author Ladutenko Konstantin <kostyfisik at gmail (.) com>
/// @date   Tue Sep  3 00:38:27 2013
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
#include "ucomplex.h"
#include "nmie-wrapper.h"
#include "nmie.h"
#include <cstdio>
#include <cstdlib>
#include <stdexcept>
#include <vector>
namespace nmie {  
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  void MultiLayerMie::AddTargetLayer(double thickness, complex layer_index) {
    if (thickness <= 0)
      throw std::invalid_argument("Layer thickness should be positive!");
    target_thickness_.push_back(thickness);
    target_index_.push_back(layer_index);
    printf("thickness=%g\n",target_thickness_.back());
  }  // end of void  MultiLayerMie::AddTargetLayer(...)
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  void MultiLayerMie::SetCoatingThickness(std::vector<double> thickness) {
    coating_thickness_.clear();
    for (auto w : thickness)
      if (w <= 0)
        throw std::invalid_argument("Coating thickness should be positive!");
      else coating_thickness_.push_back(w);
  }
  // end of void MultiLayerMie::SetCoatingThickness(...);
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  void MultiLayerMie::SetCoatingIndex(std::vector<complex> index) {
    coating_index_.clear();
    for (auto value : index) coating_index_.push_back(value);
  }  // end of void MultiLayerMie::SetCoatingIndex(std::vector<complex> index);
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  /// nMie starts layer numeration from 1 (no separation of target
  /// and coating). So the first elment (zero-indexed) is not used
  /// and has some unused value.
  void MultiLayerMie::GenerateSizeParameter() {
    size_parameter_.clear();
    size_parameter_.push_back(0.0);
    printf("sp=%g\n",size_parameter_.back());
    double radius = 0.0;
    for (auto width : target_thickness_) {
      radius += width;
      size_parameter_.push_back(2*PI*radius / wavelength_);
      printf("sp=%g\n",size_parameter_.back());            
    }
    for (auto width : coating_thickness_) {
      radius += width;
      size_parameter_.push_back(2*PI*radius / wavelength_);
    }
  }  // end of void MultiLayerMie::GenerateSizeParameter();
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  void MultiLayerMie::GenerateIndex() {
    index_.clear();
    index_.push_back({0.0, 0.0});
    for (auto index : target_index_) index_.push_back(index);
    for (auto index : coating_index_) index_.push_back(index);
  }  // end of void MultiLayerMie::GenerateIndex();
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  void MultiLayerMie::RunMie(double *Qext_out, double *Qsca_out,
                             double *Qabs_out, double *Qbk_out) {
    GenerateSizeParameter();
    GenerateIndex();
    if (size_parameter_.size() != index_.size())
      throw std::invalid_argument("Each size parameter should have only one index!");
    int L = static_cast<int>(size_parameter_.size()) - 1;
    double Qext, Qsca, Qabs, Qbk, Qpr, g, Albedo;
    int nt = 0;
    double *Theta = (double *)malloc(nt * sizeof(double));
    complex *S1 = (complex *)malloc(nt * sizeof(complex *));
    complex *S2 = (complex *)malloc(nt * sizeof(complex *));
    //Nmax in nmie.c has for (i = 1; i <= L; i++) {... x[i] ...}
    const int max=100;
    double x[max];
    complex m[max];
    x[0] = 0.0;
    m[0] = {0.0,0.0};
    for (int i = 1; i<=L; ++i) {
      x[i] = size_parameter_[i];
      m[i] = index_[i];
    }
    // double *x = &(size_parameter_.front());
    // complex *m = &(index_.front());
    // for (int i = 1; i<=L; ++i)
    //   printf("l=%i x=%g re(m)=%g im(m)=%g\n",i,x[i], m[i].r, m[i].i);
    int terms = 0;
    terms = nMie(L, x, m, nt, Theta,
                 &Qext, &Qsca, &Qabs, &Qbk, &Qpr, &g, &Albedo,
                 S1,S2);
    free(Theta);
    free(S1);
    free(S2);
    if (terms == 0) {
      *Qext_out = Qfaild_;
      *Qsca_out = Qfaild_;
      *Qabs_out = Qfaild_;
      *Qbk_out = Qfaild_;
      throw std::invalid_argument("Failed to evaluate Q!");
    }
    //printf("%g\t%g\t%g\t%g\t%g\n", wavelength_, Qext,Qsca,Qabs,Qbk);        
    *Qext_out = Qext;
    *Qsca_out = Qsca;
    *Qabs_out = Qabs;
    *Qbk_out = Qbk;
    
    // double *x = (double *)malloc((L+1) * sizeof(double));
    // complex *m = (complex *)malloc((L+1) * sizeof(complex));
    // free(x);
    // free(m);
    
    //printf ("L %i, ", L);
  }  // end of void MultiLayerMie::RunMie();
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  void MultiLayerMie::RunMieDebug(int LL, double xx[], complex mm[], int nThetaa, double Thetaa[], double *Qext, double *Qsca, double *Qabs, double *Qbk, double *Qpr, double *g, double *Albedo, complex S1a[], complex S2a[]) {
    GenerateSizeParameter();
    GenerateIndex();
    if (size_parameter_.size() != index_.size())
      throw std::invalid_argument("Each size parameter should have only one index!");
    int L = static_cast<int>(size_parameter_.size()) - 1;
    //double Qext, Qsca, Qabs, Qbk, Qpr, g, Albedo;
    int nt = 0;
    double *Theta = (double *)malloc(nt * sizeof(double));
    complex *S1 = (complex *)malloc(nt * sizeof(complex *));
    complex *S2 = (complex *)malloc(nt * sizeof(complex *));
    //Nmax in nmie.c has for (i = 1; i <= L; i++) {... x[i] ...}
    double *x = &(size_parameter_.front());
    complex *m = &(index_.front());
    // for (int i = 1; i<=L; ++i)
    //   printf("l=%i x=%g re(m)=%g im(m)=%g\n",i,x[i], m[i].r, m[i].i);
    int terms = 0;
    terms = nMie(L, x, m, nt, Theta,
                 Qext, Qsca, Qabs, Qbk, Qpr, g, Albedo,
                 S1,S2);
    // free(Theta);
    // free(S1);
    // free(S2);

    GenerateSizeParameter();
    GenerateIndex();
    if (size_parameter_.size() != index_.size())
      throw std::invalid_argument("Each size parameter should have only one index!");
    L = static_cast<int>(size_parameter_.size()) - 1;
    // double Qext, Qsca, Qabs, Qbk, Qpr,
    nt = 0;
    terms = 0;
    for (int i = 0; i<=L; ++i) {
      printf("x[i]=%g vs sp[i]=%g diff=%g\n", xx[i],size_parameter_[i],
             xx[i]-size_parameter_[i]);
      m[i] = index_[i];
    }
    //    x[1] = 4.71239;
    //x[1] = size_parameter_[1];
    //x[2] = 9.42478;
    // printf("sp=%g\n",size_parameter_.back()); 
    // printf("sp[2]=%g\n",size_parameter_[2]);
    // const double sp = size_parameter_.back();
    //size_parameter_[2] = 10.0;
    //x[2] = size_parameter_.back();
    //x[2] = 2.0;
    terms = nMie(L, x, m, nt, Theta,
                 Qext, Qsca, Qabs, Qbk, Qpr, g, Albedo,
                 S1,S2);
    //printf("%g\t%g\t%g\t%g\t%g\n", wavelength_, Qext,Qsca,Qabs,Qbk);        
    
    // double *x = (double *)malloc((L+1) * sizeof(double));
    // complex *m = (complex *)malloc((L+1) * sizeof(complex));
    // free(x);
    // free(m);
    
    //printf ("L %i, ", L);
    if (terms == 0) {
      // *Qext_out = Qfaild_;
      // *Qsca_out = Qfaild_;
      // *Qabs_out = Qfaild_;
      // *Qbk_out = Qfaild_;
      throw std::invalid_argument("Failed to evaluate Q!");
    }
  }  // end of void MultiLayerMie::RunMie();
  ///MultiLayerMie::
}  // end of namespace nmie
