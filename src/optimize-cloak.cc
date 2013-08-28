///
/// @file   core-bi-shell.c
/// @author Ladutenko Konstantin <kostyfisik at gmail (.) com>
/// @date   Thu Jul 18 14:27:28 2013
/// @copyright 2013 Ladutenko Konstantin
///
/// core-bi-shell is free software: you can redistribute it and/or modify
/// it under the terms of the GNU General Public License as published by
/// the Free Software Foundation, either version 3 of the License, or
/// (at your option) any later version.
///
/// core-bi-shell is distributed in the hope that it will be useful,
/// but WITHOUT ANY WARRANTY; without even the implied warranty of
/// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
/// GNU General Public License for more details.
///
/// You should have received a copy of the GNU General Public License
/// along with core-bi-shell.  If not, see <http://www.gnu.org/licenses/>.
///
/// core-bi-shell uses nmie.c from scattnlay by Ovidio Pena
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
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include "./nmie/ucomplex.h"
#include "./nmie/nmie.h"
#include "./nmie/Au-dispersion.h"
const double PI = 3.141592653589793238463;
 
int main(int argc, char *argv[]) {
  int L = 3; // Core bi-shell ball.
  // All real sizes should be in nanometers!
  double WL;
  double core_radius = 70; // Ball radius.
  double shell_width1 = 20;
  double shell_width2 = 40;
  double shell_radius1 = core_radius + shell_width1;
  double shell_radius2 = core_radius + shell_width1 + shell_width2;

  //double dielectric_index_r = 4, dielectric_index_i = 0.0;
  double dielectric_index_r = 4, dielectric_index_i = 0.0000000001;
  complex core_refractive_index, shell_refractive_index1,   shell_refractive_index2;
  core_refractive_index.r = dielectric_index_r;
  core_refractive_index.i = dielectric_index_i;
  shell_refractive_index2 = core_refractive_index;

  double Qext, Qsca, Qabs, Qbk, Qpr, g, Albedo;
  int nt = 0;
  double *Theta = (double *)malloc(nt * sizeof(double));
  complex *S1 = (complex *)malloc(nt * sizeof(complex *));
  complex *S2 = (complex *)malloc(nt * sizeof(complex *));
  //Nmax in nmie.c has for (i = 1; i <= L; i++) {... x[i] ...} 
  double *x = (double *)malloc((L+1) * sizeof(core_radius));
  complex *m = (complex *)malloc((L+1) * sizeof(core_refractive_index));
  // To fit notaion from Wen Yang paper we do not use zero index.
  x[0] = 0;
  m[0].r = 0;
  m[0].i = 0;
  int sample;
  double core_size_parameter, shell_size_parameter1, shell_size_parameter2;
  FILE *fp;
  char fname [300];
  snprintf(fname, 300, "dielectric(n%g_k%g)-inner-shell(Au)-r-total(%gnm).spectra",
           dielectric_index_r, dielectric_index_i, shell_radius2);
  fp = fopen(fname, "w");
  printf("%s\n", fname);
  int terms = 0;
  long unstable = 0, iterations = 0;

  fprintf(fp,"#WL,\tQext,\tQsca,\tQabs,\tQbk\n");
  for (sample = 0; sample < dispersion_samples; sample++) {
  /* for (sample = 0; sample < 1; sample++) { */
    WL = Au_dispersion[sample][d_wl]; // wavelength.
    shell_refractive_index1.r = Au_dispersion[sample][d_n];
    shell_refractive_index1.i = Au_dispersion[sample][d_k];
    core_size_parameter = 2*PI*core_radius/WL;
    shell_size_parameter1 = 2*PI*shell_radius1/WL;
    shell_size_parameter2 = 2*PI*shell_radius2/WL;
    x[1] = core_size_parameter;
    m[1] = core_refractive_index;
    x[2] = shell_size_parameter1;
    m[2] = shell_refractive_index1;
    x[3] = shell_size_parameter2;
    m[3] = shell_refractive_index2;
    terms = nMie(L, x, m, nt, Theta,
         &Qext, &Qsca, &Qabs, &Qbk, &Qpr, &g, &Albedo,
         S1,S2);
    iterations++;
    if (terms == 0) {
      Qext = 0; Qabs = 0; Qsca = 0; Qbk = 0;
      unstable++; 
      if (unstable%2000 == 0) {printf("*"); fflush(stdout);}
      printf(".");
      //continue;
    }
    fprintf(fp,"%g\t%g\t%g\t%g\t%g\n", WL, Qext,Qsca,Qabs,Qbk);        
  }
  fclose(fp);
  if (unstable) {
    printf("%g: - ",x[1]);      
    printf("Sum is not stable (%li of %li)! Change intput parameters!\n",unstable, iterations);
  }

  free(x);
  free(m);
  free(Theta);
  free(S1);
  free(S2);
  return 0;
}
