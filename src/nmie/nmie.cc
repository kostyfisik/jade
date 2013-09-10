//**********************************************************************************//
//    Copyright (C) 2009  Ovidio Pena <ovidio@bytesfall.com>                        //
//                                                                                  //
//    This file is part of scattnlay                                                //
//                                                                                  //
//    This program is free software: you can redistribute it and/or modify          //
//    it under the terms of the GNU General Public License as published by          //
//    the Free Software Foundation, either version 3 of the License, or             //
//    (at your option) any later version.                                           //
//                                                                                  //
//    This program is distributed in the hope that it will be useful,               //
//    but WITHOUT ANY WARRANTY; without even the implied warranty of                //
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                 //
//    GNU General Public License for more details.                                  //
//                                                                                  //
//    The only additional condition is that we expect that all publications         //
//    describing  work using this software , or all commercial products             //
//    using it, cite the following reference:                                       //
//    [1] O. Pena and U. Pal, "Scattering of electromagnetic radiation by           //
//        a multilayered sphere," Computer Physics Communications,                  //
//        vol. 180, Nov. 2009, pp. 2348-2354.                                       //
//                                                                                  //
//    You should have received a copy of the GNU General Public License             //
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.         //
//**********************************************************************************//
/// @copyright 2013 Ladutenko Konstantin
/// @author Ladutenko Konstantin <kostyfisik at gmail (.) com>
/// Changes:
///  2013.07.12
///   Added auto-check of result stability. Fast and dirty. nMieFast for
///   draft simulations. Renamed original nMie to nMieBase
///  
//**********************************************************************************//
// This library implements the algorithm for a multilayered sphere described by:    //
//    [1] W. Yang, "Improved recursive algorithm for light scattering by a          //
//        multilayered sphere,” Applied Optics,  vol. 42, Mar. 2003, pp. 1710-1720. //
//                                                                                  //
// You can find the description of all the used equations in:                       //
//    [2] O. Pena and U. Pal, "Scattering of electromagnetic radiation by           //
//        a multilayered sphere," Computer Physics Communications,                  //
//        vol. 180, Nov. 2009, pp. 2348-2354.                                       //
//                                                                                  //
// Hereinafter all equations numbers refer to [2]                                   //
//**********************************************************************************//
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "ucomplex.h"
#include "nmie.h"

#define round(x) ((x) >= 0 ? (int)((x) + 0.5):(int)((x) - 0.5))

// Calculate Nstop - equation (17)
int Nstop(double xL) {
  int result;

  if (xL <= 8) {
    result = round(xL + 4*pow(xL, 1/3) + 1);
  } else if (xL <= 4200) {
    result = round(xL + 4.05*pow(xL, 1/3) + 2);
  } else {
    result = round(xL + 4*pow(xL, 1/3) + 2);
  }

  return result;
}

int Nmax(int L, double x[], complex m[]) {
  int i, result;

  result = Nstop(x[L + 1]);
  for (i = 1; i <= L; i++) {
    if (result < Cabs(RCmul(x[i], m[i]))) {
      result = round(Cabs(RCmul(x[i], m[i])));
    }
    if (result < Cabs(RCmul(x[i - 1], m[i]))) {
      result = round(Cabs(RCmul(x[i - 1], m[i])));
    }
  }

  //return result + 300; 
  return result + 5;  // Tig: May be it save to ommit +5 also.
  //return result + 15;
}

// Calculate an - equation (5)
complex calc_an(int n, double XL, complex Ha, complex mL, complex PsiXL, complex ZetaXL, complex PsiXLM1, complex ZetaXLM1) {
  complex Num = Csub(Cmul(Cadd(Cdiv(Ha, mL), Complex(n/XL, 0)), PsiXL), PsiXLM1);
  complex Denom = Csub(Cmul(Cadd(Cdiv(Ha, mL), Complex(n/XL, 0)), ZetaXL), ZetaXLM1);

  return Cdiv(Num, Denom);
}

// Calculate bn - equation (6)
complex calc_bn(int n, double XL, complex Hb, complex mL, complex PsiXL, complex ZetaXL, complex PsiXLM1, complex ZetaXLM1) {
  complex Num = Csub(Cmul(Cadd(Cmul(Hb, mL), Complex(n/XL, 0)), PsiXL), PsiXLM1);
  complex Denom = Csub(Cmul(Cadd(Cmul(Hb, mL), Complex(n/XL, 0)), ZetaXL), ZetaXLM1);

  return Cdiv(Num, Denom);
}

// Calculates S1_n - equation (25a)
complex calc_S1_n(int n, complex an, complex bn, double Pin, double Taun) {
  return RCmul((double)(n + n + 1)/(double)(n*n + n), Cadd(RCmul(Pin, an), RCmul(Taun, bn)));
}

// Calculates S2_n - equation (25b) (it's the same as (25a), just switches Pin and Taun)
complex calc_S2_n(int n, complex an, complex bn, double Pin, double Taun) {
  return calc_S1_n(n, an, bn, Taun, Pin);
}

//**********************************************************************************//
// This function calculates the actual scattering parameters and amplitudes         //
//                                                                                  //
// Input parameters:                                                                //
//   L: Number of layers                                                            //
//   x: Array containing the size parameters of the layers [1..L]                   //
//   m: Array containing the relative refractive indexes of the layers [1..L]       //
//   nTheta: Number of scattering angles                                            //
//   Theta: Array containing all the scattering angles where the scattering         //
//          amplitudes will be calculated                                           //
//                                                                                  //
// Output parameters:                                                               //
//   Qext: Efficiency factor for extinction                                         //
//   Qsca: Efficiency factor for scattering                                         //
//   Qabs: Efficiency factor for absorption (Qabs = Qext - Qsca)                    //
//   Qbk: Efficiency factor for backscattering                                      //
//   Qpr: Efficiency factor for the radiation pressure                              //
//   g: Asymmetry factor (g = (Qext-Qpr)/Qsca)                                      //
//   Albedo: Single scattering albedo (Albedo = Qsca/Qext)                          //
//   S1, S2: Complex scattering amplitudes                                          //
//                                                                                  //
// Return value:                                                                    //
//   Number of multipolar expansion terms used for the calculations                 //
//
// double expansion_corrector -- try to correct stability.
//**********************************************************************************//

int nMieBase(int L, double x[], complex m[], int nTheta, double Theta[], double *Qext, double *Qsca, double *Qabs, double *Qbk, double *Qpr, double *g, double *Albedo, complex S1[], complex S2[], int expansion_corrector) {
  int n_max = Nmax(L, x, m) + expansion_corrector;

  complex an, bn, anP1, bnP1, Qbktmp;

//  complex D1_lmlx[n_max + 2][L + 1], D1_lmlxM1[n_max + 2][L + 1];
//  complex **D1_lmlx, **D1_lmlxM1;
//  complex D3_lmlx[n_max + 1][L + 1], D3_lmlxM1[n_max + 1][L + 1];
//  complex D1XL[n_max + 2], D3XL[n_max + 1];
//  complex PsiZeta_lmlx[n_max + 1][L + 1], PsiZeta_lmlxM1[n_max + 1][L + 1];
//  complex *PsiXL, *ZetaXL, *PsiZetaXL;
//  complex Q[n_max + 1][L + 1];
//  complex Hb[n_max + 1][L + 1];//, Hb[n_max + 1][L + 1];
//  double Pi[n_max + 1][nTheta], Tau[n_max + 1][nTheta];
  complex z1, z2;
  complex Num, Denom;
  complex G1, G2;
  complex Temp;
  double Tmp, x2;

  int n, l, t;

  // Allocate memory to the arrays
  complex **D1_lmlx = (complex **) malloc((n_max + 2)*sizeof(complex *));
  complex **D1_lmlxM1 = (complex **) malloc((n_max + 2)*sizeof(complex *));

  complex **D3_lmlx = (complex **) malloc((n_max + 1)*sizeof(complex *));
  complex **D3_lmlxM1 = (complex **) malloc((n_max + 1)*sizeof(complex *));

  complex **PsiZeta_lmlx = (complex **) malloc((n_max + 1)*sizeof(complex *));
  complex **PsiZeta_lmlxM1 = (complex **) malloc((n_max + 1)*sizeof(complex *));

  complex **Q = (complex **) malloc((n_max + 1)*sizeof(complex *));

  complex **Ha = (complex **) malloc((n_max + 1)*sizeof(complex *));
  complex **Hb = (complex **) malloc((n_max + 1)*sizeof(complex *));

  double **Pi = (double **) malloc((n_max + 1)*sizeof(double *));
  double **Tau = (double **) malloc((n_max + 1)*sizeof(double *));

  for (n = 0; n < (n_max + 2); n++) {
    D1_lmlx[n] = (complex *) malloc((L + 1)*sizeof(complex));
    D1_lmlxM1[n] = (complex *) malloc((L + 1)*sizeof(complex));
  }

  for (n = 0; n < (n_max + 1); n++) {
    D3_lmlx[n] = (complex *) malloc((L + 1)*sizeof(complex));
    D3_lmlxM1[n] = (complex *) malloc((L + 1)*sizeof(complex));

    PsiZeta_lmlx[n] = (complex *) malloc((L + 1)*sizeof(complex));
    PsiZeta_lmlxM1[n] = (complex *) malloc((L + 1)*sizeof(complex));

    Q[n] = (complex *) malloc((L + 1)*sizeof(complex));

    Ha[n] = (complex *) malloc((L + 1)*sizeof(complex));
    Hb[n] = (complex *) malloc((L + 1)*sizeof(complex));

    Pi[n] = (double *) malloc(nTheta*sizeof(double));
    Tau[n] = (double *) malloc(nTheta*sizeof(double));
  }

  complex *D1XL = (complex *) malloc((n_max + 2)*sizeof(complex));
  complex *D3XL = (complex *) malloc((n_max + 1)*sizeof(complex));

  complex *PsiXL = (complex *) malloc((n_max + 1)*sizeof(complex));
  complex *ZetaXL = (complex *) malloc((n_max + 1)*sizeof(complex));
  complex *PsiZetaXL = (complex *) malloc((n_max + 1)*sizeof(complex));

  // Initialize the scattering parameters
  *Qext = 0;
  *Qsca = 0;
  *Qabs = 0;
  *Qbk = 0;
  Qbktmp = Complex(0, 0);
  *Qpr = 0;
  *g = 0;
  *Albedo = 0;

  // Initialize Pi, Tau and the scattering amplitudes
  for (t = 0; t < nTheta; t++) {
    Pi[0][t] = 0.0;
    Tau[0][t] = 0.0;
    S1[t] = Complex(0, 0);
    S2[t] = Complex(0, 0);
  }

  //********************************************************//
  // Calculate D1, D3 and PsiZeta for z1 in the first layer //
  //********************************************************//
  z1 = RCmul(x[1], m[1]);

  // Downward recurrence for D1 - equations (16a) and (16b)
  D1_lmlx[n_max + 1][1] = Complex(0, 0);
  for (n = n_max + 1; n > 0; n--) {
    D1_lmlx[n - 1][1] = Csub(Cdiv(Complex(n, 0), z1), Cdiv(Complex(1, 0), Cadd(D1_lmlx[n][1], Cdiv(Complex(n, 0), z1))));
  }

  // Upward recurrence for PsiZeta and D3 - equations (18a) - (18d)
  PsiZeta_lmlx[0][1] = RCmul(0.5, Csub(Complex(1, 0), Cmul(Complex(cos(2*z1.r), sin(2*z1.r)), Complex(exp(-2*z1.i), 0))));
  D3_lmlx[0][1] = Complex(0, 1);
  for (n = 1; n <= n_max; n++) {
    PsiZeta_lmlx[n][1] = Cmul(PsiZeta_lmlx[n - 1][1], Cmul(Csub(Cdiv(Complex(n, 0), z1), D1_lmlx[n - 1][1]), Csub(Cdiv(Complex(n, 0), z1), D3_lmlx[n - 1][1])));

    D3_lmlx[n][1] = Cadd(D1_lmlx[n][1], Cdiv(Complex(0, 1), PsiZeta_lmlx[n][1]));
  }

  //******************************************************************//
  // Calculate Ha and Hb in the first layer - equations (7a) and (8a) //
  //******************************************************************//
  for (n = 1; n <= n_max; n++) {
    Ha[n][1] = D1_lmlx[n][1];
    Hb[n][1] = D1_lmlx[n][1];
  }

  //*******************************************//
  // Iteration from the layer 2 to the layer L //
  //*******************************************//
  for (l = 2; l <= L; l++) {
    //**************************************************************//
    //Calculate D1, D3 and PsiZeta for z1 and z2 in the layers 2..L //
    //**************************************************************//
    z1 = RCmul(x[l], m[l]);
    z2 = RCmul(x[l - 1], m[l]);

    // Downward recurrence for D1 - equations (16a) and (16b)
    D1_lmlx[n_max + 1][l] = Complex(0, 0);
    D1_lmlxM1[n_max + 1][l] = Complex(0, 0);
    for (n = n_max + 1; n > 0; n--) {
      D1_lmlx[n - 1][l] = Csub(Cdiv(Complex(n, 0), z1), Cdiv(Complex(1, 0), Cadd(D1_lmlx[n][l], Cdiv(Complex(n, 0), z1))));
      D1_lmlxM1[n - 1][l] = Csub(Cdiv(Complex(n, 0), z2), Cdiv(Complex(1, 0), Cadd(D1_lmlxM1[n][l], Cdiv(Complex(n, 0), z2))));
    }

    // Upward recurrence for PsiZeta and D3 - equations (18a) - (18d)
    PsiZeta_lmlx[0][l] = RCmul(0.5, Csub(Complex(1, 0), Cmul(Complex(cos(2*z1.r), sin(2*z1.r)), Complex(exp(-2*z1.i), 0))));
    PsiZeta_lmlxM1[0][l] = RCmul(0.5, Csub(Complex(1, 0), Cmul(Complex(cos(2*z2.r), sin(2*z2.r)), Complex(exp(-2*z2.i), 0))));

    D3_lmlx[0][l] = Complex(0, 1);
    D3_lmlxM1[0][l] = Complex(0, 1);

    for (n = 1; n <= n_max; n++) {
      PsiZeta_lmlx[n][l] = Cmul(PsiZeta_lmlx[n - 1][l], Cmul(Csub(Cdiv(Complex(n, 0), z1), D1_lmlx[n - 1][l]), Csub(Cdiv(Complex(n, 0), z1), D3_lmlx[n - 1][l])));
      PsiZeta_lmlxM1[n][l] = Cmul(PsiZeta_lmlxM1[n - 1][l], Cmul(Csub(Cdiv(Complex(n, 0), z2), D1_lmlxM1[n - 1][l]), Csub(Cdiv(Complex(n, 0), z2), D3_lmlxM1[n - 1][l])));

      D3_lmlx[n][l] = Cadd(D1_lmlx[n][l], Cdiv(Complex(0, 1), PsiZeta_lmlx[n][l]));
      D3_lmlxM1[n][l] = Cadd(D1_lmlxM1[n][l], Cdiv(Complex(0, 1), PsiZeta_lmlxM1[n][l]));
    }

    //******************************************//
    //Calculate Q, Ha and Hb in the layers 2..L //
    //******************************************//

    // Upward recurrence for Q - equations (19a) and (19b)
    Num = RCmul(exp(-2*(z1.i - z2.i)), Complex(cos(-2*z2.r) - exp(-2*z2.i), sin(-2*z2.r)));
    Denom = Complex(cos(-2*z1.r) - exp(-2*z1.i), sin(-2*z1.r));
    Q[0][l] = Cdiv(Num, Denom);

    for (n = 1; n <= n_max; n++) {
      Num = Cmul(Cadd(Cmul(z1, D1_lmlx[n][l]), Complex(n, 0)), Csub(Complex(n, 0), Cmul(z1, D3_lmlx[n - 1][l])));
      Denom = Cmul(Cadd(Cmul(z2, D1_lmlxM1[n][l]), Complex(n, 0)), Csub(Complex(n, 0), Cmul(z2, D3_lmlxM1[n - 1][l])));

      Tmp = (x[l - 1]*x[l - 1])/(x[l]*x[l]);

      Q[n][l] = Cdiv(Cmul(RCmul(Tmp, Q[n - 1][l]), Num), Denom);
    }

    // Upward recurrence for Ha and Hb - equations (7b), (8b) and (12) - (15)
    for (n = 1; n <= n_max; n++) {
      //Ha
      G1 = Csub(Cmul(m[l], Ha[n][l - 1]), Cmul(m[l - 1], D1_lmlxM1[n][l]));
      G2 = Csub(Cmul(m[l], Ha[n][l - 1]), Cmul(m[l - 1], D3_lmlxM1[n][l]));

      Temp = Cmul(Q[n][l], G1);

      Num = Csub(Cmul(G2, D1_lmlx[n][l]), Cmul(Temp, D3_lmlx[n][l]));
      Denom = Csub(G2, Temp);

      Ha[n][l] = Cdiv(Num, Denom);

      //Hb
      G1 = Csub(Cmul(m[l - 1], Hb[n][l - 1]), Cmul(m[l], D1_lmlxM1[n][l]));
      G2 = Csub(Cmul(m[l - 1], Hb[n][l - 1]), Cmul(m[l], D3_lmlxM1[n][l]));

      Temp = Cmul(Q[n][l], G1);

      Num = Csub(Cmul(G2, D1_lmlx[n][l]), Cmul(Temp, D3_lmlx[n][l]));
      Denom = Csub(G2, Temp);

      Hb[n][l] = Cdiv(Num, Denom);
    }
  }

  //************************************//
  //Calculate D1, D3 and PsiZeta for XL //
  //************************************//
  z1 = Complex(x[L], 0);

  // Downward recurrence for D1XL - equations (16a) and (16b)
  D1XL[n_max + 1] = Complex(0, 0);
  for (n = n_max + 1; n > 0; n--) {
    D1XL[n - 1] = Csub(Cdiv(Complex(n, 0), z1), Cdiv(Complex(1, 0), Cadd(D1XL[n], Cdiv(Complex(n, 0), z1))));
  }

  //Upward recurrence for PsiXL, ZetaXL and D3XL - equations (18b), (18d) and (20a) - (21b)
  PsiXL[0] = Complex(sin(z1.r), 0);
  ZetaXL[0] = Complex(sin(z1.r), -cos(z1.r));

  PsiZetaXL[0] = RCmul(0.5, Csub(Complex(1, 0), Cmul(Complex(cos(2*z1.r), sin(2*z1.r)), Complex(exp(-2*z1.i), 0))));
  D3XL[0] = Complex(0, 1);
  for (n = 1; n <= n_max; n++) {
    PsiXL[n] = Cmul(PsiXL[n - 1], Csub(Cdiv(Complex(n, 0), z1), D1XL[n - 1]));
    ZetaXL[n] = Cmul(ZetaXL[n - 1], Csub(Cdiv(Complex(n, 0), z1), D3XL[n - 1]));

    PsiZetaXL[n] = Cmul(PsiZetaXL[n - 1], Cmul(Csub(Cdiv(Complex(n, 0), z1), D1XL[n - 1]), Csub(Cdiv(Complex(n, 0), z1), D3XL[n - 1])));
    D3XL[n] = Cadd(D1XL[n], Cdiv(Complex(0, 1), PsiZetaXL[n]));
  }

  //*********************************************************************//
  // Finally, we calculate the an and bn coefficients and the resulting  //
  // scattering parameters                                               //
  //*********************************************************************//
  x2 = x[L]*x[L];

  anP1 = calc_an(1, x[L], Ha[1][L], m[L], PsiXL[1], ZetaXL[1], PsiXL[0], ZetaXL[0]);
  bnP1 = calc_bn(1, x[L], Hb[1][L], m[L], PsiXL[1], ZetaXL[1], PsiXL[0], ZetaXL[0]);
  for (n = 1; n < n_max; n++) {
    an = anP1;
    bn = bnP1;

    anP1 = calc_an(n + 1, x[L], Ha[n + 1][L], m[L], PsiXL[n + 1], ZetaXL[n + 1], PsiXL[n], ZetaXL[n]);
    bnP1 = calc_bn(n + 1, x[L], Hb[n + 1][L], m[L], PsiXL[n + 1], ZetaXL[n + 1], PsiXL[n], ZetaXL[n]);

    // Equation (27)
    *Qext = *Qext + (double)(n + n + 1)*(an.r + bn.r);
    // Equation (28)
    *Qsca = *Qsca + (double)(n + n + 1)*(an.r*an.r + an.i*an.i + bn.r*bn.r + bn.i*bn.i);
    // Equation (29)
    *Qpr = *Qpr + ((n*(n + 2)/(n + 1))*((Cadd(Cmul(an, Conjg(anP1)), Cmul(bn, Conjg(bnP1)))).r) + ((double)(n + n + 1)/(n*(n + 1)))*(Cmul(an, Conjg(bn)).r));

    // Equation (33)
    Qbktmp = Cadd(Qbktmp, RCmul((double)((n + n + 1)*(1 - 2*(n % 2))), Csub(an, bn)));

    //****************************************************//
    // Calculate Pi_n and Tau_n for all values of Theta   //
    // Equations (26a) - (26c)                            //
    //****************************************************//
    for (t = 0; t < nTheta; t++) {
      Pi[n][t] = ((n == 1) ? 1.0 : (((double)(n + n - 1)*cos(Theta[t])*Pi[n - 1][t] - (double)n*Pi[n - 2][t])/((double)(n - 1))));
      Tau[n][t] = (double)n*cos(Theta[t])*Pi[n][t] - (double)(n + 1)*Pi[n - 1][t];

      S1[t] = Cadd(S1[t], calc_S1_n(n, an, bn, Pi[n][t], Tau[n][t]));
      S2[t] = Cadd(S2[t], calc_S2_n(n, an, bn, Pi[n][t], Tau[n][t]));
    }
  }

  *Qext = 2*(*Qext)/x2;                                 // Equation (27)
  printf("Ssca=%g\n",*Qsca);
  *Qsca = 2*(*Qsca)/x2;                                 // Equation (28)
  *Qpr = *Qext - 4*(*Qpr)/x2;                           // Equation (29)

  *Qabs = *Qext - *Qsca;                                // Equation (30)
  *Albedo = *Qsca / *Qext;                              // Equation (31)
  *g = (*Qext - *Qpr) / *Qsca;                          // Equation (32)

  *Qbk = (Qbktmp.r*Qbktmp.r + Qbktmp.i*Qbktmp.i)/x2;    // Equation (33)

  // Free the memory used for the arrays
  for (n = 0; n < (n_max + 2); n++) {
    free(D1_lmlx[n]);
    free(D1_lmlxM1[n]);
  }

  for (n = 0; n < (n_max + 1); n++) {
    free(D3_lmlx[n]);
    free(D3_lmlxM1[n]);

    free(PsiZeta_lmlx[n]);
    free(PsiZeta_lmlxM1[n]);

    free(Q[n]);

    free(Ha[n]);
    free(Hb[n]);

    free(Pi[n]);
    free(Tau[n]);
  }

  free(D1_lmlx);
  free(D1_lmlxM1);

  free(D3_lmlx);
  free(D3_lmlxM1);

  free(PsiZeta_lmlx);
  free(PsiZeta_lmlxM1);

  free(Q);

  free(Ha);
  free(Hb);

  free(Pi);
  free(Tau);

  free(D1XL);
  free(D3XL);

  free(PsiXL);
  free(ZetaXL);
  free(PsiZetaXL);

  return n_max;
}

int isClose (double A, double B) {
  double epsilon = 1e-13;
  if (B == 0 && A == 0) return 1;
  if (A*(1+epsilon) > B && B*(1+epsilon) > A) return 1;
  return 0;
}
int nMie(int L, double x[], complex m[], int nTheta, double Theta[],
         double *Qext, double *Qsca, double *Qabs, double *Qbk, double *Qpr,
                  double *g, double *Albedo, complex S1[], complex S2[]) {
  double uQext, uQsca, uQabs, uQbk;
  double uQext_old, uQsca_old, uQabs_old, uQbk_old;
  double corrector = 15.0;
  int terms;
  terms = nMieBase(L, x, m, nTheta, Theta, &uQext, &uQsca, &uQabs, &uQbk, Qpr, g, Albedo, S1,S2,corrector);
  uQext_old = uQext;  uQsca_old = uQsca;   uQabs_old = uQabs;   uQbk_old = uQbk;      
  int unstable = 1;
  while (unstable) {
    corrector +=1.0;
    if (corrector > 400.0) {
      terms = 0;
      break;
    }    
    terms = nMieBase(L, x, m, nTheta, Theta, &uQext, &uQsca, &uQabs, &uQbk, Qpr, g, Albedo, S1,S2,corrector);
    if (isClose(uQext_old,uQext) &&
        isClose(uQabs_old,uQabs) &&
        isClose(uQsca_old,uQsca) &&
        isClose(uQbk_old,uQbk)  ) unstable = 0;
    uQext_old = uQext;  uQsca_old = uQsca;   uQabs_old = uQabs;   uQbk_old = uQbk;      
    
    /* else printf("+"); */
  }
  if (terms) {
    *Qext = uQext;
    *Qsca = uQsca;
    *Qabs = uQabs;
    *Qbk = uQbk;    
  /* } else { */
  /*   printf("\n####### Unstable ! #####################\n"); */
  }
  return terms;
}
int isCloseFast (double A, double B) {
  double epsilon = 1e-6;
  if (B == 0 && A == 0) return 1;
  if (A*(1+epsilon) > B && B*(1+epsilon) > A) return 1;
  return 0;
}

int nMieFast(int L, double x[], complex m[], int nTheta, double Theta[],
         double *Qext, double *Qsca, double *Qabs, double *Qbk, double *Qpr,
                  double *g, double *Albedo, complex S1[], complex S2[]) {
  double uQext, uQsca, uQabs, uQbk;
  double uQext_old, uQsca_old, uQabs_old, uQbk_old;
  double corrector = 1.0;
  nMieBase(L, x, m, nTheta, Theta, &uQext, &uQsca, &uQabs, &uQbk, Qpr, g, Albedo, S1,S2,corrector);
  uQext_old = uQext;  uQsca_old = uQsca;   uQabs_old = uQabs;   uQbk_old = uQbk;      
  int unstable = 1;
  int terms;
  while (unstable) {
    corrector +=10.0;
    if (corrector > 120.0) {
      terms = 0;
      break;
    }    
    terms = nMieBase(L, x, m, nTheta, Theta, &uQext, &uQsca, &uQabs, &uQbk, Qpr, g, Albedo, S1,S2,corrector);
    if (isClose(uQext_old,uQext) &&
        isClose(uQabs_old,uQabs) &&
        isClose(uQsca_old,uQsca) &&
        isClose(uQbk_old,uQbk)  ) unstable = 0;
    /* else printf("+"); */
    uQext_old = uQext;  uQsca_old = uQsca;   uQabs_old = uQabs;   uQbk_old = uQbk;     
  }
  if (terms) {
    *Qext = uQext;
    *Qsca = uQsca;
    *Qabs = uQabs;
    *Qbk = uQbk;    
  /* } else { */
  /*   printf("\n####### Unstable ! #####################\n"); */
  }
  return terms;
}
