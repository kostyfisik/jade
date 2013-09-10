#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include "./nmie/ucomplex.h"
#include "./nmie/nmie.h"
#include "./nmie/nmie-wrapper.h"
#define MAXLAYERS 1100
#define MAXTHETA 800
#define pi 3.14159265358979323846

//***********************************************************************************//
// This is the main function of 'scattnlay', here we read the parameters as          //
// arguments passed to the program which should be executed with the following       //
// syntaxis:                                                                         //
// ./scattnlay -l Layers x1 m1.r m1.i [x2 m2.r m2.i ...] [-t ti tf nt] [-c comment]  //
//                                                                                   //
// When all the parameters were correctly passed we setup the integer L (the         //
// number of layers) and the arrays x and m, containing the size parameters and      //
// refractive indexes of the layers, respectively and call the function nMie.        //
// If the calculation is successful the results are printed with the following       //
// format:                                                                           //
//                                                                                   //
//    * If no comment was passed:                                                    //
//        'Qext, Qsca, Qabs, Qbk, Qpr, g, Albedo'                                    //
//                                                                                   //
//    * If a comment was passed:                                                     //
//        'comment, Qext, Qsca, Qabs, Qbk, Qpr, g, Albedo'                           //
//***********************************************************************************//
int main(int argc, char *argv[]) {
  char comment[200];
  //int has_comment = 0;
  int i, j, L = 0;
  double x[MAXLAYERS];
  complex m[MAXLAYERS];
  double Theta[MAXTHETA];
  complex S1[MAXTHETA], S2[MAXTHETA];
  double Qext, Qabs, Qsca, Qbk, Qpr, g, Albedo;
  double ti = 0.0, tf = 90.0;
  int nt = 0;

  if (argc < 5) {
    printf("Insufficient parameters.\nUsage: %s -l Layers x1 m1.r m1.i [x2 m2.r m2.i ...] [-t ti tf nt] [-c comment]\n", argv[0]);
    return -1;
  }

  x[0] = 0;
  m[0].r = 0;
  m[0].i = 0;
  strcpy(comment, "");
  for (i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-l") == 0) {
      i++;
      L = atoi(argv[i]);
      if (L > MAXLAYERS) {
        printf("The number of layers (%i) is higher tan the maximum allowed (%i)\n", L, MAXLAYERS);
        return -1;
      }
      if (argc < 3*(L + 1)) {
        printf("Insufficient parameters.\nUsage: %s -l Layers x1 m1.r m1.i [x2 m2.r m2.i ...] [-t ti tf nt] [-c comment]\n", argv[0]);
        return -1;
      } else {
        for (j = 1; j <= L; j++) {
          i++;
          x[j] = atof(argv[i]);
          i++;
          m[j].r = atof(argv[i]);
          i++;
          m[j].i = atof(argv[i]);
        }
      }
    } else if (strcmp(argv[i], "-t") == 0) {
      i++;
      ti = atof(argv[i]);
      i++;
      tf = atof(argv[i]);
      i++;
      nt = atoi(argv[i]);
      if (nt > MAXTHETA) {
        printf("The number of theta values (%i) is higher tan the maximum allowed (%i)\n", nt, MAXTHETA);
        return -1;
      }
    } else if (strcmp(argv[i], "-c") == 0) {
      i++;
      strcpy(comment, argv[i]);
      //has_comment = 1;
    } else { i++; }
  }

  if (nt < 0) {
    printf("Error reading Theta.\n");
    return -1;
  } else if (nt == 1) {
    Theta[0] = ti*pi/180.0;
  } else {
    for (i = 0; i < nt; i++) {
      Theta[i] = (ti + (double)i*(tf - ti)/(nt - 1))*pi/180.0;
    }
  }
  nmie::MultiLayerMie multi_layer_mie;
  double lambda_work = 3.75; // cm
  const double a_thickness = 0.75*lambda_work;  // 2.8125 cm
  printf("x[1]=%g, x[2]=%g, a=%g\n", x[1]*lambda_work/2/pi, x[2]*lambda_work/2/pi,
         a_thickness);
  multi_layer_mie.AddTargetLayer(a_thickness, {2.0, 0.0001});
  //multi_layer_mie.AddTargetLayer(a_thickness, {1.0, 0.0});
  multi_layer_mie.AddTargetLayer((x[2]-x[1])*lambda_work/2.0/pi, {1.0, 0.0});
  multi_layer_mie.SetWavelength(lambda_work);
  multi_layer_mie.RunMie(&Qext, &Qsca, &Qabs, &Qbk);
  printf("RUN   MLM %g\t%g\t%g\t%g\n", Qext, Qsca,Qabs,Qbk);

  multi_layer_mie.RunMieDebug(L, x, m, nt, Theta, &Qext, &Qsca, &Qabs, &Qbk, &Qpr, &g, &Albedo, S1, S2);
  printf("RUN debug %g\t%g\t%g\t%g\n", Qext, Qsca,Qabs,Qbk);
  nMie(L, x, m, nt, Theta, &Qext, &Qsca, &Qabs, &Qbk, &Qpr, &g, &Albedo, S1, S2);
  printf("RUN  orig %g\t%g\t%g\t%g\n", Qext, Qsca,Qabs,Qbk);

}


