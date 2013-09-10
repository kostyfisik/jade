#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include "./nmie/ucomplex.h"
#include "./nmie/nmie.h"
#define MAXLAYERS 1100
#define MAXTHETA 800
#define PI 3.14159

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
  int has_comment = 0;
  int i, j, L;
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
      has_comment = 1;
    } else { i++; }
  }

  if (nt < 0) {
    printf("Error reading Theta.\n");
    return -1;
  } else if (nt == 1) {
    Theta[0] = ti*PI/180.0;
  } else {
    for (i = 0; i < nt; i++) {
      Theta[i] = (ti + (double)i*(tf - ti)/(nt - 1))*PI/180.0;
    }
  }

  nMie(L, x, m, nt, Theta, &Qext, &Qsca, &Qabs, &Qbk, &Qpr, &g, &Albedo, S1, S2);

  if (has_comment) {
    printf("%6s, %+.5e, %+.5e, %+.5e, %+.5e, %+.5e, %+.5e, %+.5e\n", comment, Qext, Qsca, Qabs, Qbk, Qpr, g, Albedo);
  } else {
    printf("%+.5e, %+.5e, %+.5e, %+.5e, %+.5e, %+.5e, %+.5e\n", Qext, Qsca, Qabs, Qbk, Qpr, g, Albedo);
  }

  if (nt > 0) {
    printf(" Theta,         S1.r,         S1.i,         S2.r,         S2.i\n");

    for (i = 0; i < nt; i++) {
      printf("%6.2f, %+.5e, %+.5e, %+.5e, %+.5e\n", Theta[i]*180.0/PI, S1[i].r, S1[i].i, S2[i].r, S2[i].i);
    }
  }
}


