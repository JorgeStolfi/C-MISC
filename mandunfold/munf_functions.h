/* h-functions for {mandunfold.c} */
/* Last edited on 2008-05-23 16:36:31 by stolfi */

#ifndef munf_functions_H
#define munf_functions_H

#include <math.h>
#include <complex.h>

typedef double complex complex_map_t(double complex c);
  /* Type of a mandelbrot unrolling function. */
  
void munf_get_function
  ( int func, 
    complex_map_t **h, 
    char **h_desc, 
    double complex *ctr, 
    double complex *rad
  );
  /* Gets the procedure, description, and default plot window for 
    the unrolling function number {func}. */

/* FUNCTIONS AND THEIR ATTRIBUTES */

/* Identity function - {W} same as {M}. */
extern char *h00_desc;
extern double complex h00_ctr;
extern double complex h00_rad;
double complex h00(double complex c);

/* Quadratic; unfolds the cardioid of {M#0} to a disk on the real axis. */
extern char *h01_desc;
extern double complex h01_ctr;
extern double complex h01_rad;
double complex h01(double complex c);

/* Quartic; unfolds the cardioids of {M#0} and {M#1}. */
extern char *h02_desc;
extern double complex h02_ctr;
extern double complex h02_rad;
double complex h02(double complex c);

/* Cubic; unfolds the cardioids of {M#0} and {M#1}. */
extern char *h03_desc;
extern double complex h03_ctr;
extern double complex h03_rad;
double complex h03(double complex c);

/* Cosine; unfolds the cardioid of {M#0}. */
extern char *h04_desc;
extern double complex h04_ctr;
extern double complex h04_rad;
double complex h04(double complex c);
  
/* Gaussian bell; unfolds the cardioid of {M#0} and moves wtip to {oo}. */
extern char *h05_desc;
extern double complex h05_ctr;
extern double complex h05_rad;
double complex h05(double complex c);

/* Rational bell; unfolds the cardioid of {M#0} and moves wtip to {oo}. */
extern char *h06_desc;
extern double complex h06_ctr;
extern double complex h06_rad;
double complex h06(double complex c);

/* Cosine; unfolds the cardioid of {M#0} and mirrors about wtip. */
extern char *h07_desc;
extern double complex h07_ctr;
extern double complex h07_rad;
double complex h07(double complex c);

/* Quadratic; unfolds the cardioid of {M#0} and puts figure upright. */
extern char *h08_desc;
extern double complex h08_ctr;
extern double complex h08_rad;
double complex h08(double complex c);

#endif
