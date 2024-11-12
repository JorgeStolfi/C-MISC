/* Last edited on 2008-05-23 17:13:56 by stolfi */

#include <munf_functions.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

#define W_TIP_RE (-1.40115519)
  /* The western tip of the circles on the real axis of {M#0}. */

char *h00_desc = "c";
double complex h00_ctr = 00.00000 + 00.00000*I;
double complex h00_rad = 02.01000 + 02.01000*I;
double complex h00(double complex c)
  { return c; } 

char *h01_desc = "1/4 - c^2";
double complex h01_ctr = 00.00000 + 00.00000*I;
double complex h01_rad = 01.70000 + 01.70000*I;
double complex h01(double complex c)
  { return 0.25 - c*c; } 

char *h02_desc = "2(1 - c^2)^2 - 7/4";
double complex h02_ctr = 00.00000 + 00.00000*I;
double complex h02_rad = 01.50000 + 01.50000*I;
double complex h02(double complex c)
  { double complex t = 1 - c*c;
    return 2*t*t - 1.75;
  } 

char *h03_desc = "c(c^2 - 3)/2 - 3/4";
double complex h03_ctr = 00.00000 + 00.00000*I;
double complex h03_rad = 02.50000 + 02.50000*I;
double complex h03(double complex c)
  { double complex t = c*c - 3;
    return c*t/2 - 0.75;
  } 

char *h04_desc = "cos(c) - 3/4";
double complex h04_ctr = 00.00000 + 00.00000*I;
double complex h04_rad = 09.00000 + 03.00000*I;
double complex h04(double complex c)
  { return ccos(c) - 0.75; } 
  
char *h05_desc = "(1/4 + g)*exp(-c^2/2) - g";
double complex h05_ctr = 00.00000 + 00.00000*I;
double complex h05_rad = 09.00000 + 09.00000*I;
double complex h05(double complex c)
  { double complex n = cexp(-c*c/2);
    return 0.25*n + W_TIP_RE*(1-n);
  } 

char *h06_desc = "(1/4 + g)*(1/(1 + c^2)) - g";
double complex h06_ctr = 00.00000 + 00.00000*I;
double complex h06_rad = 12.00000 + 12.00000*I;
double complex h06(double complex c)
  { double complex n = 1/(1 + c*c);
    return 0.25*n + W_TIP_RE*(1-n);
  } 

char *h07_desc = "(1/4 + g)*(1 + cos(c))/2 - g";
double complex h07_ctr = 00.00000 + 00.00000*I;
double complex h07_rad = 09.00000 + 09.00000*I;
double complex h07(double complex c)
  { double complex n = (1 + ccos(c))/2;
    return 0.25*n + W_TIP_RE*(1-n);
  } 

char *h08_desc = "c^2 + 1/4";
double complex h08_ctr = 00.00000 + 00.00000*I;
double complex h08_rad = 01.70000 + 01.70000*I;
double complex h08(double complex c)
{ return c*c + 0.25; }

void munf_get_function
  ( int func, 
    complex_map_t **h, 
    char **h_desc, 
    double complex *ctr, 
    double complex *rad
  )
  {
    switch (func)
      {
        case 0: (*h) = &h00; (*h_desc) = h00_desc; (*ctr) = h00_ctr; (*rad) = h00_rad; break;
        case 1: (*h) = &h01; (*h_desc) = h01_desc; (*ctr) = h01_ctr; (*rad) = h01_rad; break;
        case 2: (*h) = &h02; (*h_desc) = h02_desc; (*ctr) = h02_ctr; (*rad) = h02_rad; break;
        case 3: (*h) = &h03; (*h_desc) = h03_desc; (*ctr) = h03_ctr; (*rad) = h03_rad; break;
        case 4: (*h) = &h04; (*h_desc) = h04_desc; (*ctr) = h04_ctr; (*rad) = h04_rad; break;
        case 5: (*h) = &h05; (*h_desc) = h05_desc; (*ctr) = h05_ctr; (*rad) = h05_rad; break;
        case 6: (*h) = &h06; (*h_desc) = h06_desc; (*ctr) = h06_ctr; (*rad) = h06_rad; break;
        case 7: (*h) = &h07; (*h_desc) = h07_desc; (*ctr) = h07_ctr; (*rad) = h07_rad; break;
        case 8: (*h) = &h08; (*h_desc) = h08_desc; (*ctr) = h08_ctr; (*rad) = h08_rad; break;
        default: fprintf(stderr, "invalid function %d\n", func); exit (1);
      }
  }
