/* See SP1DSpline.h */
/* Last edited on 2005-06-06 11:50:58 by stolfi */

#include <SP1DSpline.h>

#include <math.h>
#include <affirm.h>
#include <nat.h>
#include <bool.h>
#include <stdio.h>
#include <stdlib.h>
#include <values.h>

/* INTERNAL PROTOTYPES */

void SP1DSpline_Coeffs_Plus(int c, int i, double *a);
  /* Returns in {a[0..d]} the coeffs of the polynomial that matches
    the Hermite-like mother spline element of index {i} in the range
    {[0 _ 1]}. */

void SP1DSpline_Coeffs_Minus(int c, int i, double *a);
  /* Returns in {a[0..d]} the coeffs of the polynomial that matches
    the Hermite-like mother spline element of index {i} in the range
    {[-1 _ 0]}. */

double SP1DSpline_Eval_Plus(int c, int i, double z);
  /* Evaluates at {z} the polynomial that matches the
    Hermite-like mother spline element of index {i} 
    in the range {[0 _ 1]}. */

/* EVALUATION */

/* The formulas below provide elements that have unit L1 norm: */

void SP1DSpline_Coeffs_Plus(int c, int i, double *a)
  {
    affirm((i >= 0) && (i <= c), "invalid index");
    switch(c)
      {
        case 0: 
          /* Unit triangle {1 - z}. */
          affirm(i == 0, "bad basis element index");
          a[0] =  +1.0; a[1] =  -1.0;
          break;
          
        case 1: 
          switch(i)
            { case 0:
                /* Cubic bell {(1 + 2·z)·(1 - z)^2}. */
                a[0] =  +1.0; a[1] =   0.0; a[2] =  -3.0; a[3] =  +2.0;
                break;
              case 1:
                /* Cubic swing {6·z·(1-z)^2}. */
                a[0] =   0.0; a[1] =  +6.0; a[2] = -12.0; a[3] =  +6.0;
                break;
              default:
                affirm(FALSE, "bad basis element index");
            }
          break;
         
        case 2: 
          switch(i)
            { case 0:
                /* Quintic bell {(1 + 3·z + 6·z^2)·(1-z)^3}. */
                a[0] =  +1.0; a[1] =   0.0; a[2] =   0.0; 
                a[3] = -10.0; a[4] = +15.0; a[5] =  -6.0;
                break;
              case 1:
                /* Quintic swing {5·z·(1 + 3·z)·(1-z)^3}. */
                a[0] =   0.0; a[1] =  +5.0; a[2] =   0.0; 
                a[3] = -30.0; a[4] = +40.0; a[5] = -15.0; 
                break;
              case 2:
                /* Quintic camel {30·z^2·(1-z)^3}. */
                a[0] =   0.0; a[1] =   0.0; a[2] = +30.0; 
                a[3] = -90.0; a[4] = +90.0; a[5] = -30.0; 
                break;
              default:
                affirm(FALSE, "bad basis element index");
            }
          break;
         
        default:
          affirm(FALSE, "bad time continuity order");
      }
  }

void SP1DSpline_Coeffs_Minus(int c, int i, double *a)
  {
    affirm((i >= 0) && (i <= c), "invalid index");
    switch(c)
      {
        case 0: 
          /* Unit triangle {z}. */
          affirm(i == 0, "bad basis element index");
          a[0] =   0.0; a[1] =  +1.0;
          break;
          
        case 1: 
          switch(i)
            { case 0:
                /* Cubic bell {(3 - 2·z)·z^2}. */
                a[0] =   0.0; a[1] =   0.0; a[2] =  +3.0; a[3] =  -2.0;
                break;
              case 1:
                /* Cubic swing {6·(z - 1)·z^2}. */
                a[0] =   0.0; a[1] =  0.0; a[2] =   -6.0; a[3] =  +6.0;
                break;
              default:
                affirm(FALSE, "bad basis element index");
            }
          break;
         
        case 2: 
          switch(i)
            { case 0:
                /* Quintic bell {(10 -15·z + 6·z^2)·z^3}. */
                a[0] =   0.0; a[1] =   0.0; a[2] =   0.0; 
                a[3] = +10.0; a[4] = -15.0; a[5] =  +6.0;
                break;
              case 1:
                /* Quintic swing {5·(z-1)·(1 + 3·z)·z^3}. */
                a[0] =   0.0; a[1] =   0.0; a[2] =   0.0; 
                a[3] = -20.0; a[4] = +35.0; a[5] = -15.0; 
                break;
              case 2:
                /* Quintic camel {30·(1-z)^2·z^3}. */
                a[0] =   0.0; a[1] =   0.0; a[2] =   0.0; 
                a[3] = +30.0; a[4] = -60.0; a[5] = +30.0; 
                break;
              default:
                affirm(FALSE, "bad basis element index");
            }
          break;
         
        default:
          affirm(FALSE, "bad time continuity order");
      }
  }

double SP1DSpline_Eval_Plus(int c, int i, double z)
  { affirm((i >= 0) && (i <= c), "invalid index");
    double u = 1.0 - z;
    if (c == 0)
      { if (i == 0)
          { return u; }
        else
          { affirm(FALSE, "duh?"); }
      }
    else if (c == 1)
      { if (i == 0)
          { return (2*z + 1)*u*u; }
        else if (i == 1)
          { return 6*z*u*u; }
        else 
          { affirm(FALSE, "duh?"); }
      }
    else if (c == 2)
      { if (i == 0)
          { return ((6*z + 3)*z + 1)*u*u*u; }
        else if (i == 1)
          { return 5*z*(3*z + 1)*u*u*u; }
        else if (i == 2)
          { return 30*(z*z*u*u*u); }
        else 
          { affirm(FALSE, "duh?"); }
      }
    else
      { affirm(FALSE, "not implemented yet"); }
    return 0.0;
  }

void SP1DSpline_Coeffs(int c, int i, int j, double *a)
  {
    affirm(c >= 0, "bad continuity");
    int d = 2*c + 1;
    if (j == 0) 
      { SP1DSpline_Coeffs_Plus(c, i, a); }
    else if (j == -1)
      { SP1DSpline_Coeffs_Minus(c, i, a); }
    else
      { int k;
        for(k = 0; k <= d; k++) { a[k] = 0.0; }
      }
  }

double SP1DSpline_Eval(int c, int i, double t)
  {
    if ((t <= -1.0) || (t >= +1.0))
      { return 0.0; }
    else
      { if (t >=  0.0) 
          { return SP1DSpline_Eval_Plus(c, i, t); }
        else
          { double s = (i % 2 == 0 ? 1.0 : -1.0);
            return s * SP1DSpline_Eval_Plus(c, i, -t);
          }
      }
  }

double SP1DSpline_Eval_0(int c, int i)
  {
    return (i == 0 ? 1.0 : 0.0);
  }
