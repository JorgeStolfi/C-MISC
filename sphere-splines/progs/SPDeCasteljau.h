/* SPDeCasteljau.h -- The DeCasteljau algorithm for Bezier evaluation */
/* Last edited on 2005-10-27 07:57:31 by stolfi */

#ifndef SPDeCasteljau_H
#define SPDeCasteljau_H

#include <r3.h>
#include <r6.h>
#include <vec.h>
#include <js.h>
#include <nat.h>
   
typedef double_vec_t  BezCoeff_vec_t;
  /* Bezier Control values. Their order is defined by 
    {SPDeCasteljau_BezLabelToIndex} (below); for degree 6 it is
       C600, C510, C501, C420, C411, C402, ... , C006. */

#define BezCoeff_vec_new double_vec_new

void SPDeCasteljau_EvalGen(BezCoeff_vec_t coeff, int deg, r3_t *a, double *f);
  /* Computes the value {f} of a spline with degree {deg} and 
    coefficients {coeff}, at  a point with barycentric coords {a[0..2]} */

void SPDeCasteljau_GradGen(BezCoeff_vec_t coeff, int deg, r3_t *a, double *f, r3_t *df);
  /* Computes the value {f} and the gradient {df} (in R^3) of a spline
    with degree {deg} and coefficients {coeff}, at a point with
    barycentric coords {a[0..2]} */

void SPDeCasteljau_HessGen(BezCoeff_vec_t coeff, int deg, r3_t *a, double *f, r3_t *df, r6_t *ddf);
  /* Computes the value {f}, the gradient {df}, and the Hessian {ddf}
    (in R^3) of a spline with degree {deg} and coefficients {coeff},
    at a point with barycentric coords {a[0..2]} */

/* Hand-optimized versions: */

void SPDeCasteljau_Eval0(BezCoeff_vec_t coeff, r3_t *a, double *f);
void SPDeCasteljau_Eval1(BezCoeff_vec_t coeff, r3_t *a, double *f);
void SPDeCasteljau_Eval2(BezCoeff_vec_t coeff, r3_t *a, double *f);
void SPDeCasteljau_Eval3(BezCoeff_vec_t coeff, r3_t *a, double *f);
void SPDeCasteljau_Eval4(BezCoeff_vec_t coeff, r3_t *a, double *f); 
void SPDeCasteljau_Eval5(BezCoeff_vec_t coeff, r3_t *a, double *f);
void SPDeCasteljau_Eval6(BezCoeff_vec_t coeff, r3_t *a, double *f);
void SPDeCasteljau_Eval7(BezCoeff_vec_t coeff, r3_t *a, double *f);

void SPDeCasteljau_Grad0(BezCoeff_vec_t coeff, r3_t *a, double *f, r3_t *df);
void SPDeCasteljau_Grad1(BezCoeff_vec_t coeff, r3_t *a, double *f, r3_t *df);
void SPDeCasteljau_Grad2(BezCoeff_vec_t coeff, r3_t *a, double *f, r3_t *df);
void SPDeCasteljau_Grad3(BezCoeff_vec_t coeff, r3_t *a, double *f, r3_t *df);
void SPDeCasteljau_Grad4(BezCoeff_vec_t coeff, r3_t *a, double *f, r3_t *df); 
void SPDeCasteljau_Grad5(BezCoeff_vec_t coeff, r3_t *a, double *f, r3_t *df);
void SPDeCasteljau_Grad6(BezCoeff_vec_t coeff, r3_t *a, double *f, r3_t *df);
void SPDeCasteljau_Grad7(BezCoeff_vec_t coeff, r3_t *a, double *f, r3_t *df);

/* INDEXING BEZIER COEFFICIENTS */

typedef struct BezLabel { unsigned char e[3]; } BezLabel;
  /* Standard label for a Bezier Control value or coefficient 
    of a homogeneous polynomial, eg. (6,0,0) or (3,2,1).  

    The term order for degree {d} is given by the algorithm

      poly = 0; 
      index = 0;
      for (i = d; i >=0; i--)
        { for (j = d - i; j >= 0; j--)
            { k = d - i - j;
              poly = poly + c[index] * B_d{i,j,k}
              index = index + 1;
            }
        }
  */

BezLabel SPDeCasteljau_IndexToBezLabel(nat_t i, nat_t deg);

nat_t SPDeCasteljau_BezLabelToIndex(BezLabel e, nat_t deg);
  /* Maps standard labels to position in the ControlValue array; e.q. 
     (6,0,0) <--> 0, (5,1,0) <--> 1, ... , (0,0,6) <--> 27 */

#endif
