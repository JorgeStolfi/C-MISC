/* SPMatrixGSL.h -- interface to the GSL sparse matrix lib */
/* Last edited on 2011-07-30 14:17:48 by stolfilocal */

#ifndef SPMatrixGSL_H
#define SPMatrixGSL_H

#include <SPMatrix.h>
#include <vec.h>
#include <bool.h>
#include <nat.h>

#include <gsl/gsl_matrix_double.h>
#include <stdio.h>

gsl_matrix *SPMatrix_To_GSLMatrix(SPMatrix *A);
  /* Unpacks matrix {A} into a GSL matrix. */

SPMatrix SPMatrix_From_GSLMatrix(gsl_matrix *GA);
  /* Packs the GSL matrix {GA} into an {SPMatrix}. */

#endif
 
