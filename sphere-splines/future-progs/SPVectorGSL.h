/* SPVectorGSL.h -- interface to the GSL sparse vector lib */
/* Last edited on 2009-02-10 11:19:44 by stolfi */

#ifndef SPVectorGSL_H
#define SPVectorGSL_H

#include <vec.h>

#include <gsl/gsl_vector_double.h>
#include <SPVector.h>

#include <stdio.h>

gsl_vector *SPVector_To_GSLVector(SPVector *v);
  /* Unpacks vector {v} into a GSL vector. */

SPVector SPVector_From_GSLVector(gsl_vector *Gv);
  /* Packs the GSL vector {Gv} into an {SPVector}. */

#endif
 
