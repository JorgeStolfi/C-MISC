/* SPVector.h -- sparse matrices */
/* Last edited on 2009-01-09 22:14:29 by stolfi */

#ifndef SPVector_H
#define SPVector_H

#include <vec.h>

#include <stdio.h>

typedef double_vec_t SPVector;

void SPVector_Write(FILE *wr, SPVector v);
  /* Write  the  vector {v} to {wr}, in a format that can be read back. */
  
SPVector SPVector_Read(FILE *rd);
  /* Reads a vector from {rd}, assuming it was written by {SPVector_Write}. */

#endif
 
