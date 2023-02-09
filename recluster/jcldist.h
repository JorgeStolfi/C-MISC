/* Last edited on 2007-08-15 22:36:06 by stolfi */

#ifndef jcldist_h
#define jcldist_h

#include <jclbasic.h>

/* DISTANCE MATRIX */

void WriteMatrix(FILE *f, double **d, long *it, long N, byte **L, long skip);
  /*
    Prints the distance matrix "d[it[0..N-1]][it[0..N-1]]" to file "f".
    Uses "L[it[0..N-1]][0..skip-1]" as the row and column labels. */
    
int FormatDist(double d);
  /* 
    Converts a distance value from the range "[0 _ 1]" to 
    an integer in the range "[0..99]". */
    
#endif
