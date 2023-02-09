/* Last edited on 1997-11-24 00:51:47 by stolfi */
/* See the authorship and copyright notice at the end of this file. */

/* IMAGE OUTPUT */

#ifndef jclimage_h
#define jclimage_h

#include <stdio.h>
#include <jclbasic.h>

typedef struct { byte R, G, B; } pixel;

void SelectColors(pixel *map, long NC);
  /*
    Stores in "map[0..NC-1]" a palette of colors suitable for displaying
    values in the range "[0_1]".  Entry "map[0]" is black, "map[NC-1}" is
    white, and intermediate entries have intermediate intensities
    (with varying hues).
  */
  
void WritePseudoColorImage(FILE* f, byte *g, long M, long N, pixel *map, long NC);
  /*
    Writes a pseudocolor image with "M" rows and "N" columns.
    The pixels are "map[g[0..M*N-1]]", in row order,
    where "map[0..NC-1]" is the color palette.
  */
    
byte *MatrixToImage(
    double **d, long *row, long M, long *col, long N, 
    double dMin, double dMax, long NC
  );
  /*
    Converts the distances "d[row[0..M-1]][col[0..N-1]]" to integers
    in the range "[0..NC-1]", where "NC" is 255 or less.  Maps a value
    0.0 or less to "0", "dMin" to 1, "dMax" or more to NC-1,
    and intermediate values to intermediate values. Returns a pointer 
    to a newly allocated array of "M*N" bytes.
  */
    
#endif
