/* SPOverlapTable.h -- Triangle overlap among two triangulations */
/* Last edited on 2002-12-09 10:08:35 by anamaria */

#ifndef SPOverlapTable_H
#define SPOverlapTable_H

#include <basic.h>
#include <arrays.h>
#include <stdio.h>

typedef struct natArrayArray { nat nel; natArray *el; } natArrayArray;
  /* A two-dimensional irregular array of unsigned ingegers. */

typedef natArrayArray OverlapTable; 
  /* An {OverlapTable} {tb} lists the pairs of overlapping
    triangles between two triangulations {S} and {T}.
    Specifically, the triangle {i} of {S} overlaps triangles
    {tb.e[i].e[0..ni-1]} of {T}, where {ni = tb.e[i].ne}. */
      
void SPOverlapTable_Write(FILE *wr, OverlapTable tb);
  /* Writes {tb} to {wr}, in a format that can be read back. */

OverlapTable SPOverlapTable_Read(FILE *rd);
  /* Reads from {rd} a table that was written by {Write} above. */

/* Arrays of arrays of unsigned ints ({natArray}s) */

natArray natArrayArray_New(nat nel);
void natArrayArray_Expand(natArrayArray *arrp, nat i);
void natArrayArray_Trim(natArrayArray *arrp, nat nel);

#endif
