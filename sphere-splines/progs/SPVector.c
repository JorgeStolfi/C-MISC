/* See SPVector.h */
/* Last edited on 2023-02-12 07:50:04 by stolfi */

#include <SPVector.h>

#include <SPBasic.h>

#include <vec.h>
#include <filefmt.h>
#include <fget.h>
#include <nget.h>
#include <affirm.h>
#include <nat.h>
#include <bool.h>

// #include <gsl/gsl_vector_double.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>

#define SPVector_FileFormat "2005-10-23"
    
void SPVector_Write(FILE *wr, SPVector v)
  { int nv = v.ne;
    int i;
    /* Count non-zero elements: */
    int nw = 0;
    for (i = 0; i < nv; i++)
      { if (v.e[i] != 0.0) { nw++; } }
    /* Write the nonzero elements: */
    filefmt_write_header(wr, "SPVector", SPVector_FileFormat);
    fprintf(wr, "dim = %d\n", nv);
    fprintf(wr, "elems = %d\n", nw);
    for (i = 0; i < nv; i++)
      { if (v.e[i] != 0.0) 
          { fprintf(wr, "%d %24.16e\n", i, v.e[i]); }
      }
    filefmt_write_footer(wr, "SPVector");
    fflush(wr);
  }

SPVector SPVector_Read(FILE *rd)
  { SPVector v;
    filefmt_read_header(rd, "SPVector", SPVector_FileFormat);
    int nv = nget_int32(rd, "dim"); fget_eol(rd); 
    int nr = nget_int32(rd, "elems"); fget_eol(rd); 
    { v = double_vec_new(nv);
      int iAnt = -1;
      int j;
      for (j = 0; j < nr; j++)
        { int i = fget_int32(rd);
          affirm((i > iAnt) && (i < nv), "bad index");
          v.e[i] = fget_double(rd);
          fget_eol(rd);
        }
    }  
    filefmt_read_footer(rd, "SPVector");
    return v;
  }
