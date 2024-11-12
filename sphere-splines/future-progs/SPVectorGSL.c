/* See SPVectorGSL.h */
/* Last edited on 2007-11-10 11:01:38 by anamaria */

#include <SPVectorGSL.h>
#include <SPVector.h>

#include <SPBasic.h>

#include <vec.h>
#include <filefmt.h>
#include <fget.h>
#include <nget.h>
#include <affirm.h>
#include <nat.h>
#include <bool.h>

#include <gsl/gsl_vector_double.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>

gsl_vector *SPVector_To_GSLVector(SPVector *v)
  { int nv = v->ne;
    gsl_vector *Gv = gsl_vector_alloc(nv);
    int i;
    for (i = 0; i < nv; i++) { gsl_vector_set(Gv, i, v->e[i]); }
    return Gv;
  }

SPVector SPVector_From_GSLVector(gsl_vector *Gv)
  { int nv = Gv->size;
    SPVector v = double_vec_new(nv);
    int i;
    for (i = 0; i < nv; i++) { v.e[i] = gsl_vector_get(Gv, i);  }
    return v;
  }


