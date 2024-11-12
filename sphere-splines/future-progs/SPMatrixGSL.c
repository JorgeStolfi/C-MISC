/* See SPMatrixGSL.h */
/* Last edited on 2009-01-09 22:12:56 by stolfi */

#include <SPMatrixGSL.h>
#include <SPVectorGSL.h>

#include <SPMatrix.h>
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
#include <gsl/gsl_matrix_double.h>
#include <gsl/gsl_linalg.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#include <limits.h>
  
gsl_matrix *SPMatrix_To_GSLMatrix(SPMatrix *A)
  {
    gsl_matrix *GA = gsl_matrix_calloc(A->rows, A->cols);
    int nA = A->ents.ne;
    int kA;
    for (kA = 0; kA < nA; kA++)
      { MatEntry *Ak = &(A->ents.e[kA]);
        gsl_matrix_set(GA, Ak->row, Ak->col, Ak->va);
      }
    return GA;
  }

SPMatrix SPMatrix_From_GSLMatrix(gsl_matrix *GA)
  {
    MatEntry_vec_t Aents = MatEntry_vec_new(GA->size1);
    int nA = 0;
    int i;
    for (i = 0; i < GA->size1; i++)
      { int j;
        for (j = 0; j < GA->size2; j++)
          { double GAij = gsl_matrix_get(GA, i, j);
            if (GAij != 0.0)
              { MatEntry_vec_expand(&Aents, nA);
                MatEntry *Aij = &(Aents.e[nA]);
                Aij->row = i; Aij->col = j;
                Aij->va = GAij;
                nA++;
              }
          }
      }
    MatEntry_vec_trim(&Aents, nA);
    return (SPMatrix){GA->size1, GA->size2, Aents};
  }
