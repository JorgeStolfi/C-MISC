/* See {gch_vector.h} */
/* Last edited on 2014-07-22 02:33:55 by stolfilocal */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>

#include <gch_vector.h>

gch_vector_t *gch_vector_new(int N)
  { gch_vector_t *vec = (gch_vector_t*)malloc(sizeof(gch_vector_t));
    int i;
    vec->N = N;
    vec->coord = (mpz_t*)malloc(sizeof(mpz_t)*N);
    for (i = 0; i < N; i++) { mpz_init(vec->coord[i]); }
    return vec;
  }

void gch_vector_free(gch_vector_t *vec)
  { int i;
    for (i = 0; i < vec->N; i++) { mpz_clear(vec->coord[i]); }
    free(vec->coord);
    free(vec);
  }

gch_vector_t* gch_vector_copy(gch_vector_t *vec)
  { gch_vector_t *copy = gch_vector_new(vec->N);
    int i;
    for (i = 0; i < vec->N; i++) { mpz_set(copy->coord[i],vec->coord[i]); }
    return copy;
  }

void gch_vector_print(FILE *wr, gch_vector_t *vec)
  { int i;
    fprintf(wr, "( ");
    for (i = 0; i < vec->N; i++)
      { mpz_out_str (wr, 10, vec->coord[i]);
        fprintf(wr, " ");
      }
    fprintf(wr, ")\n");
  }
