/* Last edited on 2014-07-22 02:34:14 by stolfilocal */

#ifndef gch_vector_H
#define gch_vector_H

#define _GNU_SOURCE
#include <stdio.h>
#include <gmp.h>

typedef struct gch_vector_t 
  { int N;
    mpz_t *coord;
  } gch_vector_t;

gch_vector_t *gch_vector_new(int N);

void gch_vector_free(gch_vector_t *vec);

gch_vector_t* gch_vector_copy(gch_vector_t *vec);

void gch_vector_print(FILE *wr, gch_vector_t *vec);

#endif
