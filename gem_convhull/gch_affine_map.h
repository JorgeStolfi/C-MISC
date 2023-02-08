/* Last edited on 2014-07-22 03:38:05 by stolfilocal */
#ifndef gch_affine_map_H
#define gch_affine_map_H

#define _GNU_SOURCE

#include <gch_vector.h>
#include <gch_matrix.h>

typedef struct gch_affine_map_t
  { int M, N;
    gch_vector_t *v0;
    gch_matrix_t *mat; /* N vectors of dimension M */
  } gch_affine_map_t;
  /* An affine map that takes the column {M}-vector {v} 
    to the column {N}-vector {u = mat*(v-v0)}. */

gch_affine_map_t* gch_affine_map_new(int M, int N, gch_vector_t *basis[]);
  /* Alocates a new affine map. */

void gch_affine_map_free(gch_affine_map_t *bt);
  /* Reclaims the space of affine map {bt}. */

gch_vector_t *gch_affine_map_map(gch_affine_map_t *bt, gch_vector_t *v);
  /* Applies the affine transformation {bt} to the vector {v}. */

#endif
