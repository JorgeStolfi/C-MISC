/* Last edited on 2014-07-22 03:36:59 by stolfilocal */
#ifndef gch_hyperplane_H
#define gch_hyperplane_H

#define _GNU_SOURCE
#include <stdio.h>

#include <gmp.h>

#include <gch_matrix.h>
#include <gch_vector.h>

typedef struct gch_hyperplane_t
  { int maxN;
    int N;
    mpz_t *coef;
  } gch_hyperplane_t;

gch_hyperplane_t *gch_hyperplane_new(int maxN);

void gch_hyperplane_free(gch_hyperplane_t *plane);

void gch_hyperplane_from_matrix(gch_hyperplane_t *plane, int N, gch_matrix_t *mat); 
  /* The gch_matrix mat must have N+1 cols and N rows */
  
void gch_hyperplane_from_points(gch_hyperplane_t *plane, int N, gch_vector_t* base[]); 
  /* Array base must have N vectors of N coordinates */
  
void gch_hyperplane_point_test(gch_hyperplane_t *plane, mpz_t rop, gch_vector_t *p);

int gch_hyperplane_point_test_sign(gch_hyperplane_t *plane, gch_vector_t *p);

void gch_hyperplane_neg(gch_hyperplane_t *plane);

void gch_hyperplane_print(FILE *wr, gch_hyperplane_t *plane);

int  gch_hyperplane_is_null(gch_hyperplane_t *plane);

#endif
