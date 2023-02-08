/* Last edited on 2014-07-22 03:10:30 by stolfilocal */
#ifndef gch_matrix_H
#define gch_matrix_H

#define _GNU_SOURCE
#include <stdio.h>

#include <gmp.h>
#include <gch_vector.h>

typedef struct gch_matrix_t
  { int rows, cols;
    mpz_t *data;
  } gch_matrix_t;

#define gch_matrix_elem(mat,row,col) (mat->data[(row)*(mat->cols)+(col)])

gch_matrix_t* gch_matrix_new(int rows, int cols);

void gch_matrix_free(gch_matrix_t *mat);

void gch_matrix_set_row(gch_matrix_t *mat, int row, gch_vector_t *vec); 
  /* The vector must have dimension {mat.cols-1} ??? Check ??? */
  
void gch_matrix_sub_det(gch_matrix_t *mat, mpz_t rop, int level, int *rowList, int *colList);

void gch_matrix_det(gch_matrix_t *mat, mpz_t rop);

gch_vector_t* gch_matrix_map_col(gch_matrix_t *mat, gch_vector_t *vec); 
  /* Inner product. */

void gch_matrix_print(FILE *wr, gch_matrix_t *mat);

#endif
