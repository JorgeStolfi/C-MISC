/* See {gch_matrix.h} */
/* Last edited on 2014-07-22 03:10:47 by stolfilocal */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>

#include <gmp.h>

#include <gch_vector.h>

#include <gch_matrix.h>

gch_matrix_t* gch_matrix_new(int rows, int cols)
  { gch_matrix_t *mat = (gch_matrix_t*) malloc(sizeof(gch_matrix_t));
    int i;
    mat->rows = rows;
    mat->cols = cols;
    mat->data = (mpz_t*) malloc(sizeof(mpz_t)*rows*cols);
    for (i = 0; i < rows*cols; i++)
      { mpz_init(mat->data[i]); }
    return mat;
  }

void gch_matrix_free(gch_matrix_t *mat)
  { int i;
    for (i = 0; i < (mat->rows)*(mat->cols); i++)
      { mpz_clear(mat->data[i]); }
    free(mat->data);
    free(mat);
  }

void gch_matrix_set_row(gch_matrix_t *mat, int row, gch_vector_t *vec)
  { // a gch_vector of dimension cols-1
    int i;
    mpz_set_si(gch_matrix_elem(mat,row,0),1);
    for (i = 1; i < mat->cols; i++)
      { mpz_set(gch_matrix_elem(mat,row,i), vec->coord[i-1]); }
  }

void gch_matrix_sub_det(gch_matrix_t *mat, mpz_t rop, int level, int *rowList, int *colList)
  {
    if (level == 1)
      { mpz_set(rop, gch_matrix_elem(mat,rowList[0],colList[0])); }
    else
      { int i,j,k;
        int *nextRowList = (int*)malloc(sizeof(int)*(level-1));
        int *nextColList = (int*)malloc(sizeof(int)*(level-1));
        mpz_t subVal;

        mpz_init(subVal);

        mpz_set_si(rop,0);

        /* Build row lists: */
        for (j = 0; j < level-1; j++) { nextRowList[j] = rowList[j+1]; }

        for (i = 0; i < level; i++)
          { /* Build column lists: */
            k = 0;
            for (j = 0; j < level; j++)
              { if (j != i)
                  {  nextColList[k++] = colList[j]; }
              }
            // computando parcela da somatória
            gch_matrix_sub_det(mat,subVal, level-1,nextRowList,nextColList);
            if (i%2 == 0)
              { mpz_addmul(rop,gch_matrix_elem(mat,rowList[0],colList[i]),subVal); }
            else
              { mpz_submul(rop,gch_matrix_elem(mat,rowList[0],colList[i]),subVal); }
          }
        mpz_clear(subVal);
        free(nextRowList);
        free(nextColList);
      }
  }

void gch_matrix_det(gch_matrix_t *mat, mpz_t rop)
  { int i;
    int *completeList = (int*)malloc(sizeof(int)*(mat->rows));
    for (i = 0; i < mat->rows; i++)  { completeList[i] = i; }
    gch_matrix_sub_det(mat, rop, mat->rows, completeList, completeList);
    free(completeList);
  }

gch_vector_t* gch_matrix_map_col(gch_matrix_t *mat, gch_vector_t *vec)
  { gch_vector_t *product = gch_vector_new(mat->rows);
    int i,j;
    for (i = 0; i < mat->rows; i++)
      { for (j = 0; j < mat->cols; j++)
          { mpz_addmul(product->coord[i],gch_matrix_elem(mat,i,j),vec->coord[j]); }
      }
    return product;
  }

void gch_matrix_print(FILE * wr, gch_matrix_t *mat)
  { int i,j;
    for (i = 0; i < mat->rows; i++)
      { for (j = 0; j < mat->cols; j++)
          { mpz_out_str (wr, 10, gch_matrix_elem(mat,i,j));
            fprintf(wr, "   ");
          }
        fprintf(wr, "\n");
      }
  }
