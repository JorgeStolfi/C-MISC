/* See {gch_hyperplane.h} */
/* Last edited on 2014-07-22 03:37:49 by stolfilocal */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>

#include <gch_hyperplane.h>

gch_hyperplane_t *gch_hyperplane_new(int maxN)
  {
    gch_hyperplane_t *plane = (gch_hyperplane_t*) malloc(sizeof(gch_hyperplane_t));
    int i;
    plane->maxN = maxN;
    plane->coef = (mpz_t*)malloc(sizeof(mpz_t)*(maxN+1));
    for (i = 0; i < maxN+1; i++) { mpz_init(plane->coef[i]); }
    return plane;
  }

void gch_hyperplane_free(gch_hyperplane_t *plane)
  {
    int i;
    for (i = 0; i < plane->maxN+1; i++) { mpz_clear(plane->coef[i]); }
    free(plane->coef);
    free(plane);
  }

void gch_hyperplane_from_matrix(gch_hyperplane_t *plane, int N, gch_matrix_t *mat)
  {
    /* The matrix {mat} must have N+1 cols and N rows */
    int i;
    int *matRows = (int*)malloc(sizeof(int)*N);
    int *matCols = (int*)malloc(sizeof(int)*N);

    plane->N = N;

    for (i = 0; i < N; i++)
      { matRows[i] = i;
        matCols[i] = i+1;
      }

    for (i = 0; i <= N; i++)
      { if (i > 0) { matCols[i-1]--; }
        gch_matrix_sub_det(mat,plane->coef[i],N,matRows,matCols);
        if ((N+i)%2 == 1) { mpz_neg(plane->coef[i],plane->coef[i]); }
      }

    free(matRows);
    free(matCols);
  }

void gch_hyperplane_from_points(gch_hyperplane_t *plane, int N, gch_vector_t* base[])
  { /* Array {base} must have N vectors of N coordinates */
    int i,j;
    gch_matrix_t *mat = gch_matrix_new(N,N+1);
    for (i = 0; i < N; i++)
      { mpz_set_si(gch_matrix_elem(mat,i,0),1); }
    for (i = 0; i < N; i++)
      { for (j = 0; j < N; j++)
          { mpz_set(gch_matrix_elem(mat,i,j+1),base[i]->coord[j]); }
      }
    gch_hyperplane_from_matrix(plane,N,mat);
    gch_matrix_free(mat);
  }
  
void gch_hyperplane_point_test(gch_hyperplane_t *plane, mpz_t rop, gch_vector_t *p)
  { int i;
    mpz_set(rop,plane->coef[0]);
    for (i = 0; i < plane->N; i++)
      { mpz_addmul(rop, plane->coef[i+1], p->coord[i]); }
  }

int gch_hyperplane_point_test_sign(gch_hyperplane_t *plane, gch_vector_t *p)
  { int sign;
    mpz_t test;
    mpz_init(test);
    gch_hyperplane_point_test(plane,test,p);
    sign = mpz_sgn(test);
    mpz_clear(test);
    return sign;
  }

void gch_hyperplane_neg(gch_hyperplane_t *plane)
  { int i;
    for (i = 0; i <= plane->N; i++)
      { mpz_neg(plane->coef[i],plane->coef[i]); }
  }

void gch_hyperplane_print(FILE *wr, gch_hyperplane_t *plane)
  {
    int i;
    fprintf(wr, " < ");
    for (i = 0; i <= plane->N; i++)
      { if (i > 0) { fprintf(wr, ","); }
        mpz_out_str (wr, 10, plane->coef[i]);
      }
    fprintf(wr, " > \n");
  }

int  gch_hyperplane_is_null(gch_hyperplane_t *plane)
  { int i;
    for (i = 0; i <= plane->N; i++)
      { if (mpz_cmp_si(plane->coef[i],0) != 0) { return 0; } }
    return 1; 
  }

