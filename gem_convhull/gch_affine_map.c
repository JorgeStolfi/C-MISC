/* See {gch_affine_map.h} */
/* Last edited on 2014-07-22 02:16:35 by stolfilocal */
#include <stdlib.h>
#include <gmp.h>

#include <gch_affine_map.h>

gch_affine_map_t* gch_affine_map_new(int M, int N, gch_vector_t *basis[])
  { gch_affine_map_t *bt = (gch_affine_map_t*)malloc(sizeof(gch_affine_map_t));
    int i,j;
    bt->M = M;
    bt->N = N;
    bt->v0 = gch_vector_copy(basis[0]);
    bt->mat = gch_matrix_new(N,M);
    for (i = 1; i < N+1; i++)
      { for (j = 0; j < M; j++)
          { mpz_sub(gch_matrix_elem(bt->mat,i-1,j),basis[i]->coord[j],bt->v0->coord[j]); }
      }
    return bt;
  }

void gch_affine_map_free(gch_affine_map_t *bt)
  { gch_matrix_free(bt->mat);
    gch_vector_free(bt->v0);
    free(bt);
  }

gch_vector_t* gch_affine_map_map(gch_affine_map_t *bt, gch_vector_t *v)
  { int i;
    gch_vector_t *u = gch_vector_new(v->N);
    for (i = 0; i < bt->M; i++)
      { mpz_sub(u->coord[i],v->coord[i],bt->v0->coord[i]); }
    gch_vector_t *w = gch_matrix_map_col(bt->mat, u);
    gch_vector_free(u);
    return w;
  }
