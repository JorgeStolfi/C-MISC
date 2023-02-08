/* See {gch_geometry.h}  */
/* Last edited on 2014-07-22 03:41:10 by stolfilocal */

#include <gmp.h>
#include <gch_geometry.h>
#include <gch_matrix.h>

int gch_geometry_affine_hull_basis_get(int dim, gch_vector_list_sub_t *S, int *basisIndex)
  {
    int affineHullDim;

    /* Find first facet's hyperplane H: */
    gch_hyperplane_t *H = gch_geometry_find_first_facet(dim, S, basisIndex);

    if (H->N < dim)
      { affineHullDim = H->N-1; }
    else 
      { /* Find a point q of S not contained in H: */
        int q = 0;
        while((q < S->count) && (gch_hyperplane_point_test_sign(H,S->vectors[q]) == 0)) { q++; }
        if (q == S->count)
          { /* If there is NO q in S that is not contained in H: */
            affineHullDim = dim-1;
          }
        else 
          { /* If there IS a q in S that is not contained in H: */
            basisIndex[dim] = q;
            affineHullDim = dim;
          }
      }
    gch_hyperplane_free(H);
    return affineHullDim;
  }

int gch_geometry_find_adj_facet(gch_vector_list_sub_t *S, gch_hyperplane_t *G, gch_hyperplane_t *K)
  {
    int p = -1; /* Index of the chosen point */
    int q;
    mpz_t G_p, K_p, G_q, K_q;
    mpz_t P1, P2;

    mpz_init(G_p);
    mpz_init(G_q);
    mpz_init(K_p);
    mpz_init(K_q);
    mpz_init(P1);
    mpz_init(P2);

    for (q = 0; q < S->count; q++)
      { gch_hyperplane_point_test(G,G_q,S->vectors[q]);
        if (mpz_cmp_si(G_q, 0) != 0)
          { gch_hyperplane_point_test(K, K_q, S->vectors[q]);
            if (p == -1)
              { p = q;
                mpz_set(G_p,G_q);
                mpz_set(K_p,K_q);
              }
            else 
              { mpz_mul(P1,K_q,G_p);
                mpz_mul(P2,K_p,G_q);
                if (mpz_cmp(P1,P2) < 0)
                  { p = q;
                    mpz_set(G_p,G_q);
                    mpz_set(K_p,K_q);
                  }
              }
          }
      }

    mpz_clear(G_p);
    mpz_clear(G_q);
    mpz_clear(K_p);
    mpz_clear(K_q);
    mpz_clear(P1);
    mpz_clear(P2);

    return p;
  }

gch_hyperplane_t* gch_geometry_find_first_facet(int dim, gch_vector_list_sub_t *S, int *basisIndex)
  {
    int i;
    int test;
    int subspaceDim;
    gch_matrix_t *mat = gch_matrix_new(dim, dim+1);
    /* Get the point {p0} of {S} with least {X[0]} coordinate: */
    int p0 = 0;
    for (i = 1; i < S->count; i++)
      { if (mpz_cmp(S->vectors[i]->coord[0],S->vectors[p0]->coord[0]) < 0) { p0 = i; } }

    /* Store {p0} in the array {basisIndex}:  */
    basisIndex[0] = p0;

    /* Copy point {p0} to row 0 of the matrix: */
    gch_matrix_set_row(mat,0,S->vectors[p0]);

    gch_hyperplane_t *H = gch_hyperplane_new(dim);
    gch_hyperplane_t *Haux = gch_hyperplane_new(dim);
    gch_hyperplane_t *K = gch_hyperplane_new(dim);

    /* We have a 0-hyperplane of a 1-dimensional space, defined by {points[0]}: */
    subspaceDim = 1;
    gch_hyperplane_from_matrix(H,subspaceDim,mat);
    while(subspaceDim < dim)
      { /* Increase the dimension of hyperplane H by adding a new coefficient of value = 0: */
        subspaceDim++;
        H->N = subspaceDim;

        /* Check if H is already a supporting hyperplane containing a (subspaceDim-1)-face: */
        int j;
        for (j = 0; j < S->count; j++)
          { if (gch_hyperplane_point_test_sign(H,S->vectors[j]) == 0)
              { /* Copy point j to row subspaceDim-1 of the matrix */
                gch_matrix_set_row(mat,subspaceDim-1,S->vectors[j]);
                gch_hyperplane_from_matrix(Haux,subspaceDim,mat);
                if (!gch_hyperplane_is_null(Haux)) { break; }
              }
          }

        if (j < S->count)
          { basisIndex[subspaceDim-1] = j; }
        else
          { /* Find a point q of S not contained in H: */
            int q = 0;
            while((q < S->count)&&((test = gch_hyperplane_point_test_sign(H,S->vectors[q])) == 0)) { q++; }

            if (q == S->count)
              { H->N = subspaceDim-1; /* S is contained in a subspace of dimension < dim */
                break; /* Leave the while(subspaceDim < dim) */
              }
            else
              { /* Negate H if necessary: */
                if (test < 0) { gch_hyperplane_neg(H); }
                /* Copy point q to row subspaceDim-1 of the matrix: */
                gch_matrix_set_row(mat,subspaceDim-1,S->vectors[q]);
                /* Find the hyperplane K: */
                gch_hyperplane_from_matrix(K,subspaceDim,mat);
                /* Get a point p (not contained in H) from the new supporting hyperplane */
                int p = gch_geometry_find_adj_facet(S, H, K);
                /* Copy point p to row subspaceDim-1 of the gch_matrix: */
                gch_matrix_set_row(mat,subspaceDim-1,S->vectors[p]);
                /* Get the new supporting gch_hyperplane: */
                gch_hyperplane_from_matrix(H,subspaceDim,mat);
                /* Store p in the array basisIndex: */
                basisIndex[subspaceDim-1] = p;
              }
          }
      }

    gch_matrix_free(mat);
    gch_hyperplane_free(K);
    gch_hyperplane_free(Haux);

    return H;
  }
