/* See {gch_convexhull.h} */
/* Last edited on 2015-12-01 16:25:01 by stolfilocal */

#define _GNU_SOURCE
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <gmp.h>

#include <affirm.h>

#include <gem.h>
#include <gem_bary.h>

#include <gch_util.h>
#include <gch_vector.h>
#include <gch_vector_list.h>
#include <gch_list.h>
#include <gch_tree.h>
#include <gch_hyperplane.h>
#include <gch_geometry.h>
#include <gch_computed_solutions.h>

#include <gch_convexhull.h>
    
int gch_node_compare(gem_ref_t a, gem_ref_t b, int d);
 /* Compares the nodes {a} and {b} by lexicographic order
   of their vertex numbers ???.  Result is in {-1,0,+1}. */

void gch_get_vertex_tuple(gem_ref_t p, int d, int vdata[]);
  /* ??? */


gem_ref_t gch_convexhull_trivial(gch_vector_list_sub_t *S, int odim, gch_list_t* facets);
  /* Builds the convex hull of a set {S} of points of {R^1} whose affine span is {R^1};
    that is, a list of two or more numbers, not all equal.
    Returns the result as a barycentric gem with dimension {odim}, consisting
    of two nodes (the min and max of those numbers) sharing the wall with wall-color 0.
    All the other walls will be unattached. 
    
    The data fields of each gem nodes are set (with {gem_set_data}) to the original 
    indices {S->realindex[i]} of the extremal points.
    
    Returns the cell that corresponds to the point with smallest orginal index.  */

gem_ref_t gch_convexhull_find_rec
  ( int dim, 
    gch_vector_list_sub_t *S, 
    int odim, 
    gch_list_t* facets, 
    gch_computed_solutions_t *cs
  );
  /* Computes the convex hull of the point set {S}, as a barycentric gem.
    Assumes that each point is given by {dim} Cartesian coordinates, 
    and the affine span of {S} is {R^dim}. 
    The gem will have dimension {odim}; requires {dim<=odim}.
    The gem cell walls with wall-colors {dim..odim} will be border walls.
    
    The facets of the hull will be appended to the {facets} list.
    
    ??? What about {cs}? */

gem_ref_t gch_convexhull_trivial(gch_vector_list_sub_t *S, int odim, gch_list_t* facets)
  {
    assert(S->count >= 2);
    int verbose = 1;
    int ind = 6 + 2*(odim - 1); /* Indentation for debug messages. */
    int rMin, rMax; /* Original indices of extremal points, sorted by index. */
    
    { 
      int iMin = 0;
      int iMax = 0;
      gch_vector_t *vMin = S->vectors[0];
      gch_vector_t *vMax = S->vectors[0];
      int i;
      for (i = 1; i < S->count; i++)
        { gch_vector_t *vi = S->vectors[i];
          if (mpz_cmp(vi->coord[0], vMin->coord[0]) < 0) { iMin = i; vMin = vi; }
          if (mpz_cmp(vi->coord[0], vMax->coord[0]) > 0) { iMax = i; vMax = vi; }
        }
      /* Get original indices of points: */
      rMin = S->realIndex[iMin];
      rMax = S->realIndex[iMax];
      /* Sort the indices in incerasing order: */
      assert(rMin != rMax);
      if (rMin > rMax) { int temp = rMin; rMin = rMax; rMax = temp; }
    }
    
    /* !!! Check whether dimension is 1 or 0 !!!  */
    if (verbose) { fprintf(stderr, "%*s  trivial hull corners: %d %d\n", ind, "", rMin, rMax); }
    gem_ref_t nMin = gem_node_new(odim); gem_set_data(nMin, rMin);
    gem_ref_t nMax = gem_node_new(odim); gem_set_data(nMax, rMax);
    gem_splice(nMin, nMax, 0);
    if (facets != NULL)
      { gch_list_add(facets,(void*) nMin);
        gch_list_add(facets,(void*) nMax);
      }
    return nMin;
  }

gem_ref_t gch_convexhull_find_rec
  ( int dim, 
    gch_vector_list_sub_t *S, 
    int odim, 
    gch_list_t* facets, 
    gch_computed_solutions_t *cs
  )
  {
    assert(cs->dim == odim);
    int verbose = 1;
    int ind = 6 + 2*(odim - dim); /* Indentation for debug messages. */

    if (verbose) { fprintf(stderr, "%*sENTER gch_convexhull_find_rec dim = %d\n", ind, "", dim); }
    gch_problem_input_t* pi = NULL;
    if (dim <= odim-2)
      { if (verbose) { fprintf(stderr, "%*s  ???\n", ind, ""); }
        pi = gch_problem_input_new(S);
        gem_ref_t gem = gch_solution_find(cs, dim, pi, facets);
        if (gem != NULL)
          { gch_problem_input_free(pi);
            return gem;
          }
      }
    
    gem_ref_t convexHull;  /* The gem node to return. */
    if (dim == 1)
      { if (verbose) { fprintf(stderr, "%*s  computing trivial hull\n", ind, ""); }
        convexHull = gch_convexhull_trivial(S, odim, facets);
      }
    else
      { if (verbose) { fprintf(stderr, "%*s  computing non-trivial hull\n", ind, ""); }
        int i;

        gch_vector_t** Basis = (gch_vector_t**)malloc(sizeof(gch_vector_t*)*dim);
        
        /* Get the hyperplane {H} of some facet: */
        if (verbose) { fprintf(stderr, "%*s  obtaining the first facet:", ind, ""); }
        int *basisIndex = (int*)malloc(sizeof(int)*dim);
        gch_hyperplane_t *H = gch_geometry_find_first_facet(dim, S, basisIndex);

        /* Store into array {points} the points that lie on {H}: */
        int *points = (int*)malloc(sizeof(int)*(S->count));
        int numPoints = 0;
        for (i = 0; i < S->count; i++)
          { if (gch_hyperplane_point_test_sign(H, S->vectors[i]) == 0)
              { if (verbose) { fprintf(stderr, " %d", S->realIndex[i]); }
                points[numPoints] = i;
                numPoints++;
              }
          }
        if (verbose) { fprintf(stderr, "\n"); }
        
        /* Project the points that lie on {H} onto {R^(dim-1)}: */
        if (verbose) { fprintf(stderr, "%*s  projecting the points of the first facet\n", ind, ""); }
        for (i = 0; i < dim; i++) { Basis[i] = S->vectors[basisIndex[i]]; }
        gch_affine_map_t *bt = gch_affine_map_new(dim, dim-1, Basis);
        gch_vector_list_sub_t *subS = gch_vector_list_sub_proj_new(S, points, numPoints, bt);
        
        /* Build recursively the gem of the first facet: */
        if (verbose) { fprintf(stderr, "%*s  building the gem of the first facet\n", ind, ""); }
        gch_list_t *Q = gch_list_new();  /* List of {dim-2}-faces that comprise the horizon. */
        gem_ref_t F = gch_convexhull_find_rec(dim-1, subS, odim, Q, cs);
        gch_vector_list_sub_free(subS);

        convexHull = F;

        /* Add the first facet to the list of facets of the hull:  */
        if (facets != NULL) { gch_list_add(facets, (void*)F); }

        /* The gift-wrapping loop: */
        gch_hyperplane_t *K = gch_hyperplane_new(dim);
        int *vertexTupleIndex = (int*)malloc(sizeof(int)*dim);
               
        while(!gch_list_empty(Q))
          { 
            if (verbose) { fprintf(stderr, "%*s    obtaining another ridge face\n", ind, ""); }
            /* Get some {dim-2}-face {G} from the horizon: */
            gem_ref_t G = (gem_ref_t )gch_list_first_content(Q);

            /* Get an affine basis of the {dim-2}-face {G}: */
            gch_get_vertex_tuple(G, dim-1, vertexTupleIndex);

            for (i = 0; i < dim; i++)
              { Basis[i] = gch_vector_list_real_index(S,vertexTupleIndex[i]); }

            /* Set {H} to an hyperplane containing {G}: */
            gch_hyperplane_from_points(H, dim, Basis);

            /* Find a point {q} of {S} not contained in {H}: */
            if (verbose) { fprintf(stderr, "%*s    finding the point {q}\n", ind, ""); }
            int q = 0;
            long test;
            while ((test = gch_hyperplane_point_test_sign(H, S->vectors[q])) == 0) { q++; }
            /* Negate {H} if necessary: */
            if (test < 0) { gch_hyperplane_neg(H); }

            /* Find the hyperplane {K}:  */
            Basis[dim-1] = S->vectors[q];
            gch_hyperplane_from_points(K, dim, Basis);

            if (gch_hyperplane_point_test_sign(K,gch_vector_list_real_index(S,vertexTupleIndex[dim-1])) < 0)
              { gch_hyperplane_neg(K); }
              

            /* Get a point {p} (not contained in {H}) from the new supporting hyperplane: */
            if (verbose) { fprintf(stderr, "%*s    finding the point {p}\n", ind, ""); }
            int p = gch_geometry_find_adj_facet(S, H, K);
        
            /* Get the new supporting hyperplane: */
            Basis[dim-1] = S->vectors[p];
            gch_hyperplane_from_points(H, dim, Basis);

            /* Store into {points} the points contained in {H}: */
            if (verbose) { fprintf(stderr, "%*s    collecting points of the new facet:", ind, ""); }
            numPoints = 0;
            for (i = 0; i < S->count; i++)
              { if (gch_hyperplane_point_test_sign(H, S->vectors[i]) == 0)
                  { if (verbose) { fprintf(stderr, " %d", S->realIndex[i]); }
                    points[numPoints] = i;
                    numPoints++;
                  }
              }
            if (verbose) { fprintf(stderr, "\n"); }
        
            /* Recursive call to find the new facet ({dim-1}-face) {F1} and its {dim-2}-faces {Q1}: */
            subS = gch_vector_list_sub_proj_new(S, points, numPoints, bt);
            gch_list_t *Q1 = gch_list_new();
            gem_ref_t F1 = gch_convexhull_find_rec(dim-1, subS, odim, Q1, cs);
            gch_vector_list_sub_free(subS);

            /* Make sure that {convexHull} is the `lowest' gem node: */
            if (gch_node_compare(F1,convexHull,dim-1) < 0) { convexHull = F1; }

            /* Report a ridge to the previous level: */
            if (facets != NULL) { gch_list_add(facets,F1); }

            /* Add to the horizon {Q} the {dim-2}-faces of {F1} that are not in the ridge, splice the rest: */
            if (verbose) { fprintf(stderr, "%*s    updating the horizon:", ind, ""); }
            while(!gch_list_empty(Q1))
              { gem_ref_t G1 = (gem_ref_t)gch_list_first_extract(Q1);

                /* Scan the list {Q} of {dim-2}-faces in the horizon: */
                gch_list_iterator_t iterQ;
                gch_list_reset(Q, &iterQ);
                while(!gch_list_at_end(iterQ))
                  { if (gch_node_compare((gem_ref_t)(gch_list_content(iterQ)),G1,dim-2) == 0) break;
                    gch_list_advance(&iterQ);
                  }
                if (!gch_list_at_end(iterQ))
                  { if (verbose) { fprintf(stderr, " S"); }
                    gem_ref_t G0 = (gem_ref_t)gch_list_content(iterQ);
                    gem_bary_splice(G1,G0,dim-2);
                    gch_list_remove(Q,iterQ);
                  }
                else
                  { if (verbose) { fprintf(stderr, " A"); }
                    gch_list_add(Q,(void*)G1); 
                  }
              }
            if (verbose) { fprintf(stderr, "\n"); }
            gch_list_free(Q1);
          }
        gch_affine_map_free(bt);
        free(vertexTupleIndex);
            
        gch_list_free(Q);
        gch_hyperplane_free(H);
        gch_hyperplane_free(K);

        free(basisIndex);
        free(Basis);
        free(points);
      }
    if (dim <= odim-2)
      { gch_computed_solutions_add(cs, dim, pi, convexHull, facets); }
    if (verbose) { fprintf(stderr, "%*sEXIT gch_convexhull_find_rec dim = %d\n", ind, "", dim); }
    return convexHull;
  }

gem_ref_t gch_convexhull_find(gch_vector_list_t *vl, int *dimP)
  {
    int verbose = 1;

    gch_vector_list_sub_t *S1 = gch_vector_list_sub_complete_new(vl);
    
    int odim = vl->dim; /* Nominal dimension of the gem.  The actual dimension depends of {vl}. */

    /* Get an affine basis for the spanned space, and its dimension {dim}: */
    if (verbose) { fprintf(stderr, "  obtaining a basis for the spanned subspace...\n"); }
    int basisIndex[vl->dim+1]; /* An affine basis for the space spanned by {vl}. */
    int dim = gch_geometry_affine_hull_basis_get(vl->dim, S1, basisIndex);

    /* Projects point set {S1} to a set {S} that spans {R^dim}: */
    gch_vector_list_sub_t *S;
    if (dim == vl->dim)
      { if (verbose) { fprintf(stderr, "    points span the whole R^%d\n", vl->dim); }
        S = S1;
      }
    else
      { if (verbose) { fprintf(stderr, "    projecting the points from R^%d to R^%d\n", vl->dim, dim); }
        int i;
        gch_vector_t **subspaceBasis = (gch_vector_t**)malloc(sizeof(gch_vector_t*)*(dim+1));
        for (i = 0; i < dim+1; i++) { subspaceBasis[i] = S1->vectors[basisIndex[i]]; }
        gch_affine_map_t *bt = gch_affine_map_new(vl->dim, dim, subspaceBasis);
        S = gch_vector_list_sub_complete_proj_new(vl, bt);
        gch_vector_list_sub_free(S1);
        gch_affine_map_free(bt);
        free(subspaceBasis);
      }

    if (verbose) { fprintf(stderr, "    finding the convex hull recursively\n"); }

    gch_computed_solutions_t *cs = gch_computed_solutions_new(dim);
    gem_ref_t gem = gch_convexhull_find_rec(dim, S, odim, NULL, cs);
    
    gch_computed_solutions_free(cs);
    gch_vector_list_sub_free(S);
    
    (*dimP) = dim;
    return gem;
  }

int gch_node_compare(gem_ref_t a, gem_ref_t b, int d)
  { demand((d >= 0) && (d <= gem_DIM_MAX), "invalid dimension");
    demand((gem_node_dim(a) >= d) && (gem_node_dim(b) >= d), "invalid node dimensions");
    int ad = gem_get_data(a);
    int bd = gem_get_data(b);
    int i,j;
    if (ad < bd) return -1;
    if (ad > bd) return +1;

    for(i = 0; i < d; i++) 
      { gem_ref_t at = a;
        gem_ref_t bt = b;
        for(j = i; j >= 0; j--) 
          { at = gem_step(at,j);
            bt = gem_step(bt,j);
          }
        ad = gem_get_data(at);
        bd = gem_get_data(bt);
        if (ad < bd) return -1;
        if (ad > bd) return +1;
      }
    return 0;
  }

void gch_get_vertex_tuple(gem_ref_t p, int d, int vdata[])
  {
    int i,j;
    gem_ref_t a =  p;
    vdata[0] = gem_get_data(a);
    for (i = 0; i < d; i++) 
      { a = p;
        for (j = i; j >= 0; j--) { a = gem_step(a,j); }
        vdata[i+1] = gem_get_data(a);
      }
  }
