/* Cost functions for TSP/Cycle-cover problems. */
/* Last edited on 2023-03-31 04:19:55 by stolfi */

#ifndef tsp_cost_H
#define tsp_cost_H

/* Must define _GNU_SOURCE in order to get {asprintf} */
#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <bool.h>

#include <tsp_lib.h>

void gen_arc_cost_matrix(int32_t nv, Distr_t distr, int32_t d, double *c);
/* Stores in {c[0..nv*nv-1]} a symmetric cost matrix randomly selected
  from the distribution {distr}. When the latter specifies Euclidean
  costs, the vertex positions are chosen in {d}-dimensional space. The
  matrix is stored by rows, i.e. the cost of arc {(i,j)} is stored in
  {c[i*nv+j]}. */

double tour_cost(int32_t nv, vtx_t *v, double *c);
/* Interprets {v} as a tour, i.e. the set of arcs 
  {v[i]->v[(i+1)%nv]} for {i} in {0..nv-1}. Returns the 
  total cost, based on the arc cost matrix {c}. */

double perm_cost(int32_t nv, vtx_t *v, double *c);
/* Interprets {v} as a perm, i.e. the set of arcs 
  {i->v[i]} for {i} in {0..nv-1}. Returns the the 
  total cost, based on the arc cost matrix {c}. */
  
void compute_perm_costs(int32_t nv, bool_t tour, double *c, int32_t np, vtx_t *p, double *Z);
/* Computes costs {Z[0..np-1]} of the {np} permutations of {0..nv-1}
  stored in {p[0..np*nv-1]}, using the arc cost matrix {c[0..nv*nv-1]}.
  The perms are interpreted as tours if {tour=TRUE}, or as spins if
  {tour=FALSE}. */

void sort_costs(int32_t np, int32_t *ix, double *Z);
/* Performs index-sorting of the costs {Z[0..np-1]}.  Namely, 
  fills {ix[0..np-1]} with a permutation of the integers {0..np-1}, 
  such that {Z[ix[i]]} is sorted in increasing order. */

void gen_point_coords(Distr_t distr, int32_t d, double *x);
/* Stores in {x[0..d-1]} a point of {R^d}, randomly drawn from 
  the distribution specified by {distr}. */

double gen_arc_cost(Distr_t distr, int32_t d);
/* Returns an arc cost randomly drawn from the distribution
  specified by {distr} and {d}. */

void gen_gauss_point_coords(int32_t d, double *x);
/* Stores in {x[0..d-1]} a random point with symmetric Gaussian 
  distribution, centered at the origin, with 
  variance 1 in each coordinate (total variance {d}). */
  
void gen_uball_point_coords(int32_t d, double *x);
/* Stores in {x[0..d-1]} a random point with uniform 
  distribution in the ball centered at the origin with 
  unit radius.  Rather slow for large {d}. */
  
void gen_ucube_point_coords(int32_t d, double *x);
/* Stores in {x[0..d-1]} a random point with uniform 
  distribution in the cube centered at the origin with 
  side 2. */
  

#endif
