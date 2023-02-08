/* Last edited on 2014-07-22 03:19:22 by stolfilocal */
#ifndef gch_geometry_H
#define gch_geometry_H

#include <gch_vector.h>
#include <gch_vector_list.h>
#include <gch_hyperplane.h>

int gch_geometry_find_adj_facet(gch_vector_list_sub_t *S, gch_hyperplane_t *G, gch_hyperplane_t *K);
   /* Returns the list index of the chosen point. */

gch_hyperplane_t* gch_geometry_find_first_facet(int dim, gch_vector_list_sub_t *S, int *basisIndex); 
  /* Expects that {dim} is the dimension of the{S} points. On return, {basisIndex} will have
    the listIndexes of {dim} points of S that form a basis of the resulting hyperplane.
    The array {basisIndex} must have size {dim}. */
  
int gch_geometry_affine_hull_basis_get(int dim, gch_vector_list_sub_t *S, int *basisIndex); 
  /* The array {basisIndex} must have size {dim+1}. */

#endif
