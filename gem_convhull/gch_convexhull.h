/* Last edited on 2014-07-22 13:47:27 by stolfilocal */
#ifndef gch_convexhull_H
#define gch_convexhull_H

#define _GNU_SOURCE

#include <gem.h>

#include <gch_vector_list.h>

gem_ref_t gch_convexhull_find(gch_vector_list_t *vl, int *dimP);
  /* Finds the convex hull of the points whose ???homo/cart??? coordinates are 
    the elements of {vl}. Returns some node of the barucentric subdivision of the hull.
    Also stores in {*dimP} the dimension of the space spanned by the points,
    which is the dimension of the gem. */

#endif
