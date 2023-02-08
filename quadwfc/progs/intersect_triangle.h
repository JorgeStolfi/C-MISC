/* Ray-triangle intersection. */
/* Last edited on 2005-09-05 22:18:28 by stolfi */

#ifndef intersect_triangle_H
#define intersect_triangle_H

#include <basic.h>

void intersect_triangle 
  ( r3_t *u, r3_t *v, 
    r3_t *a, r3_t *b, r3_t *c, 
    /*out:*/ 
    r3_t *q, r3_t *nq, 
    double *alpha
  );
  /* Returns in {*q} the intersection point of segment {u,v} with
    triangle {a,b,c}. Also returns the triangle's normal in {*nq}, and
    the fraction until the intersection in {*alpha}. If there is no
    intersection, returns {*alpha = INF}. */

#endif
