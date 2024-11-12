/* SPRange.h -- find the range of values of a spherical function. */
/* Last edited on 2005-06-06 11:22:13 by stolfi */

#ifndef SPRange_H
#define SPRange_H

#include <SPTriang.h>
#include <SPBasic.h>
#include <js.h>
 
double SPRange_OnTriangle
  ( S2Point *p, S2Point *q, S2Point *r, 
    ScalarField func, 
    int smpOrder
  );
  /* Computes the maximum absolute function value {fabs(func(p))} 
    over the spherical triangle with corners {p,q,r}, using 
    {(N+1)(N+2)/2} sample points, where {N = smpOrder}. */

double SPRange_OnTriangulation(ScalarField func, Triangulation *tri, int smpOrder);
  /* Calls {OnTriangle(func,p,q,r,smpOrder)} for for every triangle of {tri},
    and returns the max of all results. */

double SPRange_OnSphere(ScalarField func, int smpOrder);
  /* Calls {OnTriangle(func,p,q,r,smpOrder)} for 8 triangles that cover
    the sphere, and returns the max of all results. */

#endif
