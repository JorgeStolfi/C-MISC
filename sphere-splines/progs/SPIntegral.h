/* SPIntegral.h -- integration of functions on the sphere */
/* Last edited on 2005-06-06 11:24:18 by stolfi */

#ifndef SPIntegral_H
#define SPIntegral_H

#include <SPTriang.h>
#include <vec.h>
#include <js.h>

/* SIMPLE INTEGRALS */

/* These procedures accumulate their result on the argument {*sum},
  and in the correction term {*corr}, using Kahan's summation
  formula. The client must initialize both {sum} and {corr} with 0. */

void SPIntegral_OnFlatTriangle
  ( R3Point *u, R3Point *v, R3Point *w,
    double func(R3Point *p),
    double *sum, 
    double *corr
  );
  /* Computes the integral of function {func} over the flat (planar)
    triangle with corners {p,q,r}, using  Gaussian quadrature of order 5
    (error {O(h^6)}). */

void SPIntegral_OnSphericalTriangle
  ( S2Point *u, S2Point *v, S2Point *w,
    double func(S2Point *p),
    double *sum, 
    double *corr
  );
  /* Computes the integral of function {func} over the spherical
    triangle with corners {p,q,r}, using Gaussian quadrature of 
    order 5 (error {O(h^6)}) on each part. */

void SPIntegral_BySamples
  ( double func(S2Point *p),
    S2Point_vec_t sp,
    double_vec_t sw,
    double *sum,
    double *corr
  );
  /* Computes {SUM { sw[k]*func(sp[k]) : k = 0..N-1 }}, where {sp} is
    a client-specified set of sample points, {sw} their
    corresponding weights (measures of the domain elements),
    and {N = sp.ne}. */
 
/* INTEGRAL ON THE SPHERE */

double SPIntegral_OnSphere
  ( double func(S2Point *p),
    Triangulation *tri
  );
  /* Computes the integral of {func(p)} over the whole sphere, using
    the sample points and weights stored in {tri}. If {tri} is NULL,
    the procedure uses the default triangulation 
    {SPIntegral_GetDefaultTriangulation()}. */

/* GENERATING SAMPLE POINTS IN TRIANGLES */

/* For the tools in this section, the parameter {*ns} is the number of
  samples already stored in {sp}. The new samples are stored starting
  at {sp.e[*ns]}, and the value of {*ns} is incremented accordingly.
  The integration weight corresponding to {sp.e[i]} is stored in
  {sw.e[i]}. The arrays {sp} and {sw} are allocated and/or expanded
  as needed. */

void SPIntegral_GaussSampleTriangle
  ( S2Point *u, S2Point *v, S2Point *w, 
    S2Point_vec_t *sp, 
    double_vec_t *sw, 
    int *ns
  );
  /* Stores in the {sp[*ns..*ns+12]} and {sw[*ns..*ns+12]} the 13 
    integration knots and weights used by {SPIntegral_OnSphericalTriangle}
    for the spherical triangle with corners {u,v,w}. */

void SPIntegral_SuperSampleTriangle
  ( S2Point *u, S2Point *v, S2Point *w, 
    int smpOrder, 
    S2Point_vec_t *sp, 
    double_vec_t *sw, 
    int *ns
  );
  /* Stores in the vector {sp}, starting at {sp.e[*ns]}, a certain
    number of integration sample points distributed in the interior of
    the spherical triangle with corners {u,v,w}.
    
    Currently, the points are obtained by subdividing the triangle
    into a triangular mesh of {smpOrder^2} smaller triangles, and
    then picking 13 Gaussian knots in each triangle. */

void SPIntegral_RecursiveSampleTriangle
  ( S2Point *u, S2Point *v, S2Point *w, 
    int nMin,
    S2Point_vec_t *sp, 
    double_vec_t *sw, 
    int *ns
  );
  /* Stores in the vector {sp}, starting at {sp.e[*ns]}, at least
    {nMin} integration sample points distributed in the interior of
    the spherical triangle with corners {u,v,w}.
    
    If {nMin == 1}, picks the center of the triangle. Else, if {nMin
    <= 3}, picks 3 points. Else subdivides the triangle into 4 smaller
    triangles, and recursively picks at least {ceil(nMin/4)} samples
    in each part. */

void SPIntegral_SampleOctants
  ( int smpOrder, 
    S2Point_vec_t *sp, 
    double_vec_t *sw, 
    int *ns
  );
  /* Calls {SPIntegral_SuperSampleTriangle} for the eight octants
     of the sphere (the faces of the regular octahedral triangulation). */

void SPIntegral_SampleTriangulation
  ( Triangulation *tri, 
    int smpOrder,
    S2Point_vec_t *sp, 
    double_vec_t *sw, 
    int *ns
  );
  /* Calls {SPIntegral_SuperSampleTriangle} for each triangle of {tri}. */
  
/* DEFAULT SAMPLING ORDER AND DEFAULT TRIANGULATION */

void SPIntegral_SetDefaultSamplingOrder(int smpOrder);
int SPIntegral_GetDefaultSamplingOrder(void);
  /* Sets/gets the integration sampling order for triangles in
     triangulations. The default sampling order should not
     be changed once set. */

Triangulation *SPIntegral_GetDefaultTriangulation(void);
  /* Obtains the default triangulation (the regular octahedron). 
    Clients must call {SPIntegral_SetDefaultSamplingOrder} before
    calling this procedure. */

/* DEFAULT TRIANGULATION */

/* The following procedures define the default triangulation used
  when computing integrals for functions that don't have their
  own triangulation: */

void SPTriang_SetDefaultTriangulation(Triangulation *tri);
  /* Saves {tri} as the default triangulation. */

/* OBSOLETE */

double SPIntegral_OnTriangleXXX
  ( S2Point *u, S2Point *v, S2Point *w,
    double func(S2Point *p),
    int depth
  );
  /* Computes the integral of function {func} over the spherical
    triangle with corners {p,q,r}, using subdivision into {4^depth}
    subtriangles and Gaussian quadrature of order 5 (error {O(h^6)})
    on each part. */

double SPIntegral_OnTwoTrianglesXXX
  ( S2Point *pa, S2Point *qa, S2Point *ra,
    S2Point *pb, S2Point *qb, S2Point *rb,
    double func(S2Point *p),
    int depth
  );  
  /* Decomposes the intersection of triangles {pa,qa,ra}
    and {pb,qb,rb} into zero or more triangles,
    and calls {OnTriangle} on each of them. */

#endif
