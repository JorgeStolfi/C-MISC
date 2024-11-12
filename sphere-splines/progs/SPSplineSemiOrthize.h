/* SPSplineSemiOrthize.h -- Semi-orthogonalization of spline bases */
/* Last edited on 2005-10-06 21:49:39 by stolfi */

#ifndef SPSplineSemiOrthize_H
#define SPSplineSemiOrthize_H

#include <SPTriang.h>
#include <SPFunction.h>
#include <SPSpline.h>
#include <SPH3.h>
#include <SPBasic.h>

#include <r3.h>
#include <vec.h>
#include <bool.h>
#include <nat.h>

#include <stdio.h>
#include <stdlib.h>

void SPSplineSemiOrthize_Basis
  ( Basis bas, 
    Triangulation *tri,
    bool_t weirdToo,
    Metric dot,
    bool_t precise,
    int gramSchmidt,
    int eigenfuncs,
    Metric mdot,
    bool_t verbose
  );
  /* Applies a partial orthonormalization procedure to the elements of
    basis {bas}.
    
    The elements must have been previously sorted in an order
    compatible with domain containment.
    
    The elements are processed in increasing order; each normal
    element {bas[i]}, with support {D}, is made {dot}-orthogonal to all
    elements {bas[j]} whose support is strictly contained in {D} (and
    which therefore must have been already processed). The {dot}
    product is assumed to be symmetric.
    
    Then, the elements in each support class are made {dot}-orthogonal to
    each other:
    
      If {gramSchmidt > 0}, the elements in the class are made
      {dot}-orthogonal by the Gram-Schmidt method. (This step is
      repeated {gramSchmidt} times for each class.)
    
      Them if {eigenfuncs > 0}, the elements in the class are replaced
      by an orthonormal set of stationary functions of the metric
      induced by {mdot}. (This step is repeated {eigenfuncs} times for
      each class.)
    
    Any weird elements in the basis must be at the end of {bas}. If
    {weirdToo} is true, they are orthogonalized among themselves, as
    described above. They are never combined with normal ones.
    
    In any case, all elements (including the weird ones) are
    made `mostly positive' (with {SPFunction_GenMakePositive(f,dot)})
    and normalized with {SPFunction_GenNormalize(f,dot)}.
    
    All elements in {bas} are modified by the procedure via the
    {scale} and {add}, and thus must be {add}-compatible with
    each other. 
  
    The {precise} parameter is used when building the sub-matrix of dot
    products used to make a class orthogonal to previously porcessed
    classes. If {precise=FALSE}, assumes that, after {bas[i]} and
    {bas[j]} have been processed, {dot(bas[i],bas[j]) is exactly 1
    when {i==j}, and zero when {i != j} and the supports are the same
    or nested.  If {precise=TRUE} makes no such assumptions and computes
    all dot products.  */

void SPSplineSemiOrthize_Check
  ( Basis bas, 
    int ini, int lim, 
    Triangulation *tri,
    bool_t weirdToo,
    Metric dot,
    double *maxOrthoError,
    double *maxNormError,
    double maxError,
    bool_t verbose
  );
  /* Checks whether the elements {bas[ini..lim-1]} are
    semi-orthonormalized. 
    
    More precisely, sets {maxOrthoError} to the maximum absolute value of
    all dot products {dot(bas[i],bas[j])} whose supports are nested;
    and sets {maxNormError} to the maximum absolute difference between
    {dot(bas[i],bas[i])} and 1. Stops when any of those errors becomes
    greater than {maxError}. */

#endif
