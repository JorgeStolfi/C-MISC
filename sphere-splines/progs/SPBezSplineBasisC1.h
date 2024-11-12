/* SPBezSplineBasisC1.h -- Construction of C1 spline bases. */
/* Last edited on 2005-06-06 11:44:42 by stolfi */

#ifndef SPBezSplineBasisC1_H
#define SPBezSplineBasisC1_H

#include <SPTriang.h>
#include <SPFunction.h>
#include <SPSpline.h>
#include <SPHBezFunction.h>
#include <SPNHBezFunction.h>
#include <SPDeCasteljau.h>
#include <SPBasic.h>
#include <vec.h>
#include <bool.h>

void BuildEdgeBasisC1
  ( int deg,
    bool_t newStyle,
    Arc e,
    Basis *bas,
    int *nb,
    char *triFile,
    Triangulation *tri
  );
  /* Builds the homogeneous ANS basis elements of continuity 1 and
    degree {deg} associated to arc {e} of {tri}. Stores them in {bas}
    starting at index {*nb}, and increments {*nb} accordingly. Only
    works for {deg >= 5}. */

void BuildVertBasisC1
  ( int deg,
    bool_t newStyle,
    Arc e,
    Basis *bas,
    int *nb,
    char *triFile,
    Triangulation *tri
  );
  /* Builds the homogeneous ANS basis elements of continuity 1 and
    degree {deg} associated to vertex {v == Org(e)}. Stores them in
    {bas} starting at index {*nb}, and increments {*nb} accordingly.
    Only works for {deg >= 5}.
    
    The {newStyle} parameter selects between the original ANS basis
    (FALSE, with {m+3} star elements at each vertex of order {m}) and
    the improved version (TRUE, with {m} boat elements + 3 star
    elements). */

r3_t ComputeEdgeContinuityCoeffs(Arc e, Triangulation *tri);
  /* Computes the `diamond rule' continuity coefficients for the
    arc {e}; namely, the numbers {bn}, {bo}, {bd}, such that
    {vp == bn*vn + bo*vo + bd*vd}, where
      {vp = Dest(Oprev(e)).pos}
      {vn = Dest(Onext(e)).pos}
      {vo = Org(e).pos}
      {vd = Dest(e).pos}
  */
  
void CheckPiecesForC1Continuity
  ( Triangulation *tri, 
    Arc a, 
    PieceData *pa, 
    PieceData *pb
  );
  /* Checks whether the Bezier spline pieces {pa} (for face {Left(a)})
    and {pb} (for face {Right(a)}) are C1-continuous across the arc
    {a}. 
    
    At present, {pa} and {pb} must be be spherical polynomials of the
    same degree, both homogeneous or both non-homogeneous, in
    barycentric Bezier form. */

void CheckCoeffsForC1Continuity
  ( Triangulation *tri, 
    Arc a, 
    int deg,
    BezCoeff_vec_t ca, 
    BezCoeff_vec_t cb
    );
    
#endif



