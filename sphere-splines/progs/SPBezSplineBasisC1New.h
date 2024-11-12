/* SPBezSplineBasisC1New.h -- New construction of C1 spline bases. */
/* Last edited on 2005-06-06 11:49:02 by stolfi */

#ifndef SPBezSplineBasisC1New_H
#define SPBezSplineBasisC1New_H

#include <SPTriang.h>
#include <SPFunction.h>
#include <SPSpline.h>
#include <SPHBezFunction.h>
#include <SPNHBezFunction.h>
#include <SPDeCasteljau.h>
#include <SPBezSplineBasisC1.h>
#include <vec.h>
#include <js.h>
#include <SPBasic.h>

void BuildEdgeBasisC1New
  ( int deg,
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

void BuildVertBasisC1New
  ( int deg,
    Arc e,
    Basis *bas,
    int *nb,
    char *triFile,
    Triangulation *tri
  );
  /* Builds the homogeneous ANS basis elements of continuity 1 and
    degree {deg} associated to vertex {v == Org(e)}. Stores them in
    {bas} starting at index {*nb}, and increments {*nb} accordingly.
    Only works for {deg >= 5}. */

r3x3_t VertexFrame(S2Point *u, S2Point *v);
  /* Computes an orthonormal frame of {R^3} where the first vector
    is collinear with {u}, and the second one is in the plane
    spanned by {u} and {v}. The frame axes are returned as COLUMNS of a
    3x3 matrix. */

r3x3_t EdgeFrame(S2Point *u, S2Point *v);
  /* Computes a frame of {R^3} where the first two vectors lie in the
    plane spanned by {u} and {v}, and the third vector is
    orthogonal to both; moreover, the first vector is orthogonal to {v}
    and has scalar 1 with {u}, the second vector is orthogonal to {u}
    and has scalar 1 with {v}, and the third vector has unit norm.
    The frame axes are returned as COLUMNS of a 3x3 matrix. */

#endif
