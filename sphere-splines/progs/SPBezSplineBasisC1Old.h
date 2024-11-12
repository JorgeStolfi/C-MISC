/* SPBezSplineBasisC1Old.h -- Original ANS construction of C1 spline bases. */
/* Last edited on 2005-06-06 11:48:53 by stolfi */

#ifndef SPBezSplineBasisC1Old_H
#define SPBezSplineBasisC1Old_H

#include <SPTriang.h>
#include <SPFunction.h>
#include <SPSpline.h>
#include <SPHBezFunction.h>
#include <SPNHBezFunction.h>
#include <SPDeCasteljau.h>
#include <SPBezSplineBasisC1.h>
#include <SPBasic.h>
#include <vec.h>
#include <js.h>

void BuildEdgeBasisC1Old
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

void BuildVertBasisC1Old
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
  
#endif



