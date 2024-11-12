/* SPBezSplineBasisC0.h --  Internal procs for building C0 spline bases. */
/* Last edited on 2005-06-06 11:44:33 by stolfi */

#ifndef SPBezSplineBasisC0_H
#define SPBezSplineBasisC0_H

#include <SPTriang.h>
#include <SPFunction.h>
#include <SPSpline.h>
#include <SPHBezFunction.h>
#include <SPNHBezFunction.h>
#include <SPDeCasteljau.h>
#include <SPBasic.h>
#include <vec.h>
#include <js.h>

void BuildEdgeBasisC0
  ( int deg,
    Arc e,
    Basis *bas,
    int *nb,
    char *triFile,
    Triangulation *tri
  );
  /* Builds the homogeneous ANS basis elements of continuity 0 and 
    degree {deg} associated to edge {Edge(e)} of {tri}.  Stores them in {bas}
    starting at index {*nb}, and increments {*nb} accordingly. 
    
    Only works for {deg >= 1}. Barring degeneracies, there are {deg-1}
    C0 edge elements for each edge {e}, all of them spanning the two
    triangles adjecent to {e}. */

void BuildVertBasisC0
  ( int deg,
    Arc e,
    Basis *bas,
    int *nb,
    char *triFile,
    Triangulation *tri
  );
  /* Builds the homogeneous ANS basis elements of continuity 0 and degree
    {deg} associated to vertex {v == Org(e)}. Stores them in {bas}
    starting at index {*nb}, and increments {*nb} accordingly. 
    
    Only works for {deg >= 1}. Actually, there is only one C0 vertex
    element for each vertex. */

#endif
