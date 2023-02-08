/* Delaunay triangulation of points on the plane. */
/* Last edited on 2005-09-15 18:48:12 by stolfi */

#ifndef triangulate_H
#define triangulate_H

#include <basic.h>
#include <quad.h>
#include <wavefront.h>

qarc_t triangulate(int ns, sref_t st[]);
  /* The divide-and-conquer algorithm for computing the Delaunay
    triangulation of a set of sites {st[0..ns-1]} (points of the
    plane). Returns an arc {e} on the triangulation, such that
    {RIGHT(e)} is the outermost (unbounded) face.

    Considers only the X and Y coordinates of the {curr} field of each
    sample, ignoring the Z coordinate. The {ORG} and {DEST} fields of
    every arc will point to a sample record from {st}. The {LEFT} and
    {RIGHT} fields will point to newly created {face_t} records. The
    outermost face will have {omot=TRUE}, all the others will have
    {omit=FALSE}.

    For details, see "Primitives for the Manipulation of General
    Subdivisions and the Computation of Voronoi Diagrams", by L. Guibas,
    J. Stolfi, ACM Transactions on Graphics, April 1985. */

/* Geometric tools: */

sign_t orient_xy (r3_t *a, r3_t *b, r3_t *c);
  /* Orientation of triangle {a,b,c} projected on the XY plane: 
    -1 clockwise, 0 degenerate (collinear), +1 counterclockwise. */

sign_t plane_side(r3_t *a, r3_t *b, r3_t *c, r3_t *d);
  /* Orientation of triangle {a,b,c} as seen from {d}:  
    -1 clockwise, 0 degenerate (the four points are coplanar),
    +1 counterclockwise. */

#endif

