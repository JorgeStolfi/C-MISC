/* SPDelaunay.h -- Delaunay triangulation on a sphere. */
/* Last edited on 2005-06-06 11:25:43 by stolfi */ 

#ifndef SPDelaunay_H
#define SPDelaunay_H

/*
  Delaunay triangulation of a set of points on the sphere.
  See {Primitives for the Manipulation of General Subdivisions and the
  Computation of Voronoi Diagrams}, by L. Guibas and J. Stolfi, ACM
  TOG, April 1985.

  The procedures in this module use only the quad-edge structure and
  the {Site} coordinates. The structure can be interpreted either as a
  Delaunay triangulation of the sphere, or as a convex polyhedron.
  These procedures do not use or create the {Face} records, and,
  except where stated, do not expect edges and sites to have
  meaningful {num} fields. When the structure is finished, use
  {SPTriang_BuildTables} to create these missing fields.
*/

#include <SPTriang.h>
#include <SPQuad.h>
#include <SPBasic.h>
#include <vec.h>
#include <bool.h>
#include <stdio.h>

Arc SPDelaunay_BuildInc(SiteRef_vec_t site);
  /* Builds the convex hull of the given sites incrementally, by
    creating a tetrahedron with the first four sites and inserting the
    others one by one.  Returns an arc of the hull. */

bool_t SPDelaunay_EdgeIsDelaunay(Arc e);    
  /* TRUE iff the edge {e} has the local Delaunay property, namely if 
    it is not degenerate and is a convex edge of the polyhedron. */

void SPDelaunay_CheckDelaunay(Arc_vec_t arc);
  /* Checks consistency of a spherical Delaunay triangulation, namely
    whether the corresponding polyhedron is convex. Aborts if
    inconsistent. Will accept a polyhedron that wraps more than once
    around the origin. */

Arc SPDelaunay_InsertSite(Site *s, Arc e);
  /* Adds a new site {s} to the convex hull {e}.
    Returns an arc of the new hull. */

Arc SPDelaunay_BuildTetrahedron (Site *a, Site *b, Site *c, Site *d); 
  /* Builds a convex tetrahedron with vertices {a,b,c,d}.
    Returns one arc of the tetrahedron. The topology is such
    that {Onext} and {Lnext} move counterclockwise when seen
    from outside the hull of those four points. */
    
void SPDelaunay_Retriangulate(Arc e, int deg, Arc_vec_t *tr);
  /* Retriangulates the face Left(e), resulting from the elimination of
    a vertex, which must have {deg} sides.
    Also returns an array {tr[0..nTr]} containing one side arc out 
    of each new triangle. */
  
bool_t SPDelaunay_FixDelaunay(Triangulation *tri);
  /* Attempts to swap edges in a spherical triangulation so as to make
    it into a Delaunay triangulation (i.e. so as to make the
    corresponding polyhedron convex). Vertex records and position are
    preserved; edge and face records are reused. Requires edge numbers
    to be consistent with their indices, i.e. {EdgeNum(arc.e[i]) ==
    i/2} for all {i}. Recomputes the coordinate matrices.
    
    Returns TRUE on success, FALSE on failure (namely when the
    triangulation is so messed up that edge swapping doesn't work.) */
  
/* DEFAULT TRIANGULATION */

Triangulation *SPDelaunay_RegularTetrahedron(int smpOrder);
Triangulation *SPDelaunay_RegularOctahedron(int smpOrder);
Triangulation *SPDelaunay_RegularIcosahedron(int smpOrder);
  /* Constructs the named regular convex polyhedra. Each octant will
    be sampled with order {smpOrder} (see {SPIntegral_SampleTriangle}). */

#endif

