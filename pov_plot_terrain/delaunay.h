#ifndef delaunay_H
#define delaunay_H

/* Delaunay triangulation. */
/* Last edited on 2007-12-24 00:34:59 by stolfi */

/*
** The divide-and-conquer algorithm for computing the Delaunay
** triangulation of a set of points (sites) on the plane. 
** For details, see 
**
**   "Primitives for the Manipulation of General Subdivisions 
**   and the Computation of Voronoi Diagrams"
**
**   L. Guibas, J. Stolfi, ACM Transactions on Graphics, April 1985
**
** Implemented by Jim Roth (DEC CADM Advanced Group) on May 1986.
** Adapted by J. Stolfi on April 1993.
** See the copyright notice at the end of this file.
*/

#include <quad.h>
#include <bool.h>
#include <sign.h>
#include <r2.h>


typedef struct delaunay_site_t{
	r2_t p;
	int index;
} delaunay_site_t;

/* The Delaunay triangulation: */

quad_arc_t delaunay_build(delaunay_site_t sites[], int nsites);
  /* Builds the Delaunay trinagulation of {site[0..nsites-1]}.
    Returns an arc {e} on the perimeter of the triangulation,
    oriented so that {LEFT(e)} is the outer (unbounded) face. */

/* Quad-edge data pointers: */

#define ORG(e) ((delaunay_site_t *) quad_odata(e))
#define DEST(e) ((delaunay_site_t *) quad_ddata(e))

#define SET_ORG(e,v) quad_odata(e) = (void *)(v)
#define SET_DEST(e,v) quad_ddata(e) = (void *)(v)

/* Topology tools: */

quad_arc_t connect(quad_arc_t a, quad_arc_t b);

/* Geometric tools: */

sign_t orient(delaunay_site_t *a, delaunay_site_t *b, delaunay_site_t *c);
  /* Orientation of the three sites {a,b,c}: positive means CCW,
    negative means CW, zero means they are collinear. */
   
bool_t leftof (delaunay_site_t *s, quad_arc_t e);
  /* TRUE iff the site {s} lies to the left of the line {ORG(e)-->DEST(e)}. */
    
bool_t rightof (delaunay_site_t *s, quad_arc_t e);
  /* TRUE iff the site {s} lies to the right of the line {ORG(e)-->DEST(e)}. */
    
bool_t incircle (delaunay_site_t *a, delaunay_site_t *b, delaunay_site_t *c, delaunay_site_t *d);
  /* TRUE iff the site {d} lies inside the circle {a,b,c}. */

#endif

/*
** Copyright notice:
**
** Copyright 1996 Institute of Computing, Unicamp.
**
** Permission to use this software for any purpose is hereby granted,
** provided that any substantial copy or mechanically derived version
** of this file that is made available to other parties is accompanied
** by this copyright notice in full, and is distributed under these same
** terms. 
**
** DISCLAIMER: This software is provided "as is" with no explicit or
** implicit warranty of any kind.  Neither the authors nor their
** employers can be held responsible for any losses or damages
** that might be attributed to its use.
**
** End of copyright notice.
*/
