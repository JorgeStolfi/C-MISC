/* Delaunay triangulation by straightline divide-and-conquer. */
/* Last edited on 2007-12-24 00:35:52 by stolfi */ 

/* 
** Written by J. Stolfi on april 1993, based on an original
** implementation by Jim Roth (DEC CADM Advanced Group, May 1986).  
** See the copyright notice at the end of this file.
*/ 

#include <delaunay.h>

#include <quad.h>
#include <bool.h>
#include <sign.h>
#include <sign_get.h>

/* Internal prototypes: */

void sort_sites(delaunay_site_t sites[], int nsites);

void rec_delaunay(
    delaunay_site_t sites[],  /* The sites */
    int sl, int sh,       /* Consider only sites[sl..sh-1] */
    quad_arc_t *le,       /* Output: leftmost and */
    quad_arc_t *re        /*   rightmost edges of traingulation. */
  );

quad_arc_t delaunay_build(delaunay_site_t sites[], int nsites)
  {
    quad_arc_t le, re;
    sort_sites(sites, nsites);
    rec_delaunay(sites, 0, nsites, &le, &re);
    return (le);
  }

/* Shell-sort the sites into x order, breaking ties by y: */

void sort_sites(delaunay_site_t sites[], int nsites)
  {
    int gap, i, j;
    delaunay_site_t tmp;

    for (gap = nsites/2; gap > 0; gap /= 2)
      for (i = gap; i < nsites; i++)
	for (
	    j = i-gap; 
	    j >= 0 && 
	      ( sites[j].p.c[0] != sites[j+gap].p.c[0] ? 
                (sites[j].p.c[0] > sites[j+gap].p.c[0]) : 
                (sites[j].p.c[1] > sites[j+gap].p.c[1])
              );
	    j -= gap
	  ) 
	  {
	    tmp = sites[j]; sites[j] = sites[j+gap]; sites[j+gap] = tmp;
	  }
  }

/* Connect two vertices with a new edge: */

quad_arc_t connect(quad_arc_t a, quad_arc_t b)
  {
    quad_arc_t e;

    e = quad_make_edge();
    SET_quad_org(e, quad_dst(a));
    SET_quad_dst(e, quad_org(b));
    quad_splice(e, quad_lnext(a));
    quad_splice(quad_sym(e), b);
    return e;
  }

/* Recursively create the Delaunay triangulation of a sorted set of sites. */

void rec_delaunay(
    delaunay_site_t sites[],
    int sl, int sh,
    quad_arc_t *le, quad_arc_t *re
  )
  {
    if (sh == sl+2) 
      {
	quad_arc_t a = quad_make_edge();
	SET_quad_org(a, &sites[sl]); 
        SET_quad_dst(a, &sites[sl+1]);
	*le = a; *re = quad_sym(a);
      }
    else if (sh == sl+3) 
      {
	quad_arc_t a = quad_make_edge();
	quad_arc_t b = quad_make_edge();
	int ct = orient(&sites[sl], &sites[sl+1], &sites[sl+2]);
	quad_splice(quad_sym(a), b);
	SET_quad_org(a, &sites[sl]); 
        SET_quad_dst(a, &sites[sl+1]);
	SET_quad_org(b, &sites[sl+1]);  
        SET_quad_dst(b, &sites[sl+2]);
	if (ct == 0.0) 
	  { *le = a; *re = quad_sym(b); }
	else 
	  { quad_arc_t c = connect(b, a);
	    if (ct > 0.0) 
	      { *le = a; *re = quad_sym(b); }
	    else 
	      { *le = quad_sym(c); *re = c; }
	  }
      }
    else
      {
	quad_arc_t ldo, ldi, rdi, rdo;
	quad_arc_t basel, lcand, rcand;

        int sm = (sl+sh)/2;

        rec_delaunay(sites, sl, sm, &ldo, &ldi);
	rec_delaunay(sites, sm, sh, &rdi, &rdo);

	while (1) 
          {
	    if (leftof(quad_org(rdi), ldi)) ldi = quad_lnext(ldi);
	    else if (rightof(quad_org(ldi), rdi)) rdi = quad_onext(quad_sym(rdi));
	    else break;
	  }

	basel = connect(quad_sym(rdi), ldi);
	if (quad_org(ldi) == quad_org(ldo)) ldo = quad_sym(basel);
	if (quad_org(rdi) == quad_org(rdo)) rdo = basel;

	while (1) 
          {

	    lcand = quad_onext(quad_sym(basel));
	    if (rightof(quad_dst(lcand), basel))
	      while (incircle(quad_dst(basel), quad_org(basel), quad_dst(lcand), quad_dst(quad_onext(lcand)))) 
                { quad_arc_t t = quad_onext(lcand); quad_destroy_edge(lcand); lcand = t; }

	    rcand = quad_oprev(basel);
	    if (rightof(quad_dst(rcand), basel))
	      while (incircle(quad_dst(basel), quad_org(basel), quad_dst(rcand), quad_dst(quad_oprev(rcand)))) 
                { quad_arc_t t = quad_oprev(rcand); quad_destroy_edge(rcand); rcand = t; }

	    if (!rightof(quad_dst(lcand), basel) && !rightof(quad_dst(rcand), basel)) break;

	    if ( !rightof(quad_dst(lcand), basel) ||
		 ( rightof(quad_dst(rcand), basel) && 
                   incircle(quad_dst(lcand), quad_org(lcand), quad_org(rcand), quad_dst(rcand))
                 )
               )
	      basel = connect(rcand, quad_sym(basel));
	    else
	      basel = connect(quad_sym(basel), quad_sym(lcand));
	  }
	*le = ldo; *re = rdo;
      }
  }

/* Test if point to right of given edge: */

bool_t rightof(delaunay_site_t *s, quad_arc_t e)
  {
    return orient(s, quad_dst(e), quad_org(e)) > 0;
  }

/* Test if point to left of given edge: */

bool_t leftof(delaunay_site_t *s, quad_arc_t e)
  {
    return orient(s, quad_org(e), quad_dst(e)) > 0;
  }

/* Counterclockwise triangle predicate: */

sign_t orient(delaunay_site_t *a, delaunay_site_t *b, delaunay_site_t *c)
  {
    double x1 = a->p.c[0], y1 = a->p.c[1];
    double x2 = b->p.c[0], y2 = b->p.c[1];
    double x3 = c->p.c[0], y3 = c->p.c[1];
    double det = (x2*y3-y2*x3) - (x1*y3-y1*x3) + (x1*y2-y1*x2);
    return sign_double(det);
  }

/* InCircle predicate: */

bool_t incircle(delaunay_site_t *a, delaunay_site_t *b, delaunay_site_t *c, delaunay_site_t *d)
  {
    if ((a == b) || (a == c) || (a == d) || (b == c) || (b == d) || (c == d)) 
      { return FALSE; }

    double x1 = a->p.c[0], y1 = a->p.c[1];
    double x2 = b->p.c[0], y2 = b->p.c[1];
    double x3 = c->p.c[0], y3 = c->p.c[1];
    double x4 = d->p.c[0], y4 = d->p.c[1];
    
    double da = ((y4-y1)*(x2-x3)+(x4-x1)*(y2-y3))*((x4-x3)*(x2-x1)-(y4-y3)*(y2-y1));
    double db = ((y4-y3)*(x2-x1)+(x4-x3)*(y2-y1))*((x4-x1)*(x2-x3)-(y4-y1)*(y2-y3));

    return da > db;
  }

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
** NOTE: this copyright notice does not claim to supersede any copyrights
** that may apply to the original DEC implementation of the quad-edge
** data structure.
**
** DISCLAIMER: This software is provided "as is" with no explicit or
** implicit warranty of any kind.  Neither the authors nor their
** employers can be held responsible for any losses or damages
** that might be attributed to its use.
**
** End of copyright notice.
*/
 
