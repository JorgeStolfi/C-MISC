/* See {triangulate.h}. */
/* Last edited on 2009-02-10 10:06:06 by stolfi */

#include <basic.h>
#include <triangulate.h>
#include <mesh.h>

#include <wavefront.h>

#include <r3.h>
#include <hr3.h>
#include <r3x3.h>

/* Delaunay triangulation by straightline divide-and-conquer. */

/* 
** Written by J. Stolfi on april 1993, based on an original
** implementation by Jim Roth (DEC CADM Advanced Group, May 1986).  
** See the copyright notice at the end of this file.
*/ 

/* #include <quad.h> */

/* Internal prototypes: */

bool_t rightof (r3_t *s, qarc_t e);
bool_t leftof (r3_t *s, qarc_t e);
bool_t incircle (r3_t *a, r3_t *b, r3_t *c, r3_t *d);
bool_t goes_after (r3_t *a, r3_t *b);

qarc_t connect(qarc_t a, qarc_t b);

void sort_sites(int ns, sref_t st[]);

void rec_delaunay(
    sref_t st[],           /* The sites */
    int sl, int sh,        /* Consider only sites[sl..sh] */
    qarc_t *le,            /* Output: leftmost and */
    qarc_t *re,            /*   rightmost edges of traingulation. */
    int level              /* Recursion level (for debugging). */
  );

qarc_t triangulate(int ns, sref_t st[])
  {
    qarc_t le, re;
    sort_sites(ns, st);
    rec_delaunay(st, 0, ns-1, &le, &re, 0);
    create_face_records(le);
    /* Make the outermost face into a hole: */
    face_t *f = RIGHT(le); f->omit = TRUE;
    return (le);
  }

/* Sort the sites into x order, breaking ties by y: */

void sort_sites(int ns, sref_t st[])
  {
    int i, j;
    for (i = 1; i < ns; i++)
      { for (j = i - 1; j >= 0 && goes_after(&(st[j]->curr.pos), &(st[j+1]->curr.pos)); j--) 
          { segment_t *tmp = st[j]; st[j] = st[j+1]; st[j+1] = tmp; }
      }
  }

bool_t goes_after(r3_t *a, r3_t *b)
  {
    return ( a->c[0] != b->c[0] ? (a->c[0] > b->c[0]) : (a->c[1] > b->c[1]) );
  }

/* Connect two vertices with a new edge: */

qarc_t connect(qarc_t a, qarc_t b)
  {
    qarc_t e;

    e = quad_make_edge();
    SET_quad_org(e, quad_dst(a));
    SET_quad_dst(e, quad_org(b));
    quad_splice(e, quad_lnext(a));
    quad_splice(quad_sym(e), b);
    return e;
  }

/* Recursively create the Delaunay triangulation of a sorted set of sites. */

void rec_delaunay(
    sref_t st[],
    int sl, int sh,
    qarc_t *le, 
    qarc_t *re,
    int level
  )
  {
    /* fprintf(stderr, "%*s  triangulating [%d..%d]", 2*level, "", sl, sh); */
    /* fprintf(stderr, "  X = [%5g _ %5g]", st->e[sl]->curr.pos.c[0], st->e[sh]->curr.pos.c[0]); */
    /* fprintf(stderr, "\n"); */
    
    if (sh == sl) 
      {	/* Only one sample. */
        affirm(FALSE, "cannot triangulate a single sample");
      }
    else if (sh == sl+1) 
      {
	/* Only two samples. */
        qarc_t a = quad_make_edge();
	SET_quad_org(a, st[sl]); 
        SET_quad_dst(a, st[sl+1]);
	*le = a; *re = quad_sym(a);
      }
    else if (sh == sl+2) 
      {
	/* Only three samples. */
        qarc_t a = quad_make_edge();
	qarc_t b = quad_make_edge();
        segment_t *u = st[sl];
	segment_t *v = st[sl+1];
	segment_t *w = st[sl+2];
	double ct = orient_xy(&(u->curr.pos), &(v->curr.pos), &(w->curr.pos));
	quad_splice(quad_sym(a), b);
	SET_quad_org(a, u); SET_quad_dst(a, v);
	SET_quad_org(b, v); SET_quad_dst(b, w);
	if (ct == 0) 
	  { *le = a; *re = quad_sym(b); }
	else 
	  { qarc_t c = connect(b, a);
	    if (ct > 0) 
	      { *le = a; *re = quad_sym(b); }
	    else 
	      { *le = quad_sym(c); *re = c; }
	  }
      }
    else
      {
	qarc_t ldo, ldi, rdi, rdo;
	qarc_t basel, lcand, rcand;

        int sm = (sl+sh)/2;

        rec_delaunay(st, sl, sm, &ldo, &ldi, level+1);
	rec_delaunay(st, sm+1, sh, &rdi, &rdo, level+1);

	while (1) 
          {
	    if (leftof(&ORGP(rdi), ldi))
              { ldi = quad_lnext(ldi); }
	    else if (rightof(&ORGP(ldi), rdi)) 
              { rdi = quad_onext(quad_sym(rdi)); }
	    else break;
	  }

	basel = connect(quad_sym(rdi), ldi);
	if (quad_org(ldi) == quad_org(ldo)) { ldo = quad_sym(basel); }
	if (quad_org(rdi) == quad_org(rdo)) { rdo = basel; }
        
	while (1) 
          {
	    /* if (debug) { debug_arc("basel =", basel); } */
            lcand = quad_onext(quad_sym(basel));
	    if (rightof(&DESTP(lcand), basel))
              { 
                while (incircle(&DESTP(basel), &ORGP(basel), &DESTP(lcand), &DESTP(quad_onext(lcand)))) 
                  { qarc_t t = quad_onext(lcand); 
                    /* if (debug) { debug_arc("killed left", lcand); } */
                    quad_destroy_edge(lcand); 
                    lcand = t;
                  }
              }

	    rcand = quad_oprev(basel);
	    if (rightof(&DESTP(rcand), basel))
	      { while (incircle(&DESTP(basel), &ORGP(basel), &DESTP(rcand), &DESTP(quad_oprev(rcand)))) 
                  { qarc_t t = quad_oprev(rcand); 
                    /* if (debug) { debug_arc("killed right", rcand); } */
                    quad_destroy_edge(rcand);
                    rcand = t;
                  }
              }

	    if (!rightof(&DESTP(lcand), basel) && !rightof(&DESTP(rcand), basel)) break;

	    if ( !rightof(&DESTP(lcand), basel) ||
		 ( rightof(&DESTP(rcand), basel) && 
                   incircle(&DESTP(lcand), &ORGP(lcand), &ORGP(rcand), &DESTP(rcand))
                 )
               )
	      { basel = connect(rcand, quad_sym(basel)); }
	    else
	      { basel = connect(quad_sym(basel), quad_sym(lcand)); }
	  }
	*le = ldo; *re = rdo;
      }
  }

/* Test if point to right of given edge: */

bool_t rightof(r3_t *s, qarc_t e)
  {
    return orient_xy(s, &DESTP(e), &ORGP(e)) > 0;
  }

/* Test if point to left of given edge: */

bool_t leftof(r3_t *s, qarc_t e)
  {
    return orient_xy(s, &ORGP(e), &DESTP(e)) > 0;
  }

/* Counterclockwise triangle predicate: */

sign_t orient_xy(r3_t *a, r3_t *b, r3_t *c)
  {
    double x1 = a->c[0], y1 = a->c[1];
    double x2 = b->c[0], y2 = b->c[1];
    double x3 = c->c[0], y3 = c->c[1];

    double det = (x2*y3-y2*x3) - (x1*y3-y1*x3) + (x1*y2-y1*x2);
    
    if (det < 0)
      { return -1; }
    else if (det > 0)
      { return +1; }
    else
      { return 0; }
  }

sign_t plane_side(r3_t *a, r3_t *b, r3_t *c, r3_t *d)
  {
    hr3_point_t ah = (hr3_point_t){{{ 1.0, a->c[0], a->c[1], a->c[2] }}};
    hr3_point_t bh = (hr3_point_t){{{ 1.0, b->c[0], b->c[1], b->c[2] }}};
    hr3_point_t ch = (hr3_point_t){{{ 1.0, c->c[0], c->c[1], c->c[2] }}};
    hr3_point_t dh = (hr3_point_t){{{ 1.0, d->c[0], d->c[1], d->c[2] }}};
    return hr3_orient(&ah, &bh, &ch, &dh);
  }

/* InCircle predicate: */

bool_t incircle(r3_t *a, r3_t *b, r3_t *c, r3_t *d)
  {
    if ((a == b) || (a == c) || (a == d) || (b == c) || (b == d) || (c == d)) { return FALSE; }
    
    double x1 = a->c[0], y1 = a->c[1];
    double x2 = b->c[0], y2 = b->c[1];
    double x3 = c->c[0], y3 = c->c[1];
    double x4 = d->c[0], y4 = d->c[1];

    return ((y4-y1)*(x2-x3)+(x4-x1)*(y2-y3))*((x4-x3)*(x2-x1)-(y4-y3)*(y2-y1)) >
	   ((y4-y3)*(x2-x1)+(x4-x3)*(y2-y1))*((x4-x1)*(x2-x3)-(y4-y1)*(y2-y3));
  }

