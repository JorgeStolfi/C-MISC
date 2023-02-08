/* See {delaunay_plot.h}. */
/* Last edited on 2009-01-06 04:22:27 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#include <jsstring.h>
#include <quad.h>
#include <pswr.h>
#include <bool.h>

#include <delaunay.h>
#include <delaunay_plot.h>

PSStream *open_ps_stream(double r, char *prefix, bool_t eps);
void close_ps_stream(PSStream *fps, bool_t eps);
delaunay_site_t far_left(delaunay_site_t *ap, delaunay_site_t *bp);
delaunay_site_t left_voronoi_vertex(quad_arc_t e);
double max_site_radius(delaunay_site_t *st, int nsites);
delaunay_site_t circumcenter(delaunay_site_t *ap, delaunay_site_t *bp, delaunay_site_t *cp);

#define RBIG 10.0

void draw_voronoi_edges(PSStream *fps, quad_arc_t e);
void draw_delaunay_edges(PSStream *fps, quad_arc_t e);
void draw_sites(PSStream *fps, delaunay_site_t *st, int nsites);

void plot_delaunay (quad_arc_t e, delaunay_site_t *st, int nsites, char *prefix, bool_t eps)
  { 
    double r = max_site_radius(st, nsites);
  
    /* Create Postscript document or EPS figure stream. */
    PSStream *fps = open_ps_stream(r, prefix, eps);
    
    /* Start a new picture: */
    double wm = 2.4*r;
    double xmin = -wm/2; double xmax = +wm/2;
    double ymin = -wm/2; double ymax = +wm/2;
    pswr_new_picture(fps, xmin,xmax, ymin, ymax);
    
    /* Plot Voronoi edges, dashed: */    
    float penwd = (eps ? 0.20 : 0.10);
    pswr_set_pen(fps, 0,0,0, penwd, 1.0, 1.0);
    draw_voronoi_edges(fps, e);
    
    /* Plot Delaunay edges, solid: */
    pswr_set_pen(fps, 0,0,0, penwd, 0.0, 0.0);
    draw_delaunay_edges(fps, e);
    
    /* Plot sites: */
    pswr_set_pen(fps, 0,0,0, penwd, 0.0, 0.0);
    pswr_set_fill_color(fps, 1.00,0.00,0.75);
    draw_sites(fps, st, nsites);

    if (! eps)
      { /* Add caption and frame: */
        pswr_set_pen(fps, 0,0,0, 0.10, 0.0, 0.0);
        pswr_add_caption(fps, "Voronoi/Delaunay diagram", 0.5);
        pswr_set_pen(fps, 0,0,0, 0.20, 0.0, 0.0);
        pswr_frame(fps);
      }
    /* We are done: */
    pswr_close_stream(fps);
  }

void draw_voronoi_edges(PSStream *fps, quad_arc_t e)
  {
    auto void draw_voronoi_edge(quad_arc_t e, void *closure);
    
    quad_enum(e, draw_voronoi_edge, NULL);

    void draw_voronoi_edge(quad_arc_t e, void *closure)
      { delaunay_site_t lp = left_voronoi_vertex(e);
        delaunay_site_t rp = left_voronoi_vertex(quad_sym(e));
        pswr_segment(fps, lp.p.c[0], lp.p.c[1], rp.p.c[0], rp.p.c[1]);
      }
   }

void draw_delaunay_edges(PSStream *fps, quad_arc_t e)
  {
    auto void draw_delaunay_edge(quad_arc_t e, void *closure);
    
    quad_enum(e, draw_delaunay_edge, NULL);
    
    void draw_delaunay_edge(quad_arc_t e, void *closure)
      { delaunay_site_t *op = quad_org(e);
        delaunay_site_t *dp = quad_dst(e);
        pswr_segment(fps, op->p.c[0], op->p.c[1], dp->p.c[0], dp->p.c[1]);
      }
 }

void draw_sites(PSStream *fps, delaunay_site_t *st, int nsites)
  {
    int i;
    for (i = 0; i < nsites; i++) 
      { double x = st[i].p.c[0];
        double y = st[i].p.c[1];
        pswr_dot(fps, x, y, 0.5,  TRUE, TRUE); 
      }
  }

  
delaunay_site_t left_voronoi_vertex(quad_arc_t e)
  { delaunay_site_t *ap =  quad_org(e);
    delaunay_site_t *bp =  quad_dst(e);
    delaunay_site_t *cp =  quad_dst(quad_lnext(e));
    if (orient(ap, bp, cp) > 0)
      { /* Internal face, compute circumcenter: */
        return circumcenter(ap, bp, cp);
      }
    else
      { /* External face, compute an "almost infinte" point on the left: */
        return far_left(ap, bp);
      }
  }

delaunay_site_t circumcenter(delaunay_site_t *ap, delaunay_site_t *bp, delaunay_site_t *cp)
  { double xa = ap->p.c[0], ya = ap->p.c[1];
    double xb = bp->p.c[0], yb = bp->p.c[1];
    double xc = cp->p.c[0], yc = cp->p.c[1];
    double A00 = 2*(xa-xb), A01 = 2*(ya-yb), B0 = (xa-xb)*(xa+xb)+(ya-yb)*(ya+yb);
    double A10 = 2*(xa-xc), A11 = 2*(ya-yc), B1 = (xa-xc)*(xa+xc)+(ya-yc)*(ya+yc);
    double det = A00*A11 - A01*A10;
    delaunay_site_t v;
    v.p.c[0] = (B0*A11 - B1*A01)/det;
    v.p.c[1] = (A00*B1 - A10*B0)/det;
    return v;
  }
  
double max_site_radius(delaunay_site_t *st, int nsites)
  { double r2max = 0.0;
    int i;
    for (i = 0; i < nsites; i++) 
      { double x = st[i].p.c[0];
        double y = st[i].p.c[1];
        double r2i = x*x + y*y;
        if (r2i > r2max) { r2max = r2i; }
      }
    return sqrt(r2max);
  }

delaunay_site_t far_left(delaunay_site_t *ap, delaunay_site_t *bp)
  { double xa = ap->p.c[0], ya = ap->p.c[1];
    double xb = bp->p.c[0], yb = bp->p.c[1];
    double xm = (xa + xb)/2, ym = (ya + yb)/2;
    double xd = (xb - xa), yd = (yb - ya);
    double d = sqrt(xd*xd + yd*yd);
    delaunay_site_t v;
    v.p.c[0] = xm - RBIG * yd/d;
    v.p.c[1] = ym + RBIG * xd/d;
    return v;
  }
  
PSStream *open_ps_stream(double r, char *prefix, bool_t eps)
  { 
    double mm = (72.0/25.4); /* One mm in pt. */
    double xfigsz = 150.00*mm; /* Figure X size excluding margin (pt). */
    double yfigsz = 150.00*mm; /* Figure Y size excluding margin (pt). */
    double fmrg = 3.0; /* Figure margin width (pt). */
    double pmrg = 2.0; /* Picture margin width (pt). */

    /* Add caption only if there is a user caption, or it is not EPS. */
    /* Select a good figure size: */
    PSStream *fps = pswr_new_stream
      ( /* prefix */                txtcat(prefix, "-"),
        /* file */                  NULL,
        /* eps */                   eps,
        /* docName */               "doc",
        /* paperSize */             "letter",
        /* landscape */             FALSE,
        /* hPageSize, vPageSize */  xfigsz + 2*fmrg, yfigsz + 2*fmrg
      );
    pswr_set_canvas_layout
      ( fps,
        /* hPicSize, vPicSize */     xfigsz, yfigsz,
        /* adjustPicSize */          FALSE,
        /* hPicMargin,vPicMargin */  pmrg, pmrg,
        /* captionLines */           (eps ? 0 : 1),  
        /* vCount, hCount */         0, 0  /* Let {pswr} choose it. */
      ); 
    return fps;
  }  
