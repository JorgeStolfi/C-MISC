/* See {geomodel_plot.h}. */
/* Last edited on 2023-02-12 11:45:59 by stolfi */

#include <geomodel_plot.h>

#include <basic.h>
#include <geomodel.h>
#include <plot_utils.h>

#include <r2.h>
#include <r3.h>
#include <r4.h>
#include <r3x3.h>
#include <r4x4.h>
#include <rmxn.h>
#include <hr3.h>
#include <epswr.h>
#include <affirm.h>
#include <bool.h>
#include <sign.h>
#include <frgb.h>

#include <values.h>
#include <stdio.h>
#include <math.h>

/* INTERNAL PROTOTYPES */

void paint_reflector(
    epswr_figure_t *eps,
    hr3_pmap_t *map,        /* Perspective projection matrix. */
    reflector_t *rf,        /* Reflector to plot. */
    frgb_t *color,          /* Base color of reflector. */
    r3_t *dLight,           /* Direction towards main light source. */
    double shadow           /* Amount of darkening by shadow. */
  );
  /* Paints the reflector {rfl} to file {eps}. */

/* IMPLEMENTATIONS */
    
void plot_geomodel(
    epswr_figure_t *eps,
    hr3_pmap_t *map,         /* Perspective projection matrix. */
    geomodel_t *geo,         /* Mesh to plot. */
    frgb_t *color,           /* Base color of mesh. */
    r3_t *dLight,            /* Direction towards main light source. */
    double shadow            /* Amount of darkening by shadow. */
  )  
  {
    int irf;
    for (irf = 0; irf < geo->nrf; irf++)
      { reflector_t *rf = &(geo->rf[irf]);
        paint_reflector(eps, map, rf, color, dLight, shadow);
      }
  }

void paint_reflector(
    epswr_figure_t *eps,
    hr3_pmap_t *map,        /* Perspective projection matrix. */
    reflector_t *rf,        /* Reflector to plot. */
    frgb_t *color,          /* Base color of reflector. */
    r3_t *dLight,           /* Direction towards main light source. */
    double shadow           /* Amount of darkening by shadow. */
  )  
  {
    auto void PaintTriangle(r3_t *P, r3_t *Q, r3_t *R);
    auto void DrawSegment(r3_t *P, r3_t *Q);
    
    void PaintTriangle(r3_t *P, r3_t *Q, r3_t *R)
      { paint_triangle(eps, P, Q, R, map, color, dLight, shadow); }
    
    void DrawSegment(r3_t *P, r3_t *Q)
      { draw_edge(eps, P, Q, map); }
      
    int nx = rf->nv[0];
    
    int xinf = rf->bb[0].end[0], xsup = rf->bb[0].end[1];
    double xstep = (xsup - xinf)/(nx-1);
    
    int ny = rf->nv[1];
    int yinf = rf->bb[1].end[0], ysup = rf->bb[1].end[1];
    double ystep = (ysup - yinf)/(ny-1);
    
    int ix, iy;
    int pass;
    for (pass = 0; pass < 2; pass++)
      {
        for (iy = 0; iy < ny; iy++)
          {
            r3_t S; /* Previous point on the previous row. */
            r3_t Q; /* Previous point on the same row. */
            for (ix = 0; ix < nx; ix++)
              { 
                r3_t P, R;
                P = (r3_t){{ xinf + ix*xstep, yinf + iy*ystep, rf->z[ix + nx*iy] }};
                if (iy > 0)
                  { R = (r3_t){{ P.c[0], P.c[1]-ystep, rf->z[ix + nx*(iy-1)] }};
                    if (ix > 0)
                      {
                        if (pass == 0)
                          {
                            PaintTriangle(&S, &Q, &R);
                            PaintTriangle(&P, &Q, &R);
                          }
                        else if (pass == 1)
                          { DrawSegment(&Q, &R); }
                      }
                    if (pass == 1) { DrawSegment(&R, &P); }
                    S = R;
                  }
                if (ix > 0) { if (pass == 1) { DrawSegment(&P, &Q); } }
                Q = P;
              }
          }
      }
  }
