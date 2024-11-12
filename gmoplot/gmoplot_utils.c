/* See {gmoplot_utils.h} */

/* Copyright © 2005 by the State University of Campinas (UNICAMP). */
/* See the copyright, authorship, and warranty notice at end of file. */
/* Last edited on 2011-12-29 20:31:49 by stolfilocal */

#include <gmoplot_utils.h>

#include <r3.h>
#include <r4.h>
#include <hr3.h>
#include <r3x3.h>
#include <r4x4.h>
#include <rmxn.h>
#include <hr3.h>
#include <gmo.h>
#include <gmo_ps_plot.h>
#include <argparser_geo.h>

#include <pswr.h>
#include <affirm.h>
#include <argparser.h>
#include <bool.h>
#include <sign.h>
#include <interval.h>

#include <values.h>
#include <stdio.h>
#include <math.h>

gmo_t *gmoplot_read_named_mesh(char *gmoName)
  { if (gmoName == NULL)
      { return NULL; }
    else
      { 
        FILE *rd = open_read(txtcat(gmoName, ".tri"), TRUE);
        gmo_t *P = read_mesh(rd);
        fclose(rd);
        return P;
      }
  }
    
void gmoplot_plot_everything
  ( PSStream *fps,           /* Poststcript stream, ready to plot. */
    gmo_t *P,                /* GEometric model to plot, or NULL */
    int N,                   /* Mesh subdivision parameter. */
    hr3_point_t *obs,        /* Observer's position. */
    hr3_point_t *ctr,        /* Observer's position. */
    interval_t win[],        /* XY window on image plane. */
    hr3_point_t *upp,        /* Camera vertical reference. */
    r3_t *dLight,            /* Direction towards main light source. */
    double edgeWidth,        /* Nominal line width in mm. */
    bool_t verbose           /* TRUE to print messages along the way. */
  )
  { double rad = xy_window_radius(win);
    hr3_pmap_t map = hr3_pmap_persp(obs, ctr, rad, upp);

    mumble("painting model...\n");;
    paint_faces(fps, &map, P, N, wcolor, dLight, 0.1); 

    mumble("Front axes...\n");;
    pswr_set_pen(fps, 0.0,0.0,0.0, 2.0 * edgeWidth, 0,0);
    draw_all_axes(fps, &map);
  }

void gmoplot_plot_both_sides
  ( PSStream *fps,          /* Postscript stream. */
    char *fileTag,          /* Prefix for for page names. */
    gmo_t *P,               /* Geometric model to plot, or NULL */
    int N,                  /* Mesh subdivision parameter. */
    hr3_point_t *obs,       /* Observer's position */
    hr3_point_t *ctr,       /* Observer's position */
    interval_t win[],       /* XY window on image plane. */
    hr3_point_t *upp,       /* Camera vertical reference */
    r3_t *dLight,           /* Direction towards main light source. */
    double edgeWidth,       /* Nominal line width in mm. */
    string_vec_t *caption,  /* Figure caption, or NULL. */
    bool_t verbose          /* TRUE to print messages along the way. */
  )
  { hr3_point_t obsx = (*obs);  /* Observer (to be flipped) */
    int side, k;

    for (side = +1; side >= -1; side -= 2)
      { if (fps->eps)
          { char *figTag = txtcat(fileTag, (side == 1 ? "-f" : "-v"));
            pswr_new_canvas(fps, figTag);
            free(figTag);
          }
        pswr_new_picture(fps, -1.0,1.0, -1.0,1.0);
        plot_everything(
          fps,
          P, N,
          &obsx, ctr, rad, upp,
          rcolor, wcolor, dLight, 
          edgeWidth,
          verbose
        );
        /* Add caption, etc: */
        if (caption != NULL)
          { int k; 
            for (k = 0; k < caption->nel; k++)
              { pswr_add_caption(fps, caption->el[k], 0.0); }
          }
        if (! fps->eps)
          { pswr_add_caption(fps, (side == 1 ? "Front view" : "Back view"), 0.0); }
        /* Reverse position of observer relative to ctr: */
        double obsw = obsx.c.c[0];
        for (k = 1; k <= 3; k++) 
          { obsx.c.c[k] = 2*obsw*ctr->c.c[k] - obsx.c.c[k]; }
      }
  }
