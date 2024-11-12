/* See {grd_draw.h} */
/* Last edited on 2023-06-27 12:47:35 by stolfi */

#define gdr_draw_C_COPYRIGHT \
  "Duh?"

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#include <bool.h>
#include <cmp.h>
#include <jsmath.h>
#include <jsfile.h>
#include <jsrandom.h>
#include <argparser.h>
#include <in.h>
#include <rn.h>
#include <vec.h> 
#include <drtree.h> 
#include <drtree_compact.h> 
#include <drtree_planar.h> 
#include <drtree_lineage.h> 
#include <drtree_plot.h> 
#include <epswr.h> 

#include <gdr_demo.h>
#include <gdr_sim.h>

#include <gdr_draw.h>

/* INTERNAL PROTOTYPES */

drtree_node_t *gdr_draw_get_nodes_from_sim_state(gdr_sim_state_t *st, int32_t y0, int32_t y1);
  /* Creates a list of {drtree_node_t} from the simulation state {st},
    clipped to the interval between year {y0} and {y1}. */

/* IMPLEMENTATIONS */

void gdr_draw_state_plot
  ( char *prefix,
    gdr_sim_state_t *st,
    int32_t y0,
    int32_t y1,
    int32_t yRef
  )
  { bool_t debug = FALSE;
    if (debug) { fprintf(stderr, "  > %s\n", __FUNCTION__); }
    
    int32_t ni = st->ni; /* Number of individuals created. */
    int32_t ny = st->ny; /* Number of years in range of interest. */
    fprintf(stderr, "    found %d individuals in %d years of interest\n", ni, ny);
    
    /* Convert {st} to {drtree} format: */
    drtree_node_t *dt = gdr_draw_get_nodes_from_sim_state(st, y0, y1);

    /* Assign rows {rdr[0..ni-1]} to individuals: */
    int32_t rdr[ni]; 
    int32_t ncols; /* Number of columns in plot. */
    int32_t nrows; /* Number of rows in plot. */
    if ((yRef >= 0) && (yRef < ny))
      { drtree_planar_arrange(ni, dt, y0, y1, rdr, &ncols, &nrows); }
    else
      { drtree_compact_arrange(ni, dt, y0, y1, rdr, &ncols, &nrows); }
    fprintf(stderr, "    cell grid size = %d x %d\n", ncols, nrows);
    assert(ncols == y1 - y0 + 1);

    double Xstep = 1.0; /* Column width, mm per year. */
    double Ystep = 2.0; /* Row height, mm per individual. */
    epswr_figure_t *eps = drtree_plot_create_eps_figure(prefix, ncols, nrows, Xstep, Ystep);
    
    fprintf(stderr, "!! Surviving lineages not highlighted\n");
    
    frgb_t clr_y0 = (frgb_t){{ 1.000f, 0.300f, 0.300f }};
    frgb_t clr_y1 = (frgb_t){{ 0.450f, 0.450f, 1.000f }};
    frgb_t clr_yRef = (frgb_t){{ 0.150f, 0.800f, 0.150f }};
    drtree_plot_time_line(eps, y0, ncols, nrows, Xstep, Ystep, &clr_y0, 0);
    drtree_plot_time_line(eps, y0, ncols, nrows, Xstep, Ystep, &clr_y1, ny-1);
    
    int32_t *fnd = NULL; /* Founder of each surviving lineage at {yRef}, or {NULL}. */
    if ((yRef >= 0) && (yRef < ny)) 
      { /* Identify lineages at {yRef}: */
        drtree_plot_time_line(eps, y0, ncols, nrows, Xstep, Ystep, &clr_yRef, yRef);
        fnd = drtree_lineage_collect_surviving(ni, dt, yRef, ny-1);
      }

    /* Plot individuals: */
    drtree_plot_individuals(eps, y0, ncols, nrows, Xstep, Ystep, ni, dt, rdr, NULL, fnd);
    
    epswr_end_figure(eps);
    
    if (debug) { fprintf(stderr, "  < %s\n", __FUNCTION__); }
    return;
  }

drtree_node_t *gdr_draw_get_nodes_from_sim_state(gdr_sim_state_t *st, int32_t y0, int32_t y1)
  { bool_t debug = TRUE;
    if (debug) { fprintf(stderr, "    > %s\n", __FUNCTION__); }

    demand(y1 > y0, "invalid year range");
    int32_t ni = st->ni; /* Number of individuals created. */
    drtree_node_t *dt = (drtree_node_t*)notnull(malloc(ni*sizeof(drtree_node_t)), "no mem");
    for (int32_t iq = 0; iq < ni; iq++)
      { drtree_node_t *q = &(dt[iq]);
        int32_t ybr = st->ybr.e[iq];
        int32_t ydt = st->ydt.e[iq];
        int32_t par = st->par.e[iq];
        if (debug) { fprintf(stderr, "      iq= %d life span = {%d..%d} par = %d\n", iq, ybr, ydt, par); }
        if (ydt == -1)
          { /* Unexpanded node: set its life as trivial. */
            assert(st->nch.e[iq] == -1);
            assert(st->lch.e[iq] == -1);
            ydt = ybr;
          }
        assert((par == -1) || ((par >= 0) && (par < iq))); 
        q->tbr = ybr;
        q->tdt = ydt;
        q->par = par;
      }
    dt = drtree_clip_time_range(ni, dt, y0, y1);
    if (debug) { fprintf(stderr, "    < %s\n", __FUNCTION__); }
    return dt;
  }
