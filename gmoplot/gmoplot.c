/* gmoplot.c -- Plots a geometric model. */
/* Last edited on 2007-12-27 01:36:54 by stolfi */

#include <r3.h>
#include <r4.h>
#include <hr3.h>
#include <r3x3.h>
#include <r4x4.h>
#include <rmxn.h>
#include <hr3.h>
#include <gmo.h>
#include <gmo_ps_plot.h>
#include <gmoplot_utils.h>
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

typedef struct options_t /* Parsed command line options */
  { char *gmoName;              /* Geometric model file (minus ".gmo" extension) or NULL */
    char *outName;              /* Output file (minus ".ps"/".eps" extension). */
    gmoplot_plot_options_t plt; /* Plot format and style options. */
  } options_t;

options_t *gmoplot_get_options(int argn, char **argc);

int main(int argn, char **argc)
  { options_t *o = gmoplot_get_options(argn, argc);

    gmoplot_plot_options_t *po = &(o->plt);
    fix_zenith(&(po->upp), &(po->obs), &(po->ctr), TRUE);

    /* Compute the projected plot window: */
    interval_t win[2];
    gmoplot_compute_window(&(po->obs), &(po->ctr), po->radius, &(po->upp), win);
    
    gmo_t *P = gmoplot_read_named_mesh(o->gmoName);
      
    /* Open the figure stream: */
    PSStream *fps = new_ps_stream
      (po->eps, o->outName, po->paperSize, po->figSize, po->caption.nel);

    /* Plot the triangulation: */
    pswr_new_picture(fps, -1.0,1.0, -1.0,1.0);
    
    gmoplot_plot_everything
      ( fps, P, 0,
        &(po->obs), &(po->ctr), win, &(po->upp),
        &(o->rcolor), &(o->wcolor), &(po->light),
        po->edgeWidth,
        TRUE
      );
    
    /* Add caption, if it is the case: */
    if(! po->eps) { pswr_add_caption(fps, o->outName, 0.0); }
    int k; 
    for (k = 0; k < po->caption.nel; k++)
      { pswr_add_caption(fps, po->caption.el[k], 0.0); }
    pswr_close_stream(fps);
    return 0;
  }
  
#define PPUSAGE argparser_SetUsage

options_t *gmoplot_get_options(int argn, char **argc)
  { options_t *o = notnull(malloc(sizeof(options_t)), "no mem");
    argparser_t *pp = argparser_new(stderr, argn, argc);
    argparser_set_usage(pp, 
      PROG_NAME " \\\n"
      "  -gmoName <name>  \\\n"
      "  -outName <name> \\\n" \
      PLOT_PARAMS_HELP " \n"
    );

    argparser_get_keyword(pp, "-gmoName");                           
    o->gmoName = argparser_get_next(pp);  

    argparser_get_keyword(pp, "-outName");                               
    o->outName = argparser_get_next(pp);  
    
    o->plt = parse_plot_options(pp);

    argparser_finish(pp);
    return o;
  }

