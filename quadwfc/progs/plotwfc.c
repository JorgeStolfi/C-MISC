#define PROG_NAME "plotwfc"
#define PROG_DESC "generate a Postscript rendering of a geomodel and wavefront"
#define PROG_VERS "1.1"

/* Copyright © 2005 by the State University of Campinas (UNICAMP). */
/* See the copyright, authorship, and warranty notice at end of file. */
/* Last edited on 2011-12-29 20:32:23 by stolfilocal */

#define PROG_HELP \
  "  " PROG_NAME " \\\n" \
  "    -triName <name>  \\\n" \
  "    [ -geoName <name> ] \\\n" \
  "    -outName <name> \\\n" \
  "    [ -rcolor R.RRR G.GGG B.BBB ] \\\n" \
  "    [ -wcolor R.RRR G.GGG B.BBB ] \\\n" \
  PLOT_PARAMS_HELP

#define PROG_INFO \
  "NAME\n" \
  "  " PROG_NAME " - " PROG_DESC "\n" \
  "\n" \
  "SYNOPSIS\n" \
  PROG_HELP "\n" \
  "\n" \
  "DESCRIPTION\n" \
  "  ???\n" \
  "\n" \
  "AUTHORS\n" \
  "  Created in aug/2005 by Lucas Freitas and Jorge Stolfi at IC-UNICAMP."

#include <basic.h>
#include <mesh.h>
#include <geomodel.h>
#include <plot_utils.h>
#include <mesh_plot.h>
#include <geomodel_plot.h>

#include <jsfile.h>
#include <jsstring.h>
#include <r3.h>
#include <r4.h>
#include <hr3.h>
#include <r3x3.h>
#include <r4x4.h>
#include <rmxn.h>
#include <hr3.h>
//#include <plst_data.h>
//#include <plst_ps_plot.h>
#include <argparser_geo.h>

#include <frgb.h>
#include <pswr.h>
#include <affirm.h>
#include <argparser.h>
#include <bool.h>
#include <sign.h>

#include <values.h>
#include <stdio.h>
#include <math.h>

typedef struct options_t /* Parsed command line options */
  { char *triName;        /* Mesh file (minus ".tri" extension) */
    char *geoName;        /* Geophysical model file (minus ".geo" extension) or NULL */
    char *outName;        /* Output file (minus ".ps"/".eps" extension). */
    frgb_t rcolor;        /* Color of reflectors. */
    frgb_t wcolor;        /* Color of wavefront. */
    plot_options_t plt;   /* Plot format and style options. */
  } options_t;

void plot_everything
  ( PSStream *fps,           /* Poststcript stream, ready to plot. */
    mesh_t *tri,             /* mesh_t to draw, or NULL */
    geomodel_t *geo,         /* geophysical model to draw, or NULL */
    int N,                   /* Mesh subdivision parameter. */
    hr3_point_t *obs,        /* Observer's position. */
    hr3_point_t *ctr,        /* Center of interest. */
    double rad,              /* Radius of region of interest. */
    hr3_point_t *upp,        /* Camera vertical reference. */
    frgb_t *rcolor,          /* Reflector color. */
    frgb_t *wcolor,          /* Wavefront color. */
    r3_t *dLight,            /* Direction towards main light source. */
    double edgeWidth,        /* Nominal line width in mm. */
    bool_t verbose           /* TRUE to print messages along the way. */
  );
  /* Plots `everything' about the mesh {tri} and model {geo}. */

void plot_both_sides
  ( PSStream *fps,           /* Postscript stream. */
    char *fileTag,           /* Prefix for for page names. */
    mesh_t *tri,             /* mesh_t to draw, or NULL */
    geomodel_t *geo,         /* geophysical model to draw, or NULL */
    int N,                   /* Mesh subdivision parameter. */
    hr3_point_t *obs,        /* Observer's position */
    hr3_point_t *ctr,        /* Center of interest. */
    double rad,              /* Radius of region of interest. */
    hr3_point_t *upp,        /* Camera vertical reference */
    frgb_t *rcolor,          /* Reflector color. */
    frgb_t *wcolor,          /* Wavefront color. */
    r3_t *dLight,            /* Direction towards main light source. */
    double edgeWidth,        /* Nominal line width in mm. */
    string_vec_t *caption,   /* Figure caption, or NULL. */
    bool_t verbose           /* TRUE to print messages along the way. */
  );
  /* Plots two views of the the mesh {tri} and model {geo}, front and back.
    Each view is produced by {plot_everything} with the given
    parameters.
    
    If {fps} is an EPS stream, calls {pswr_new_canvas} before each view.
    The {figName} parameter is set to {fileTag} plus {-f} for front,
    or {-v} for back.
    
    The observer's position is {obs} for the front view, and its
    antipode for the back view.
    
    If the given {caption} is not NULL, it will be printed under each
    view. If the file is not EPS format, there will be an extra
    caption line with words {Front view} and {Back view}, as
    appropriate. However, captions are supressed if
    {pgs->captionLines} is zero. */

options_t *pwfc_parse_options(int argn, char **argc);

mesh_t *pwfc_read_named_mesh(char *triName);

geomodel_t *pwfc_read_named_geomodel(char *geoName);

int main(int argn, char **argc)
  { options_t *o = pwfc_parse_options(argn, argc);

    plot_options_t *po = &(o->plt);
    fix_zenith(&(po->upp), &(po->obs), &(po->ctr), TRUE);

    mesh_t *tri = pwfc_read_named_mesh(o->triName);
    geomodel_t *geo = pwfc_read_named_geomodel(o->geoName);
      
    /* Open the figure stream: */
    PSStream *fps = new_ps_stream
      (po->eps, o->outName, po->paperSize, FALSE, po->figSize, po->caption.ne);

    /* Plot the triangulation: */
    pswr_new_picture(fps, -1.0,1.0, -1.0,1.0);
    plot_everything
      ( fps, tri, geo, 0,
        &(po->obs), &(po->ctr), po->radius, &(po->upp),
        &(o->rcolor), &(o->wcolor), &(po->light),
        po->edgeWidth,
        TRUE
      );
    
    /* Add caption, if it is the case: */
    if(! po->eps) { pswr_add_caption(fps, o->outName, 0.0); }
    int k; 
    for (k = 0; k < po->caption.ne; k++)
      { pswr_add_caption(fps, po->caption.e[k], 0.0); }
    pswr_close_stream(fps);
    return 0;
  }
  
mesh_t *pwfc_read_named_mesh(char *triName)
  { FILE *rd = open_read(txtcat(triName, ".tri"), TRUE);
    mesh_t *tri = read_mesh(rd);
    fclose(rd);
    return tri;
  }
    
geomodel_t *pwfc_read_named_geomodel(char *geoName)
  { 
    if (geoName == NULL)
      { return NULL; }
    else
      { 
        FILE *rd = open_read(txtcat(geoName, ".geo"), TRUE);
        geomodel_t *geo = notnull(malloc(sizeof(geomodel_t)), "no mem");
        (*geo) = read_geomodel(rd);
        fclose(rd);
        return geo;
      }
  }
    
void plot_everything
  ( PSStream *fps,           /* Poststcript stream, ready to plot. */
    mesh_t *tri,             /* mesh_t to draw, or NULL */
    geomodel_t *geo,         /* geophysical model to draw, or NULL */
    int N,                   /* Mesh subdivision parameter. */
    hr3_point_t *obs,        /* Observer's position. */
    hr3_point_t *ctr,        /* Observer's position. */
    double rad,              /* Radius of region of interest. */
    hr3_point_t *upp,        /* Camera vertical reference. */
    frgb_t *rcolor,          /* Reflector color. */
    frgb_t *wcolor,          /* Wavefront color. */
    r3_t *dLight,            /* Direction towards main light source. */
    double edgeWidth,        /* Nominal line width in mm. */
    bool_t verbose           /* TRUE to print messages along the way. */
  )
  { hr3_pmap_t map = hr3_pmap_persp(obs, ctr, rad, upp);

    mumble("geophysical model...\n");;
    plot_geomodel(fps, &map, geo, rcolor, dLight, 0.1); 

    mumble("mesh...\n");;
    plot_mesh(fps, &map, tri, N, wcolor, dLight, 0.1); 

    mumble("Front axes...\n");;
    pswr_set_pen(fps, 0.0,0.0,0.0, 2.0 * edgeWidth, 0,0);
    draw_all_axes(fps, &map);
  }

void plot_both_sides
  ( PSStream *fps,          /* Postscript stream. */
    char *fileTag,          /* Prefix for for canvas names. */
    mesh_t *tri,            /* mesh_t to draw, or NULL */
    geomodel_t *geo,        /* geophysical model to draw, or NULL */
    int N,                  /* Mesh subdivision parameter. */
    hr3_point_t *obs,       /* Observer's position */
    hr3_point_t *ctr,       /* Observer's position */
    double rad,             /* Radius of region of interest. */
    hr3_point_t *upp,       /* Camera vertical reference */
    frgb_t *rcolor,          /* Reflector color. */
    frgb_t *wcolor,          /* Wavefront color. */
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
          tri, geo, N,
          &obsx, ctr, rad, upp,
          rcolor, wcolor, dLight, 
          edgeWidth,
          verbose
        );
        /* Add caption, etc: */
        if (caption != NULL)
          { int k; 
            for (k = 0; k < caption->ne; k++)
              { pswr_add_caption(fps, caption->e[k], 0.0); }
          }
        if (! fps->eps)
          { pswr_add_caption(fps, (side == 1 ? "Front view" : "Back view"), 0.0); }
        /* Reverse position of observer relative to ctr: */
        double obsw = obsx.c.c[0];
        for (k = 1; k <= 3; k++) 
          { obsx.c.c[k] = 2*obsw*ctr->c.c[k] - obsx.c.c[k]; }
      }
  }

options_t *pwfc_parse_options(int argc, char **argv)
  {
    argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, PROG_NAME " version " PROG_VERS ", usage:\n" PROG_HELP);
    argparser_set_info(pp, PROG_INFO);
    argparser_process_help_info_options(pp);
    
    options_t *o = notnull(malloc(sizeof(options_t)), "no mem");

    argparser_get_keyword(pp, "-triName");                           
    o->triName = argparser_get_next(pp);  

    if (argparser_keyword_present(pp, "-geoName"))                           
      { o->geoName = argparser_get_next(pp); }
    else
      { o->geoName = NULL; }

    if (argparser_keyword_present(pp, "-rColor"))
      { r3_t color = argparser_get_next_r3(pp, 0.0, 1.0);
        o->rcolor.c[0] = color.c[0];
        o->rcolor.c[1] = color.c[1];
        o->rcolor.c[2] = color.c[2];
      }
    else
      { o->rcolor = (frgb_t){{1.0, 0.9, 0.8}}; }

    if (argparser_keyword_present(pp, "-wColor"))
      { r3_t color = argparser_get_next_r3(pp, 0.0, 1.0);
        o->wcolor.c[0] = color.c[0];
        o->wcolor.c[1] = color.c[1];
        o->wcolor.c[2] = color.c[2];
      }
    else
      { o->wcolor = (frgb_t){{1.0, 0.9, 0.8}}; }

    argparser_get_keyword(pp, "-outName");                               
    o->outName = argparser_get_next(pp);  
    
    o->plt = parse_plot_options(pp);

    argparser_finish(pp);
    return o;
  }

