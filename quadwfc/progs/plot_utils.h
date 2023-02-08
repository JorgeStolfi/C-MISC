/* Basic plotting operations. */
/* Last edited on 2005-08-26 00:15:23 by stolfi */

#ifndef plot_utils_H
#define plot_utils_H

#include <basic.h>

#include <bool.h>
#include <sign.h>
#include <pswr.h>
#include <frgb.h>
#include <r3.h>
#include <r4.h>
#include <r4x4.h>
#include <hr3.h>
#include <argparser.h>

/* USEFUL POINTS */

#define Origin  (hr3_point_t){{{1.0, 0.0, 0.0, 0.0}}}
#define Zenith  (hr3_point_t){{{0.0, 0.0, 0.0, 1.0}}}
#define NoPoint (hr3_point_t){{{0.0, 0.0, 0.0, 0.0}}}

/* PLOTTING */

#define Invisible ((frgb_t){{-1,-1,-1}})
/* Use this color to omit the faces. */

PSStream *new_ps_stream
  ( bool_t eps,
    char *name,
    char *paperSize,
    bool_t landscape,
    double figSize,
    int nCap
  );
  /* Creates a new Postscript stream as specified by the given
    options: either a standalone document, called "{name}.ps", or a
    bunch of EPS files, called "{name}-{NNNNNN}.eps".
    The {figSize} is in millimeters. */

double default_fig_size
  ( bool_t eps, 
    char *paperSize, 
    int nRows, 
    int nCols,
    int captionLines
  );
  /* A convenient default figure size (in mm) for the given output options.
    If {eps = FALSE}, chooses the size so as to fit the given number of rows
    and columns in the page, with 1 inch margins.  If {eps = TRUE}, ignores
    the {paperSize} and returns 150.0 mm divided by the *minimum* of 
    {nCols,nRows}. */
    
/* MESH PLOTTING */

void paint_triangle
  ( PSStream *fps,
    r3_t *P, 
    r3_t *Q, 
    r3_t *R,
    hr3_pmap_t *map,         /* Perspective projection matrix. */
    frgb_t *color,           /* Color to use. */
    r3_t *dLight,            /* Direction towards main light source. */
    double shadow            /* Amount of darkening by shadow. */
  );
  /* Paints a colored triangle in perspective, given the Cartesian 
    coordinates of the original (unprojected) corner points. */

void draw_edge
  ( PSStream *fps,
    r3_t *P, 
    r3_t *Q, 
    hr3_pmap_t *map          /* Perspective projection matrix. */
  );
  /* Draws a line segment in perspective, given 
    the original (unprojected) endpoints. */

void draw_vertex
  ( PSStream *fps,
    r3_t *P, 
    hr3_pmap_t *map          /* Perspective projection matrix. */
  );
  /* Plots a dot at the point {P}, given the 
    original (unprojected) position. */

/* UTILITIES */

void draw_axis
  ( PSStream *fps, 
    hr3_pmap_t *map, 
    r3_t *dir,
    double length
  );
  /* Draws an arrow from the origin, with the specified length
    and direction. */
    
void draw_all_axes(PSStream *fps, hr3_pmap_t *map);
  /* Draws the three Cartesian coordinate axes X, Y, Z. */

/* FINDING PERSPECTIVE PROJECTION */

void fix_zenith(hr3_point_t *upp, hr3_point_t *obs, hr3_point_t *ctr, bool_t verbose);
  /* Replaces {*upp}, if necessary, by a point is a valid zenith
    refrence point for the observer {*obs} and center of interest
    {*ctr}: namely, a point that is not {[0,0,0,0]} and is not
    collinear with {*obs} and the origin.
     
    If the given {*upp} is already valid, leaves it unchanged. 
    Otherwise, replaces it either {SPPlot_DefaultZenith()}, or some other
    cardinal point at infinity. If {upp == NULL}, does nothing. */

hr3_point_t default_zenith(void);
  /* A default zenith reference point: currently the {Z}-infinity
    point {[0,0,0,1]}. */

/* PLOTTING OPTIONS */

typedef struct plot_options_t
  { /* General plotting options: */
    bool_t eps;           /* TRUE to generate ".eps" instead of ".ps". */
    char *paperSize;      /* Document paper size (when {eps == FALSE}). */
    double figSize;       /* Figure size in mm. */
    double meshSize;      /* Plotting resolution in mm. */
    string_vec_t caption; /* Figure caption. */
    /* View parameters: */
    hr3_point_t obs;      /* Position of observer. */
    hr3_point_t ctr;      /* Center of interest. */
    double radius;        /* Radius of interest. */
    hr3_point_t upp;      /* Approximate {up} reference point, or [0,0,0,0]. */
    /* Parameters for PaintValues: */
    r3_t light;           /* Light source direction */
    double shadow;        /* Relative shadow darkening factor */
    /* Parameters for PaintEdges: */
    double edgeWidth;     /* Line width in millimeters */
  } plot_options_t;  

#define PLOT_PARAMS_HELP \
  "  [ -eps | -ps ] [ -paperSize STRING ] \\\n" \
  "  [ -caption TEXT ]... \\\n" \
  "  [ -obs OW OX OY OZ ] \\\n" \
  "  [ -ctr CW CX CY CZ ] \\\n" \
  "  [ -radius NUM ] \\\n" \
  "  [ -up UW UX UY UZ ] \\\n" \
  "  [ -meshSize NUM ] \\\n" \
  "  [ -edgeWidth MILLIMETERS ] \\\n" \
  "  [ -light DX DY DZ ] [ -shadow NUM ]"

plot_options_t parse_plot_options(argparser_t *pp);
  /* Parses from the {pp} object a set of {SPplot_options_t}
    suitable for drawing triangulations only (without functions
    on them), as described by {PLOT_PARAMS_HELP}. */

#endif
