#define PROG_NAME "plotmap"
#define PROG_DESC "generate a Postscript plot of a street map"
#define PROG_VERS "2023-02-21"

#define PROG_COPYRIGHT "Copyright © 2023 Universidade Estadual de Campinas (UNICAMP)"
/* Last edited on 2024-06-22 17:39:34 by stolfi */

#define PROG_HELP \
  "  PROG_NAME \\\n" \
  "    [ -size XSIZE YSIZE ] \\\n" \
  "    [ -grid ]\\\n" \
  "    [ -clip XMIN XMAX YMIN YMAX ] \\\n" \
  "    [ -showRect XMIN XMAX YMIN YMAX ]... \\\n" \
  "    [ -showCircle XCTR YCTR RADIUS ]... \\\n" \
  "    [ -showBall XCTR YCTR DMAX ]... \\\n" \
  "    < MAP.rnt \\\n" \
  "    > PLOT.eps"
  
#define PROG_INFO \
  "SYNOPSIS\n" \
  "" PROG_HELP "\n" \
  "\n" \
  "DESCRIPTION\n" \
  "  Reads a street map description from {stdin}, produces an" \
  " Encapsulated Postscript plot of the same on {stdout}.  See" \
  " the map file format in the {stmap.h} interface.\n" \
  "\n" \
  "OPTIONS\n" \
  "  !!! TO BE WRITTEN !!!\n" \
  "\n" \
  "DOCUMENTATION OPTIONS\n" \
  argparser_help_info_HELP_INFO "\n" \
  "\n" \
  "SEE ALSO\n" \
  "  gnuplot(1).\n" \
  "\n" \
  "AUTHOR\n" \
  "  Created in ~2002 by Jorge Stolfi, IC-UNICAMP.\n" \
  "\n" \
  "MODIFICATION HISTORY\n" \
  "  2023-02-21 J.Stolfi: Converted to use {epswr.h} instead of {pswr}\n" \
  "  2023-02-21 J.Stolfi: Converted to use {argparser.h}\n" \
  "  2023-02-21 J.Stolfi: Converted {int} to {int32_t} etc.\n" \
  "\n" \
  "WARRANTY\n" \
  argparser_help_info_NO_WARRANTY "\n" \
  "\n" \
  "RIGHTS\n" \
  "  " PROG_COPYRIGHT ".\n" \
  "\n" \
  argparser_help_info_STANDARD_RIGHTS

#define _GNU_SOURCE
#include <values.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>

#include <epswr.h>
#include <argparser.h>
#include <rn.h>

#include <stmap.h>

#define pltmap_MAXMARKS  1000

typedef struct options_t
  { double xSize, ySize;     /* Figure dimensions (mm). */
    char *paper;             /* Paper size for non-encapsulated output. */
    double grid;             /* Coordinate grid spacing, or 0 for no grid. */
    Interval xClip, yClip;   /* Coordinate clipping ranges (m). */
    
    int32_t nRect;           /* Number of Euclidean circles to draw. */
    Interval *xRect;         /* X-ranges of rectangles to draw (m). */
    Interval *yRect;         /* Y-ranges of rectangles to draw (m). */
    
    int32_t nCircle;         /* Number of Euclidean circles to draw. */
    double *rCircle;         /* Radii of Euclidean circles to draw (m). */
    Point *ctrCircle;        /* Centers of Euclidean circles to draw (m). */
    
    int32_t nBall;           /* Number of street-cost balls to draw. */
    double *rBall;           /* Radii of street-cost balls to draw (m). */
    Point *ctrBall;          /* Centers of street-cost balls to draw (m). */
    
  } options_t;
  /* Command line options passed to the program. */
    
/* PROTOTYPES */

int32_t main(int32_t argc, char **argv);
options_t *parse_options(int32_t argc, char **argv);
void get_next_arg_string(char **varp, int32_t *argnp, int32_t argc, char **argv, char *usage);
void get_next_arg_double(double *varp, int32_t *argnp, int32_t argc, char **argv, char *usage);
void arg_error(char *msg, char *arg, char *pname, char *usage);
void plot_map(epswr_figure_t *eps, Map *m, options_t *o);

/* IMPLEMENTATIONS */

int32_t main(int32_t argc, char **argv)
  { 
    options_t *o = parse_options(argc, argv);

    Map *m = st_map_read(stdin);
    
    double mrg = 4.0; /* Margin width (pt). */
    bool_t eps_verbose = FALSE;
    epswr_figure_t *eps = epswr_new_named_figure
      ( NULL, NULL, NULL, -1, NULL,
        o->xSize*epswr_pt_per_mm, o->ySize*epswr_pt_per_mm, 
        mrg, mrg, mrg, mrg, 
        eps_verbose
      );
      
    plot_map(eps, m, o);
    epswr_end_figure(eps);
    return 0;
  }
  
void plot_map(epswr_figure_t *eps, Map *m, options_t *o)
  /* 
    Plots the map {m} on the Postscript file {f}.
    If the X and Y ranges in {o} are not empty, clips the 
    map to those ranges. */
  {
    Interval xPlot, yPlot;      /* Coordinate ranges of plot area */
    Interval xClip, yClip;      /* Clipping ranges */
    int32_t vntups[m->nv];      /* Number of TUPs installed at each vertex. */
    int32_t vcover[m->nv];      /* Vertex coverage by TUP street-cost balls. */
    int32_t ecover[m->ne];      /* Edge coverage by TUP street-cost balls. */
    float vwidth[m->nv];        /* Vertex diameters (mm). */
    float ewidth[m->ne];        /* Edge diameters (mm). */
    RGBColor vcolor[m->nv];     /* Vertex colors. */
    RGBColor ecolor[m->ne];     /* Edge colors. */
    int32_t vi, ei, i; 

    if ((o->xClip.lo < o->xClip.hi) && (o->yClip.lo < o->yClip.hi))
      { double tol;       /* Tolerance for vertex and edge clipping */
        /* Plot area is exactly what the user asked: */
        xPlot = o->xClip; yPlot = o->yClip;
        /* Clipping rectangle is a bit wider, just in case: */
        xClip = o->xClip; yClip = o->yClip;
        tol = 0.01 * ((xPlot.hi - xPlot.lo) + (yPlot.hi - yPlot.lo));
        st_interval_widen(&xClip, tol);
        st_interval_widen(&yClip, tol);
      }
    else
      { /* Plot area is the map's bounding box, with some extra room: */
        st_map_get_bbox(m, &xPlot, &yPlot);
        st_interval_widen(&xPlot, 120.0);
        st_interval_widen(&yPlot, 120.0);
        st_adjust_rect_shape(&xPlot, &yPlot, o->xSize, o->ySize);
        /* No clipping: */
        xClip.lo = INFINITY; xClip.hi = -INFINITY; /* Empty interval. */
        yClip.lo = INFINITY; yClip.hi = -INFINITY; /* Empty interval. */
      }

    /* Set client window: */
    epswr_set_client_window(eps, xPlot.lo, xPlot.hi, yPlot.lo, yPlot.hi);
    
    /* Plot coordinate grid: */
    if (o->grid > 0.0)
      { epswr_set_pen(eps, 0.0,1.0,1.0,  0.10, 0.0, 0.0);
        epswr_coord_lines(eps, epswr_axis_HOR, 0.0, o->grid);
        epswr_coord_lines(eps, epswr_axis_VER, 0.0, o->grid);
      }
      
    /* Draw rectangles: */
    epswr_set_pen(eps, 0.0,0.5,0.0, 0.25, 0.0, 0.0);
    for (i = 0; i < o->nRect; i++)
      { Interval *xr = &(o->xRect[i]);
        Interval *yr = &(o->yRect[i]);
        epswr_rectangle(eps, xr->lo,xr->hi, yr->lo,yr->hi, FALSE, TRUE);
      }

    /* Draw circles: */
    epswr_set_pen(eps, 0.0,0.25,0.0, 0.10, 0.0, 0.0);
    for (i = 0; i < o->nCircle; i++)
      { double x = o->ctrCircle[i].c[0];
        double y = o->ctrCircle[i].c[1];
        epswr_circle(eps, x, y, o->rCircle[i], FALSE, TRUE);
        epswr_segment(eps, x - 15.0, y, x + 15, y);
        epswr_segment(eps, x, y - 15.0, x, y + 15);
      }

    for (vi = 0; vi < m->nv; vi++) { vntups[vi] = 0; }
    int32_t uBall[o->nBall];
    float rBall[o->nBall];
    for (i = 0; i < o->nBall; i++)
      { /* Find vertex nearest to specified center: */
        int32_t ui = st_map_nearest_vertex(m, o->ctrBall[i]);
        uBall[i] = ui;
        rBall[i] = (float)o->rBall[i]; 
        vntups[ui]++;
      }
    st_compute_coverage(m, uBall, rBall, o->nBall, vcover, ecover);  

    /* Colors to use for map elems covered by 0, 1, or 2+ TUPs: */
    RGBColor uncColor = (RGBColor){0.000f, 0.000f, 1.000f}; /* Y=0.10 */
    RGBColor sngColor = (RGBColor){0.750f, 0.000f, 0.600f}; /* Y=0.26 */ 
    RGBColor mulColor = (RGBColor){1.000f, 0.670f, 0.000f}; /* Y=0.70 */
    
    /* Color for TUP-containing vertices: */
    RGBColor tupColor = (RGBColor){0.500f, 0.000f, 0.300f};  /* Y=0.18 */

    /* Color for significant intersections (degree != 2): */
    RGBColor vtxColor = (RGBColor){0.00, 0.00, 0.00};  /* Y=0.00 */

    /* Define vertex styles according to TUP count, degree, and coverage: */
    for (vi = 0; vi < m->nv; vi++) 
      { VertexData *vd = m->vd[vi];
        float *vwd = &(vwidth[vi]); 
        RGBColor *vcl = &(vcolor[vi]);
        /* Vertex color depends on coverage: */
        if (vntups[vi] > 0) 
          { *vcl = tupColor; }
        else if (vcover[vi] == 0)
          { *vcl = vtxColor; /* Irrelevant if vertex has degree 2 */ }
        else if (vcover[vi] == 1)
          { *vcl = sngColor; }
        else 
          { *vcl = mulColor; }
        
        /* TUP locations are big, trivial vertices are invisible: */
        if (vntups[vi] > 0) 
          { *vwd = 1.440f; }
        else if (vd->deg == 2) 
          { *vwd = 0.000f; }
        else
          { *vwd = 0.720f; }
      }

    /* Define edge styles according to coverage: */
    for (ei = 0; ei < m->ne; ei++)
      { float *ewd = &(ewidth[ei]); 
        RGBColor *ecl = &(ecolor[ei]);
        *ewd = (ecover[ei] == 0 ? 0.240f : 0.480f);
        if (ecover[ei] == 0)
          { *ecl = uncColor; }
        else if (ecover[ei] == 1)
          { *ecl = sngColor; }
        else
          { *ecl = mulColor; }
      }

    /* Plot map proper: */
    st_map_plot(eps, m, xClip, yClip, vwidth, vcolor, ewidth, ecolor);
  }

/* OPTION PARSING */
  
#define plotmap_MIN_GRID (2.0)
  /* Min grid line spacing. */

#define plotmap_MAX_MAP_SIZE (90000.0)
  /* Max map coordinates/size (metres) */

#define plotmap_MIN_PLOT_SIZE (30.0)
#define plotmap_MAX_PLOT_SIZE (3000.0)
  /* Range of valid values for the "-size" option (mm). */

options_t *parse_options(int32_t argc, char **argv)
   /*
     Parses the command line options, returns a record with their options. */
  {
   
    options_t *o = (options_t *)malloc(sizeof(options_t));

    /* Defaults: */
    o->xSize = 140.0; /* 14 cm square. */
    o->ySize = 140.0; /* 14 cm square. */
    o->grid = 0.0; /* No grid */
    
    o->xClip.lo = INFINITY; o->xClip.hi = -INFINITY; /* No clipping. */
    o->yClip.lo = INFINITY; o->yClip.hi = -INFINITY; /* No clipping. */
    
    o->nRect = 0;  /* No rectangles */
    o->xRect = (Interval *)malloc(pltmap_MAXMARKS * sizeof(Interval));
    o->yRect = (Interval *)malloc(pltmap_MAXMARKS * sizeof(Interval));
    
    o->nCircle = 0; /* No Euclidean circles */
    o->ctrCircle = (Point *)malloc(pltmap_MAXMARKS * sizeof(Point));
    o->rCircle = rn_alloc(pltmap_MAXMARKS);
    
    o->nBall = 0;   /* No street-cost balls. */
    o->ctrBall = (Point *)malloc(pltmap_MAXMARKS * sizeof(Point));
    o->rBall = rn_alloc(pltmap_MAXMARKS);

    /* Initialize argument parser: */
    argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, PROG_NAME " version " PROG_VERS ", usage:\n" PROG_HELP);
    argparser_set_info(pp, PROG_INFO);
    argparser_process_help_info_options(pp);

    /* Parse the options: */
    
    if (argparser_keyword_present(pp, "-size"))
      { o->xSize = argparser_get_next_double(pp, plotmap_MIN_PLOT_SIZE, plotmap_MAX_PLOT_SIZE);
        o->ySize = argparser_get_next_double(pp, plotmap_MIN_PLOT_SIZE, plotmap_MAX_PLOT_SIZE);
      }
    
    if (argparser_keyword_present(pp, "-grid"))
      { o->grid = argparser_get_next_double(pp, plotmap_MIN_GRID, plotmap_MAX_MAP_SIZE);
      }
    
    if (argparser_keyword_present(pp, "-clip"))
      { o->xClip.lo = argparser_get_next_double(pp, -INF, +INF);
        o->xClip.hi = argparser_get_next_double(pp, o->xClip.lo + 0.1, +INF);
        o->yClip.lo = argparser_get_next_double(pp, -INF, +INF);
        o->yClip.hi = argparser_get_next_double(pp, o->yClip.lo + 0.1, +INF);
      }
    
    while (argparser_keyword_present(pp, "-showRect"))
      { int32_t i = o->nRect;
        if (o->nRect >= pltmap_MAXMARKS) { argparser_error(pp, "too many rectangles"); }
        o->xRect[i].lo = argparser_get_next_double(pp, -plotmap_MAX_MAP_SIZE, +plotmap_MAX_MAP_SIZE);            
        o->xRect[i].hi = argparser_get_next_double(pp, o->xRect[i].lo + 0.1, plotmap_MAX_MAP_SIZE);
        o->yRect[i].lo = argparser_get_next_double(pp, -plotmap_MAX_MAP_SIZE, +plotmap_MAX_MAP_SIZE);            
        o->yRect[i].hi = argparser_get_next_double(pp, o->yRect[i].lo + 0.1, plotmap_MAX_MAP_SIZE);
        o->nRect++;
      }
    
    while (argparser_keyword_present(pp, "-showCircle"))
      { int32_t i = o->nCircle;
        if (o->nCircle >= pltmap_MAXMARKS) { argparser_error(pp, "too many circles"); }
        o->ctrCircle[i].c[0] = argparser_get_next_double(pp, -plotmap_MAX_MAP_SIZE, +plotmap_MAX_MAP_SIZE);
        o->ctrCircle[i].c[1] = argparser_get_next_double(pp, -plotmap_MAX_MAP_SIZE, +plotmap_MAX_MAP_SIZE);
        o->rCircle[i] = argparser_get_next_double(pp, 0.001, plotmap_MAX_MAP_SIZE);
        o->nCircle++;
      }
    
    while (argparser_keyword_present(pp, "-showBall"))
      { int32_t i = o->nBall;
        if (o->nBall >= pltmap_MAXMARKS) { argparser_error(pp, "too many balls"); }
        o->ctrBall[i].c[0] = argparser_get_next_double(pp, -plotmap_MAX_MAP_SIZE, +plotmap_MAX_MAP_SIZE);
        o->ctrBall[i].c[1] = argparser_get_next_double(pp, -plotmap_MAX_MAP_SIZE, +plotmap_MAX_MAP_SIZE);
        o->rBall[i] = argparser_get_next_double(pp, 0.001, plotmap_MAX_MAP_SIZE);
        o->nBall++;
      }

    argparser_finish(pp);

    return o;
  }
