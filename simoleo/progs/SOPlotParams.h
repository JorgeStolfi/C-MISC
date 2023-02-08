/* SOPlotParams.h -- plotting options from command line. */
/* Last edited on 2007-01-04 00:22:54 by stolfi */

#ifndef SOPlotParams_H
#define SOPlotParams_H

#include <SOParams.h>
#include <SOBasic.h>
#include <values.h>

typedef struct PlotOptions
  { /* General plotting options: */
    bool_t eps;             /* TRUE to generate ".eps" instead of ".ps" */
    char *paperSize;      /* Document paper size (when {eps == FALSE}). */
    double figSize;       /* Figure size in mm. */
    double meshSize;      /* Plotting resolution in mm. */
    string_vec_t caption; /* Figure caption. */
    /* Parameters for grid drawing: */
    bool_t noGrid;          /* TRUE to omit grid. */
    double gridWidth;     /* Line width for grid lines, in mm. */
    /* Parameters for value-color mapping and isolines: */
    bool_t noIsolines;      /* TRUE to omit isolines. */
    bool_t noBands;         /* TRUE to omit color bands between isolines. */
    double fRange;        /* Nominal function range is {[-fRange _ fRange]}. */
    double autoRange;     /* Degree of automatic adjustment for {fRange}. */
    double fStep;         /* Value difference between successive levels */
    double isolineWidth;  /* Line width for isolines, in millimeters */
  } PlotOptions;  

#define SOPlotParams_GridHelp \
  "  [ -eps | -ps ] [ -paperSize STRING ] \\\n" \
  "  [ -figSize NUM ] \\\n" \
  "  [ -noGrid ] [ -meshSize NUM ] \\\n" \
  "  [ -caption TEXT ]... \\\n" \
  "  [ -gridWidth MILLIMETERS ]"

PlotOptions SOPlotParams_GridParse(SOParams_T *pp);
  /* Parses from the {pp} object a set of {SOPlotOptions}
    suitable for drawing grids only (without functions
    on them), as described by {SOPlotParams_GridHelp}. */

#define SOPlotParams_FunctionHelp \
  "  [ -eps | -ps ] [ -paperSize STRING ] \\\n" \
  "  [ -caption TEXT ]... \\\n" \
  "  [ -noIsolines ] [ -noBands ] \\\n" \
  "  [ -fRange NUM ] [ -fStep NUM ] [ -autoRange NUM ] \\\n" \
  "  [ -figSize NUM ] [ -meshSize NUM ] \\\n" \
  "  [ -isolineWidth MILLIMETERS ]"
    
PlotOptions SOPlotParams_FunctionParse(SOParams_T *pp);
  /* Parses from the {pp} object a set of {SOPlotOptions}
    suitable for drawing   functions, as 
    described by {SOPlotParams_FunctionHelp}. */

#endif

