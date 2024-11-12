/* SO2DPlot.h -- General procedures for SOFunction plotting over a 2-Dimensional space */
/* Last edited on 2007-01-04 00:20:00 by stolfi */

#ifndef SO2DPlot_H
#define SO2DPlot_H

#include <SOGrid.h>
#include <SOFunction.h>
#include <SOBasic.h>
#include <SOPlot.h>

#include <dg_tree.h>
#include <dg_grid.h>

#include <vec.h>
#include <stdio.h>

void SO2DPlot_Values
  ( SOPlot_File *psf,     /* Plot file. */
    SOFunction *f,        /* Function to be plotted. */
    SOGrid_Tree *tree,    /* Reference grid for adaptive domain subdivision. */ 
    interval_t boxR[2],  /* Rectangle relative to root cell. */
    double fPlotMin,      /* Minimum value for color mapping. */
    double fPlotMax,      /* Maximum value for color mapping. */
    double fStep,         /* Value increment between isolines/bands. */
    double lineWidth,     /* Line width for isoline plotting. */
    int minDepth,         /* Minimum subdivision depth for plotting. */
    int extraDepth,       /* Extra depth of subdivision for leaves of {tree}. */
    int maxDepth,         /* Maximum subdivision depth for plotting. */
    bool_t bands,           /* TRUE paints between isolines, FALSE leaves blank. */
    bool_t isolines,        /* TRUE draws the isolines, FALSE omits them. */ 
    double *fObsMin,      /* (IN/OUT) Min value seen during plot. */ 
    double *fObsMax       /* (IN/OUT) Max value seen during plot. */ 
  );
  /* Plots values of a general SOFunction {f} on file {psf}, within
    the rectangle {boxR}, in root-cell-relative coordinates. The plot consists of
    isolines and/or discrete color steps at certain values of the
    function. The levels are spaced {fStep} apart, shifted so as to
    include the levels {±fStep/2}.
    
    Positive isolines are drawn in red, negative ones in blue. For
    color bands, zero is mapped to white, values {fPlotMax} (if positive)
    and higher are mapped to red, and values {fPlotMin} (if negative) and
    lower are mapped to blue.
   
    The curves are obtained by subdividing each cell recursively, into
    subcells, and then subdividing the leaf subcells into triangles.
    Within each triangle, the function is evaluated at the corners and
    interpolated linearly in the interior.

    The recursive subdivision of the plotting area into plot cells is
    governed by the given cell {tree}, modified by the parameters
    {extraDepth} and {maxDepth}. Specifically, any cell of {tree}
    below level {maxDepth} is ignored, and any leaf cell {C} at depth
    {k < maxDpth} is subdivided until level {min(k + extraDepth,
    maxDepth)}. A null {tree} is interpreted as the trivial tree with
    a single node (the root cell), which is therefore uniformly
    subdivided to depth {min(extraDepth,maxDepth)}.
   
    Note that this method may produce discontinuities along
    the border between two leaf cells of {tree} which have
    different depths.  Increasing the {extraDepth} parameter
    should eventually render those discontinuities invisible.
    
    The variables {fObsMin} and {fObsMax} should be initialized by 
    the caller; the procedure will update them with the 
    maximum and minimum values actually seen during the plot. */

void SO2DPlot_Tree
  ( SOPlot_File *psf,     /* Plot file. */
    SOGrid_Tree *tree,    /* The grid to plot. */
    interval_t boxR[2],  /* Low corner of region to plot. */
    double lineWidth,     /* Line width for grid cell boundaries. */
    int maxDepth          /* Omit cells below this depth. */
  );
  /* Draws the outlines of the cells of the given {tree} to file {psf}, within
    the rectangle {boxR}, in root-cell-relative coordinates. 
    Cells at levels greated than {maxDepth} are omitted. */

SOPlot_File *SO2DPlot_FunctionFigure
  ( SOPlot_File *fps,          /* Plot file (or NULL). */
    SOFunction *f,             /* Function to plot. */
    SOGrid_Tree *tree,         /* Optional grid to plot. */
    interval_t boxR[2],       /* Low corner of region to plot. */
    double fPlotMin,           /* Nominal minimum {f} value, for color scale. */
    double fPlotMax,           /* Nominal maximum {f} value, for color scale. */
    double fStep,              /* Value step for isolines and color bands. */
    double isolineWidth,       /* Line width for isoline plotting. */
    double gridWidth,          /* Line width for grid drawing. */
    int meshDepth,             /* Depth of bisection recursion for plotting. */
    bool_t bands,                /* TRUE plots color bands. */ 
    bool_t isolines,             /* TRUE draws isolines. */ 
    bool_t grid,                 /* TRUE plots the tree */
    SOPlot_PageState *pgs,     /* Page layout and state for {fps}. */
    char *docName,             /* Document name (minus extension). */
    char *figName,             /* Figure name (minus extension). */
    double *fObsMin,           /* (IN/OUT) Min value seen during plot. */ 
    double *fObsMax            /* (IN/OUT) Max value seen during plot. */ 
  );
  /* Starts a new figure (with {SOPlot_BeginFigure}) in the figure stream 
    defined by the file {psFile}, the page state {pgs}, the document
    name {docName}, and the figure name {figName}; then plots the function
    {f} in that figure with {SO2DPlot_Values} and {SO2DPlot_Tree}. */

/* FUNCTION RANGE ADJUSTMENT AND ESTIMATION */

double SO2DPlot_EstRange(SOFunction *f, int verbose);
 /* Returns an estimate of the maximum of {fabs(f(p))} over the root cell. */

void SO2DPlot_FixRange
  ( double *fRange, 
    double *fStep,
    double alpha, 
    SOFunction *f, 
    int verbose
  );
  /* Replaces {*fRange} by a suitable mixture of the given {*fRange} and
    the estimated maximum absolute value of {f(p)}, namely 
    {estRange = SO2DPlot_EstRange(f)}.
    
    The parameter {alpha} specifies the degree of mixing: {alpha = 0}
    leaves {*fRange} unchanged, {alpha = 1} replaces it by {estRange}.
    However, if the given {*fRange} is not valid (i.e., {fRange <= 0}),
    then replaces it by {estRange}, as if {alpha} was 1.
    
    If {alpha} is not zero, the computed range {*fRange} is rounded to a
    nice value.
    
    In any case (even when {alpha = 0}), checks whether the final 
    {*fStep} is reasonable for the final {*fRange}; if not, fixes 
    it appropriately. */
  
double SO2DPlot_RoundToNice(double x);
  /* Rounds the number {x}, which must be positive, to a 
    nice value (namely 0.25, 0.50, or 1.00 times a power of 10). */

#endif
