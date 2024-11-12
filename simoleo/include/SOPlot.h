/* SOPlot.h -- Basic plot routines, for arbitrary figures. */
/* Last edited on 2007-12-23 22:46:29 by stolfi */

#ifndef SOPlot_H
#define SOPlot_H

#include <SOBasic.h>

/* FILE AND PAGE CONTROL */

typedef FILE SOPlot_File;  /* A file containing Postscript drawings. */

typedef struct SOPlot_PageState  /* State of a PS file, for SOPlot_Open. */
  { /* Pagination parameters: */
    bool_t eps;           /* TRUE: encapsulated figure, FALSE: multipage doc. */
    char *paperSize;    /* Paper size ({letter}, {a3}, etc.) or NULL. */
    double hPageSize;   /* Page width (in pt). */
    double vPageSize;   /* Page height (in pt). */
    double hFigSize;    /* Figure width (in pt). */
    double vFigSize;    /* Figure height (in pt). */
    double hFigMargin;  /* Left/right margin for each figure (in pt). */
    double vFigMargin;  /* Top/bottom margin for each figure (in pt). */
    int captionLines;   /* Number of caption lines (including fig name). */
    int hCount;         /* Number of figures in each row. */
    int vCount;         /* Number of figures in each column. */
    /* Derived parameters: */
    double hSlotSize;   /* Total fig slot width, incl. margins. */
    double vSlotSize;   /* Total fig slot height, incl. margins and caption. */
    double hFirst;      /* Left X coord of first slot. */
    double vFirst;      /* Top Y coord of first slot. */
    /* Pagination state: */
    int curPage;        /* Current page number. */
    int hIndex;         /* Column index of next available figure slot. */
    int vIndex;         /* Row index of next available figure slot. */
  } SOPlot_PageState;
  /* An {SOPlot_PageState} keeps track of the pagination parameters
    and pagination state for a Postscript file which is meant to
    contain several figures. 
    
    The fields {hPageSize,vPageSize} give the total page size.
    Individual figures are arranged as a grid of {hCount} by {vCount}
    rectangular slots. Each slot measures {hFigSize} by {vFigSize},
    plus a blank margin specified by {hFigMargin}, {vFigMargin}, plus
    enough space below the figure for {captionLines} lines of caption
    text.
    
    After each call to {SOPlot_BeginFigure}, the current figure slot
    is the one with indices {hIndex,vIndex} (starting from 0, at the
    top left), in the page {curPage} (starting from 1). */

SOPlot_PageState *SOPlot_PageState_New
  ( bool_t eps, 
    char *paperSize,    /* Paper size ({letter}, {a3}, etc.). */                   
    double hFigSize,    /* Figure width (in mm). */
    double vFigSize,    /* Figure height (in mm). */
    double hFigMargin,  /* Left/right margin for each figure (in mm). */
    double vFigMargin,  /* Top/bottom margin for each figure (in mm). */
    int captionLines,   /* Number of caption lines below each figure. */
    int vCount,         /* Number of figures in each row. */   
    int hCount          /* Number of figures in each column. */
  );
  /* Allocates a new {SOPlot_PageState} record and initializes it 
    appropriately.  
    
    If {eps} is false, the page dimension fields are derived from the
    {paperSize} string (see {pswr.h}). The array of figures is
    horizontally centered on the page, with a larger margin at bottom
    than at top. If either of {vCount} or {hCount} is zero, it
    defaults to the maximum number that will fit on the page, 
    after leaving 1 inch of margin along each side.
    
    If {eps} is TRUE, the {paperSize} string is ignored, and the
    (virtual) page size is set to be {hCount*(hFigSize + 2*hFigMargin)}
    by {vCount*(vFigSize + 2*vFigMargin + captionLines*(10pt))}.
    If {vCount} and/or {hCount} is zero, they default to 1. */

SOPlot_File *SOPlot_BeginFigure
  ( SOPlot_File *fps,          /* Postscript file, or NULL. */
    SOPlot_PageState *pgs,     /* Page layout and state. */
    char *docName,             /* Document name (minus extension). */
    char *figName,             /* Figure name (minus extension). */
    double xmin, double xmax,  /* Client X plotting range. */
    double ymin, double ymax   /* Client Y plotting range. */
  );
  /* Prepares a Postscript file for plotting another figure of a
      function.  Be sure to call {SOPlot_EndFigure}
    at the end of it.
    
    The actual effect depends on the type of Postscript file, as
    defined by {pgs->eps}. For multipage, non-encapsulated postscript
    documents ({pgs->eps = FALSE}):
    
      If {fps} is NULL, a new output file is opened with name {docName}
      plus the extension ".ps", and it is initialized with
      {ps_begin_document} and {ps_begin_page}. The plotting window is
      set to the first (top left) slot on the first page.
      
      If {fps} is not NULL, the plotting window is set to the next available
      figure slot, as indicated by the page state {pgs}. A new page 
      is started (with {ps_end_page} and {ps_begin_page}) if necessary.
      
    For encapsulated figures ({pgs->eps = TRUE}): 
    
      If {fps} is NULL, a new file is opened with name {figName}
      and extension ".eps". (The {docName} is ignored).  The
      file is initialized with {ps_begin_figure}. The plotting window
      is set to the first (top left) slot on the virtual page.
      
      If {fps} is not NULL, the window is set to the next available
      figure slot, as indicated by the page state {pgs}. If all the 
      slots in the current (virtual) page are full, the file {fps} is 
      finalized (with {ps_end_figure}) and closed, and a new file is
      opened and initialized, as above. 
    
    Regardless of the {eps} flag, the plotting scale parameters are
    set up so that the rectangle  {[xmin _ xmax] \x [ymin _ ymax]} 
    in client coordinates will fit centered in the selected slot
    (excluding its margins), with equal scale factors on both axes.
    
    In any case, the resulting file is returned, open and ready for 
    plotting. */
    
void SOPlot_Close(SOPlot_File *fps, SOPlot_PageState *pgs);
  /* Terminates any plots written to {fps}, and closes it.
  
    Specifically, if {fps} is a multipage document, calls {ps_end_page}
    and {ps_end_document}; if {fps} is an encapsulated figure, calls
    {ps_end_figure}. */

double SOPlot_DefaultFigSize
  ( bool_t eps, 
    char *paperSize, 
    int nRows, 
    int nCols,
    int captionLines
  );
  /* A convenient default (square) figure size for the given output options.
    If {eps = FALSE}, chooses the size so as to fit the given number of rows
    and columns in the page, with 1 inch margins.  If {eps = TRUE}, ignores
    the {paperSize} and returns 150.0 mm divided by the *minimum* of 
    {nCols,nRows}. */

double SOPlot_Brightness(Color *C);
  /* The relative brightness (luminance) of color {C},
    namely {0.299 * R + 0.587 * G + 0.114 * B}. */

Color SOPlot_ColorScale(double s, Color *A, Color *B);
  /* Returns a color that lies {s} of the way from color {A} 
    to color {B}.  Assumes that the distance between two 
    colors is proportional to the difference in their brightnesses.  */
    
Color SOPlot_InterpolateColor
  ( double f,         /* Function value */
    double fPlotMin,  /* Minimum function value. */
    Color *CMin,      /* Color to use for {fPlotMin}. */
    bool_t clipMin,     /* TRUE maps values below {fPlotMin} to invisible. */
    double fPlotMax,  /* Maximum function value. */
    Color *CMax,      /* Color to use for {fPlotMax}. */
    bool_t clipMax      /* TRUE maps values above {fPlotMax} to invisible. */
  );
  /* Compute the color corresponding to {f}, by linear
    interpolation from {[fPlotMin _ fPlotMax]} to {[CMin _ CMax]}.
    The meaning of {clipMin} and {clipMax} is the same as
    in {SOPlot_PaintValues}. */

Color SOPlot_ClipColor(Color *C);
  /* Clips the color {C} to the unit cube, preserving 
    its brightness and hue. */

#endif
