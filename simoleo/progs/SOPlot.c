/* See SOPlot.h */
/* Last edited on 2008-07-14 22:22:35 by stolfi */

#include <SOPlot.h>

#include <SOGrid.h>
#include <SOBasic.h>
#include <vec.h>
#include <pswr.h>
#include <affirm.h>
#include <nat.h>
#include <math.h>
#include <values.h>
#include <stdio.h>

#define TwoPi (2.0 * M_PI)

#define Invisible ((Color){{-1,-1,-1}})

/* FILE AND PAGE CONTROL */

SOPlot_PageState *SOPlot_PageState_New
  ( bool_t eps, 
    char *paperSize,                       
    double hFigSize,    /* Figure width (in mm). */
    double vFigSize,    /* Figure height (in mm). */
    double hFigMargin,  /* Left/right margin for each figure (in mm). */
    double vFigMargin,  /* Top/bottom margin for each figure (in mm). */
    int captionLines,   /* Number of caption lines below each figure. */
    int hCount,         /* Number of figures in each row. */   
    int vCount          /* Number of figures in each column. */
  )
  { double mm = 72.0/25.4; 
    SOPlot_PageState *pgs = 
      (SOPlot_PageState *)notnull(malloc(sizeof(SOPlot_PageState)), "out of mem");
    double min_margin = 72.0; /* Min page margin width in pt. */
    double vCapSize = 10.0*((double)captionLines); /* Caption height (in pt). */
    
    pgs->eps = eps;
    
    /* Save figure parameters (converted to pt): */
    pgs->hFigSize = hFigSize*mm;
    pgs->vFigSize = vFigSize*mm;
    pgs->hFigMargin = hFigMargin*mm;
    pgs->vFigMargin = vFigMargin*mm;
    pgs->captionLines = captionLines;

    /* CCompute slot dimensions (in pt): */
    pgs->hSlotSize = 2*pgs->hFigMargin + pgs->hFigSize;
    pgs->vSlotSize = 2*pgs->vFigMargin + pgs->vFigSize + vCapSize;

    if (eps)
      { /* Provide default row/column counts, if not given: */
        if (hCount <= 0) { hCount = 1; }
        if (vCount <= 0) { vCount = 1; }
        /* Ignore {paperSize}, compute virtual page dimensions: */
        pgs->paperSize = NULL;
        pgs->hPageSize = hCount*pgs->hSlotSize;
        pgs->vPageSize = vCount*pgs->vSlotSize;
        /* Compute the top left corner of top left slot: */
        pgs->hFirst = 0;
        pgs->vFirst = vCount*pgs->vSlotSize;
      }
    else
      { /* Get actual page dimensions: */
        pgs->paperSize = paperSize;
        ps_get_paper_dimensions(paperSize, FALSE, &(pgs->hPageSize), &(pgs->vPageSize));
        /* Compute row/column counts from page size, if not given: */
        if (hCount <= 0) 
          { hCount = (pgs->hPageSize - 2*min_margin) / pgs->hSlotSize; }
        if (vCount <= 0)
          { vCount = (pgs->vPageSize - 2*min_margin) / pgs->vSlotSize; }
        /* Just in case the computed counts are negative: */
        if (hCount <= 0) { hCount = 1; }
        if (vCount <= 0) { vCount = 1; }
        /* Compute top left corner of top left slot: */
        pgs->hFirst = (pgs->hPageSize - hCount*pgs->hSlotSize)/2;
        pgs->vFirst = (pgs->vPageSize + vCount*pgs->vSlotSize)/2;
      }
    pgs->hCount = hCount;
    pgs->vCount = vCount;
    /* Prepare for first figure: */
    pgs->curPage = -1; /* Indicates that the file hasn't been opened yet */
    pgs->hIndex = 0;
    pgs->vIndex = pgs->vCount;
    return pgs;
  }

SOPlot_File *SOPlot_BeginFigure
  ( SOPlot_File *fps,
    SOPlot_PageState *pgs, 
    char *docName,
    char *figName, 
    double xmin, double xmax, 
    double ymin, double ymax
  )
  { /* If a file was provided, advance to the next slot in same page: */
    if (fps != NULL) 
      { /* Assumes that {SOPlot_BeginFigure} was called at least 
          once on this figure. Adjusts {pgs} to point to the next 
          figure slot on the same page. If there is no such
          slot (the page is already full), leaves {vIndex >= vCount}. */
        affirm(pgs->curPage > 0, "page state inconsistent");
        pgs->hIndex++;
        if (pgs->hIndex >= pgs->hCount) 
          { pgs->hIndex = 0; pgs->vIndex++; }
      }
    /* Make sure that the slot exists: */
    if (pgs->eps)
      { if ((fps == NULL) || (pgs->curPage <= 0) || (pgs->vIndex >= pgs->vCount))
          { /* No file, or uninitialized file, or virtual page exausted. */
            if ((fps != NULL) && (pgs->curPage > 0))
              { /* Open file with stuff in it: */
                ps_end_figure(fps); fclose(fps); fps = NULL;
              }
            if (fps == NULL) { fps = open_write(txtcat(figName, ".eps"), TRUE); }
            /* Now {fps} is an open file with no stuff in it: */
            ps_begin_figure(fps, 0.0, pgs->hPageSize, 0.0, pgs->vPageSize);
            pgs->curPage = 1;
            pgs->hIndex = 0;
            pgs->vIndex = 0;
          }
        /* Now {fps} is an open file, with {ps_begin_figure} called: */
        ps_set_pen(fps, 0.0,0.0,0.0,  0.15,  0.0,0.0);
      }
    else
      { if ((fps == NULL) || (pgs->curPage <= 0))
          { /* No file, or un-initialized file; start new document: */
            if (fps == NULL) { fps = open_write(txtcat(docName, ".ps"), TRUE); }
            ps_begin_document(fps, pgs->paperSize);
            ps_begin_page(fps, NULL);
            pgs->curPage = 1;
            pgs->hIndex = 0;
            pgs->vIndex = 0;
          }
        else if (pgs->vIndex >= pgs->vCount)
          { /* Page exhausted, start new page: */
            ps_end_page(fps);
            ps_begin_page(fps, NULL);
            pgs->curPage++;
            pgs->hIndex = 0;
            pgs->vIndex = 0;
          }
        ps_set_pen(fps, 0.0,0.0,0.0,  0.15,  0.0,0.0);
      }
    /* OK, now set scales: */
    { /* Get actual plot window {[hmin_hmax]\x[vmin_vmax]} for this slot: */
      double hFigSize = pgs->hFigSize;
      double vFigSize = pgs->vFigSize;
      double hmin = pgs->hFirst + pgs->hIndex*pgs->hSlotSize + pgs->hFigMargin;
      double vmax = pgs->vFirst - pgs->vIndex*pgs->vSlotSize - pgs->vFigMargin;
      double hmax = hmin + hFigSize;
      double vmin = vmax - vFigSize;
      /* Scale so as to fit the client rectangle exactly: */
      double hscale = hFigSize/(xmax - xmin);
      double vscale = vFigSize/(ymax - ymin);
      double scale = (hscale < vscale ? hscale : vscale);
      double rx = hFigSize/2/scale;
      double ry = vFigSize/2/scale;
      double xctr = (xmax + xmin)/2;
      double yctr = (ymax + ymin)/2;
      ps_set_window
        ( fps, 
          xctr-rx, xctr+rx, yctr-ry, yctr+ry, 
          hmin,    hmax,    vmin,    vmax,
          1,1
        );
      fprintf(stderr, "plot ranges: \n");
      fprintf
        ( stderr, "  client [%8.4f _ %8.4f] x [%8.4f _ %8.4f]\n",
          xctr-rx, xctr+rx, yctr-ry, yctr+ry
        );
      fprintf
        ( stderr, "  device [%8.4f _ %8.4f] x [%8.4f _ %8.4f]\n",
          hmin,    hmax,    vmin,    vmax
        );
    }     
    return fps;
  }

void SOPlot_Close(SOPlot_File *fps, SOPlot_PageState *pgs)
  { if (fps == NULL) { return; }
    if (pgs->curPage <= 0)
      { fprintf(stderr, "** empty Postscript file?\n"); }
    else if (pgs->eps)
      { ps_end_figure(fps); }
    else
      { ps_end_page(fps); ps_end_document(fps); }
    fclose(fps);
    pgs->curPage = -1;
  }
  
double SOPlot_DefaultFigSize
  ( bool_t eps, 
    char *paperSize, 
    int nRows, 
    int nCols,
    int captionLines
  )
  { if (nRows <= 0) { nRows = 1; }
    if (nCols <= 0) { nCols = 1; }
    if (eps)
      { return 150.0 / (nRows < nCols ? nRows : nCols); }
    else
      { double hsize, vsize;
        ps_get_paper_dimensions(paperSize, FALSE, &hsize, &vsize);
        hsize = (hsize/72.0 - 2.0)*25.0/nCols; 
        vsize = ((vsize - 10.0*nRows*captionLines)/72.0 - 2.0)*25.0/nRows;
        if (hsize < 0.0) { hsize = 10.0; }
        if (vsize < 0.0) { vsize = 10.0; }
        return (hsize < vsize ? hsize : vsize);
      }
  }

double SOPlot_Brightness(Color *C)
  { 
    return 0.299 * C->c[0] + 0.587 * C->c[1] + 0.114 * C->c[2];
  }

Color SOPlot_ColorScale(double s, Color *A, Color *B)
  { double yA = SOPlot_Brightness(A) + 0.02;
    double yB = SOPlot_Brightness(B) + 0.02;
    double u, v;
    if (fabs(yA - yB) < 0.00001)
      { /* Same brightness, use linear interpolation: */
        v = s;
      }
    else
      { /* Interpolate brigtnesses in log scale: */
        double y = yA*exp(s*log(yB/yA));
        v = (y - yA)/(yB - yA);
        if (v < 0.0) { v = 0.0; }
        if (v > 1.0) { v = 1.0; }
      }
    u = 1.0 - v;
    return (Color)
      {{
        u*A->c[0] + v*B->c[0],
        u*A->c[1] + v*B->c[1],
        u*A->c[2] + v*B->c[2]
      }};
  }

Color SOPlot_InterpolateColor
  ( double f,         /* Function value */
    double fPlotMin,  /* Minimum function value. */
    Color *CMin,      /* Color to use for {fPlotMin}. */
    bool_t clipMin,     /* TRUE maps values below {fPlotMin} to invisible. */
    double fPlotMax,  /* Maximum function VALUE. */
    Color *CMax,      /* Color to use for {fPlotMax}. */
    bool_t clipMax      /* TRUE maps values above {fPlotMax} to invisible. */
  )
  { if (f < fPlotMin)
      { if (clipMin) { return Invisible; } else { return *CMin; } }
    else if (f > fPlotMax)
      { if (clipMax) { return Invisible; } else { return *CMax; } }
    else
      { double s = (f - fPlotMin)/(fPlotMax - fPlotMin);
        double t = 1.0 - s;
        return (Color)
          {{
            s*CMax->c[0] + t*CMin->c[0],
            s*CMax->c[1] + t*CMin->c[1],
            s*CMax->c[2] + t*CMin->c[2]
          }};
      }
  }
  
Color SOPlot_ClipColor(Color *C)
  { 
    double m = 0; int i;
    for (i = 0; i < 3; i++) 
      { double ci = fabs(C->c[i]); if (ci > m) { m = ci; } }
    if (m <= 1.0)
      { return *C; } 
    else
      { Color A = (Color){{C->c[0]/m, C->c[1]/m, C->c[2]/m}};
        double yC = SOPlot_Brightness(C);
        double yA = yC/m;
        double s = (yC < 1.0 ? yC : 1.0) - yA;
        return (Color)
          {{
            A.c[0] + s*(1.0 - A.c[0]),
            A.c[1] + s*(1.0 - A.c[1]),
            A.c[2] + s*(1.0 - A.c[2])
          }};
      }
  }
