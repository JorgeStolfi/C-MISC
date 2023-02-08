/* See SOPlotParms.h */
/* Last edited on 2007-01-04 00:22:51 by stolfi */

#include <SOPlotParams.h>
#include <SOParams.h>
#include <SOBasic.h>
#include <SOPlot.h>
#include <values.h>

/* DEFAULT OPTION VALUES */

/* Maximum cell/subcell size in mm: */
#define DEFAULT_MESH_SIZE 1.0

PlotOptions SOPlotParams_GridParse(SOParams_T *pp)
  {
    PlotOptions o;

    /* General plotting options: */
    
    if (SOParams_KeywordPresent(pp, "-eps"))
      { o.eps = TRUE; } 
    else if (SOParams_KeywordPresent(pp, "-ps"))
      { o.eps = FALSE; }
    else
      { o.eps = FALSE; }

    if (SOParams_KeywordPresent(pp, "-paperSize"))
      { o.paperSize = SOParams_GetNext(pp); }
    else
      { o.paperSize = "letter"; }  

    o.noGrid = SOParams_KeywordPresent(pp, "-noGrid");

    if (SOParams_KeywordPresent(pp, "-meshSize"))
      { o.meshSize = SOParams_GetNextDouble(pp, 0.0, 10000.0); }
    else
      { o.meshSize = DEFAULT_MESH_SIZE; }

    o.caption = string_vec_new(10);
    { int nCapLines = 0;
      while (SOParams_KeywordPresent(pp, "-caption"))
        { string_vec_expand(&(o.caption), nCapLines);
          o.caption.el[nCapLines] = SOParams_GetNext(pp);
          nCapLines++;
        }
      string_vec_trim(&(o.caption), nCapLines);
    }

    if (SOParams_KeywordPresent(pp, "-figSize"))
      { o.figSize = SOParams_GetNextDouble(pp, 5.0, 2000.0); }
    else
      { /* Guess: user caption, plus 2 lines of figure name if not {eps}: */
        int nCap = o.caption.nel + (o.eps ? 0 : 2);
        /* Default figure size -- fits in page, allowing caption: */
        o.figSize = SOPlot_DefaultFigSize(o.eps, o.paperSize, 1, 1, nCap);
      }
    
    if (SOParams_KeywordPresent(pp, "-gridWidth"))
      { o.gridWidth = SOParams_GetNextDouble(pp, 0.0, INFTY); }
    else
      { o.gridWidth = 0.15; }

    /* Parameters for value-color mapping and isolines: */

    o.fStep = 0.1;
    o.fRange = -1.0;  /* Invalid range, to force computation. */ 
    o.autoRange = 1.0;
    o.isolineWidth = 0.15;
    o.noBands = FALSE;
    o.noIsolines = FALSE;

    return o;
  }

PlotOptions SOPlotParams_FunctionParse(SOParams_T *pp)
  {
    PlotOptions o = SOPlotParams_GridParse(pp);

    /* Parameters for value-color mapping and isolines: */
    
    o.noIsolines = SOParams_KeywordPresent(pp, "-noIsolines");
    o.noBands = SOParams_KeywordPresent(pp, "-noBands");

    if (SOParams_KeywordPresent(pp, "-fRange"))
      { o.fRange = SOParams_GetNextDouble(pp, 0.0, INFTY); }
    else
      { o.fRange = -1.0; /* Empty range, to force computation: */ }

    if (SOParams_KeywordPresent(pp, "-autoRange"))
      { o.autoRange = SOParams_GetNextDouble(pp, 0.0, 2.0); }
    else
      { if (o.fRange > 0.0)
          { /* If {-fRange} but no {-autoRange}, don't adjust given range: */
            o.autoRange = 0.0;
          }
        else
          { /* If neither {-fRange} nor {-autoRange}, must adjust: */
            o.autoRange = 1.0;
          }
      }
        
    if (SOParams_KeywordPresent(pp, "-fStep"))
      { o.fStep = SOParams_GetNextDouble(pp, 0.0, INFTY); }
    else
      { o.fStep = -1.0; /* Invalid step, to force computation: */ }

    if (SOParams_KeywordPresent(pp, "-isolineWidth"))
      { o.isolineWidth = SOParams_GetNextDouble(pp, 0.0, INFTY); }
    else
      { o.isolineWidth = 0.15; }

    return o;
  }
