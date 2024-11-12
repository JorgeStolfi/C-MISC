/* See SPPlotParms.h */
/* Last edited on 2008-05-24 12:29:01 by stolfi */

#include <SPPlotOptions.h>
#include <SPPlot.h>
#include <SPOptions.h>
#include <SPBasic.h>
#include <r3.h>
#include <r4.h>
#include <SPH3.h>
#include <js.h>
#include <values.h>

SPPlotOptions_t SPPlotOptions_TriangParse(SPOptions_Parser_t *pp)
  {
    SPPlotOptions_t o;

    /* General plotting options: */
    
    if (SPOptions_TestKeyword(pp, "-eps"))
      { o.eps = TRUE; } 
    else if (SPOptions_TestKeyword(pp, "-ps"))
      { o.eps = FALSE; }
    else
      { o.eps = FALSE; }

    if (SPOptions_TestKeyword(pp, "-paperSize"))
      { o.paperSize = SPOptions_GetNext(pp); }
    else
      { o.paperSize = "letter"; }  

    if (SPOptions_TestKeyword(pp, "-meshSize"))
      { o.meshSize = SPOptions_GetNextDouble(pp, 0.0, 10000.0); }
    else
      { o.meshSize = 2.0; }

    o.caption = string_vec_new(10);
    { int nCapLines = 0;
      while (SPOptions_TestKeyword(pp, "-caption"))
        { string_vec_expand(&(o.caption), nCapLines);
          o.caption.e[nCapLines] = SPOptions_GetNext(pp);
          nCapLines++;
        }
      string_vec_trim(&(o.caption), nCapLines);
    }

    if (SPOptions_TestKeyword(pp, "-figSize"))
      { o.figSize = SPOptions_GetNextDouble(pp, 5.0, 2000.0); }
    else
      { /* Guess: user caption, plus 2 lines of figure name if not {eps}: */
        int nCap = o.caption.ne + (o.eps ? 0 : 2);
        /* Default figure size -- fits in page, allowing caption: */
        o.figSize = SPPlot_DefaultFigSize
          (o.eps, o.paperSize, 1, 1, o.projLonLat, nCap);
      }
    
    if (SPOptions_TestKeyword(pp, "-lineWidth"))
      { o.lineWidth = SPOptions_GetNextDouble(pp, 0.0, DBL_MAX); }
    else
      { o.lineWidth = 0.15; }

    /* View parameters: */
    
    o.projLonLat = SPOptions_TestKeyword(pp, "-projLonLat");

    if (SPOptions_TestKeyword(pp, "-obs"))
      { o.obs.h = SPOptions_GetNextR4(pp, -DBL_MAX, DBL_MAX); }
    else
      { /* INvalid observer, to force computation: */
        o.obs = (SPH3_Point){{{0.0,0.0,0.0,0.0}}};
      }

    if (SPOptions_TestKeyword(pp, "-upp"))
      { o.upp.h = SPOptions_GetNextR4(pp, -DBL_MAX, DBL_MAX); }
    else
      { o.upp = (SPH3_Point){{{0.0,0.0,0.0,0.0}}}; }
      
    if (SPOptions_TestKeyword(pp, "-radius"))
      { o.radius = SPOptions_GetNextDouble(pp, 0.0, 2.0); }
    else
      { /* Invalid radius to force computation: */
        o.radius = -1.0;
      }
    
    /* Lighting parameters: */

    if (SPOptions_TestKeyword(pp, "-light"))
      { o.light = SPOptions_GetNextDir(pp); }
    else
      { o.light = (r3_t){{0.1, 0.2, 1.0}};
        r3_dir(&(o.light), &(o.light));
      }

    if (SPOptions_TestKeyword(pp, "-shadow"))
      { o.shadow = SPOptions_GetNextDouble(pp, 0.0, 1.0); }
    else
      { o.shadow = 0.1; }

    /* Just in case: */
    
    o.autoObs = 0.0;
    o.autoRadius = 0.0;

    /* Parameters for IsoLines/PaintValues: */

    o.fStep = 0.1;
    o.fRange = -1.0;  /* Invalid range, to force computation. */ 
    o.autoRange = 1.0;

    return o;
  }

SPPlotOptions_t SPPlotOptions_FunctionParse(SPOptions_Parser_t *pp)
  {
    SPPlotOptions_t o = SPPlotOptions_TriangParse(pp);

    /* Automatic view selection: */
    
    if (SPOptions_TestKeyword(pp, "-autoObs"))
      { o.autoObs = SPOptions_GetNextDouble(pp, 0.0, 1.0); }
    else
      { /* If "-obs" was given, default {autoObs} is 0, else near 1: */
        if (r4_norm(&(o.obs.h)) != 0.0)
          { /* If "-obs" but no "-autoObs",  stick with given obs: */
            o.autoObs = 0.0;
          }
        else
          { /* If neither "-obs" nor "-autoObs", adjust obs: */
            o.autoObs = 0.95;
          }
      }
    
    if (SPOptions_TestKeyword(pp, "-autoRadius"))
      { o.autoRadius = SPOptions_GetNextDouble(pp, 0.0, 2.0); }
    else
      { if (o.radius > 0.0)
          { /* If "-radius" but no "-autoRadius", don't adjust given radius: */
            o.autoRadius = 0.0; 
          }
        else
          { /* If neither "-radius" nor "-autoRadius", must adjust: */
            o.autoRadius = 1.0;
          }
      }
    
    /* Parameters for PaintValues: */
    
    if (SPOptions_TestKeyword(pp, "-fRange"))
      { o.fRange = SPOptions_GetNextDouble(pp, 0.0, DBL_MAX); }
    else
      { o.fRange = -1.0; /* Empty range, to force computation: */ }

    if (SPOptions_TestKeyword(pp, "-autoRange"))
      { o.autoRange = SPOptions_GetNextDouble(pp, 0.0, 2.0); }
    else
      { if (o.fRange > 0.0)
          { /* If "-fRange" but no "-autoRange", don't adjust given range: */
            o.autoRange = 0.0;
          }
        else
          { /* If neither "-fRange" nor "-autoRange", must adjust: */
            o.autoRange = 1.0;
          }
      }
        
    /* Parameters for IsoLines: */

    if (SPOptions_TestKeyword(pp, "-fStep"))
      { o.fStep = SPOptions_GetNextDouble(pp, 0.0, DBL_MAX); }
    else
      { o.fStep = -1.0; /* Invalid step, to force computation: */ }

    return o;
  }
