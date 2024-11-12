/* SPPlotFunction.c -- Plots level curves of a spherical function */
/* Last edited on 2011-12-12 11:58:07 by stolfi */

#define PROG_NAME "SPPlotFunction"

#include <SPFunction.h>
#include <SPSpline.h>
#include <SPTriang.h>
#include <SPPlot.h>
#include <SPIntegral.h>
#include <SPOptions.h>
#include <SPPlotOptions.h>
#include <SPH3.h>
#include <r3.h>
#include <r4.h>
#include <pswr.h>
#include <SPBasic.h>
#include <bool.h>
#include <stdio.h>
#include <values.h>
#include <string.h>

#define PROG_USAGE \
  PROG_NAME " \\\n" \
  "  -funcName NAME \\\n" \
  "  [ -bothSides ] \\\n" \
  "  [ {-gridDots | -gridLines} NLON NLAT ] \\\n" \
  "  [ -showTriang ] [ -triName NAME ] \\\n" \
  "  [ -showSupp [ -supp SW SX SY SZ ] ] \\\n" \
  "  -outName NAME \\\n" \
  SPPlotOptions_FunctionHelp

typedef struct Options /* Parsed command line options */
  { char *funcName;       /* Function file (minus ".sfn" extension) or "-". */
    char *triName;        /* Optional triangulation file name (minus ".tri").  */
    bool_t bothSides;     /* TRUE plots both sides of sphere. */
    bool_t showTriang;    /* If func is an {SPSpline}, show its triang.  */
    SPH3_Plane supp;      /* Optional supporting plane, or <0,0,0,0>. */
    bool_t showSupp;      /* If func is an {SPSpline}, show supp. plane.  */
    /* Lon-Lat grid plotting: */
    int gridNLon;         /* Number of steps in longitude. */
    int gridNLat;         /* Number of steps in latitude. */
    bool_t gridDots;      /* TRUE shows grid as face-centered dots, FALSE as rectangles. */
    char *outName;        /* Output file (minus ".ps"/".eps" extension) */
    SPPlotOptions_t plt;  /* Plot format and style options. */
  } Options;

Options GetOptions(int argn, char **argc);
SPFunction *ReadFunction(char *name);
Triangulation *GetTriangulation(SPFunction *f, char *name, int smpOrder);
SPH3_Plane GetSupportingPlane(SPFunction *f, SPH3_Plane *supp);

int main(int argn, char **argc)
  { Options o = GetOptions(argn, argc);
    SPPlotOptions_t *po = &(o.plt);
    int smpOrder = 8; /* Rough integrals are OK here. */
    SPIntegral_SetDefaultSamplingOrder(smpOrder);
    SPFunction *f = ReadFunction(o.funcName);
    Triangulation *tri = GetTriangulation(f, o.triName, smpOrder);
    SPH3_Plane supp = (o.showSupp ? GetSupportingPlane(f, &(o.supp)) : Omega);
    double fMinObs = INFINITY, fMaxObs = -INFINITY;
    SPH3_Point *obs; /* Given observer, or NULL */

    /* Fix perspective parameters, if needed: */
    if (po->projLonLat)
      { obs = NULL; }
    else
      { obs = &(po->obs); 
        SPPlot_FixView
          ( &(po->obs), po->autoObs, 
            &(po->upp), 
            &(po->radius), po->autoRadius, 
            f, TRUE
          );
      }
      
    /* Plot full radius if plotting both sides: */
    if (o.bothSides)
      { if (po->radius != SPPlot_FullRadius) 
          { fprintf(stderr, "** \"-radius\" ignored, using %5.3f\n", SPPlot_FullRadius); }
        po->radius = SPPlot_FullRadius;
      }
      
    /* Fix range, if needed: */
    SPPlot_FixRange(&(po->fRange), &(po->fStep), po->autoRange, f, TRUE);
    
    /* Compute the mesh size in relative scale: */
    double relMeshSize = po->radius * po->meshSize/(po->figSize/2);
      
    /* Open the figure stream: */
    SPPlot_Stream *fps = SPPlot_NewStream
      ( po->eps, o.outName, po->paperSize, po->figSize, po->projLonLat, po->caption.ne );

    auto double evalFunction(S2Point *p);
    
    double evalFunction(S2Point *p)
      { return f->m->eval(f, p); }

    /* Note that {SPPlot_MultipleViews} calls {SPPlot_BeginFigure} internally: */
    SPPlot_MultipleViews
      ( /* {fps} */          fps, 
        /* {funcTag} */      "fn",
        /* {func} */         evalFunction, 
        /* {tri} */          tri, 
        /* {relMeshSize} */  relMeshSize,
        /* {showTriang} */   o.showTriang,
        /* {supp} */         &supp, 
        /* {fRange} */       po->fRange,
        /* {fStep} */        po->fStep,
        /* {obs} */          obs,
        /* {upp} */          &(po->upp),
        /* {rad} */          po->radius,
        /* {dLight} */       &(po->light),
        /* {lineWidth} */    po->lineWidth,
        /* {gridNLon} */     o.gridNLon,
        /* {gridNLat} */     o.gridNLat,
        /* {gridDots} */     o.gridDots,
        /* {aSide} */        +1, 
        /* {bSide} */        (o.bothSides ? -1 : 0), 
        /* {caption} */      &(po->caption), 
        /* {index} */        0, 
        /* {time} */         0.0,
        /* {error} */        0.0,
        /* {capAlign} */     0.5,
        /* {verbose} */      TRUE,
        /* {fMinObs} */      &fMinObs,
        /* {fMaxObs} */      &fMaxObs
      );

    pswr_close_stream(fps);

    fprintf(stderr, "observed range = [ %g _ %g ]\n", fMinObs, fMaxObs);
    return 0;
  }
  
SPFunction *ReadFunction(char *name)
  { char *fileName = (strcmp(name, "-") == 0 ? "-" : txtcat(name, ".sfn"));
    FILE *rd = open_read(fileName, TRUE);
    SPFunction *f = SPFunction_Read(rd);
    if (rd != stdin) { fclose(rd); }
    return f;
  }

Triangulation *GetTriangulation(SPFunction *f, char *name, int smpOrder)
  { SPSpline *fpw;
    if (strcmp(name,"") != 0)
      { char *fileName = (strcmp(name, "-") == 0 ? "-" : txtcat(name, ".tri"));
        Triangulation *tri = SPTriang_ReadCached(fileName, smpOrder);
        return tri;
      }
    else if ((fpw = SPSpline_Cast(f)) != NULL) 
      { return fpw->d->tri; }
    else
      { fprintf(stderr, "** warning - no triangulation - type = %s\n", f->type);
        return NULL;
      }
  }

SPH3_Plane GetSupportingPlane(SPFunction *f, SPH3_Plane *supp)
  { SPSpline *fpw;
    if (r4_norm(&(supp->f)) != 0.0)
      { return (*supp); }
    else if ((fpw = SPSpline_Cast(f)) != NULL)
      { return fpw->d->supp; }
    else
      { fprintf(stderr, "** warning - no supporting plane - type = %s\n", f->type);
        return Omega;
      }
  }

Options GetOptions(int argn, char **argc)
  { Options o;
    SPOptions_Parser_t *pp = SPOptions_NewParser(stderr, argn, argc);

    SPOptions_SetUsage(pp, PROG_USAGE);

    if (SPOptions_TestKeyword(pp, "-funcName"))
      { o.funcName = SPOptions_GetNext(pp); }
    else
      { o.funcName = "-"; }  

    o.showTriang = SPOptions_TestKeyword(pp, "-showTriang");
    
    if (SPOptions_TestKeyword(pp, "-triName"))
      { o.triName = SPOptions_GetNext(pp); }
    else
      { o.triName = ""; }  

    o.showSupp = SPOptions_TestKeyword(pp, "-showSupp");
    
    if (SPOptions_TestKeyword(pp, "-supp"))
      { o.supp.f = SPOptions_GetNextR4(pp, -DBL_MAX, DBL_MAX); }
    else
      { o.supp = NoPlane; }  

    bool_t getGrid = FALSE;
    o.gridNLon = o.gridNLat = 0;
    if (SPOptions_TestKeyword(pp, "-gridDots"))
      { getGrid = TRUE; o.gridDots = TRUE; }
    else if (SPOptions_TestKeyword(pp, "-gridLines"))
      { getGrid = TRUE; o.gridDots = FALSE; }
    else
      { getGrid = FALSE; o.gridDots = FALSE; } 
    if (getGrid)
      { o.gridNLon = SPOptions_GetNextInt(pp, 1, INT_MAX);
        o.gridNLat = SPOptions_GetNextInt(pp, 1, INT_MAX);
      } 

    SPOptions_GetKeyword(pp, "-outName");                               
    o.outName = SPOptions_GetNext(pp);  

    o.plt = SPPlotOptions_FunctionParse(pp);

    o.bothSides = SPOptions_TestKeyword(pp, "-bothSides");
    
    SPOptions_Finish(pp);
    return o;
  }
