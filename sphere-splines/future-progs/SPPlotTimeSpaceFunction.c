/* SPPlotTimeSpaceFunction.c -- Plots level curves of a spherical function */
/* Last edited on 2008-05-24 12:29:41 by stolfi */

#include <SPTimeSpaceFunction.h>
/* #include <SPTimeSpaceSpline.h> */
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

typedef struct Options /* Parsed command line options */
  { char *funcName;       /* Function file (minus ".tsf" extension) or "-". */
    char *triName;        /* Optional triangulation file name (minus ".tri").  */
    bool_t showTriang;    /* Show the triangulation, if any.  */
    SPH3_Plane supp;      /* Optional supporting plane, or <0,0,0,0>. */
    bool_t showSupp;      /* Show the supp. plane, if any.  */
    bool_t bothSides;     /* TRUE plots both sides of sphere. */
    char *outName;        /* Output file (minus ".ps"/".eps" extension) */
    SPPlotOptions_t plt;  /* Plot format and style options. */
    /* Time stepping options: */
    double timeStart;     /* Time of first plot. */
    double timeStep;      /* Time step. */
    int nTimes;           /* Number of epochs to plot. */
  } Options;

Options GetOptions(int argn, char **argc);
SPTimeSpaceFunction *ReadFunction(char *name);
Triangulation *GetTriangulation(SPTimeSpaceFunction *f, char *name, int smpOrder);
SPH3_Plane GetSupportingPlane(SPTimeSpaceFunction *f, SPH3_Plane *supp);

int main(int argn, char **argc)
  { Options o = GetOptions(argn, argc);
    SPPlotOptions_t *po = &(o.plt);
    int smpOrder = 8; /* Rough integrals are OK here. */
    SPIntegral_SetDefaultSamplingOrder(smpOrder);
    SPTimeSpaceFunction *f = ReadFunction(o.funcName);
    Triangulation *tri = GetTriangulation(f, o.triName, smpOrder);
    SPH3_Plane supp = (o.showSupp ? GetSupportingPlane(f, &(o.supp)) : Omega);
    
    double fGlobMinObs = INFINITY, fGlobMaxObs = -INFINITY;

    /* Fix perspective parameters, if needed: */
    SPPlot_FixView
      ( &(po->obs), po->autoObs, 
        &(po->upp), 
        &(po->radius), po->autoRadius, 
        NULL, TRUE
      );
      
    /* Plot full radius if plotting both sides: */
    if (o.bothSides)
      { if (po->radius != SPPlot_FullRadius) 
          { fprintf(stderr, "** \"-radius\" ignored, using %5.3f\n", SPPlot_FullRadius); }
        po->radius = SPPlot_FullRadius;
      }
      
    /* Fix range, if needed: */
    SPPlot_FixRange(&(po->fRange), &(po->fStep), po->autoRange, NULL, TRUE);
    
    /* Compute the mesh size in relative scale: */
    double relMeshSize = po->radius * po->meshSize/(po->figSize/2);
      
    /* Open the figure stream: */
    SPPlot_Stream *fps = SPPlot_NewStream
      ( po->eps, o.outName, po->paperSize, po->figSize, po->projLonLat, po->caption.ne );

    int ie; /* Epoch index. */
    for (ie = 0; ie < o.nTimes; ie++)
      { 
        double t = o.timeStart + ie*o.timeStep;
        
        auto double evalFunction(S2Point *p);

        double evalFunction(S2Point *p)
          { return f->m->eval(f, p, t); }
          
        /* Observed function range in this epoch: */
        double fMinObs = INFINITY, fMaxObs = -INFINITY;

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
            /* {obs} */          &(po->obs),
            /* {upp} */          &(po->upp),
            /* {rad} */          po->radius,
            /* {dLight} */       &(po->light),
            /* {lineWidth} */    po->lineWidth,
            /* {gridNLon} */     0,
            /* {gridNLat} */     0,
            /* {gridDots} */     FALSE,
            /* {aSide} */        +1, 
            /* {bSide} */        (o.bothSides ? -1 : 0), 
            /* {caption} */      &(po->caption), 
            /* {index} */        ie, 
            /* {time} */         t,
            /* {error} */        0.0,
            /* {capAlign} */     0.5,
            /* {verbose} */      TRUE,
            /* {fMinObs} */      &fMinObs,
            /* {fMaxObs} */      &fMaxObs
          );

       fprintf(stderr, "observed range in epoch = [ %g _ %g ]\n", fMinObs, fMaxObs);
       if (fMinObs < fGlobMinObs) { fGlobMinObs = fMinObs; }
       if (fMaxObs > fGlobMaxObs) { fGlobMaxObs = fMaxObs; }
     }

    pswr_close_stream(fps);

    fprintf
      ( stderr, "observed range in all epochs = [ %g _ %g ]\n", 
        fGlobMinObs, fGlobMaxObs
      );
    return 0;
  }
  
SPTimeSpaceFunction *ReadFunction(char *name)
  { char *fileName = (strcmp(name, "-") == 0 ? "-" : txtcat(name, ".tsf"));
    FILE *rd = open_read(fileName, TRUE);
    SPTimeSpaceFunction *f = SPTimeSpaceFunction_Read(rd);
    fclose(rd);
    return f;
  }

Triangulation *GetTriangulation(SPTimeSpaceFunction *f, char *name, int smpOrder)
  { /* SPSpline *fpw; */
    if ((name != NULL) && (strcmp(name,"") != 0))
      { char *fileName = (strcmp(name, "-") == 0 ? "-" : txtcat(name, ".tri"));
        Triangulation *tri = SPTriang_ReadCached(fileName, smpOrder);
        return tri;
      }
    /* else if ((fpw = SPSpline_Cast(f)) != NULL)  */
    /*   { return fpw->d->tri; } */
    else
      { fprintf(stderr, "** warning - no triangulation - type = %s\n", f->type);
        return NULL;
      }
  }

SPH3_Plane GetSupportingPlane(SPTimeSpaceFunction *f, SPH3_Plane *supp)
  { /* SPSpline *fpw; */
    if (r4_norm(&(supp->f)) != 0.0)
      { return (*supp); }
    /* else if ((fpw = SPSpline_Cast(f)) != NULL) */
    /*   { return fpw->d->supp; } */
    else
      { fprintf(stderr, "** warning - no supporting plane - type = %s\n", f->type);
        return Omega;
      }
  }


Options GetOptions(int argn, char **argc)
  { Options o;
    SPOptions_Parser_t *pp = SPOptions_NewParser(stderr, argn, argc);

    SPOptions_SetUsage(pp, 
      "SPPlotTimeSpaceFunction \\\n"
      "  -funcName NAME \\\n"
      "  [ -showTriang ] [ -triName NAME ] \\\n"
      "  [ -showSupp [ -supp SW SX SY SZ ] ] \\\n"
      "  -outName NAME \\\n"
      SPPlotOptions_FunctionHelp " \n"
    );

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
    
    if (o.showSupp && SPOptions_TestKeyword(pp, "-supp"))
      { o.supp.f = SPOptions_GetNextR4(pp, -DBL_MAX, DBL_MAX); }
    else
      { o.supp = NoPlane; }  

    SPOptions_GetKeyword(pp, "-outName");                               
    o.outName = SPOptions_GetNext(pp);  

    o.plt = SPPlotOptions_FunctionParse(pp);

    o.bothSides = SPOptions_TestKeyword(pp, "-bothSides");
    
    SPOptions_GetKeyword(pp, "-timeStart");                               
    o.timeStart = SPOptions_GetNextDouble(pp, -DBL_MAX, +DBL_MAX);
    
    SPOptions_GetKeyword(pp, "-timeStep");                               
    o.timeStep = SPOptions_GetNextDouble(pp, -DBL_MAX, +DBL_MAX);
    
    SPOptions_GetKeyword(pp, "-nTimes");                               
    o.nTimes = SPOptions_GetNextInt(pp, 0, 1000);

    SPOptions_Finish(pp);
    return o;
  }
