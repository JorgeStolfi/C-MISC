/* SOPlotFunction - Plots an arbitrary function from R^2 to R. */
/* Last edited on 2004-06-20 10:09:46 by stolfi */

/*
  This program plots an arbitrary function from R^2 to R,
  with graded color and/or level curves.
*/

#include <SOGrid.h>
#include <SOFunction.h>
#include <SOIntegral.h>
#include <SOParams.h>
#include <SOPlotParams.h>
#include <SOBasic.h>
#include <SO2DPlot.h>
#include <SOApprox.h>

#include <affirm.h>
#include <nat.h>
#include <vec.h>

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <values.h>
#include <limits.h>

#define MAXSTR 50
#define MAX_PLOT_STEPS 10

typedef struct Options
  { char *funcName;   /* Filename of function to approximate (minus ".fun"). */
    char *treeName;   /* Filename of grid tree (minus ".tree"), or "". */
    char *outName;    /* Prefix for output file names. */
    /* Plotting options: */
    PlotOptions plt;  /* Plot format and style options. */
  } Options;

Options *GetOptions(int argn, char **argc);
SOGrid_Tree *ReadTree(char *treeName);

int main(int argn, char **argc)
  { 
    Options *o = GetOptions(argn, argc);
    SOFunction *f = SOApprox_ReadFunction(addext(o->funcName, ".fun")); 
    SOGrid_Tree *tree = SOApprox_ReadTree(addext(o->treeName, ".tree"));

    double fMin = INFTY, fMax = -INFTY;

    affirm(f->d->pDim == 2, "function domain must be two-dimensional");
    affirm(f->d->fDim == 1, "function must return scalar result");
    
    SOIntegral_GaussOrder = 6; 

    SOApprox_PlotFunction(f, tree, o->outName, &(o->plt), &fMin, &fMax);

    fprintf(stderr, "observed function extrema:");
    fprintf(stderr, " min = %16.12f max = %16.12f\n", fMin, fMax);
    fprintf(stderr, "\n");
    
    return 0;
  }

#define PPUSAGE SOParams_SetUsage

Options *GetOptions(int argn, char **argc)
  {
    Options *o = (Options *)notnull(malloc(sizeof(Options)), "no mem");
    SOParams_T *pp = SOParams_NewT(stderr, argn, argc);

    PPUSAGE(pp, "SOPlotFunction \\\n");
    PPUSAGE(pp, "  -funcName NAME \\\n");
    PPUSAGE(pp, "  [ -treeName NAME ] \\\n");
    PPUSAGE(pp, "  -outName NAME \\\n");
    PPUSAGE(pp, SOPlotParams_FunctionHelp " \n");

    SOParams_GetKeyword(pp, "-funcName");                               
    o->funcName = SOParams_GetNext(pp);  
       
    if (SOParams_KeywordPresent(pp, "-treeName"))
      { o->treeName = SOParams_GetNext(pp); }
    else
      { o->treeName = ""; }
       
    SOParams_GetKeyword(pp, "-outName");                               
    o->outName = SOParams_GetNext(pp);  

    /* Plotting options: */

    o->plt = SOPlotParams_FunctionParse(pp);

    SOParams_Finish(pp);
    return o;
  }

