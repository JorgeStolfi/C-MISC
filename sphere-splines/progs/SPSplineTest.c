/* SPSplineTest -- Test of piecewise functions. */
/* Last edited on 2008-05-24 12:30:14 by stolfi */

#include <SPTriang.h>
#include <SPFunction.h>
#include <SPSpline.h>
#include <SPIntegral.h>
#include <SPProcFunction.h>
#include <SPOptions.h>
#include <SPH3.h>
#include <r4.h>
#include <SPBasic.h>
#include <bool.h>
#include <stdio.h>
#include <values.h>

/* This program reads a triangulation {T} and generates a piecewise
  function {f} that, in each triangle of {T}, is a clone of a specific 
  procedurally-defined function {func}.  This is useful to check, 
  for instance, whether point location and barycentric coordinate
  conversion are working correctly. */

typedef struct Options
  { char *triName;    /* Name of triangulation file (minus the ".tri"). */
    char *funcName;   /* Name of function to insert in each piece. */
    char *outName;    /* Name of basis file (minus the ".bas" extension) */
    bool_t bary;      /* TRUE to pass barycentric coords to each piece. */
    SPH3_Plane supp;  /* Supporting plane. */
  } Options;
  
Options GetOptions(int argn, char **argc);
r4_t ParseR4Option(SPOptions_Parser_t *pp);

int main(int argn, char **argc)
  { Options o = GetOptions(argn, argc);
    int smpOrder = 4; /* Rough integrals are OK here. */
    SPIntegral_SetDefaultSamplingOrder(smpOrder);
    char *triFile = txtcat(o.triName, ".tri");
    FILE *rd = open_read(triFile, TRUE);
    Triangulation *tri = SPTriang_Read(rd, smpOrder);
    SPFunction *func = (SPFunction *)SPProcFunction_FromName(o.funcName);
    SPSpline *fpw = SPSpline_Replicate(triFile, tri, o.bary, func);
    fpw->d->supp = o.supp;
    { SPFunction *f = (SPFunction *)fpw;
      FILE *wr = open_write(txtcat(o.outName, ".sfn"), TRUE);
      f->m->write(f, wr);
      fclose(wr);
    }
    return 0;
  }


Options GetOptions(int argn, char **argc)
  {
    Options o;
    SPOptions_Parser_t *pp = SPOptions_NewParser(stderr, argn, argc);
    SPOptions_SetUsage(pp, 
      "SPFindBezSplineBasis \\\n"
      "  -triName NAME \\\n"
      "  -funcName NAME \\\n"
      "  [ -bary ] \\\n"
      "  [ -supp W X Y Z ] \\\n"
      "  -outName NAME\n"
    );

    SPOptions_GetKeyword(pp, "-triName");                               
    o.triName = SPOptions_GetNext(pp);  

    SPOptions_GetKeyword(pp, "-funcName");                               
    o.funcName = SPOptions_GetNext(pp);  

    o.bary = SPOptions_TestKeyword(pp, "-bary");
    
    if (SPOptions_TestKeyword(pp, "-supp"))
      { o.supp.f = ParseR4Option(pp); }
    else
      { o.supp = NoPlane; }

    SPOptions_GetKeyword(pp, "-outName");                               
    o.outName = SPOptions_GetNext(pp);  

    SPOptions_Finish(pp);
    return o;
  }

r4_t ParseR4Option(SPOptions_Parser_t *pp)
  { r4_t p;
    int i;
    for (i = 0; i <= 3; i++)
      { p.c[i] = SPOptions_GetNextDouble(pp, -DBL_MAX, +DBL_MAX); }
    return p;
  } 
