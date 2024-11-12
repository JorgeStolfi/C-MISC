/* Test dot products of spherical functions. */
/* Last edited on 2008-05-24 12:25:48 by stolfi */

#define PROG_NAME "SPDotTest"

#include <SPOptions.h>
#include <SPQuad.h>
#include <SPTriang.h>
#include <SPFunction.h>
#include <SPSpline.h>
#include <SPIntegral.h>
#include <SPH3.h>
#include <SPBasic.h>
#include <vec.h>
#include <affirm.h>
#include <nat.h>
#include <bool.h>
#include <stdio.h>
#include <string.h>
#include <values.h>
#include <limits.h>

typedef struct Options
  { char *funcNameF;      /* Name of first function file (minus ".sfn"). */
    char *funcNameG;      /* Name of second function file (minus ".sfn"). */
    char *weightName;     /* Name of weight function file (minus ".sfn"). */
    char *triName;        /* Triangulation file name (minus ".tri") or "" if none. */
    int smpOrder;         /* Triangle sampling order for spherical integrals. */
    bool_t verbose;       /* TRUE to mumble while working. */
    double evExact;       /* Exact value of {<f,g>}, or INFINITY if unknown. */
    double vgExact;       /* Exact value of {<f,g>}, or INFINITY if unknown. */
    double slExact;       /* Ditto, for dot product of gradients. */
  } Options;

Options GetOptions(int argn, char **argc);
SPFunction *ReadFunction(char *name);
  
Triangulation *GetTriangulation
  ( char *triName, 
    int smpOrder, 
    SPFunction *f, 
    SPFunction *g
  );

void GetSamples
  ( Triangulation *tri,
    int smpOrder,
    S2Point_vec_t *sp,
    double_vec_t *wp
  );

void PrintValueAndError(char *dotName, double dotValue, double exact);
  /* Prints the value {dotValue} computed by the procedure {dotName},
    and its deviation from the {exact} value (if the latter is not
    {INFINITY}). */
    
int main(int argn, char **argc)
  { Options o = GetOptions(argn, argc);
    SPIntegral_SetDefaultSamplingOrder(o.smpOrder);
    SPFunction *f = ReadFunction(o.funcNameF);
    SPFunction *g = ReadFunction(o.funcNameG);
    SPFunction *w = ReadFunction(o.weightName);
    Triangulation *tri = GetTriangulation(o.triName, o.smpOrder, f, g); 
    S2Point_vec_t sp = (S2Point_vec_t){ 0, NULL };
    double_vec_t wp = (double_vec_t){ 0, NULL };
    
    GetSamples(tri, o.smpOrder, &sp, &wp);
    
    /* Plain inner products: */
    { fprintf(stderr, "Dot product of values = < f | g >:\n");
      if (o.evExact != INFINITY)
        { fprintf(stderr, "  %-22s = %+22.15e\n", "exact", o.evExact); }

      double evD = SPFunction_Dot(f, NoFMap, g, NoFMap, w, tri, o.verbose);
      PrintValueAndError("Dot", evD, o.evExact);

      double evR = SPFunction_RawDot(f, NoFMap, g, NoFMap, w, tri, o.verbose);
      PrintValueAndError("RawDot", evR, o.evExact);

      double evC = SPFunction_CustomDot(f, NoFMap, g, NoFMap, w, sp, wp);
      PrintValueAndError("CustomDot", evC, o.evExact);
    }

    /* Motion-value inner products: */
    { fprintf(stderr, "Motion dot and value = < (SGrd f) ¤ V | g >:\n");
      if (o.vgExact != INFINITY)
        { fprintf(stderr, "  %-22s = %+22.15e\n", "exact", o.vgExact); }

      double vgD = SPFunction_VelSGrdDot(f, g, w, tri, o.verbose);
      PrintValueAndError("VelSGrdDot", vgD, o.vgExact);

      double vgR = SPFunction_RawVelSGrdDot(f, g, w, tri, o.verbose);
      PrintValueAndError("RawVelSGrdDot", vgR, o.vgExact);

      double vgC = SPFunction_CustomVelSGrdDot(f, g, w, sp, wp);
      PrintValueAndError("CustomVelSGrdDot", vgC, o.vgExact);
    }

    /* Laplacian-value inner products: */
    { fprintf(stderr, "Laplacian dot value = < SLap f | g > = < SGrd f | SGrd g >:\n");
      if (o.slExact != INFINITY)
        { fprintf(stderr, "  %-22s = %+22.15e\n", "exact", o.slExact); }

      double slD = SPFunction_SLapDot(f, g, w, tri, o.verbose);
      PrintValueAndError("SLapDot", slD, o.slExact);

      double slR = SPFunction_RawSLapDot(f, g, w, tri, o.verbose);
      PrintValueAndError("RawSLapDot", slR, o.slExact);

      double slC = SPFunction_CustomSLapDot(f, g, w, sp, wp);
      PrintValueAndError("CustomDotProductGrad", slC, o.slExact);
    }
    return 0;
  }
  
void PrintValueAndError(char *dotName, double dotValue, double exact)
  {
    fprintf(stderr, "  %-22s = %+22.15e", dotName, dotValue);
    if (exact != INFINITY)
      { fprintf(stderr, " error = %+22.15e", dotValue - exact); }
    fprintf(stderr, "\n");
  }

SPFunction *ReadFunction(char *name)
  { if (strcmp(name, "") == 0)
      { return NULL; }
    else
      { FILE *rd = open_read(txtcat(name, ".sfn"), TRUE);
        SPFunction *f = SPFunction_Read(rd);
        fclose(rd);
        return f;
      }
  }

Triangulation *GetTriangulation
  ( char *name, 
    int smpOrder, 
    SPFunction *f, 
    SPFunction *g
  )
  { Triangulation *tri = NULL;
    if ((name != NULL) && (strcmp(name,"") != 0))
      { tri = SPTriang_ReadCached(txtcat(name, ".tri"), smpOrder); }
    if (tri == NULL) 
      { tri = SPSpline_CommonTriangulation(f, g); }
    if (tri == NULL) 
      { tri = SPIntegral_GetDefaultTriangulation(); }
    return tri;
  }

void GetSamples
  ( Triangulation *tri,
    int smpOrder,
    S2Point_vec_t *sp,
    double_vec_t *wp
  )
  { int ns = 0;
    fprintf(stderr, "selecting sample points for custom integrals:\n");
    affirm(tri != NULL, "missing triangulation");
    fprintf(stderr, "  SPIntegral_SampleTriangulation(smpOrder = %d)\n", smpOrder);
    SPIntegral_SampleTriangulation(tri, smpOrder, sp, wp, &ns);
    S2Point_vec_trim(sp, ns);
    double_vec_trim(wp, ns);
  }
  

Options GetOptions(int argn, char **argc)
  { Options o;
    SPOptions_Parser_t *pp = SPOptions_NewParser(stderr, argn, argc);
    
    SPOptions_SetUsage(pp, 
      PROG_NAME " \\\n"
      "  -funcNames FNAME GNAME \\\n"
      "  [ -weightName WNAME ] \\\n"
      "  [ -triName TRINAME ] \\\n"
      "  -smpOrder NUM [ -verbose ]\\\n"
      "  -exact EVDOT VGDOT SLDOT "
    );

    SPOptions_GetKeyword(pp, "-funcNames");                               
    o.funcNameF = SPOptions_GetNext(pp);  
    o.funcNameG = SPOptions_GetNext(pp);  

    if (SPOptions_TestKeyword(pp, "-weightName")) 
      { o.weightName = SPOptions_GetNext(pp); }
    else
      { o.weightName = ""; }

    if (SPOptions_TestKeyword(pp, "-triName")) 
      { o.triName = SPOptions_GetNext(pp); }
    else
      { o.triName = ""; }
      
    SPOptions_GetKeyword(pp, "-smpOrder");                               
    o.smpOrder = SPOptions_GetNextInt(pp, 1, 1000);  

    o.verbose = SPOptions_TestKeyword(pp, "-verbose");
    
    if (SPOptions_TestKeyword(pp, "-exact")) 
      { o.evExact = SPOptions_GetNextDouble(pp, -DBL_MAX, DBL_MAX); 
        o.vgExact = SPOptions_GetNextDouble(pp, -DBL_MAX, DBL_MAX);
        o.slExact = SPOptions_GetNextDouble(pp, -DBL_MAX, DBL_MAX);
      }
    else
      { o.evExact = INFINITY; 
        o.vgExact = INFINITY; 
        o.slExact = INFINITY;
      }

    SPOptions_Finish(pp);
    return o;
  }
