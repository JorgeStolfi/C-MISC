/*  Spherical spline approximation of an arbitrary function. */
/* Last edited on 2006-03-19 12:48:27 by stolfi */

#define PROG_NAME "SPComputeApprox"

/*
  This program computes a least squares approximation {g}, for a 
  spherical function {f}, of the form
  
    g = SUM { u[i]*bas[i] : i = 0..dim-1 }       (1)
    
  where {bas[0..dim-1]} is an arbitrary (but add-compatible) basis.
  The coefficients {u[i]} are found by solving the linear system
    
    G u == b    (2)
    
  where {G} is the rigidity matrix, {G[i,j] == <bas[j] | bas[i]>},
  and the vector {b} is given by {b[i] == <f | bas[i]>}.
*/
    
#include <SPTriang.h>
#include <SPMatrix.h>
#include <SPVector.h>
#include <SPSys.h>
#include <SPFunction.h>
#include <SPSpline.h>
#include <SPApprox.h>
#include <SPOptions.h>
#include <SPIntegral.h>
#include <SPPlot.h>
#include <SPPlotOptions.h>
#include <SPH3.h>

#include <pswr.h>
#include <r4.h>
#include <affirm.h>
#include <nat.h>
#include <bool.h>
#include <SPBasic.h>
#include <vec.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <values.h>
#include <limits.h>

#define PROG_USAGE \
  PROG_NAME " \\\n" \
  "  -funcName NAME \\\n" \
  "  -basisName NAME [-matName NAME ] \\\n" \
  "  [ -triName NAME ] \\\n" \
  SPSys_LinOptionsHelp " \\\n" \
  "  [ -smpOrder NUM ] \\\n" \
  "  -outName NAME \\\n" \
  "  [ -plot ] [ -verbose ] \\\n" \
  SPPlotOptions_FunctionHelp


typedef struct Options
  { char *funcName;         /* Filename of function to approximate (minus ".sfn"). */
    char *basisName;        /* Prefix of basis file names. */
    char *matName;          /* Prefix of matrix file names (default {basisName}). */
    char *triName;          /* Name of triangulation file (minus ".tri"), or "". */
    char *outName;          /* Prefix for output file names. */
    int smpOrder;           /* Triangle sampling order for spherical integrals. */
    bool_t verbose;         /* TRUE prints right-hand-side vectors, etc.. */
    SPSys_LinOptions_t lso; /* Method used to solve the linear system: */
    /* Plotting options: */
    bool_t plot;            /* TRUE to generate function and error plots. */
    SPPlotOptions_t plt;    /* Plot format and style options. */
  } Options;

Triangulation *GetTriangulation(char *name, int smpOrder, Basis F);
  /* Gets triangulation from file {triName}, or gets
    the common triangulation from {F} if {triName} is empty. */

Options *GetOptions(int argn, char **argc);

int main(int argn, char **argc)
  { Options *o = GetOptions(argn, argc);
    SPSys_LinOptions_t *lso = &(o->lso);
    char *lsoText = SPSys_LinOptionsText(lso);
    SPPlotOptions_t *po = &(o->plt);
    SPIntegral_SetDefaultSamplingOrder(o->smpOrder);
    SPFunction *f = SPApprox_ReadFunction(txtcat(o->funcName, ".sfn"));
    Basis bas = SPApprox_ReadBasis(txtcat(o->basisName, ".bas"));
    Triangulation *tri = GetTriangulation(o->triName, o->smpOrder, bas);
    SPFunction *g = NULL; /* Computed approximation. */
    double fMax; /* Max abs values of true solution {f}. */
    double gMax; /* Max abs values of computed solution {g}. */
    double eMax; /* Max abs values of {f-g}. */ 
    double eAvg; /* Root-mean-square value of {f-g}. */ 
    
    int dim = bas.ne;
    SPVector u = double_vec_new(dim); /* Coefficients of solution {g}. */
    SPVector b = double_vec_new(dim);
    SPVector v = double_vec_new(dim);
    SPVector d = double_vec_new(dim); /* Coefficients of correction {dg}. */
    
    affirm((tri == NULL) || (tri->smpOrder == o->smpOrder), "inconsistent smpOrder");
    
    /* Matrices (may be null, depending on method): */
    SPMatrix G;   /* Matrix of basis dot products. */
    SPMatrix GL;  /* Left factor of {G}. */
    SPVector GD;  /* Middle factor of {G}. */
    SPMatrix GR;  /* Right factor of {G}. */
    /* Compute needed matrices: */
    SPApprox_GetBasisMatrices(o->matName, o->lso.mth, &G, &GL, &GD, &GR); 
    
    int iter;
    /* Compute initial solution {g} (coeffs {u}): */
    { int i; for (i = 0; i < b.ne; i++) { b.e[i] = 0.0; } }
    fprintf(stderr, "=== first solution ===\n");
    fprintf(stderr, "computing right-hand side of system...\n");
    SPApprox_ComputeSystemRightHandSide(f, NoFMap, bas, NULL, tri, b, o->verbose);
    if ((o->lso.mth == SPSys_LM_ConjGrad) || (o->lso.mth == SPSys_LM_GaussSeidel))
      { fprintf(stderr, "guessing initial solution...\n");
        SPSys_GuessSol(v);
      }
    fprintf(stderr, "solving linear system...\n");
    SPSys_LinSolve(G, GL, GD, GR, b, lso, v, u, TRUE);
    fprintf(stderr, "building initial solution...\n"); 
    g = SPApprox_BuildFunction(bas, u, o->outName, f);
    fprintf(stderr, "summary for initial solution:\n"); 
    fprintf(stderr, "basis name = %s dimension = %d\n", o->basisName, b.ne);
    fprintf(stderr, "solution method = %s\n", lsoText);
    SPApprox_PrintMaxErrorValues(g, f, tri, &gMax, &fMax, &eMax, &eAvg);
    
    /* Apply residual correction, if so requested: */
    for (iter = 1; iter <= o->lso.residual; iter++)
      { 
        auto double eval_residue(double u, S2Point *p);
          /* Converts {u = f(p)} into the residual {f(p)-g(p)}. */
        
        double eval_residue(double u, S2Point *p)
          { return u - g->m->eval(g, p); }
        
        FuncMap ResidMap = (FuncMap){eval_residue, FALSE, "f(p)-g(p)"};
      
        fprintf(stderr, "=== residual correction (iter = %d) ===\n", iter);
        { int i; for (i = 0; i < b.ne; i++) { b.e[i] = 0.0; } }
        fprintf(stderr, "computing residual right-hand side...\n");
        SPApprox_ComputeSystemRightHandSide(f, ResidMap, bas, NULL, tri, b, o->verbose);
        
        fprintf(stderr, "solving residual system...\n");
        SPSys_LinSolveResidual(G, GL, GD, GR, b, lso, v, d, TRUE);

        /* Apply correction {d} to coeff vector {u}: */
        fprintf(stderr, "adding correction to solution...\n");
        SPSys_AdjustSolution(u, d, o->verbose);

        /* Recompute solution {g} and write to disk: */
        fprintf(stderr, "building corrected solution...\n"); 
        g = SPApprox_BuildFunction(bas, u, o->outName, f);
        fprintf(stderr, "summary for refinement iteration %d:\n", iter); 
        SPApprox_PrintMaxErrorValues(g, f, tri, &gMax, &fMax, &eMax, &eAvg);
      }
    
    /* Plot final solution: */
    { if (o->plot)
        { /* Fix observer and zenith if needed: */
          SPPlot_FixView(&(po->obs), 0.0, &(po->upp), NULL, 0.0, NULL, TRUE);

          /* Number of caption lines: */
          int nCap = (po->eps ? 0 : 2);
          
          /* Fix observer and zenith if needed: */
          double figSize = SPPlot_DefaultFigSize
            (po->eps, po->paperSize, 1, 2, po->projLonLat, nCap);

          /* Open the figure stream: */
          SPPlot_Stream *fps = SPPlot_NewStream
            (po->eps, o->outName, po->paperSize, figSize, po->projLonLat, nCap);

          /* Compute mesh size in sphere units: */
          double relMeshSize = po->meshSize/(figSize/2);
          
          if (po->radius < 1.0) { fprintf(stderr, "** \"-radius\" ignored\n"); }
          
          auto double evalf(S2Point *p);
          double evalf(S2Point *p) { return f->m->eval(f, p); }
          
          SPApprox_PlotFunctionAndError
            ( fps, o->outName, g, evalf,
              tri, relMeshSize, (tri != NULL), TRUE,
              (gMax > 2.0*fMax ? gMax : fMax), eMax, 
              &(po->obs), &(po->upp),
              &(po->caption),
              0, 0.0
            );
          pswr_close_stream(fps);
        }
    }
    return 0;
  }

Triangulation *GetTriangulation(char *name, int smpOrder, Basis F)
  { if ((name != NULL) && (*name != '\000'))
      { FILE *rd = open_read(txtcat(name, ".tri"), TRUE);
        Triangulation *tri = SPTriang_Read(rd, smpOrder);
        fclose(rd);
        return tri;
      }
    else
      { return SPSpline_BasisTriangulation(F); }
  }

Options *GetOptions(int argn, char **argc)
  {
    SPOptions_Parser_t *pp = SPOptions_NewParser(stderr, argn, argc);

    SPOptions_SetUsage(pp, PROG_USAGE);

    Options *o = (Options *)notnull(malloc(sizeof(Options)), "no mem");

    SPOptions_GetKeyword(pp, "-funcName");                               
    o->funcName = SPOptions_GetNext(pp); 
    
    SPOptions_GetKeyword(pp, "-basisName");                               
    o->basisName = SPOptions_GetNext(pp);  
    
    if (SPOptions_TestKeyword(pp, "-matName"))
      { o->matName = SPOptions_GetNext(pp); }
    else
      { o->matName = o->basisName; }
       
    if (SPOptions_TestKeyword(pp, "-triName"))                               
      { o->triName = SPOptions_GetNext(pp); }
    else
      { o->triName = ""; }

    SPOptions_GetKeyword(pp, "-smpOrder");
    o->smpOrder = SPOptions_GetNextInt(pp, 1, INT_MAX);

    SPOptions_GetKeyword(pp, "-outName");                               
    o->outName = SPOptions_GetNext(pp);  
    
    o->verbose = SPOptions_TestKeyword(pp, "-verbose");
                                                 
    /* Method and parameters for linear system solving: */
    o->lso = SPSys_LinOptionsParse(pp, FALSE);
    if (o->lso.mth == SPSys_LM_NONE)
      { SPOptions_Error(pp, "must specify a linear solution method"); }
      
    /* Plotting options: */

    o->plot = SPOptions_TestKeyword(pp, "-plot");
    
    o->plt = SPPlotOptions_FunctionParse(pp);

    SPOptions_Finish(pp);
    return o;
  }

