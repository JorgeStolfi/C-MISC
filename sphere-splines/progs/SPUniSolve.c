/* Solves the Helmholtz eqn by single-scale iteration. */
/* Last edited on 2023-10-15 03:25:11 by stolfi */

#define PROG_NAME "SPUniSolve"

/*
  This program solves a differential equation on the sphere
  
    { D(f)(p) == FMap(f(p), p) }         (1)
    
  where {D} is a linear differential operator, {FMap} is a given
  function, and {f} is the function to be determined. For example, 
  the Helmholz differential equation on the sphere
  
    { D(f)(p) = - SLap(f)(p) + c*f(p) }        (2)
    
  where {SLap} is the spherical Laplacian operator, and {c} is a given
  constant.
  
  It is assumed that an approximate solution {g} to (2) in the
  corresponding function space {V} can be found by solving the
  non-linear system
  
    { H a == b }         (3)
    
  where 
  
    { a[i] } is coefficient {i} of the solution 
           {g} (to be determined) in the basis {bas}; 
  
    { b[i] == <F | bas[i]> }, where {F} is the spherical function
            {p -> FMap(g(p),p)} and {<|>} is the scalar product.
           
    { H[i,j] == <SGrd(bas[j]) | SGrd(bas[i])> + c*<bas[j]|bas[i]> }
            where {SGrd} is the spherical gradient operator.
        
  Note that this system is not trivial since the right-hand side {b}
  depends on the function {g}, and therefore on the coefficient vector
  {a}, usually in a non-linear way.  So (3) must be solved
  iteratively.
*/

#define _GNU_SOURCE
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <values.h>
#include <limits.h>

#include <pswr.h>
#include <r4.h>
#include <rn.h>
#include <affirm.h>
#include <nat.h>
#include <bool.h>

#include <SPTriang.h>
#include <SPMatrix.h>
#include <SPVector.h>
#include <SPSys.h>
#include <SPTriang.h>
#include <SPFunction.h>
#include <SPSpline.h>
#include <SPIntegral.h>
#include <SPProcFunction.h>
#include <SPApprox.h>
#include <SPFuncMap.h> 
#include <SPOptions.h>
#include <SPPlot.h> 
#include <SPPlotOptions.h> 
#include <SPH3.h>

#define PROG_USAGE \
  PROG_NAME " \\\n" \
  "  -rhsName NAME -coeff NUM -solName NAME \\\n" \
  "  -basisName NAME [ -matName NAME ] \\\n" \
  SPSys_LinOptionsHelp " \\\n" \
  SPSys_GenOptionsHelp " \\\n" \
  "  [ -stopError NUM ] [ -minval NUM ] [ -smpOrder NUM ] \\\n" \
  "  -outName NAME \\\n" \
  "  [ -plotAll ] [ -plotFinal ] [ -verbose ] \\\n" \
  SPPlotOptions_FunctionHelp 

typedef struct Options
  { char *rhsName;          /* Name of right-hand-side operator {FMap}. */
    double coeff;           /* Coefficient {c} of Helmholtz equation. */
    char *solName;          /* Name of true solution. */
    char *basisName;        /* Prefix for basis filenames (minus extension). */
    char *matName;          /* Prefix for matrix filenames (default {basisName}). */
    char *outName;          /* Prefix for output filenames. */
    bool_t verbose;         /* TRUE prints right-hand-side vectors, etc.. */
    /* Parameters for PDE integration: */
    double minVal;          /* Threshold for matrix factor cleanup. */ 
    SPSys_LinOptions_t lso; /* Parameters for linear system solving */
    SPSys_GenOptions_t gso; /* Parameters for solving the non-linear system. */
    double stopError;       /* Stop when this close to true solution. */
    int smpOrder;           /* Triangle sampling order for spherical integrals. */
    /* Plotting options: */
    bool_t plotFinal;       /* TRUE to plot final solution and error. */
    bool_t plotAll;         /* TRUE to plot every iteration and error. */
    SPPlotOptions_t plt;    /* Plot format and style options. */
  } Options;

Options *GetOptions(int argn, char **argc);
SPFuncMap GetRHSOperator(char *rhsName);
SPFunction *GetTrueSolution(char *solName);

int main(int argn, char **argc)
  { Options *o = GetOptions(argn, argc);
    SPSys_LinOptions_t *lso = &(o->lso);
    SPSys_GenOptions_t *gso = &(o->gso);
    SPPlotOptions_t *po = &(o->plt);
    SPIntegral_SetDefaultSamplingOrder(o->smpOrder);
    Basis bas = SPApprox_ReadBasis(txtcat(o->basisName, ".bas"));
    int dim = bas.ne;
    Triangulation *tri = SPSpline_BasisTriangulation(bas);
    SPFuncMap RHS = SPFuncMap_FromName(o->rhsName);
    SPFunction *s = GetTrueSolution(o->solName);
    
    /* A function that evaluates the true solution: */
    auto double eval_s(S2Point *p);
    double eval_s(S2Point *p) { return s->m->eval(s, p); }
    
    SPVector a = double_vec_new(dim); /* Coeffs of current solution. */
    SPVector b = double_vec_new(dim); /* Right-hand side vector. */
    SPVector y = double_vec_new(dim);
    SPVector aNew = double_vec_new(dim);
    
    SPFunction *g;  /* Current solution. */
    
    /* Fix perspective parameters, if needed: */
    SPPlot_FixView(&(po->obs), 0.0, &(po->upp), NULL, 0.0, NULL, TRUE);

    /* Number of caption lines: */
    int nCap = (po->eps ? 0 : 2);

    /* Open the figure stream: */
    SPPlot_Stream *fps = SPPlot_NewStream
      (po->eps, o->outName, po->paperSize, po->figSize, po->projLonLat, nCap);
    
    /* Compute mesh size in sphere units: */
    double relMeshSize = po->meshSize/(po->figSize/2);

    FILE *errWr = open_write(txtcat(o->outName, ".erp"), TRUE);
    
    double gMax, sMax, eMax, eAvg;

    /* Matrices (may be null, depending on method): */
    SPMatrix H;   /* Matrix of basis dot products. */
    SPMatrix HL;  /* Left factor of {H}. */  
    SPVector HD;  /* Middle factor of {H}. */  
    SPMatrix HR;  /* Right factor of {H}. */ 
    /* Compute needed matrices: */
    SPApprox_GetHelmholtzMatrices(o->matName, o->coeff, lso->mth, o->minVal, &H, &HL, &HD, &HR);

    double aDiff = INFINITY; /* Max change in coeffs between iterations. */
    int iter = 0;
    fprintf(errWr, "# %6s %6s %16s\n",  "scale", "iter", "rmsError");
    fprintf(stderr, "guessing initial approximation...\n"); 
    SPSys_GuessSol(a); 
    while (TRUE)
      { 
        char *iterTag = fmt_int(iter,6); 
        char *iterName = txtcat3(o->outName,"-", iterTag);
        fprintf(stderr, "=== iteration %d ===\n", iter); 
        fprintf(stderr, "building approximation...\n"); 
        g = SPApprox_BuildFunction(bas, a, iterName, s);
        fprintf
          ( stderr, "basis = %s dim = %d rhs = %s sol = %s\n", 
            o->basisName, bas.ne, o->rhsName, o->solName
          );
        fprintf(stderr, "summary for iteration %d\n", iter); 
        SPApprox_PrintMaxErrorValues(g, s, tri, &gMax, &sMax, &eMax, &eAvg);
        fprintf(errWr, "  %6d %6d %16.12f %16.12f\n",  0, iter, eAvg, eMax);
        
        /* Check termination criteria: */
        double aNorm = rn_L_inf_norm(dim, a.e);
        if (SPSys_GenStopCondition(gso, o->stopError, aNorm, iter, eMax, aDiff, TRUE))
          { break; }
        
        if (o->plotAll)
          { if (po->radius < 1.0) { fprintf(stderr, "** \"-radius\" ignored\n"); }
            SPApprox_PlotFunctionAndError
              ( fps, iterTag, g, eval_s, 
                tri, relMeshSize, (tri != NULL), FALSE,
                (gMax > 2.0*sMax ? gMax : sMax), eMax,
                &(po->obs), &(po->upp),
                &(po->caption),
                iter, 0.0
              );
          }
        
        /* Do another iteration: */
        fprintf(stderr, "computing right-hand side...\n"); 
        SPApprox_ComputeSystemRightHandSide(g, RHS, bas, NULL, tri, b, o->verbose);
        fprintf(stderr, "solving system...\n"); 
        SPSys_LinSolve(H, HL, HD, HR, b, lso, y, aNew, TRUE);
        fprintf(stderr, "updating solution coefficients...\n"); 
        aDiff = SPSys_UpdateSolution(aNew, a);
        iter++;
      }
    fclose(errWr);

    if ((o->plotFinal) || (o->plotAll))
      { if (po->radius < 1.0) { fprintf(stderr, "** \"-radius\" ignored\n"); }
        char *finTag = "fin";
        SPApprox_PlotFunctionAndError
          ( fps, finTag, g, eval_s, 
            tri, relMeshSize, (tri != NULL), TRUE,
            (gMax > 2.0*sMax ? gMax : sMax), eMax,
            &(po->obs), &(po->upp),
            &(po->caption), 
            iter, 0.0
          );
      }
    pswr_close_stream(fps);
    return 0;
  }
  
SPFuncMap GetRHSOperator(char *rhsName)
  { return SPFuncMap_FromName(rhsName); }

SPFunction *GetTrueSolution(char *solName)
  { SPFunction *f = (SPFunction *)SPProcFunction_FromName(solName);
    if (f == NULL) 
      { fprintf(stderr, "Unknown solution = \"%s\"\n", solName);
        affirm(FALSE, "aborted");
      }
    return f;
  }


Options *GetOptions(int argn, char **argc)
  {
    Options *o = (Options *)notnull(malloc(sizeof(Options)), "no mem");
    SPOptions_Parser_t *pp = SPOptions_NewParser(stderr, argn, argc);

    SPOptions_SetUsage(pp, PROG_USAGE);

    SPOptions_GetKeyword(pp, "-rhsName");                               
    o->rhsName = SPOptions_GetNext(pp);  

    SPOptions_GetKeyword(pp, "-coeff");                               
    o->coeff = SPOptions_GetNextDouble(pp, -DBL_MAX, +DBL_MAX);  

    SPOptions_GetKeyword(pp, "-solName");                               
    o->solName = SPOptions_GetNext(pp);  

    SPOptions_GetKeyword(pp, "-basisName");                               
    o->basisName = SPOptions_GetNext(pp);  
       
    if (SPOptions_TestKeyword(pp, "-matName"))
      { o->matName = SPOptions_GetNext(pp); }
    else
      { o->matName = o->basisName; }

    SPOptions_GetKeyword(pp, "-outName");                               
    o->outName = SPOptions_GetNext(pp);  
    
    o->verbose = SPOptions_TestKeyword(pp, "-verbose");
                                                 
    if (SPOptions_TestKeyword(pp, "-minVal"))
      { o->minVal = SPOptions_GetNextDouble(pp, 0.0, DBL_MAX); }
    else 
      { o->minVal = 0.0; }
    
    /* Parameters for the linear equation system solver: */
    o->lso = SPSys_LinOptionsParse(pp, FALSE);
    if (o->lso.mth == SPSys_LM_NONE)
      { SPOptions_Error(pp, "must specify a linear solution method"); }

    /* Parameters for the global (non-linear) equation system solver: */
    o->gso = SPSys_GenOptionsParse(pp, FALSE);

    if (SPOptions_TestKeyword(pp, "-stopError"))
      { o->stopError = SPOptions_GetNextDouble(pp, 0, +DBL_MAX); }
    else
      { o->stopError = 0.0; }

    SPOptions_GetKeyword(pp, "-smpOrder");
    o->smpOrder = SPOptions_GetNextInt(pp, 1, INT_MAX);

    /* Plotting options: */

    o->plotAll = SPOptions_TestKeyword(pp, "-plotAll");

    o->plotFinal = SPOptions_TestKeyword(pp, "-plotFinal");

    o->plt = SPPlotOptions_FunctionParse(pp);
    
    SPOptions_Finish(pp);
    return o;
  }

