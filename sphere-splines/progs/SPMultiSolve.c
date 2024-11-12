/* Solves the Helmholtz eqn by multiscale iteration. */
/* Last edited on 2023-10-15 03:25:24 by stolfi */

#define PROG_NAME "SPMultiSolve"

/* This program solves a differential equation on the sphere
  
    D(f)(p) == FMap(f(p), p)       (1)
    
  where {D} is a linear differential operator, {FMap} is a given
  function, and {f} is the function to be determined. For example, the
  Helmholtz equation has
  
    D(f)(p) = SLap(f)(p) + c*f(p),     (2)
    
  where {SLap} is the spherical Laplacian operator, and {c} is a given
  constant.
  
  The program uses a previously built sequence {m} of approximation
  spaces {V[r] = m[r].V}, at various scales of resolution {r =
  0..maxScale} ({0} is the coarsest, {maxScale} is the finest), with
  respective bases {B[r] = m[r].bas}.
  
  It is assumed that, at each level {r} of the hierarchy, an approximate
  solution {g} to (1) in the corresponding function space {V} can be
  found by solving the non-linear system
  
    H a == b         (3)
    
  where 
  
    {a[i]} is coefficient {i} of the solution 
           {g} (to be determined) in the basis {B} of {V}; 
  
    {b[i] == <F|B[i]>}, where {F} is the spherical function
            {p -> FMap(g(p),g)} and {<|>} is the 
            integral scalar product.
           
    {H[i,j] == <D(B[j])|B[i]>}
    
  For the Helmholtz operator (2), it can be shown that 
  
    H[i,j] == <SGrd(B[j])|SGrd(B[i])> + c*<B[j]|B[i]>  (4)
  
  where {SGrd} is the spherical gradient operator.
        
  Note that the right-hand side {b} of (3) generally depends on the
  function {g}, and therefore on the coefficient vector {a}, usually
  in a non-linear way. So (3) must be solved iteratively.
  
  The multiscale method accelerates the resolution of (3) by
  performing each iteration on all meshes. After performing a few
  iterations at some scale {r}, we either raise the current solution
  {g[r]} from space {V[r]} to the coarser space {V[r-1]}, or lower it
  onto the finer space {V[r+1]}.
  
  The raising and lowering operations are defined in terms of the
  spaces {W[r]} and {Z[r]}, where {W[r]} is the orthogonal projection
  of the coarser space {V[s] == V[r-1]} onto {V[r]}, and {Z[r]} is
  the orthogonal complement of {W[r]} in {V[r]}.  It is assumed that
  {W[r]} and {V[s]} have the same dimension, and that, moreover,
  
     inf { |v'|/|v| : v in V[s] } >> 0     (3)
     
  where {v'} is the orthogonal projection {v} onto {V[r]}.
  
  The coarsening operation takes a vector {gHi[r]} in the
  high-resolution space {V[r]} and produces a closest approximation
  {gLo[r-1]} in the lower-resolution space {V[r-1]}; the approximation may
  be simply the orthogonal projection of {gHi[r]} onto {V[r-1]}. Note
  that this projection can be viewed as a projection onto the subspace
  {W[r]} of {V[r]}, followed by a dimension-preserving projection from
  {W[r]} to {V[r-1]}.
  
  The refining operation takes two vectors {gLo[r-1]} in {V[r-1]} and
  {gHiOld[r]} in {V[r]}, and produces a new vector {gHiNew[r]} in {V[r]} that
  combines all the information that can be taken from {gLo[r-1]}, whith
  watever information is present in {gHiOld[r]} but cannot be represented
  in {V[r-1]}.  That is, we decompose {gHiOld[r]} into a {fine-scale}
  component {zOld[r]} in {Z[r]} and a {coarse-scale} component {wOld[r]}
  in {W[r]}; we discard the latter, and define {gHiNew[r]} as {zOld[r]}
  plus the best possible approximation {wNew[r]} of {gLo[r-1]} in {W[r]}.
  
  For the coarsening and refining operations we need three matrices:
  {L[r]} is a bijection from {V[r-1]} to {W[r]}, {R[r]} is the
  orthogonal projection from {V[r]} onto {V[r-1]}, and {T[r]} is the
  orthogonal projection of {V[r]} onto {Z[r]}. Then coarsening a
  vector {gHi[r]} from {V[r]} to {V[r-1]} means computing {gLo[r-1] ==
  R[r]gHi[r]}; and refining a vector {gLo[r-1]} from {V[r-1]} to
  {V[r]} means computing {gHiNew[r] == T[r]gHiOld[r] + L[r]gLo[r-1]}.
  
  It turns out that (check!) 
  
    L[r] == G[r]^{-1} F[r]
    R[r] == F[r]^% G[r-1]^{-1} 
    T[r] == ???
    
  where {A^%} is the transpose of {A}, and
   
    G[r][i,j] == <bas[r][j] | bas[r][i]>
    F[r][i,j] == <bas[r-1][j] | bas[r][i]>
    
  Note that {G[r]} is square and non-singular, with {dim(V[r])} rows and
  columns.  Note also that {F[r]} is rectangular, with {dim(V[r-1])} rows
  and {dim(V[r])} columns; and has full rank because of condition (3) above.
  
  The linear matrices {F[r]} needed to transfer the solutions between meshes
  and to solve the system at each mesh can be either read from 
  files or computed on-the-fly.  
  
*/

#define _GNU_SOURCE
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <values.h>
#include <limits.h>
 
#include <pswr.h>
#include <r3.h>
#include <rn.h>
#include <affirm.h>
#include <nat.h>
#include <bool.h>
 
#include <SPTriang.h>
#include <SPMatrix.h>
#include <SPVector.h>
#include <SPSys.h>
#include <SPFunction.h>
#include <SPSpline.h>
#include <SPIntegral.h>
#include <SPProcFunction.h>
#include <SPApprox.h> 
#include <SPFuncMap.h>
#include <SPRange.h> 
#include <SPOptions.h>
#include <SPPlot.h> 
#include <SPPlotOptions.h> 
#include <SPH3.h>

typedef struct ScaleOptions_t /* Command line arguments for each scale. */
  { char *basisName;        /* Prefix for basis files of this scale. */
    char *matName;          /* Prefix for basis matrices (default {basisName}). */
    char *transfName;       /* Prefix for scale transf mats {(i-1)->i} ("" = compute). */
    char *outName;          /* Prefix for scale-specific output files. */
    double minVal;          /* Threshold for matrix factor cleanup. */ 
    SPSys_LinOptions_t lso; /* Parameters for linear system solving */
    SPSys_GenOptions_t gso; /* Parameters for solving the non-linear system. */
    double stopError;       /* Stop when this close to true solution. */
  } ScaleOptions_t;
  
vec_typedef(ScaleOptions_vec_t,ScaleOptions_vec,ScaleOptions_t);
  
#define PROG_USAGE \
  PROG_NAME " \\\n" \
  "  -rhsName NAME -coeff NUM -solName NAME \\\n" \
  "  { -scale NN \\\n" \
  "      -basisName NAME [ -matName NAME ] \\\n" \
  "      [ -transfMatName NAME ] \\\n" \
  "      SYSOPTIONS \\\n" \
  "      [ -stopError NUM ] \\\n" \
  "      -outName NAME \\\n" \
  "  }... \\\n" \
  "  [ -smpOrder NUM ] \\\n" \
  "  [ -plotAll ] [ -plotFinal ] [ -verbose ] \\\n" \
  SPPlotOptions_FunctionHelp " \n" \
  "\n" \
  "where each SYSOPTIONS is:\n" \
  "\n" \
  "  [ -minVal NUM ] \\\n" \
  SPSys_LinOptionsHelp "\\\n" \
  SPSys_GenOptionsHelp "\n"
      
typedef struct Options 
  { char *rhsName;          /* Name of right-hand-side operator {FMap}. */
    double coeff;           /* Coefficient {c} of Helmholtz equation. */
    char *solName;          /* Name of true solution. */
    ScaleOptions_vec_t so;  /* Scale-specific options. */
    int smpOrder;           /* Triangle sampling order for spherical integrals. */
    bool_t verbose;         /* TRUE prints right-hand-side vectors, etc.. */
    /* Plotting options: */
    bool_t plotFinal;       /* TRUE to plot final solution and error. */
    bool_t plotAll;         /* TRUE to plot every iteration and error. */
    SPPlotOptions_t plt;    /* Plot format and style options. */
  } Options;
  
typedef struct ScaleData_t   /* Data for each scale. */
  { Basis bas;          /* A basis for the approximation space {V} on {tri}. */
    Triangulation *tri; /* The triangular mesh, or NULL if none. */
    SPMatrix G;         /* Dot product matrix for {bas}. */
    SPMatrix GL;        /* Left factor (Cholesky or SVD) of matrix {G}, if needed. */
    SPVector GD;        /* Middle factor (SVD) of matrix {G}, if needed. */
    SPMatrix GR;        /* Right factor (SVD) of {G}, if needed. */
    SPMatrix H;         /* Differential operator matrix for {bas}, if needed. */
    SPMatrix HL;        /* Left factor (Cholesky or SVD) of {H}, if needed. */
    SPVector HD;        /* Middle factor (SVD) of matrix {G}, if needed. */
    SPMatrix HR;        /* Right factor (SVD) of {H}, if needed. */
    SPMatrix F;         /* Matrix connecting coarser {bas} with this {bas}. */
    /* Current solution: */   
    SPVector a;         /* Coordinates of current solution {g} relative to {bas}. */
    SPFunction *g;      /* Current current solution as a spherical function. */
    /* Working vectors: */  
    SPVector aNew;      /* New solution. */ 
    SPVector b;         /* Right-hand side. */
    SPVector y;       
    SPVector z;     
  } ScaleData_t;

/* Arrays of {ScaleData_t} */

vec_typedef(ScaleData_vec_t,ScaleData_vec,ScaleData_t);

Options *GetOptions(int argn, char **argc);
SPFuncMap GetRHSOperator(char *rhsName);
SPFunction *GetTrueSolution(char *solName);
ScaleData_vec_t ReadScaleData_vec_t(Options *o, double coeff);
ScaleData_t ReadScaleData(ScaleOptions_t *so, double coeff, int r);

int main(int argn, char **argc)
  { Options *o = GetOptions(argn, argc);
    SPPlotOptions_t *po = &(o->plt);
    SPIntegral_SetDefaultSamplingOrder(o->smpOrder);
    SPFuncMap RHS = SPFuncMap_FromName(o->rhsName); /* Right-hand side */
    ScaleData_vec_t m = ReadScaleData_vec_t(o, o->coeff);
    SPFunction *s = GetTrueSolution(o->solName);
    
    auto double eval_s(S2Point *p);
    double eval_s(S2Point *p) { return s->m->eval(s, p); }
    
    /* Fix perspective parameters, if needed: */
    SPPlot_FixView(&(po->obs), 0.0, &(po->upp), NULL, 0.0, NULL, TRUE);

    /* Number of caption lines: */
    int nCap = (po->eps ? 0 : 1);

    /* Compute mesh size in sphere units: */
    double relMeshSize = po->meshSize/(po->figSize/2);

    int r;
    for (r = 0; r < m.ne; r++)
      { ScaleData_t mr = m.e[r];
        ScaleOptions_t *sor = &(o->so.e[r]);
        SPSys_LinOptions_t *lso = &(sor->lso);
        SPSys_GenOptions_t *gso = &(sor->gso);
        Triangulation *tri = SPSpline_BasisTriangulation(mr.bas);
        SPVector a = mr.a;
        SPVector aNew = mr.aNew;
        SPVector b = mr.b;
        int dim = mr.bas.ne;

        double gMax, sMax, eMax, eAvg;

        /* Open the figure stream: */
        SPPlot_Stream *fps = SPPlot_NewStream
          (po->eps, sor->outName, po->paperSize, po->figSize, po->projLonLat, nCap);

        FILE *errWr = open_write(txtcat(sor->outName, ".erp"), TRUE);

        fprintf(stderr, "=== begin scale %d ===\n", r);
        affirm(tri != NULL, "no triangulation");
        fprintf(errWr, "# %6s %6s %16s %16s\n",  "scale", "iter", "rmsError", "maxError");
        if (r == 0)
          { fprintf(stderr, "guessing initial approximation...\n"); 
            SPSys_GuessSol(a);
          }
        else
          { /* Project coarser solution onto this space: */
            SPVector aLo = m.e[r-1].a;  /* Coarse solution. */
            fprintf(stderr, "projecting solution (scale %d -> %d)...\n", r-1, r); 
            SPMatrix_MulRow(aLo, mr.F, mr.z); 
            SPMatrix_DivCol(mr.GL, mr.z, mr.y); 
            SPMatrix_DivRow(mr.y, mr.GL, a); 
          }

        double aDiff = INFINITY; /* Max change in coeffs between iterations. */
        int iter = 0;
        while (TRUE)
          { char *iterTag = fmt_int(iter, 6);
            char *iterName = txtcat3(sor->outName, "-", iterTag);
            fprintf(stderr, "=== scale %d - iteration %d ===\n", r, iter); 
            fprintf(stderr, "building approximation...\n"); 
            mr.g = SPApprox_BuildFunction(mr.bas, a, iterName, s);
            fprintf(stderr, "summary for scale %d iteration %d:\n", r, iter); 
            SPApprox_PrintMaxErrorValues(mr.g, s, tri, &gMax, &sMax, &eMax, &eAvg);
            fprintf(errWr, "  %6d %6d %16.12f\n",  r, iter, eAvg);

            double aNorm = rn_L_inf_norm(dim, a.e);
            if (SPSys_GenStopCondition(gso, sor->stopError, aNorm, iter, eMax, aDiff, TRUE))
              { break; }

            if (o->plotAll)
              { if (po->radius < 1.0) { fprintf(stderr, "** \"-radius\" ignored\n"); }
                SPApprox_PlotFunctionAndError
                  ( fps, iterTag, mr.g, eval_s, 
                    mr.tri, relMeshSize, (mr.tri != NULL), FALSE,
                    (gMax > 2.0*sMax ? gMax : sMax), eMax, 
                    &(po->obs), &(po->upp),
                    &(po->caption), 
                    iter, 0.0
                  );
              }

            fprintf(stderr, "computing right-hand side...\n"); 
            SPApprox_ComputeSystemRightHandSide
              ( mr.g, RHS, mr.bas, NULL, mr.tri, b, o->verbose );

            fprintf(stderr, "solving system...\n"); 
            SPSys_LinSolve(mr.H, mr.HL, mr.HD, mr.HR, b, lso, mr.y, aNew, TRUE);
    
            fprintf(stderr, "updating solution coefficients...\n"); 
            aDiff = SPSys_UpdateSolution(aNew, a);
            free(iterTag); free(iterName);
            iter++;
          }
        fclose(errWr);
        if ((o->plotFinal) || (o->plotAll))
          { char *finTag = SPPlot_IterationTag(INT_MAX);
            if (po->radius < 1.0) { fprintf(stderr, "** \"-radius\" ignored\n"); }
            SPApprox_PlotFunctionAndError
              ( fps, finTag, mr.g, eval_s, 
                mr.tri, relMeshSize, (mr.tri != NULL), TRUE,
                (gMax > 2.0*sMax ? gMax : sMax), eMax, 
                &(po->obs), &(po->upp),
                &(po->caption), 
                iter, 0.0
              );
          }
        fprintf(stderr, "=== end scale %d ===\n\n", r);
        pswr_close_stream(fps);
      }

    return 0;
  }

ScaleData_vec_t ReadScaleData_vec_t(Options *o, double coeff)
  { int r;
    int NS = o->so.ne;
    ScaleData_vec_t m = ScaleData_vec_new(NS);
    for (r = 0; r < NS; r++)
      { m.e[r] = ReadScaleData(&(o->so.e[r]), coeff, r); }
    return m;
  }

ScaleData_t ReadScaleData(ScaleOptions_t *so, double coeff, int r)
  { ScaleData_t SD;  
    Basis bas = SPApprox_ReadBasis(so->basisName);
    int dim = bas.ne;
    SD.bas = bas;
    SPSys_LinOptions_t *lso = &(so->lso);
    SPApprox_GetBasisMatrices
      ( so->matName, lso->mth, &(SD.G), &(SD.GL), &(SD.GD), &(SD.GR) );
    SPApprox_GetHelmholtzMatrices
      ( so->matName, coeff, lso->mth, so->minVal, &(SD.H), &(SD.HL), &(SD.HD), &(SD.HR) );
    if (r > 0) { SD.F = SPApprox_GetTransferMatrix(so->transfName); }
    SD.a = double_vec_new(dim);
    SD.aNew = double_vec_new(dim);
    SD.b = double_vec_new(dim);
    SD.y = double_vec_new(dim);
    SD.z = double_vec_new(dim);
    return SD;
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
    int r;

    SPOptions_SetUsage(pp, PROG_USAGE);

    SPOptions_GetKeyword(pp, "-rhsName");                               
    o->rhsName = SPOptions_GetNext(pp);  

    SPOptions_GetKeyword(pp, "-coeff");                               
    o->coeff = SPOptions_GetNextDouble(pp, -DBL_MAX, +DBL_MAX);  

    SPOptions_GetKeyword(pp, "-solName");                               
    o->solName = SPOptions_GetNext(pp);  

    r = 0;
    while (SPOptions_TestKeywordNext(pp, "-scale"))
      { ScaleOptions_t *sor;
        ScaleOptions_vec_expand(&(o->so), r);
        sor = &(o->so.e[r]); 

        SPOptions_GetKeywordNext(pp, "-basisName"); 
        sor->basisName = SPOptions_GetNext(pp);

        sor->matName = NULL;
        sor->transfName = NULL;
        sor->outName = NULL;
        sor->lso.mth = SPSys_LM_NONE;
        sor->gso.mth = SPSys_GM_NONE;
        while (TRUE)
          { 
            if (SPOptions_TestKeywordNext(pp, "-matName"))                               
              { sor->matName = SPOptions_GetNext(pp); }
            else if (SPOptions_TestKeywordNext(pp, "-transfName"))
              { sor->transfName = SPOptions_GetNext(pp);
                if (r == 0) 
                  { SPOptions_Error(pp, "\"-transfName\" not valid for scale 0"); }
              }
            else if (SPOptions_TestKeywordNext(pp, "-outName")) 
              { sor->outName = SPOptions_GetNext(pp); }  
            else if (SPOptions_TestKeywordNext(pp, "-minVal")) 
              { sor->minVal = SPOptions_GetNextDouble(pp, 0.0, DBL_MAX); }  
            else if 
              (SPOptions_IsNext(pp, "-linearSys"))
              { sor->lso = SPSys_LinOptionsParse(pp, TRUE); }
            else if (SPOptions_IsNext(pp, "-nonLinearSys"))
              { sor->gso = SPSys_GenOptionsParse(pp, TRUE); }
            else
              { /* Parsed everything for this scale: */
                break;
              }
          }
        if (sor->matName == NULL) 
          { SPOptions_Error(pp, "missing \"-matName\"");  }
        if (sor->transfName == NULL) 
          { SPOptions_Error(pp, "missing \"-transfName\"");  }
        if (sor->outName == NULL) 
          { SPOptions_Error(pp, "missing \"-outName\"");  }
        if (sor->lso.mth == SPSys_LM_NONE) 
          { SPOptions_Error(pp, "must specify the linear solution method");  }
        r++;
      }
    ScaleOptions_vec_trim(&(o->so), r);
    
    o->verbose = SPOptions_TestKeyword(pp, "-verbose");
    
    SPOptions_GetKeyword(pp, "-smpOrder");
    o->smpOrder = SPOptions_GetNextInt(pp, 1, INT_MAX);
        
    /* Plotting options: */

    o->plotAll = SPOptions_TestKeyword(pp, "-plotAll");

    o->plotFinal = SPOptions_TestKeyword(pp, "-plotFinal");

    o->plt = SPPlotOptions_FunctionParse(pp);
    
    SPOptions_Finish(pp);
    return o;
  }

/* Arrays of ScaleData_t: */

vec_typeimpl(ScaleData_vec_t,ScaleData_vec,ScaleData_t);

/* Arrays of ScaleOptions_t: */

vec_typeimpl(ScaleOptions_vec_t,ScaleOptions_vec,ScaleOptions_t);
