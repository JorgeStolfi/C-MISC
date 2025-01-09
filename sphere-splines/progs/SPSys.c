/* See SPSys.h */
/* Last edited on 2008-05-24 12:30:25 by stolfi */

#define _GNU_SOURCE
#include <SPSys.h>

#include <SPMatrix.h>
#include <SPVector.h>
#include <SPBasic.h>

#include <vec.h>
#include <rn.h>
#include <js.h>
#include <affirm.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <values.h>

char *SPSys_IterativeOptionsText(int maxIter, double relTol, double absTol, double omega);
  /* Returns a new string showing with the given parameters from a 
    {SPSys_LinOptions_t} or {SPSys_GenOptions_t} record. */

void SPSys_ComputeNeededFactors
  ( SPMatrix A,             
    SPSys_LinMethod_t mth,  
    double minVal,          
    SPMatrix *L,            
    SPVector *D,            
    SPMatrix *R             
  )
  { 
    if (mth == SPSys_LM_Cholesky)
      { SPMatrix_Cholesky(A, minVal, L);
        (*D) = double_vec_new(0);
        (*R) = SPMatrix_Null(A.rows,A.cols);
      }
    else if (mth == SPSys_LM_SVD)
      { SPMatrix_SVD(A, minVal, L, D, R); }
    else
      { (*L) = SPMatrix_Null(A.rows,A.cols);
        (*D) = double_vec_new(0);
        (*R) = SPMatrix_Null(A.rows,A.cols);
      }
  }
  
void SPSys_LinSolve
  ( SPMatrix A,
    SPMatrix L,
    SPVector D,
    SPMatrix R,
    SPVector y,
    SPSys_LinOptions_t *o,
    SPVector t,
    SPVector x,
    bool_t verbose
  )
  { 
    switch(o->mth)
      {
        case SPSys_LM_SVD:
          { affirm(t.ne >= D.ne, "vector {t} too short");
            double_vec_t tt = (double_vec_t){ /*ne*/ D.ne, /*e*/ t.e }; 
            SPSys_SVDSolve(L, D, R, y, tt, x, verbose);
          }
          break;
        case SPSys_LM_GaussSeidel:
          { SPSys_GaussSeidelSolve
              ( A, y, 
                o->omega, o->maxIter, o->absTol, o->relTol,
                t, x, verbose
              );
	  }
          break;
        case SPSys_LM_Cholesky:
          { SPSys_CholeskySolve(L, y, t, x, verbose); }
          break;
        case SPSys_LM_ConjGrad:
          { SPSys_ConjugateGradientSolve(A, y, t, x, verbose); }
          break;
        default:
          affirm(FALSE , "no linear system solution method");
      }
  }

void SPSys_LinSolveResidual
  ( SPMatrix A,
    SPMatrix L,
    SPVector D,
    SPMatrix R,
    SPVector dy,
    SPSys_LinOptions_t *o,
    SPVector t,
    SPVector dx,
    bool_t verbose
  )
  { 
    int i;
    switch(o->mth)
      {
        case SPSys_LM_SVD:
          SPSys_SVDSolve(L, D, R, dy, t, dx, verbose);
          break;
        case SPSys_LM_GaussSeidel:
          for (i = 0; i < t.ne; i++) { t.e[i] = 0.0; } 
          SPSys_GaussSeidelSolve(A, dy, o->omega, 2, 0.0, 0.0, t, dx, verbose);
          break;
        case SPSys_LM_Cholesky:
          SPSys_CholeskySolve(L, dy, t, dx, verbose);
          break;
        case SPSys_LM_ConjGrad:
          for (i = 0; i < t.ne; i++) { t.e[i] = 0.0; } 
          SPSys_ConjugateGradientSolve(A, dy, t, dx, verbose);
          break;
        default:
          affirm(FALSE , "no linear system solution method");
      }
  }

void SPSys_AdjustSolution(SPVector x, SPVector dx, bool_t verbose)
  { int i; double dMax = 0.0;
    for (i = 0; i < dx.ne; i++) 
      { double dxi = dx.e[i];
        if (verbose) 
          { fprintf(stderr, "  x[%3d]  old = %16.12f", i, x.e[i]);
            fprintf(stderr, "  correction = %18.15f\n", dxi);
          }
        x.e[i] += dxi;
        if (verbose) 
          { fprintf(stderr, "  new = %16.12f", x.e[i]); }
        dxi = fabs(dxi);
        if (dxi > dMax) { dMax = dxi; }
      }
    fprintf(stderr, "max correction = %16.12f\n", dMax);
  }
  
double SPSys_UpdateSolution
  ( SPVector xNew, 
    SPVector x     
  )
  { int i;
    double d = rn_dist(x.ne, x.e, xNew.e);
    fprintf(stderr, "change = %.12g\n", d);
    for (i = 0; i < xNew.ne; i++) { x.e[i] = xNew.e[i]; }
    x = xNew;
    return d;
  }

void SPSys_GuessSol(SPVector x)
  { int i;
    for (i = 0; i < x.ne; i++) { x.e[i] = 0.0; }
    x.e[0] = 1.0;
  }
  
SPSys_LinOptions_t SPSys_LinOptionsParse(SPOptions_Parser_t *pp, bool_t local)
  {
    SPSys_LinOptions_t o;
    
    o.mth = SPSys_LM_NONE;
    o.omega = 1.0;
    o.maxIter = 1;
    o.relTol = 0.0;
    o.absTol = 0.0;
    o.residual = 0;
    
    if (SPOptions_TestKeywordGen(pp, "-linearSys", local))
      { 
        if (SPOptions_TestKeyword(pp, "SVD"))
          { o.mth = SPSys_LM_SVD; }
        else if (SPOptions_TestKeyword(pp, "gaussSeidel"))
          { o.mth = SPSys_LM_GaussSeidel;
            o.maxIter = SPOptions_GetNextInt(pp, 1, 1000000);
            o.relTol = SPOptions_GetNextDouble(pp, 0.0, 1.0);
            o.absTol = SPOptions_GetNextDouble(pp, 0.0, DBL_MAX);
            if (SPOptions_TestKeywordNext(pp, "-omega"))
              { o.omega = SPOptions_GetNextDouble(pp, 0.0, DBL_MAX); }
          }
        else if (SPOptions_TestKeyword(pp, "cholesky"))
          { o.mth = SPSys_LM_Cholesky; }
        else if (SPOptions_TestKeyword(pp, "conjGrad"))
          { o.mth = SPSys_LM_ConjGrad; }
        else 
          { SPOptions_Error(pp, "bad or missing solution method"); }

        /* Residual iterations: */
        if (SPOptions_TestKeywordNext(pp, "-residual"))
          { o.residual = SPOptions_GetNextInt(pp, 0, INT_MAX); }
      }
    return o;
  }

SPSys_GenOptions_t SPSys_GenOptionsParse(SPOptions_Parser_t *pp, bool_t local)
  {
    SPSys_GenOptions_t o;
    
    o.mth = SPSys_GM_NONE;
    
    o.maxIter = 1;
    o.relTol = 0.0;
    o.absTol = 0.0;

    if (SPOptions_TestKeywordGen(pp, "-nonLinearSys", local))
      { 
        if (SPOptions_TestKeyword(pp, "iterative"))
          { o.mth = SPSys_GM_Iterative;
            o.maxIter = SPOptions_GetNextInt(pp, 1, 1000000);
            o.relTol = SPOptions_GetNextDouble(pp, 0.0, 1.0);
            o.absTol = SPOptions_GetNextDouble(pp, 0.0, DBL_MAX);
          }
        else 
          { SPOptions_Error(pp, "bad or missing solution method"); }
      }
    return o;
  }

/* SPECIALIZED METHODS */

void SPSys_SVDSolve
  ( SPMatrix L,    
    SPVector D,    
    SPMatrix R,    
    SPVector y,    
    SPVector t,    
    SPVector x,    
    bool_t verbose 
  )
  { 
    if (verbose) { fprintf(stderr, "SPSys_SVDSolve\n"); }
  
    int m = L.rows; /* Rows of original matrix. */
    int p = D.ne;  /* Rank of original matrix. */
    int n = R.rows; /* Columns of original matrix. */
    
    affirm(L.cols == p, "L cols mismatch");
    affirm(R.cols == p, "R cols mismatch");

    affirm(y.ne == m, "y size mismatch");
    affirm(x.ne == n, "x size mismatch");
    affirm(t.ne == p, "t size mismatch");

    /* Solve the system {L*t = y}, assuming that {L} is orthonormal: */
    SPMatrix_TrMulCol(L, y, t);
    
    /* Solve the system {D*t' = t}, where {D} is diagonal; leave {t'} in {t}: */
    int i;
    for (i = 0; i < D.ne; i++)
      { double Di = D.e[i];
        if (Di != 0) { Di = 1.0 / Di; }
        t.e[i] *= Di;
      }

    /* Solve the system {R^T*x = t}, assuming {R} is orthonormal: */
    SPMatrix_MulCol(R, t, x);
  }

void SPSys_CholeskySolve
  ( SPMatrix L,   
    SPVector y,   
    SPVector t,   
    SPVector x,   
    bool_t verbose
  )
  { if (verbose) { fprintf(stderr, "SPSys_CholeskySolve\n"); }
    SPMatrix_DivCol(L, y, t);
    SPMatrix_DivRow(t, L, x);
  }

void SPSys_GaussSeidelIteration(
    SPMatrix A,     
    SPVector y,     
    double omega,   
    SPVector xOld,  
    SPVector xNew,  
    bool_t verbose  
  )
  { if (verbose) { fprintf(stderr, "SPSys_GaussSeidelIteration\n"); }
    { int i; for (i = 0; i < xNew.ne; i++) { xNew.e[i] = xOld.e[i]; } }
    SPMatrix_GaussSeidel(A, y, omega, xNew);
  }
  
void SPSys_GaussSeidelSolve
  ( SPMatrix A,     
    SPVector y,     
    double omega,   
    int maxIter,    
    double absTol,  
    double relTol,  
    SPVector xOld,  
    SPVector xNew,  
    bool_t verbose  
  )
  { int N = xOld.ne;
    int iter = 0;
    double d, m;
    SPVector x = double_vec_new(N);
    if (verbose) { fprintf(stderr, "SPSys_GaussSeidelSolve\n"); }
    { int i; for (i = 0; i < N; i++) { x.e[i] = xOld.e[i]; } }
    while (iter < maxIter)
      { { int i; for (i = 0; i < N; i++) { xNew.e[i] = x.e[i]; } }
        SPMatrix_GaussSeidel(A, y, omega, xNew);
        iter++;
        d = SPSys_UpdateSolution(xNew, x);
        m = rn_norm(N, x.e);
        if ((d <= absTol) || (d <= relTol*m)) { break; }
      }
    free(x.e);
  }

void SPSys_ConjugateGradientSolve
  ( SPMatrix A,    
    SPVector y,    
    SPVector xOld, 
    SPVector xNew, 
    bool_t verbose 
  )
  { double rsq;
    if (verbose) { fprintf(stderr, "SPSys_ConjugateGradientSolve\n"); }
    { int i; for (i = 0; i < xNew.ne; i++) { xNew.e[i] = xOld.e[i]; } }
    rsq = SPMatrix_ConjugateGradient(A, y, xNew);
    if (verbose) { fprintf(stderr, "  residual norm = %16.12f\n", rsq); }
  }

bool_t SPSys_GenStopCondition
  ( SPSys_GenOptions_t *gso,
    double stopError,
    double cNorm,
    int nIter,
    double error,
    double cDiff,
    bool_t verbose           /* If TRUE, prints the stopping reason to {stderr}. */
  )
  { if (gso->mth == SPSys_GM_NONE) 
      { if (verbose) 
          { fprintf(stderr, "non-linear solving was not requested\n"); }
        return TRUE;
      }
    if (nIter >= gso->maxIter) 
      { if (verbose) 
          { fprintf(stderr, "too many non-linear interations\n"); }
        return TRUE; 
      }
    if (error <= stopError)
      { if (verbose) 
          { fprintf(stderr, "average error from true sol is small enough\n"); }
        return TRUE; 
      }
    if ((nIter > 0) && ((cDiff <= gso->relTol*cNorm) && (cDiff <= gso->absTol))) 
      { if (verbose) 
          { fprintf(stderr, "non-linear iteration apparently converged\n"); }
        return TRUE; 
      }
    return FALSE;
  }

char *SPSys_LinMethodName(SPSys_LinMethod_t mth)
  {
    switch(mth)
      {
        case SPSys_LM_NONE:
          return "NONE";
        case SPSys_LM_SVD:
          return "SVD";
        case SPSys_LM_GaussSeidel:
          return "gaussSeidel";
        case SPSys_LM_Cholesky:
          return "cholesky";
        case SPSys_LM_ConjGrad:
          return "conjGrad";
        default:
          return "INVALID";
      }
  }

char *SPSys_GenMethodName(SPSys_GenMethod_t mth)
  {
    switch(mth)
      {
        case SPSys_GM_NONE:
          return "NONE";
        case SPSys_GM_Iterative:
          return "iterative";
        default:
          return "INVALID";
      }
  }

char *SPSys_LinOptionsText(SPSys_LinOptions_t *lso)
  {
    char *m = SPSys_LinMethodName(lso->mth);
    switch(lso->mth)
      {
        case SPSys_LM_NONE:
        case SPSys_LM_SVD:
        case SPSys_LM_Cholesky:
        case SPSys_LM_ConjGrad:
          break;
        case SPSys_LM_GaussSeidel:
          t = SPSys_IterativeOptionsText(lso->maxIter, lso->relTol, lso->absTol, lso->omega);
          break;
        default:
          break;
      }
    char *s = jsprintf("%s %s -residual %d", m, (t != NULL ? t : ""), lso->residual);
    if (t != NULL) { free(t); }
    return s;
  }


char *SPSys_GenOptionsText(SPSys_GenOptions_t *gso)
  {
    char *m = SPSys_GenMethodName(gso->mth);
    switch(gso->mth)
      {
        case SPSys_GM_NONE:
          break;
        case SPSys_GM_Iterative:
          t = SPSys_IterativeOptionsText(gso->maxIter, gso->relTol, gso->absTol, /*gso->omega*/ 0.0);
          break;
        default:
          break;
      }
    char *s = jsprintf("%s %s", m, (t != NULL ? t : ""));
    if (t != NULL) { free(t); }
    return s;
  }

char *SPSys_IterativeOptionsText(int maxIter, double relTol, double absTol, double omega)
  {
    if (omega != 0)
      { char *t = jsprintf(" %d %24.16e %24.16e -omega %24.16e", maxIter, relTol, absTol, omega); }
    else
      { char *t = jsprintf(" %d %24.16e %24.16e", maxIter, relTol, absTol); }
    return t;
  }

