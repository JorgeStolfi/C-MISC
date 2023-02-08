/* See SOApprox.h */
/* Last edited on 2006-03-19 12:41:40 by stolfi */

#include <SOApprox.h>

#include <SOBasic.h>
#include <SOGrid.h>
#include <SOFunction.h>
#include <SOMatrix.h>
#include <SOBasisMatrices.h>
#include <SOLinCombFunction.h>
#include <SOPlotParams.h>
#include <SO2DPlot.h>
#include <SOPlot.h>

#include <dg_grid.h>

#include <vec.h>
#include <rn.h>
#include <affirm.h>
#include <nat.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <values.h>

SOGrid_Tree *SOApprox_ReadTree(char *fileName)
  { if (strcmp(fileName, "") == 0)
      { return NULL; }
    else
      { FILE *rd = open_read(fileName, TRUE);
        SOGrid_Tree * tree = SOGrid_Tree_read(rd);
        if (rd != stdin) { fclose(rd); }
        return tree;
      }
  } 

SOFunction *SOApprox_ReadFunction(char *fileName)
  { if (strcmp(fileName, "") == 0)
      { return NULL; }
    else
      { FILE *rd = open_read(fileName, TRUE);
        SOFunction *f = SOFunction_Read(rd);
        if (rd != stdin) { fclose(rd); }
        return f;
      }
  }

Basis SOApprox_ReadBasis(char *fileName)
  { FILE *rd = open_read(fileName, TRUE);
    Basis bas = SOFunction_ReadBasis(rd);
    if (rd != stdin) { fclose(rd); }
    return bas;
  }

SOMatrix SOApprox_ReadMatrix(char *fileName)
  { FILE *rd = open_read(fileName, TRUE);
    SOMatrix M = SOMatrix_Read(rd);
    if (rd != stdin) { fclose(rd); }
    return M;
  }

SOMatrix SOApprox_GetBasisMatrix(char *matName, bool_t cholesky)
  { if (cholesky)
      { return SOApprox_ReadMatrix(txtcat(matName, "-ev-ch.mat")); }
    else
      { return SOApprox_ReadMatrix(txtcat(matName, "-ev.mat")); }
  }
  
SOMatrix SOApprox_GetTransferMatrix(char *matName)
  { return SOApprox_ReadMatrix(txtcat(matName, "-ev.mat")); }
  
SOMatrix SOApprox_GetHelmholtzMatrix(char *matName, double c, bool_t cholesky)
  { SOMatrix EM = SOApprox_ReadMatrix(txtcat(matName, "-ev.mat"));
    SOMatrix GM = SOApprox_ReadMatrix(txtcat(matName, "-gr.mat"));
    SOMatrix H;
    fprintf(stderr, "computing the differential operator matrix...\n");
    H = SOMatrix_Mix(1.0, GM, c, EM);
    if (cholesky)
      { SOMatrix HL;
        fprintf(stderr, "computing the Cholesky factorization...\n");
        HL = SOMatrix_Cholesky(H, 0.0); 
        free(H.ents.e);
        return HL;
      }
    else
      { return H; }
  }

void SOApprox_ComputeRightHandSide
  ( SOFunction *f,      /* Current solution. */
    FuncMap *FMap,      /* Right-hand side of differential eqn. */
    Basis bas,          /* Basis for the approximation space. */
    SOFunction *w,      /* Weight for integrals. */
    SOGrid_Tree *tree,  /* Grid to use for integration. */
    double_vec_t y,     /* IN/OUT: right-hand-side vector. */ 
    bool_t verbose        /* TRUE mumbles and prints the change in {y}. */
  )
  { nat_t pDim = bas.e[0]->d->pDim;
    nat_t fDim = bas.e[0]->d->fDim;
    FuncMap NoMap = IdentityMap(fDim, pDim);

    auto void fDot(SOFunction *h, SOFunction *bi, double *iv);

    void fDot(SOFunction *h, SOFunction *bi, double *iv)
      { SOFunction_Dot(h, FMap, bi, &NoMap, w, tree, iv, FALSE); }
    
    SOBasisMatrices_RecomputeVectorGen(f, bas, fDot, y, verbose);
  }

void SOApprox_CholeskySolve
  ( SOMatrix L,       /* Cholesky factor of system's matrix. */          
    int n,            /* Number of columns in soln and rhs matrices. */
    double_vec_t y,   /* Right-hand-side of system. */ 
    double_vec_t x,   /* (OUT) Coordinates of new solution relative to {bas}. */
    double_vec_t r,   /* (WORK) temporary storage area. */
    double_vec_t s,   /* (WORK) temporary storage area. */
    double_vec_t t,   /* (WORK) temporary storage area. */
    bool_t verbose      /* TRUE mumbles along the way. */
  )
  { int i, j;
    if (verbose) { fprintf(stderr, "SOApprox_CholeskySolve\n"); }
    for (j = 0; j < n; j++)
      { for (i = 0; i < L.cols; i++) { r.e[i] = y.e[i*n + j]; } 
        SOMatrix_DivCol(L, r, t);
        SOMatrix_DivRow(t, L, s);
        for (i = 0; i < L.cols; i++) { x.e[i*n + j] = s.e[i]; }
      }
  }

void SOApprox_GaussSeidelIteration(
    SOMatrix A,         /* System's matrix. */          
    int n,              /* Number of columns in soln and rhs matrices. */
    double_vec_t y,     /* Right-hand-side of system. */ 
    double omega,       /* Overrelaxation factor. */
    double_vec_t xOld,  /* (IN) Current guess at the solution. */
    double_vec_t xNew,  /* (OUT) New solution. */
    double_vec_t r,     /* (WORK) temporary storage area. */
    double_vec_t s,     /* (WORK) temporary storage area. */
    bool_t verbose        /* TRUE mumbles along the way. */
  )
  { int i, j;
    if (verbose) { fprintf(stderr, "SOApprox_GaussSeidelIteration\n"); }
    for (i = 0; i < xNew.ne; i++) { xNew.e[i] = xOld.e[i]; }
    for (j = 0; j < n; j++)
      { for (i = 0; i < A.cols; i++) 
          { r.e[i] = xNew.e[i*n + j]; 
            s.e[i] = y.e[i*n + j]; 
          } 
        SOMatrix_GaussSeidel(A, s, omega, r);
        for (i = 0; i < A.cols; i++) 
          { xNew.e[i*n + j] = r.e[i]; } 
      }
  }
  
void SOApprox_GaussSeidelSolve
  ( SOMatrix A,         /* System's matrix. */          
    int n,              /* Number of columns in soln and rhs matrices. */
    double_vec_t y,     /* Right-hand-side of system. */ 
    double omega,       /* Relaxation factor. */
    int maxIter,        /* Maximum number of iterations. */
    double absTol,      /* Stopping criterion: absolute change in solution. */
    double relTol,      /* Stopping criterion: relative change in solution. */
    double_vec_t xOld,  /* (IN) Current solution. */
    double_vec_t xNew,  /* (OUT) New solution. */
    double_vec_t r,     /* (WORK) temporary storage area. */
    double_vec_t s,     /* (WORK) temporary storage area. */
    bool_t verbose      /* TRUE mumbles along the way. */
  )
  { nat_t N = xOld.ne;
    int iter = 0;
    double d, m;
    double_vec_t x = double_vec_new(N);
    int i, j;
    if (verbose) { fprintf(stderr, "SOApprox_GaussSeidelSolve\n"); }
    for (i = 0; i < N; i++) { x.e[i] = xOld.e[i]; }
    while (iter < maxIter)
      { for (i = 0; i < N; i++) { xNew.e[i] = x.e[i]; }
        for (j = 0; j < n; j++)
          { for (i = 0; i < A.cols; i++) 
              { r.e[i] = xNew.e[i*n + j]; 
                s.e[i] = y.e[i*n + j]; 
              } 
            SOMatrix_GaussSeidel(A, s, omega, r);
            for (i = 0; i < A.cols; i++) 
              { xNew.e[i*n + j] = r.e[i]; } 
          }
        iter++;
        d = SOApprox_UpdateSolution(xNew, x);
        m = rn_norm(N, x.e);
        if ((d <= absTol) || (d <= relTol*m)) { break; }
      }
    free(x.e);
  }

void SOApprox_ConjugateGradientSolve
  ( SOMatrix A,     /* System's matrix. */          
    int n,            /* Number of columns in soln and rhs matrices. */
    double_vec_t y,     /* Right-hand-side of system. */ 
    double_vec_t xOld,  /* (IN) Current solution. */
    double_vec_t xNew,  /* (OUT) New solution. */
    double_vec_t r,     /* (WORK) temporary storage area. */
    double_vec_t s,     /* (WORK) temporary storage area. */
    bool_t verbose      /* TRUE mumbles along the way. */
  )
  { double rsq;
    int i, j;
    if (verbose) { fprintf(stderr, "SOApprox_ConjugateGradientSolve\n"); }
    for (i = 0; i < xNew.ne; i++) { xNew.e[i] = xOld.e[i]; }
    rsq = 0;
    for (j = 0; j < n; j++)
      { for (i = 0; i < A.cols; i++) 
          { r.e[i] = xNew.e[i*n + j]; 
            s.e[i] = y.e[i*n + j]; 
          } 
        rsq += SOMatrix_ConjugateGradient(A, s, r);
        for (i = 0; i < A.cols; i++) 
          { xNew.e[i*n + j] = r.e[i]; } 
      }
    if (verbose) { fprintf(stderr, "  residual norm = %16.12f\n", rsq); }
  }
  
double SOApprox_UpdateSolution
  ( double_vec_t xNew, /* (IN) New solution vector. */
    double_vec_t x     /* (I/O)Current solution vector. */
  )
  { int i;
    double d = rn_dist(x.ne, x.e, xNew.e);
    fprintf(stderr, "change = %.12g\n", d);
    for (i = 0; i < xNew.ne; i++) { x.e[i] = xNew.e[i]; }
    x = xNew;
    return d;
  }

void SOApprox_GuessSol(double_vec_t x)
  { int i;
    for (i = 0; i < x.ne; i++){ x.e[i] = 0.0; }
    x.e[0] = 1.0;
  }
  
SOFunction *SOApprox_BuildFunction
  ( Basis bas,      /* Basis for space {V[r]}. */
    double_vec_t u,     /* Coordinates of {g} relative to  {bas} */
    char *solName,    /* Name prefix for output file. */
    char *basName,    /* Name of the (Basis bas) description file */
    SOFunction *f     /* True solution (for comparisons) */
  )
  { SOLinCombFunction *g = 
      SOLinCombFunction_Make(f->d->pDim, f->d->fDim, basName, bas, u);

    FILE *wr = open_write(txtcat(solName, "-app.fun"), TRUE);
    g->m->fn.write(g, wr);  /* Writes with SOFunction encapsulation */
    fclose(wr);
    //    SOApprox_PrintSomeValues((SOFunction *)g, f);
    return (SOFunction *)g;
  }
  
void SOApprox_PrintSomeValues(SOFunction *g, SOFunction *f)
  { int m = f->d->pDim;
    int n = f->d->fDim;
    int i, j, k;
    double x[m], fx[n], gx[n];
    int N = 20;            // Sample poins per axis.
    int NP = ipow(N, m);  // Total number of sample points. 
    
    affirm(g->d->pDim == m, "incompatible domain dimensions");
    affirm(g->d->fDim == n, "incompatible range dimensions");
    
    for (i = 0; i < NP; i++)
      { // Generate coordinates of sample point number {i}:
        int ti = i;
        for (k = m-1; k >= 0; k--)
          { x[k] = ((double)(ti % N))/((double)(N-1));
            ti /= N;
          }
 	
        fprintf(stderr, "x =       (");
        for (k = 0; k < m; k++) { fprintf(stderr, " %5.3f", x[k]); }
	fprintf(stderr, " )\n");
            
        // Evaluate
    
        g->m->eval(g, x, gx);
        fprintf(stderr, " fAppr(x) =      (");
        for (j = 0; j < n; j++) { fprintf(stderr, " %16.12f", gx[j]); }
	fprintf(stderr, " )\n");
            
        f->m->eval(f, x, fx);
        fprintf(stderr, " f(x) =          (");
        for (j = 0; j < n; j++) { fprintf(stderr, " %16.12f", fx[j]); }
	fprintf(stderr, " )\n");

        fprintf(stderr, " fAppr(x)-fCorr(x) = (");
        for (j = 0; j < n; j++)
          { double ej = gx[j] - fx[j];
            fprintf(stderr, " %16.12f", ej);
          }
        fprintf(stderr, " )\n");
      }
  }
  
void SOApprox_PrintMaxErrorValues
  ( SOFunction *g,  /* Approximation. */
    SOFunction *f,  /* Target function. */
    double *gMax,     /* OUT: max abs function value. */
    double *eMax      /* OUT: max abs error VALUE. */
  )
  {  
    int m = f->d->pDim;
    int n = f->d->fDim;
    int steps = 1, divs = 10, i, j;
    double point[m], inc = (double)1 / divs;
    double gv[n], fv[n];
    bool_t greater;

    for(i = 0; i < m; i++){steps *= divs; point[i] = 0;}
    for(j = 0; j < n; j++){gMax[j] = 0; eMax[j] = 0;}

    for(i = 0; i < steps; i++) 
      {
        point[i % m] += inc;
        g->m->eval(g, point, &(gv[0]));
        f->m->eval(f, point, &(fv[0]));
        
        greater = FALSE;
        for(j = 0; j < n; j++) if(gv[j] > gMax[j]) greater = TRUE;
        if(greater) for(j = 0; j < n; j++) gMax[j] = gv[j];

        greater = FALSE;
        for(j = 0; j < n; j++) if(fabs(gv[j]-fv[j]) > eMax[j]) greater = TRUE;
        if(greater) for(j = 0; j < n; j++) eMax[j] = fabs(gv[j]-fv[j]);
      }

    for(j = 0; j < n; j++)
      {
        fprintf(stderr, "max(fAppr)[%d] = %16.12f\n", j, gMax[j]);
        fprintf(stderr, "max(fabs(fAppr-fCorr))[%d] = %16.12f\n", j, eMax[j]);
      }
  }
  
void SOApprox_PlotFunction
  ( SOFunction *f, 
    SOGrid_Tree *tree,
    char *fName, 
    PlotOptions *plt,
    double *fObsMin,    /* (IN/OUT) Min value seen during plot. */ 
    double *fObsMax     /* (IN/OUT) Max value seen during plot. */ 
  )
  { 
    double hFigSize = plt->figSize;             /* In mm. */
    double vFigSize = plt->figSize * SQRTHALF;  /* In mm. */ 
    double figMargin = (plt->eps ? 0.5 : 2.0);  /* In mm. */
    interval_t boxR[2];     /* Rectangle to plot. */
    int captionLines = (plt->eps ? 0 : 1);      /* No caption. */
    int hCount = 1, vCount = 1;  /* Figures per row and column. */
    int meshDepth = (int)ceil(2*log(plt->figSize/plt->meshSize)/log(2));
    SOPlot_PageState *pgs = SOPlot_PageState_New
      ( plt->eps, 
        plt->paperSize, hFigSize, vFigSize, figMargin, figMargin,
        captionLines, hCount, vCount
      );
    LO(boxR[X]) = 0.0;  HI(boxR[X]) = 1.0;
    LO(boxR[Y]) = 0.0;  HI(boxR[Y]) = 1.0;
    /* Perhaps here we should fix the range according to {plt->autoRange}? */
    SOPlot_File *fps = SO2DPlot_FunctionFigure
      ( NULL, f, tree,   
        boxR, 
        -plt->fRange, +plt->fRange, plt->fStep,
        plt->isolineWidth, plt->gridWidth, meshDepth,
        (! plt->noBands), (! plt->noIsolines), (! plt->noGrid), 
        pgs, fName, fName,
        fObsMin, fObsMax
      );
    SOPlot_Close(fps, pgs);
  }

/* SOPlot_File *SOApprox_PlotFunctionAndError */
/*   ( SOPlot_File *fps,        /\* Postscript file, or NULL. *\/ */
/*     SOPlot_PageState *pgs,   /\* Page layout and state for {fps}. *\/ */
/*     char *docName,           /\* Document file name, minus extension. *\/ */
/*     char *funcName,          /\* Prefix for separate figure file names. *\/ */
/*     SOFunction *g,         /\* The computed approximation. *\/ */
/*     SOFunction *f,         /\* The target function. *\/ */
/*     SOGrid_Tree *tree,      /\* The grid. *\/ */
/*     double relMeshSize,      /\* Maximum step/triangle size (radians). *\/ */
/*     bool_t showTriang,         /\* TRUE to plot the reference grid. *\/  */
/*     double gMax,             /\* Nominal maximum of {fabs(g(p))}. *\/ */
/*     double eMax,             /\* Nominal maximum of {fabs(g(p)-f(p))}. *\/ */
/*     SOH3_Point *obs,          /\* Viewpoint. *\/ */
/*     SOH3_Point *upp           /\* Camera vertical reference *\/ */
/*   ) */
/*   { auto double gValue(double *p); */
/*     auto double gfError(double *p); */
    
/*     double gValue(double *p) */
/*       { return g->m->eval(g, p); } */

/*     double gfError(double *p) */
/*       { return g->m->eval(g, p) - f->m->eval(f, p); } */

/*     double gMinObs, gMaxObs; */
/*     double eMinObs, eMaxObs; */
/*     double gStep, eStep; */
    
/*     double rad = 1.0; */
/*     SOH3_Point obsFix = (*obs), uppFix = (*upp); */
    
/*     r3_t dLight = (r3_t){{0.0, 0.0, 1.0}}; */
/*     SOH3_Plane gSupp = SOFunction_GetSupportingPlane(g); */
    
/*     char *appCapLine = "Approximation"; */
/*     char *errCapLine = "Approximation error"; */
/*     string_vec_t appCaption = (string_vec_t){ 1, &appCapLine }; */
/*     string_vec_t errCaption = (string_vec_t){ 1, &errCapLine }; */
  
/*     SOPlot_FixView(&(obsFix), 0.0, &(uppFix), &rad, 0.0, g, TRUE); */
    
/*     gStep = SOPlot_RoundToNice(gMax/DefaultIsolines); */
/*     gMax = gStep * DefaultIsolines; */
/*     fprintf(stderr, "plot gMax = %16.12f gStep = %16.12f\n", gMax, gStep); */

/*     r4_gen_print(stderr, &(gSupp.f), "%6.3f", "support = [ ", ", ", " ]\n"); */
/*     gMinObs =   INFTY; */
/*     gMaxObs = - INFTY; */
/*     fps = SOPlot_BothSides */
/*       ( fps, pgs, */
/*         /\* docName *\/ docName, */
/*         /\* funcName *\/ txtcat(funcName, "-app"), */
/*         /\* func *\/ gValue, */
/*         /\* tree *\/ tree, */
/*         /\* relMeshSize *\/ relMeshSize,  */
/*         /\* showTriang *\/ showTriang, */
/*         /\* supp *\/ &gSupp, */
/*         /\* fRange *\/ gMax, */
/*         /\* fStep *\/ gStep, */
/*         /\* obs *\/ &obsFix,  */
/*         /\* upp *\/ &uppFix, */
/*         /\* dLight *\/ &dLight, */
/*         /\* lineWidth *\/ 0.25, */
/*         /\* caption *\/ &appCaption, */
/*         /\* verbose *\/ TRUE, */
/*         /\* fMinObs *\/ &gMinObs, */
/*         /\* fMaxObs *\/ &gMaxObs */
/*       ); */
/*     fprintf(stderr, "observed function extrema:"); */
/*     fprintf(stderr, " min = %16.12f max = %16.12f\n", gMinObs, gMaxObs); */
/*     fprintf(stderr, "\n"); */
       
/*     eStep = SOPlot_RoundToNice(eMax/DefaultIsolines); */
/*     eMax = eStep * DefaultIsolines; */
/*     fprintf(stderr, "plot eMax = %16.12f eStep = %16.12f\n", eMax, eStep); */
    
/*     eMinObs = INFTY; */
/*     eMaxObs = - INFTY; */
/*     fps = SOPlot_BothSides */
/*       ( fps, pgs, */
/*         /\* docName *\/ docName, */
/*         /\* funcName *\/ txtcat(funcName, "-err"), */
/*         /\* func *\/ gfError, */
/*         /\* tree *\/ tree, */
/*         /\* relMeshSize *\/ relMeshSize,  */
/*         /\* showTriang *\/ showTriang, */
/*         /\* supp *\/ NULL, */
/*         /\* fRange *\/ eMax, */
/*         /\* fStep *\/ eStep, */
/*         /\* obs *\/ &obsFix,  */
/*         /\* upp *\/ &uppFix, */
/*         /\* dLight *\/ &dLight, */
/*         /\* lineWidth *\/ 0.25, */
/*         /\* caption *\/ &errCaption, */
/*         /\* verbose *\/ TRUE, */
/*         /\* fMinObs *\/ &eMinObs, */
/*         /\* fMaxObs *\/ &eMaxObs */
/*       ); */
/*     fprintf(stderr, "observed error extrema:"); */
/*     fprintf(stderr, " min = %16.12f max = %16.12f\n", eMinObs, eMaxObs); */
/*     fprintf(stderr, "\n"); */
/*     return fps; */
/*   } */

