/* See SPApprox.h */
/* Last edited on 2008-05-24 12:24:54 by stolfi */

#define _GNU_SOURCE
#include <SPApprox.h>

#include <SPMatrix.h>
#include <SPVector.h>
#include <SPSys.h>
#include <SPFunction.h>
#include <SPSpline.h>
#include <SPBasisMatrices.h>
#include <SPTriang.h>
#include <SPRange.h>
#include <SPIntegral.h>
#include <SPPlot.h>
#include <SPBasic.h>

#include <vec.h>
#include <rn.h>
#include <r3.h>
#include <SPH3.h>
#include <js.h>
#include <filefmt.h>
#include <nget.h>
#include <fget.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <values.h>

/* IMPLEMENTATIONS */

void SPApprox_GetBasisMatrix(char *matName, SPMatrix *E);
  /* Obtains the rigidity matrix {E[i,j] == <B[j]|B[i]>} for some
    spherical function basis {B}. The matrix is read from the file
    "{M}-ev.mat", where {M} is the given {matName}. */

void SPApprox_GetBasisMatrixCholesky(char *matName, SPMatrix *EL);
  /* Similar to {SPApprox_GetBasisMatrix}, but instead of {E} obtains
    its Cholesky lower triangular factor {EL}, which is read
    from file "{M}-ev-choL.mat". */

void SPApprox_GetBasisMatrixSVD(char *matName, SPMatrix *EL, SPVector *ED, SPMatrix *ER);
  /* Similar to {SPApprox_GetBasisMatrix}, but instead of {E} obtains
    the factors {EL,ED,ER} of its SVD decomposition {EL*ED*ER^T} which
    are read from files "{M}-ev-svdL.mat", "{M}-ev-svdD.vec", and
    "{M}-ev-svdR.mat". */

void SPApprox_GetHelmholtzMatrix(char *matName, double c, SPMatrix *H);
  /* Obtains the differential operator matrix {H[i,j]= 
    <SGrd(B[j]) | SGrd(B[i])> + c*<B[j]|B[i]>}. where {SGrd} is the
    spherical gradient operator. The matrices {G[i,j] = <SGrd(B[j]) |
    SGrd(B[i])>} and {E[i,j] = <B[j]|B[i]>} are read from files
    "{M}-gr.mat" and "{M}-ev.mat", where {M = matName}. */

/* IMPLEMENTATIONS */

SPFunction *SPApprox_ReadFunction(char *fileName)
  { 
    FILE *rd = open_read(fileName, TRUE);
    SPFunction *f = SPFunction_Read(rd);
    if (rd != stdin) { fclose(rd); }
    return f;
  }

Basis SPApprox_ReadBasis(char *fileName)
  { 
    FILE *rd = open_read(fileName, TRUE);
    Basis bas = SPFunction_ReadBasis(rd);
    if (rd != stdin) { fclose(rd); }
    return bas;
  }

SPMatrix SPApprox_ReadMatrix(char *fileName)
  { 
    FILE *rd = open_read(fileName, TRUE);
    SPMatrix M = SPMatrix_Read(rd);
    if (rd != stdin) { fclose(rd); }
    return M;
  }

SPVector SPApprox_ReadVector(char *fileName)
  { 
    FILE *rd = open_read(fileName, TRUE);
    SPVector v = SPVector_Read(rd);
    if (rd != stdin) { fclose(rd); }
    return v;
  }

void SPApprox_GetBasisMatrix(char *matName, SPMatrix *E)
  { 
    (*E) = SPApprox_ReadMatrix(txtcat(matName, "-ev.mat"));
  }
  
void SPApprox_GetBasisMatrixCholesky(char *matName, SPMatrix *EL)
  { 
    (*EL) = SPApprox_ReadMatrix(txtcat(matName, "-ev-choL.mat"));
  }

void SPApprox_GetBasisMatrixSVD(char *matName, SPMatrix *EL, SPVector *ED, SPMatrix *ER)
  { 
    (*EL) = SPApprox_ReadMatrix(txtcat(matName, "-ev-svdL.mat"));
    (*ED) = SPApprox_ReadVector(txtcat(matName, "-ev-svdD.vec"));
    (*ER) = SPApprox_ReadMatrix(txtcat(matName, "-ev-svdR.mat"));
  }

void SPApprox_GetBasisMatrices
  ( char *matName,
    SPSys_LinMethod_t mth,
    SPMatrix *E,
    SPMatrix *EL,
    SPVector *ED,
    SPMatrix *ER
  )
  { 
    SPApprox_GetBasisMatrix(matName, E);
    switch(mth)
      {
        case SPSys_LM_SVD:
          SPApprox_GetBasisMatrixSVD(matName, EL, ED, ER);
          break;
        case SPSys_LM_Cholesky:
          SPApprox_GetBasisMatrixCholesky(matName, EL);
          (*ED) = double_vec_new(0);
          (*ER) = SPMatrix_Null(EL->rows,EL->cols);
          break;
        case SPSys_LM_GaussSeidel:
        case SPSys_LM_ConjGrad:
          (*EL) = SPMatrix_Null(E->rows,E->cols);
          (*ED) = double_vec_new(0);
          (*ER) = SPMatrix_Null(E->rows,E->cols);
          break;
        default:
          affirm(FALSE , "no linear system solution method");
      }
  }

void SPApprox_GetHelmholtzMatrix(char *matName, double c, SPMatrix *H)
  { 
    SPMatrix EM = SPApprox_ReadMatrix(txtcat(matName, "-ev.mat"));
    SPMatrix GM = SPApprox_ReadMatrix(txtcat(matName, "-sl.mat"));
    fprintf(stderr, "computing the differential operator matrix...\n");
    (*H) = SPMatrix_Mix(1.0, GM, c, EM);
  }

void SPApprox_GetHelmholtzMatrices
  ( char *matName, 
    double c, 
    SPSys_LinMethod_t mth,
    double minVal,
    SPMatrix *H,
    SPMatrix *HL,
    SPVector *HD,
    SPMatrix *HR
  )
  { 
    SPApprox_GetHelmholtzMatrix(matName, c, H);
    SPSys_ComputeNeededFactors(*H, mth, minVal, HL, HD, HR);
  }

SPMatrix SPApprox_GetTransferMatrix(char *matName)
  { return SPApprox_ReadMatrix(txtcat(matName, "-ev.mat")); }
  
void SPApprox_ComputeSystemRightHandSide
  ( SPFunction *f,      /* Current solution. */
    FuncMap FMap,       /* Right-hand side of differential eqn. */
    Basis bas,          /* Basis for the approximation space. */
    SPFunction *w,      /* Weight for integrals. */
    Triangulation *tri, /* Triangulation to use for integration. */
    SPVector y,         /* IN/OUT: right-hand-side vector. */ 
    bool_t verbose      /* TRUE mumbles and prints the change in {y}. */
  )
  { 
    auto double fDot(SPFunction *h, SPFunction *bi);
    
    double fDot(SPFunction *h, SPFunction *bi)
      { return SPFunction_Dot(h, FMap, bi, NoFMap, w, tri, FALSE); }
    
    SPBasisMatrices_RecomputeVectorGen(f, bas, fDot, y, verbose);
  }

SPFunction *SPApprox_BuildFunction
  ( Basis basis,     /* Basis for space {V[r]}. */
    SPVector u,      /* Coordinates of {g} relative to  {bas} */
    char *solName,   /* Name prefix for output file. */
    SPFunction *f    /* True solution (for comparisons) */
  )
  { 
    /* Build solution: */
    SPFunction *g = SPFunction_LinComb(u.e, basis);
    
    /* Write solution to disk: */
    { char *fileName = NULL; asprintf(&fileName, "%s-app.sfn", solName);
      FILE *wr = open_write(fileName, TRUE);
      g->m->write(g, wr); 
      fclose(wr); free(fileName);
    }
    
    /* Write coefficients to disk: */
    { char *fileName = NULL; asprintf(&fileName, "%s-app.cof", solName);
      FILE *wr = open_write(fileName, TRUE);
      SPVector_Write(wr, u);
      fclose(wr); free(fileName);
    }
    
    /* Print some values to {stderr}, for checking: */
    SPApprox_PrintSomeValues(g, f);
    return g;
  }
  
void SPApprox_PrintSomeValues(SPFunction *g, SPFunction *f)
  { 
    int i, j, k;
    for (i = -2; i <= +2; i++) 
      { for (j = -2; j <= +2; j++) 
          { for (k = -2; k <= +2; k++) 
              { if ((i != 0) || (j != 0) || (k != 0 ))
                  { double x = (double)i;
                    double y = (double)j;
                    double z = (double)k;
                    r3_t p = (r3_t){{x, y, z}};
                    double pr = r3_norm(&p);
                    r3_dir(&p, &p);
                    { double gv = g->m->eval(g, &p);
                      double fv = f->m->eval(f, &p);
                      double ev = gv - fv;
                      fprintf(stderr, "fAppr((%2d,%2d,%2d)/%5.3f)", i, j, k, pr);
                      fprintf(stderr, " = %9.6f", gv);
                      fprintf(stderr, "  fAppr-fCorr = %16.12f\n", ev);
                    }
                  }
              }
          }
      }
  }
  
void SPApprox_PrintMaxAvgValues
  ( SPFunction *f,      /* A function. */
    Triangulation *tri, /* Triangulation to use for integration. */
    double *fMax,       /* OUT: max absolute value of {f}. */
    double *fAvg        /* OUT: mean value of {f}. */ 
  )
  {
    auto double fValue(S2Point *p);
    
    double fValue(S2Point *p) { return f->m->eval(f, p); }

    (*fMax) = SPRange_OnSphere(fValue, 120);
    fprintf(stderr, "max(fabs(f)) = %16.12f\n", *fMax);
    (*fAvg) = SPIntegral_OnSphere(fValue, tri); 
    fprintf(stderr, "avg(f) = %16.12f\n", *fAvg);
  }
  
void SPApprox_PrintMaxErrorValues
  ( SPFunction *g,      /* Approximation. */
    SPFunction *f,      /* Target function. */
    Triangulation *tri, /* Triangulation to use for integration. */
    double *gMax,       /* OUT: max abs value of {g}. */
    double *fMax,       /* OUT: max abs value of {f}. */
    double *eMax,       /* OUT: max abs error {g-f}. */
    double *eAvg        /* OUT: root mean square error {g-f}. */
  )
  {
    auto double fValue(S2Point *p);
    auto double gValue(S2Point *p);
    auto double gfError(S2Point *p);
    auto double gfErrorSqr(S2Point *p);
    
    double fValue(S2Point *p) { return f->m->eval(f, p); }

    double gValue(S2Point *p) { return g->m->eval(g, p); }

    double gfError(S2Point *p) { return g->m->eval(g, p) - f->m->eval(f, p); }

    double gfErrorSqr(S2Point *p) 
      { double d = g->m->eval(g, p) - f->m->eval(f, p);
        return d*d;
      }

    (*gMax) = SPRange_OnSphere(gValue, 120);
    fprintf(stderr, "max(fabs(fAppr))       = %10.2e\n", *gMax);
    (*fMax) = SPRange_OnSphere(fValue, 120);
    fprintf(stderr, "max(fabs(fCorr))       = %10.2e\n", *fMax);
    (*eMax) = SPRange_OnSphere(gfError, 120);
    fprintf(stderr, "max(fabs(fAppr-fCorr)) = %10.2e\n", *eMax);
    (*eAvg) = sqrt(fabs(SPIntegral_OnSphere(gfErrorSqr, tri))); 
    fprintf(stderr, "rms(fAppr-fCorr)       = %10.2e\n", *eAvg);
  }
  
void SPApprox_PlotFunctionAndError
  ( SPPlot_Stream *fps,      /* Postscript file, or NULL. */
    char *funcTag,           /* Prefix for page names. */
    SPFunction *g,           /* The computed approximation. */
    double f(S2Point *p),    /* The target function. */
    Triangulation *tri,      /* Reference triangulation. */
    double relMeshSize,      /* Maximum step/triangle size (radians). */
    bool_t showTriang,       /* TRUE to plot the reference triangulation. */ 
    bool_t plotTrueSol,      /* TRUE to plot the target function {f}. */ 
    double fMax,             /* Nominal maximum of {fabs(g(p))} and {fabs(f(p))}. */
    double eMax,             /* Nominal maximum of {fabs(g(p)-f(p))}. */
    SPH3_Point *obs,         /* Viewpoint. */
    SPH3_Point *upp,         /* Camera vertical reference. */
    string_vec_t *caption,   /* Caption template. */
    int index,               /* Figure `index' for caption expansion. */
    double time              /* Figure `time' for caption expansion. */
  )
  { 
    double rad = 1.0;
    SPH3_Point obsFix = (*obs), uppFix = (*upp);
    
    r3_t dLight = (r3_t){{0.0, 0.0, 1.0}};
    SPPlot_FixView(&(obsFix), 0.0, &(uppFix), &rad, 0.0, g, TRUE);

    /* Select isoline step for function plotting, adjust {fMax}: */
    double fStep = SPPlot_RoundToNice(fMax/DefaultIsolines);
    double fRange = fStep * DefaultIsolines;

    /* Copy the caption vector: */
    int nCap = caption->ne;  /* Number of original caption lines. */
    string_vec_t intcap = string_vec_new(nCap);
    { int k; for (k = 0; k < nCap; k++) { intcap.e[k] = caption->e[k]; } }
    
    /* May need one more caption line, to tell approximation from error: */
    int iSubCap; /* index of sub-caption in {intcap}, or -1 if none. */
    if ((nCap > 0) || (! fps->eps))
      { iSubCap = nCap;
        string_vec_expand(&intcap, iSubCap); 
        intcap.e[iSubCap] = NULL; /* To be replaced. */
        string_vec_trim(&intcap, iSubCap); 
      }
    else
      { iSubCap = -1; }
  
    bool_t plotApprox = TRUE;
    bool_t plotError = TRUE; 
    
    if (plotApprox)
      { 
        /* Plot function {g}: */
        pswr_fill_row(fps);

        auto double gValue(S2Point *p);
        double gValue(S2Point *p)
          { return g->m->eval(g, p); }

        /* Get and print the supporting plane of {g}: */
        SPH3_Plane gSupp = SPFunction_GetSupportingPlane(g);
        r4_gen_print(stderr, &(gSupp.f), "%6.3f", "support = [ ", ", ", " ]\n");

        fprintf(stderr, "plot fRange = %10.2e fStep = %10.2e\n", fRange, fStep);

        double gMinObs =   INFINITY;
        double gMaxObs = - INFINITY;

        char *gTag = txtcat(funcTag, "-app");

        if (iSubCap > 0) { intcap.e[iSubCap] = "Approximation"; }

        SPPlot_MultipleViews
          ( /* {fps} */          fps, 
            /* {funcTag} */      gTag,
            /* {func} */         gValue, 
            /* {tri} */          tri, 
            /* {relMeshSize} */  relMeshSize,
            /* {showTriang} */   showTriang,
            /* {supp} */         &gSupp, 
            /* {fRange} */       fRange,
            /* {fStep} */        fStep,
            /* {obs} */          &obsFix,
            /* {upp} */          &uppFix,
            /* {rad} */          SPPlot_FullRadius,
            /* {dLight} */       &dLight,
            /* {lineWidth} */    0.25,
            /* {gridNLon} */     0,
            /* {gridNLat} */     0,
            /* {gridDots} */     FALSE,
            /* {aSide} */        +1, 
            /* {bSide} */        -1, 
            /* {caption} */      &intcap,
            /* {index} */        index, 
            /* {time} */         time,
            /* {error} */        eMax,
            /* {capAlign} */     0.5,
            /* {verbose} */      FALSE,
            /* {fMinObs} */      &gMinObs,
            /* {fMaxObs} */      &gMaxObs
          );

        fprintf(stderr, "observed function extrema:");
        fprintf(stderr, " min = %+22.14e max = %+22.14e\n", gMinObs, gMaxObs);
        fprintf(stderr, "\n");
        free(gTag);
      }
    
    if (plotTrueSol)
      { 
        /* Plot target function {f}: */
         pswr_fill_row(fps);
       
        /* Select isoline step, adjust {Max}: */
        fprintf(stderr, "plot fRange = %10.2e fStep = %10.2e\n", fRange, fStep);

        double fMinObs =   INFINITY;
        double fMaxObs = - INFINITY;

        char *fTag = txtcat(funcTag, "-sol");

        if (iSubCap > 0) { intcap.e[iSubCap] = "True solution"; }

        SPPlot_MultipleViews
          ( /* {fps} */          fps, 
            /* {funcTag} */      fTag,
            /* {func} */         f, 
            /* {tri} */          tri, 
            /* {relMeshSize} */  relMeshSize,
            /* {showTriang} */   showTriang,
            /* {supp} */         NULL, 
            /* {fRange} */       fRange,
            /* {fStep} */        fStep,
            /* {obs} */          &obsFix,
            /* {upp} */          &uppFix,
            /* {rad} */          SPPlot_FullRadius,
            /* {dLight} */       &dLight,
            /* {lineWidth} */    0.25,
            /* {gridNLon} */     0,
            /* {gridNLat} */     0,
            /* {gridDots} */     FALSE,
            /* {aSide} */        +1, 
            /* {bSide} */        -1, 
            /* {caption} */      &intcap,
            /* {index} */        index, 
            /* {time} */         time,
            /* {error} */        0.0,
            /* {capAlign} */     0.5,
            /* {verbose} */      FALSE,
            /* {fMinObs} */      &fMinObs,
            /* {fMaxObs} */      &fMaxObs
          );

        fprintf(stderr, "observed function extrema:");
        fprintf(stderr, " min = %+22.14e max = %+22.14e\n", fMinObs, fMaxObs);
        fprintf(stderr, "\n");
        free(fTag);
      }
    
    if (plotError)
      { 
        /* Plot error {e = g - f}: */
        pswr_fill_row(fps);

        auto double eValue(S2Point *p);
        double eValue(S2Point *p)
          { return g->m->eval(g, p) - f(p); }

        double eStep = SPPlot_RoundToNice(eMax/DefaultIsolines);
        eMax = eStep * DefaultIsolines;
        fprintf(stderr, "plot eMax = %10.2e eStep = %10.2e\n", eMax, eStep);

        double eMinObs = INFINITY;
        double eMaxObs = - INFINITY;

        char *eTag = txtcat(funcTag, "-err");

        if (iSubCap > 0) { intcap.e[iSubCap] = "Approximation error"; }

        SPPlot_MultipleViews
          ( /* {fps} */          fps, 
            /* {funcTag} */      eTag,
            /* {func} */         eValue, 
            /* {tri} */          tri, 
            /* {relMeshSize} */  relMeshSize,
            /* {showTriang} */   showTriang,
            /* {supp} */         NULL, 
            /* {fRange} */       eMax,
            /* {fStep} */        eStep,
            /* {obs} */          &obsFix,
            /* {upp} */          &uppFix,
            /* {rad} */          SPPlot_FullRadius,
            /* {dLight} */       &dLight,
            /* {lineWidth} */    0.25,
            /* {gridNLon} */     0,
            /* {gridNLat} */     0,
            /* {gridDots} */     FALSE,
            /* {aSide} */        +1, 
            /* {bSide} */        -1, 
            /* {caption} */      &intcap,
            /* {index} */        index, 
            /* {time} */         time,
            /* {error} */        eMax,
            /* {capAlign} */     0.5,
            /* {verbose} */      FALSE,
            /* {fMinObs} */      &eMinObs,
            /* {fMaxObs} */      &eMaxObs
          );

        fprintf(stderr, "observed error extrema:");
        fprintf(stderr, " min = %+10.2e max = %+10.2e\n", eMinObs, eMaxObs);
        fprintf(stderr, "\n");
        free(eTag);
      }

  }

