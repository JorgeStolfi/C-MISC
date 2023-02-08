/* See SOBasisMatrices.h */
/* Last edited on 2007-01-04 00:21:36 by stolfi */

#include <SOBasisMatrices.h>
#include <SOFunction.h>
//#include <SOSpline.h>
#include <SOMatrix.h>
#include <SOBasic.h>
#include <vec.h>
#include <math.h>
#include <limits.h>
 
SOMatrix SOBasisMatrices_BuildMatrixGen
  ( Basis F,
    Basis G, 
    void dot(SOFunction *f, SOFunction *g, double *iv), //double dot(SOFunction *f, SOFunction *g), 
    double minVal,
    bool_t symmetric,
    bool_t verbose
  )
  { double fg;
    nat_t m = F.ne;
    nat_t n = G.ne;
    MatEntry_vec_t Ments = MatEntry_vec_new(F.ne+G.ne);
    double_vec_t fsqr = double_vec_new(F.ne);
    double_vec_t gsqr = double_vec_new(G.ne);
    int i, j; nat_t kM = 0;
    bool_t veryverbose = FALSE;
    
    if (minVal > 0.0) 
      { /* Compute the element norms squared, for small entry cleanup: */
        if (verbose) { fprintf(stderr, "computing norms...\n"); }
        for (i = 0; i < m; i++) 
          { dot(F.e[i], F.e[i], &(fsqr.e[i])); } //{ fsqr.e[i] = dot(F.e[i], F.e[i]); }
        for (j = 0; j < n; j++) 
          { //{ gsqr.e[j]=(symmetric?fsqr.e[j]:dot(G.e[j],G.e[j]));}
	    if(symmetric) gsqr.e[j] = fsqr.e[j];
            else dot(G.e[j], G.e[j], &(gsqr.e[j]));
          } 
      }
    
    /* Compute the matrix entries: */
    if (verbose) { fprintf(stderr, "computing matrix...\n"); }
    for (i = 0; i < m; i++)
      { SOFunction *f = F.e[i];
        if (veryverbose) 
          { fprintf(stderr, "[row %d:", i); }
        else if (verbose)
          { fprintf(stderr, "-"); }
        
        for (j = 0; j < n; j++)
          { if (symmetric && (j < i))
              { nat_t ksym = SOMatrix_Find(Ments.e, kM, j, i);
                if (veryverbose) { fprintf(stderr, " !"); }
                fg = (ksym >= INT_MAX ? 0.0 : Ments.e[ksym].va); 
              }
            else 
              { SOFunction *g = G.e[j];
                if (veryverbose) { fprintf(stderr, " *"); }
                dot(f, g, &fg); //fg = dot(f, g);
              }
            if ((minVal > 0.0) && (fabs(fg) < minVal*sqrt(fsqr.e[i]*gsqr.e[j])))
              { fg = 0.0; }
            if (fg != 0.0)
              { if (veryverbose) { fprintf(stderr, "<F%d,G%d>=%.7f", i, j, fg); }
                MatEntry_vec_expand(&Ments, kM);
                { MatEntry *Mij = &(Ments.e[kM]);
                  Mij->row = i; Mij->col = j;
                  Mij->va = fg; kM++;
                }
              }
          }
        if (veryverbose) { fprintf(stderr, "]\n"); }
      }
    if (veryverbose) { fprintf(stderr, "\n"); }
    MatEntry_vec_trim(&Ments, kM);
    if (minVal > 0.0) { free(fsqr.e); free(gsqr.e); }
    return (SOMatrix){F.ne, G.ne, Ments};
  }
 
void SOBasisMatrices_RecomputeVectorGen
  ( SOFunction* f,
    Basis bas,
    SOFunction_fgDot dot,
    double_vec_t y,
    bool_t verbose
  )
  { int i, j, fd = f->d->fDim;
    double sdy2[f->d->fDim], bNew[f->d->fDim], dy;
    for(j = 0; j < fd; j++) sdy2[j] = 0.0;

    if (verbose) { fprintf(stderr, "computing right-hand side \"y\"\n"); }
    for (i = 0; i < bas.ne; i++)
      { dot(f, bas.e[i], bNew);
        for(j = 0; j < fd; j++)
          {
            if (verbose) 
              { dy = bNew[j] - y.e[i*fd + j];
                sdy2[j] += dy*dy;
                fprintf( stderr, "  y[%3d][%3d] = %16.12f  change = %16.12f\n", i, j, bNew[j], dy ); 
              }
            y.e[i*fd + j] = bNew[j];
          }
      }
    if (verbose) 
      for(j = 0; j < fd; j++)
        { fprintf(stderr, "\n");
          fprintf(stderr, "  total change in \"y\"[%3d] = %16.12f\n", j, sqrt(sdy2[j]));
        }
}
   
SOMatrix SOBasisMatrices_BuildMatrixEval
  ( Basis F,
    Basis G, 
    SOFunction *w,
    bool_t ignoreTriang,
    double minVal,
    bool_t verbose
  )
  { nat_t pDim = F.e[0]->d->pDim;
    nat_t fDim = F.e[0]->d->fDim;
    FuncMap NoMap = IdentityMap(fDim, pDim);
    
    auto void dot(SOFunction *f, SOFunction *g, double *iv);
      /* Computes the dot product of {f} and {g}. */
    
    void dot(SOFunction *f, SOFunction *g, double *iv)
      { if (ignoreTriang)
          { *iv = SOFunction_RawDot(f, &NoMap, g, &NoMap, w, NULL, FALSE); return;}
        else
          { SOFunction_Dot(f, &NoMap, g, &NoMap, w, NULL, iv, FALSE);}
      }
   
    bool_t symmetric = ((F.ne == G.ne) & (F.e == G.e));
    return SOBasisMatrices_BuildMatrixGen( F, G, dot, minVal, symmetric, verbose);
  }

SOMatrix SOBasisMatrices_BuildMatrixGrad
  ( Basis F, 
    Basis G, 
    SOFunction *w,
    bool_t ignoreTriang,
    double minVal,
    bool_t verbose
  )
  { 
    auto void dot(SOFunction *f, SOFunction *g, double *iv);
      /* Computes the dot product of the gradients of {f} and {g}. */
    
    void dot(SOFunction *f, SOFunction *g, double *iv)
      { if (ignoreTriang)
          { *iv = SOFunction_RawGradDot(f, g, w, NULL, FALSE); return;}
        else
          { SOFunction_GradDot(f, g, w, NULL, iv, FALSE); return;}
      }
    
    bool_t symmetric = ((F.ne == G.ne) & (F.e == G.e));
    return SOBasisMatrices_BuildMatrixGen( F, G, dot, minVal, symmetric, verbose);
  }
 
void SOBasisMatrices_ChangeBasis
  ( double_vec_t a,
    SOMatrix FG,
    SOMatrix GGL,
    double_vec_t b 
  )
  { nat_t n = FG.cols;
    double_vec_t c = double_vec_new(n);
    double_vec_t y = double_vec_new(n);
    SOMatrix_MulRow(a, FG, c);
    SOMatrix_DivCol(GGL, c, y);
    SOMatrix_DivRow(y, GGL, b);
  }

