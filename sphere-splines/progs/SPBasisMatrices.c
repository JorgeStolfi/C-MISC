/* See SPBasisMatrices.h */
/* Last edited on 2005-10-29 05:48:28 by stolfi */

#include <SPBasisMatrices.h>
#include <SPFunction.h>
#include <SPSpline.h>
#include <SPMatrix.h>
#include <SPVector.h>
#include <SPBasic.h>

#include <vec.h>
#include <js.h>

#include <math.h>
#include <limits.h>

SPMatrix SPBasisMatrices_BuildMatrixGen
  ( Basis F,
    Basis G, 
    double dot(SPFunction *f, SPFunction *g), 
    Triangulation *tri,
    bool_t semiOrtho,
    bool_t orthogonal,
    bool_t symmetric,
    double minVal,
    bool_t verbose
  )
  { double fg;
    int nF = F.ne;
    int nG = G.ne;
    MatEntry_vec_t Ments = MatEntry_vec_new(F.ne+G.ne);
    double_vec_t fsqr = double_vec_new(F.ne);
    double_vec_t gsqr = double_vec_new(G.ne);
    int i, j; int kM = 0;
    bool_t veryverbose = FALSE;
    bool_t cleanup = (! orthogonal) && (minVal > 0.0);
    
    if (cleanup)
      { /* Compute the element norms squared, for small entry cleanup: */
        if (verbose) { fprintf(stderr, "computing norms...\n"); }
        for (j = 0; j < nF; j++) 
          { fsqr.e[j] = dot(F.e[j], F.e[j]); }
        for (i = 0; i < nG; i++) 
          { gsqr.e[i] = (symmetric ? fsqr.e[i] : dot(G.e[i], G.e[i])); }
      }
    
    /* Compute the matrix entries: */
    if (verbose) { fprintf(stderr, "computing matrix...\n"); }
    for (i = 0; i < nG; i++)
      { SPFunction *g = G.e[i];
        
        if (veryverbose) 
          { fprintf(stderr, "[row %d:", i); }
        else if (verbose)
          { fprintf(stderr, "-"); }
        
         for (j = 0; j < nF; j++)
          { 
            SPFunction *f = F.e[j];
            
            /* Compute/fetch/assume {M[i,j] = fg = dot(F[j],G[i])}: */
            if (orthogonal && (i != j))
              { if (veryverbose) { fprintf(stderr, " +"); }
                fg = 0.0; 
              }
            else if (semiOrtho && SPSpline_NestedSupports(f, g, tri))
              { if (veryverbose) { fprintf(stderr, " ×"); }
                fg = (double)(i == j); 
              }
            else if (symmetric && (j < i))
              { int ksym = SPMatrix_Find(Ments.e, kM, j, i);
                if (veryverbose) { fprintf(stderr, " !"); }
                fg = (ksym >= INT_MAX ? 0.0 : Ments.e[ksym].va); 
              }
            else 
              { 
                if (veryverbose) { fprintf(stderr, " *"); }
                fg = dot(f, g);
              }
            
            /* Cleanup small off-diagonal entries: */
            if ((i != j) && cleanup)
              { if (fabs(fg) < minVal*sqrt(fsqr.e[j]*gsqr.e[i])) { fg = 0.0; } }

            /* Store in matrix if nonzero: */
            if (fg != 0.0)
              { if (veryverbose) { fprintf(stderr, "<F%d,G%d>=%.7f", j, i, fg); }
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
    if (cleanup) { free(fsqr.e); free(gsqr.e); }
    return (SPMatrix){F.ne, G.ne, Ments};
  }

void SPBasisMatrices_RecomputeVectorGen
  ( SPFunction* f,
    Basis bas,
    Metric dot,
    SPVector y,
    bool_t verbose
  )
  { int i;
    double sdy2 = 0.0;
    if (verbose) { fprintf(stderr, "computing right-hand side \"y\"\n"); }
    for (i = 0; i < bas.ne; i++)
      { double bNew = dot(f, bas.e[i]);
        fprintf(stderr, "-"); 
        if (verbose) 
          { double dy = bNew - y.e[i];
            sdy2 += dy*dy;
            fprintf
              ( stderr, "  y[%3d] = %16.12f  change = %16.12f\n", i, bNew, dy ); 
          }
        y.e[i] = bNew;
      }
    if (verbose) 
      { fprintf(stderr, "\n");
        fprintf(stderr, "  total change in \"y\" = %16.12f\n", sqrt(sdy2));
      }
  }
   
SPMatrix SPBasisMatrices_BuildMatrixEval
  ( Basis F,
    Basis G, 
    SPFunction *w,
    Triangulation *tri,
    bool_t semiOrtho,
    bool_t orthogonal,
    bool_t ignoreTriang,
    double minVal,
    bool_t verbose
  )
  { 
    auto double dot(SPFunction *f, SPFunction *g);
      /* Computes the dot product of {f} and {g}. */
    
    double dot(SPFunction *f, SPFunction *g)
      { if (ignoreTriang)
          { return SPFunction_RawDot(f, NoFMap, g, NoFMap, w, tri, FALSE); }
        else
          { return SPFunction_Dot(f, NoFMap, g, NoFMap, w, tri, FALSE); }
      }
   
    bool_t symmetric = ((F.ne == G.ne) & (F.e == G.e));
    return SPBasisMatrices_BuildMatrixGen
      (F, G, dot, tri, semiOrtho, orthogonal, symmetric, minVal, verbose);
  }

SPMatrix SPBasisMatrices_BuildMatrixVelSGrd
  ( Basis F, 
    Basis G, 
    SPFunction *w,
    bool_t orthogonal,
    bool_t ignoreTriang,
    double minVal,
    bool_t verbose
  )
  { 
    auto double dot(SPFunction *f, SPFunction *g);
      /* Computes the dot product of the gradients of
        {w(p)*r3_dot(v(p) ¤ SGrd(f)(p))*g(p)}. */
    
    double dot(SPFunction *f, SPFunction *g)
      { if (ignoreTriang)
          { return SPFunction_RawVelSGrdDot(f, g, w, NULL, FALSE); }
        else
          { return SPFunction_VelSGrdDot(f, g, w, NULL, FALSE); }
      }
    
    bool_t symmetric = FALSE;   /* In general. */
    return SPBasisMatrices_BuildMatrixGen
      (F, G, dot, NULL, FALSE, orthogonal, symmetric, minVal, verbose);
  }

SPMatrix SPBasisMatrices_BuildMatrixSLap
  ( Basis F, 
    Basis G, 
    SPFunction *w,
    bool_t orthogonal,
    bool_t ignoreTriang,
    double minVal,
    bool_t verbose
  )
  { 
    auto double dot(SPFunction *f, SPFunction *g);
      /* Computes the dot product of {SLap(f)} and {g}. */
    
    double dot(SPFunction *f, SPFunction *g)
      { if (ignoreTriang)
          { return SPFunction_RawSLapDot(f, g, w, NULL, FALSE); }
        else
          { return SPFunction_SLapDot(f, g, w, NULL, FALSE); }
      }
    
    bool_t symmetric = ((F.ne == G.ne) & (F.e == G.e));
    return SPBasisMatrices_BuildMatrixGen
      ( F, G, dot, NULL, FALSE, orthogonal, symmetric, minVal, verbose);
  }
 
void SPBasisMatrices_ChangeBasis
  ( SPVector a,
    SPMatrix FG,
    SPMatrix GGL,
    SPVector b 
  )
  { int n = FG.cols;
    SPVector c = double_vec_new(n);
    SPVector y = double_vec_new(n);
    SPMatrix_MulRow(a, FG, c);
    SPMatrix_DivCol(GGL, c, y);
    SPMatrix_DivRow(y, GGL, b);
  }

