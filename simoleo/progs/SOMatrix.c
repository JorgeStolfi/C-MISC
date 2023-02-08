/* See SOMatrix.h */
/* Last edited on 2006-02-28 12:06:22 by stolfi */

#include <SOMatrix.h>

#include <SOBasic.h>
#include <vec.h>
#include <filefmt.h>
#include <fget.h>
#include <nget.h>
#include <affirm.h>
#include <nat.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
  
void SOMatrix_UnpackRow(MatEntry_vec_t *Aents, nat_t row, double_vec_t Ai, nat_t *kA);
  /* Unpacks row {i} of the matrix {A} as an ordinary vector of double
    {Ai}. On input, the variable {*kA} must be the index into {Aents ==}
{    A.ents} of the first element of that row. On output, {*kA} will be
    the index of the first element of row {i+1}. */

void SOMatrix_UnpackHalfRow(MatEntry_vec_t *Lents, nat_t row, double_vec_t Li, nat_t *kL);
  /* Same as {UnpackRow}, but expects and checks that {L} is lower triangular ---
    meaning that there are no elements {L[i,j]} with {j > i}. */

void SOMatrix_PackRow(double_vec_t Ci, nat_t row, MatEntry_vec_t *Cents, nat_t *nC);
  /* Appends the non-zero entries of {Ci} as one more row of elements
    of a matrix {C}. Will reallocate the element list {Cents == C.ents} if
    necessary. On input, the variable {*nC} must be the count of
    elements already stored in {Cents^}; this count will be updated on
    exit. Don't forget to trim the element list after adding the last
    row. */

SOMatrix SOMatrix_Null(nat_t rows, nat_t cols)
  { MatEntry_vec_t ents = MatEntry_vec_new(0);
    return (SOMatrix){rows, rows, ents};
  }

SOMatrix SOMatrix_BinaryOp
  ( SOMatrix A, SOMatrix B,
    double func(int row, int col, double va, double vb),
    bool_t skip00
  )
  { MatEntry *Ael = A.ents.e; int nA = A.ents.ne;
    MatEntry *Bel = B.ents.e; int nB = B.ents.ne;
    MatEntry_vec_t Cents = MatEntry_vec_new(nA + nB);
    nat_t kA = 0, kB = 0, nC = 0;
    nat_t row, col;
    double va, vb, vc;
    affirm(A.rows == B.rows, "inconsistent rows");
    affirm(A.cols == B.cols, "inconsistent cols");
    row = -1; col = A.cols-1;
    while (TRUE)
      { /* Determine indices {[row,col]} of next element to be computed: */
        if (! skip00)
          { /* Must evaluate even elements which are missing in {A} and {B}: */
            col++; 
            if (col >= A.cols) { row++; col = 0; }
          }
        else if ((kA >= nA) && (kB >= nB))
          { row = A.rows; col = 0; }
        else if (kA >= nA)
          { row = Bel[kB].row; col = Bel[kB].col; }
        else if (kB >= nB)
          { row = Ael[kA].row; col = Ael[kA].col; }
        else if (Ael[kA].row == Bel[kB].row)
          { if (Ael[kA].col <= Bel[kB].col)
              { row = Ael[kA].row; col = Ael[kA].col; }
            else
              { row = Bel[kB].row; col = Bel[kB].col; }
          }
        else if (Ael[kA].row < Bel[kB].row)
          { row = Ael[kA].row; col = Ael[kA].col; }
        else
          { row = Bel[kB].row; col = Bel[kB].col; }
        
        /* Are we done already? */
        if (row >= A.rows) { break; }
        
        /* Get element values from {A} and {B}, and skip over them: */
        if ((kA < nA) && (row == Ael[kA].row) && (col == Ael[kA].col))
          { va = Ael[kA].va; kA++; }
        else
          { va = 0.0; }
        if ((kB < nB) && (row == Bel[kB].row) && (col == Bel[kB].col))
          { vb = Bel[kB].va; kB++; }
        else
          { vb = 0.0; }
        
        /* Compute element of result: */
        if (skip00 && (va == 0.0) && (vb == 0.0))
          { vc = 0.0; } 
        else
          { vc = func(row, col, va, vb); }

        /* Append element to {C}: */
        if (vc != 0.0)
          { MatEntry_vec_expand(&Cents, nC);
            { MatEntry *Cij = &(Cents.e[nC]);
              Cij->row = row; Cij->col = col;
              Cij->va = vc; nC++;
            }
          }
      }
    MatEntry_vec_trim(&Cents, nC);
    return (SOMatrix){A.rows, B.rows, Cents};
  }

SOMatrix SOMatrix_Mix(double alpha, SOMatrix A, double beta, SOMatrix B)
  { auto double mix(int row, int col, double va, double vb);
    
    double mix(int row, int col, double va, double vb)
      { return alpha*va + beta*vb; }
      
    return SOMatrix_BinaryOp(A, B, mix, TRUE);
  }
  
SOMatrix SOMatrix_RelDiff(SOMatrix A, SOMatrix B)
  { auto double rdiff(int row, int col, double va, double vb);
    
    double rdiff(int row, int col, double va, double vb)
      { double ma = fabs(va), mb = fabs(vb);
        double m = (ma > mb ? ma : mb);
        if (m == 0.0) 
          { return 0.0; }
        else
          { double ra = va/m, rb = vb/m;
            return 0.5*(rb - ra)/sqrt(ra*ra + rb*rb);
          }
      }
      
    return SOMatrix_BinaryOp(A, B, rdiff, TRUE);
  }

SOMatrix SOMatrix_Scale(double alpha, SOMatrix A)
  { MatEntry *Ael = A.ents.e; int nA = A.ents.ne;
    nat_t nC = 0, kA;
    MatEntry_vec_t Cents = MatEntry_vec_new((alpha == 0.0 ? 0 : nA));
    for (kA = 0; kA < nA; kA++)
      { MatEntry *Aij = &(Ael[kA]);
        MatEntry *Cij = &(Cents.e[nC]);
        double va = alpha * Aij->va;
        if (va != 0.0) { (*Cij) = (*Aij); Cij->va = va; nC++; }
      }
    MatEntry_vec_trim(&Cents, nC);
    return (SOMatrix){A.rows, A.cols, Cents};
  }

void SOMatrix_GetDiagonal(SOMatrix A, double_vec_t d)
  { MatEntry *Ael = A.ents.e; int nA = A.ents.ne;
    int nd = (A.rows < A.cols ? A.rows : A.cols);
    int kd = 0, kA;
    affirm(d.ne == nd, "d has wrong size");
    for (kA = 0; kA < nA; kA++)
      { MatEntry *Aij = &(Ael[kA]);
        if (Aij->row == Aij->col) 
          { int i = Aij->row;
            affirm(i <= kd, "matrix entries out of order");
            while (kd < i) { d.e[kd] = 0.0; kd++; }
            d.e[kd] = Aij->va;
            kd++;
          }
      }
    while (kd < nd) { d.e[kd] = 0.0; kd++; }
  }

SOMatrix SOMatrix_Mul(SOMatrix A, SOMatrix B)
  { int nA = A.ents.ne;
    int nB = B.ents.ne;
    nat_t maxnel = (nA > nB ? nA : nB);
    MatEntry_vec_t Cents = MatEntry_vec_new(maxnel);
    double_vec_t a = double_vec_new(A.rows);
    double_vec_t c = double_vec_new(B.cols);
    int i; nat_t kA = 0, nC = 0;

    affirm(A.cols == B.rows , "incompatible sizes");
    for (i = 0; i < A.rows; i++)
      { SOMatrix_UnpackRow(&A.ents, i, a, &kA);
        SOMatrix_MulRow(a, B, c);
        SOMatrix_PackRow(c, i, &Cents, &nC);
      }
    MatEntry_vec_trim(&Cents, nC);
    free(a.e); free(c.e);
    return (SOMatrix){A.rows, B.cols, Cents};
  }

void SOMatrix_MulRow(double_vec_t x, SOMatrix A, double_vec_t y)
  { MatEntry *Ael = A.ents.e; int nA = A.ents.ne;
    int j, kA;
    affirm(x.ne == A.rows, "vector has wrong size");
    affirm(y.ne == A.cols, "vector has wrong size");
    for (j = 0; j < A.cols; j++) { y.e[j] = 0.0; }
    for (kA = 0; kA < nA; kA++)
      { MatEntry *Aij = &(Ael[kA]);
        nat_t i = Aij->row, j = Aij->col;
        y.e[j] += x.e[i] * Aij->va;
      }
  }

void SOMatrix_MulCol(SOMatrix A, double_vec_t x, double_vec_t y)
  { MatEntry *Ael = A.ents.e; int nA = A.ents.ne;
    int i, kA;
    affirm(x.ne == A.cols, "vector has wrong size");
    affirm(y.ne == A.rows, "vector has wrong size");
    for (i = 0; i < A.cols; i++) { y.e[i] = 0.0; }
    for (kA = 0; kA < nA; kA++)
      { MatEntry *Aij = &(Ael[kA]);
        nat_t i = Aij->row, j = Aij->col;
        y.e[i] += Aij->va * x.e[j];
      }
  }

void SOMatrix_DivRow(double_vec_t y, SOMatrix L, double_vec_t x)
  { MatEntry *Lel = L.ents.e; int nL = L.ents.ne;
    int i, kL;

    affirm(L.rows == L.cols, "matrix is not square");
    affirm(y.ne == L.cols, "vector has wrong size");
    affirm(x.ne == L.cols, "vector has wrong size");

    for (i = 0; i < L.cols; i++) { x.e[i] = y.e[i]; }
    for (kL = nL-1; kL >= 0; kL--)
      { MatEntry *Lij = &(Lel[kL]);
          nat_t i = Lij->row, j = Lij->col;
          affirm(i >= j, "not lower triangular");
          if (i > j)
            { x.e[j] -= Lij->va * x.e[i]; }
          else
            { x.e[j] /= Lij->va; }
        }
    }

void SOMatrix_DivCol(SOMatrix L, double_vec_t y, double_vec_t x)
  { int i, j, kL;
    double_vec_t Li = double_vec_new(L.cols);
    
    affirm(L.rows == L.cols, "matrix is not square");
    affirm(y.ne == L.cols, "vector has wrong size");
    affirm(x.ne == L.cols, "vector has wrong size");
    
    kL = 0;
    for (i = 0; i < L.rows; i++)
      { SOMatrix_UnpackHalfRow(&L.ents, i, Li, &kL);
        x.e[i] = y.e[i];
        for (j = 0; j < i; j++)
          { x.e[i] -= Li.e[j]*x.e[j]; }
        x.e[i] /= Li.e[i];
      }
    free(Li.e);
  }

SOMatrix SOMatrix_PosDiv(SOMatrix L, SOMatrix A)
  { int i, kA, nX;
    MatEntry_vec_t Xents = MatEntry_vec_new(A.ents.ne);
    double_vec_t Ai = double_vec_new(A.cols);
    double_vec_t Xi = double_vec_new(A.cols);
    
    affirm(L.rows == L.cols, "L matrix is not square");
    affirm(A.cols == L.rows, "incompatible matrix sizes");
    
    kA = 0; nX = 0;
    for (i = 0; i < A.rows; i++)
      { SOMatrix_UnpackRow(&(A.ents), i, Ai, &kA);
        SOMatrix_DivRow(Ai, L, Xi);
        SOMatrix_PackRow(Xi, i, &Xents, &nX);
      }
    free(Ai.e); free(Xi.e);
    MatEntry_vec_trim(&Xents, nX);
    return (SOMatrix){A.rows, A.cols, Xents}; 
  }

SOMatrix SOMatrix_PreDiv(SOMatrix L, SOMatrix A)
  { int i, kL, kA, kX, nX;
    MatEntry_vec_t Xents = MatEntry_vec_new(A.ents.ne);
    double_vec_t Li = double_vec_new(L.cols);
    double_vec_t Xi = double_vec_new(A.cols);
    
    affirm(L.rows == L.cols, "L matrix is not square");
    affirm(A.rows == L.cols, "A matrix has wrong rows");
    
    kL = 0; kA = 0; nX = 0;
    for (i = 0; i < L.rows; i++)
      { SOMatrix_UnpackRow(&(L.ents), i, Li, &kL);
        SOMatrix_UnpackRow(&(A.ents), i, Xi, &kA);
        kX = 0;
        while(kX < nX)
          { MatEntry *Xjk = &(Xents.e[kX]);
            int j = Xjk->row, k = Xjk->col;
            Xi.e[k] -= Li.e[j] * Xjk->va;
          }
        { int k; for (k = 0; k < A.cols; k++) { Xi.e[k] /= Li.e[i]; } }
        SOMatrix_PackRow(Xi, i, &Xents, &nX);
      }
    free(Li.e); free(Xi.e);
    MatEntry_vec_trim(&Xents, nX);
    return (SOMatrix){A.rows, A.cols, Xents}; 
  }

SOMatrix SOMatrix_Cholesky(SOMatrix A, double minVal)
  { 
    /* Andre-Louis Cholesky (spelled with a 'y'), born in France in
      1875, was a geodesist in the French military. He developed his
      computational procedure to solve geodetic problems, among which
      was, in the words of his eulogy, the /problem of levelling in
      Morocco./ He died in battle in 1918. Benoit published the
      computational method in /Bulletin geodesique/ in 1924.

      The most likely French pronunciation of Cholesky is that heard
      in 'Chopin' -- the English 'Sh' sound. 
        -- David Pattison, Washington, DC */
  
    MatEntry_vec_t Lents = MatEntry_vec_new(A.rows);
    int i, j, nL = 0, kL, kA = 0;
    double Lii;
    
    double_vec_t Ai = double_vec_new(A.rows);
    double_vec_t Li = double_vec_new(A.rows);
    double_vec_t dg = double_vec_new(A.rows); /* Diagonal elements */
    
    affirm(A.rows == A.cols , "matrix must be square");
    for (i = 0; i < A.rows; i++)
      { SOMatrix_UnpackRow(&A.ents, i, Ai, &kA);
        /* Compute {L[i,j] = X[i,j] - SUM{ L[i,k]*A[j,k] : k = 0..j-1 })/L[j,j]} */
        /* where {X[i,j] = (j <= i ? A[i,j] : 0)} */
        kL = 0;
        for (j = 0; j < i; j++)
          { double sum = 0.0, corr = 0.0;
            while (kL < nL)
              { MatEntry *Lrs = &(Lents.e[kL]); 
                int k = Lrs->col;
                affirm(Lrs->row == j, "matrix indices buggy");
                kL++;
                if (k == j)
                  { double Ljj = Lrs->va;
                    affirm(Ljj != 0.0, "zero element in diagonal");
                    Li.e[j] = (Ai.e[j] - sum)/Ljj;
                    break;
                  } 
                else
                  { double Ljk = Lrs->va;
                    double term = Li.e[k]*Ljk;
                    /* Kahan's summation formula: */
                    double tcorr = term - corr;
                    double newSum = sum + tcorr;
                    corr = (newSum - sum) - tcorr;
                    sum = newSum;
                    affirm (k < j, "Cholesky's upper half nonzero");
                  }
              }
          }
        
        /* Compute {L[i,i] = sqrt(L[i,i] - SUM{ L[i,j]^2 : j = 0..i-1 })} */
        { double sum = 0.0, corr = 0.0;
          for (j = 0; j < i; j++)
            { double w = Li.e[j], term = w*w; 
              /* Kahan's summation formula: */
              double tcorr = term - corr;
              double newSum = sum + tcorr;
              corr = (newSum - sum) - tcorr;
              sum = newSum;
            }
          Lii = Ai.e[i] - sum;
        }
        if (Lii < 0.0)
          { fprintf(stderr, "matrix is not positive definite?\n");
            fprintf(stderr, "row %04d Lii = %24.15e\n", i, Lii);
            for (j = 0; j <= i; j++)
              { if (Li.e[j] != 0.0) 
                  { fprintf(stderr, "  %04d %24.15e\n", j, Li.e[j]); }
              }
            fprintf(stderr, "\n");
            affirm(FALSE, "aborted");
          }
        Li.e[i] = sqrt(Lii);
        /* Save diagonal elements: */
        dg.e[i] = Li.e[i];
        /* Cleanup small entries: */
        for (j = 0; j < i; j++)
          { if ((minVal > 0.0) && (fabs(Li.e[j]) < minVal*sqrt(dg.e[i]*dg.e[j])))
              { Li.e[j] = 0.0; }
          }
        for (j = i+1; j < A.cols; j++) { Li.e[j] = 0.0; }
        SOMatrix_PackRow(Li, i, &Lents, &nL);
      }
    free(Ai.e); free(Li.e); free(dg.e);
    MatEntry_vec_trim(&Lents, nL);
    return (SOMatrix){A.rows, A.cols, Lents};
  }
  
void SOMatrix_GaussSeidel(SOMatrix A, double_vec_t b, double omega, double_vec_t x)
  { MatEntry *Ael = A.ents.e; int nA = A.ents.ne;
    
    int i, kA = 0;
    double Aii;

    affirm(A.rows == A.cols, "matrix must be square");
    affirm(b.ne == A.cols, "vector has wrong size");
    affirm(x.ne == A.cols, "vector has wrong size");
    for (i = 0; i < A.rows; i++)
      { double s = 0.0;
        while ((kA < nA) && (Ael[kA].row  == i))
          { MatEntry *Aij = &(Ael[kA]);
            nat_t j = Aij->col;
            if (j == i)
              { Aii = Aij->va; }
            else
              { s += Aij->va * x.e[j]; }
            kA++;
          }
        x.e[i] = (1.0 - omega) * x.e[i] + omega*(b.e[i] - s)/Aii;
      }
  }

double SOMatrix_ConjugateGradient(SOMatrix A, double_vec_t b, double_vec_t x)
  {
    double rsq;
    double_vec_t xi = double_vec_new(A.cols);
    double_vec_t xj = double_vec_new(A.cols);
    double_vec_t g = double_vec_new(A.cols);
    double_vec_t h = double_vec_new(A.cols);
    double eps = 1.0e-6;

    auto bool_t ConjugateGradientIteration(void);
      /* Performs one `macro' iteration of the Conjugate Gradient
        algorithm for solving the linear system {A x == b}.
        Returns TRUE iff the iteration converged. */  
    
    bool_t ConjugateGradientIteration(void) 
      {
        double rp = 0.0, bsq = 0.0;
        double anum, aden, gam, dgg, gg;
        int iter, j;
        SOMatrix_MulCol(A, x, xi);
        for (j = 0; j < A.cols; j++)
          { bsq += b.e[j]*b.e[j];
            xi.e[j] -= b.e[j];
            rp += xi.e[j]*xi.e[j];
          }
        SOMatrix_MulRow(xi, A, g);
        for (j = 0; j < A.cols; j++)
          { g.e[j] = -g.e[j];
            h.e[j] = g.e[j];
          }
        for (iter = 1; iter <= 10*A.cols; iter++)
          { fprintf(stderr, "*");
            SOMatrix_MulCol(A, h, xi);
            anum = 0.0;
            aden = 0.0;
            for (j = 0; j < A.cols; j++)
              { anum = anum + g.e[j]*h.e[j];
                aden = aden + xi.e[j]*xi.e[j];
              }
            affirm(aden != 0.0, "divide by zero"); 
            anum = anum/aden;
            for (j = 0; j < A.cols; j++)
              { xi.e[j] = x.e[j];
                x.e[j] += anum * h.e[j];
              }
            SOMatrix_MulCol(A, x, xj);
            rsq = 0.0;
            for (j = 0; j < A.cols; j++)
              { xj.e[j] -= b.e[j];
                rsq += xj.e[j]*xj.e[j];
              }
            if ((rsq == rp) || (rsq <= bsq*((double)A.cols)*eps*eps))
              { fprintf(stderr, "1");
                return TRUE ;
              }
            if (rsq > rp)
              { x = xi;
                fprintf(stderr, "2");
                return FALSE ;
              }
            rp = rsq;
            SOMatrix_MulRow(xj, A, xi);
            gg = 0.0;
            dgg = 0.0;
            for (j = 0; j < A.cols; j++)
              { gg = gg + g.e[j]*g.e[j];
                dgg = dgg + (xi.e[j] + g.e[j])*xi.e[j];
              }
            if (gg == 0.0)
              { fprintf(stderr, "3");
                return TRUE;
              }
            gam = dgg/gg;
            for (j = 0; j < A.cols; j++)
              { g.e[j] = -xi.e[j];
                h.e[j] = g.e[j] + gam*h.e[j];
              }
          }
        fprintf(stderr, "4");
        return TRUE;
      }

    int i;
    
    affirm(A.rows == A.cols , "matrix must be square");
    affirm(b.ne == A.cols, "wrong dimensions");
    affirm(x.ne == A.cols, "wrong dimensions");
    for (i = 1; i <= 3; i++)
      { fprintf(stderr, "!");
        if (ConjugateGradientIteration()) { return rsq; }
      }
    free(xi.e); free(xj.e);
    free(g.e); free(h.e);
    return rsq;
  }

nat_t SOMatrix_Find(MatEntry *r, nat_t n, nat_t row, nat_t col)
  { /* Binary search in the element list. */
    int p = 0, q = n-1;
    while (p <= q)
      { int m = (p + q)/2;
        MatEntry *rm = &(r[m]);
        if ((rm->row < row) || ((rm->row == row) && (rm->col < col)))
          { p = m + 1; }
        else if ((rm->row > row) || ((rm->row == row) && (rm->col > col)))
          { q = m - 1; }
        else
          { return m; }
      }
    return INT_MAX;
  }

double SOMatrix_GetValue(SOMatrix A, nat_t row, nat_t col)
  {
    nat_t k = SOMatrix_Find(A.ents.e, A.ents.ne, row, col);
    return (k == INT_MAX ? 0.0 : A.ents.e[k].va);
  }

#define SOMatrix_FileFormat "2002-11-20"
    
void SOMatrix_Write(FILE *wr, SOMatrix A)
  { MatEntry *Ael = A.ents.e; int nA = A.ents.ne;
    int kA;
    filefmt_write_header(wr, "SOMatrix", SOMatrix_FileFormat);
    fprintf(wr, "rows = %d\n", A.rows);
    fprintf(wr, "cols = %d\n", A.cols);
    fprintf(wr, "elems = %d\n", nA);
    for (kA = 0; kA < nA; kA++) 
      { MatEntry *Aij = &(Ael[kA]);
        fprintf(wr, "%d %d %24.16e\n", Aij->row, Aij->col, Aij->va);
      }
    filefmt_write_footer(wr, "SOMatrix");
    fflush(wr);
  }

SOMatrix SOMatrix_Read(FILE *rd)
  { SOMatrix A;
    int nA;
    filefmt_read_header(rd, "SOMatrix", SOMatrix_FileFormat);
    A.rows = nget_int(rd, "rows"); fget_eol(rd); 
    A.cols = nget_int(rd, "cols"); fget_eol(rd);
    nA = nget_int(rd, "elems"); fget_eol(rd);
    { MatEntry_vec_t Aents = MatEntry_vec_new(nA);
      int kA;
      for (kA = 0; kA < nA; kA++)
        { Aents.e[kA].row = fget_int(rd); 
          Aents.e[kA].col = fget_int(rd); 
          Aents.e[kA].va = fget_double(rd);
          fget_eol(rd);
        }
      A.ents = Aents;
    }  
    filefmt_read_footer(rd, "SOMatrix");
    return A;
  }

void SOMatrix_UnpackRow(MatEntry_vec_t *Aents, nat_t row, double_vec_t Ai, nat_t *kA)
  { int j;
    MatEntry *Aij;
    for (j = 0; j < Ai.ne; j++) { Ai.e[j] = 0.0; }
    while (((*kA) < Aents->ne) && ((Aij = &(Aents->e[(*kA)]))->row == row))
      { Ai.e[Aij->col] = Aij->va; (*kA)++; }
  }

void SOMatrix_UnpackHalfRow(MatEntry_vec_t *Lents, nat_t row, double_vec_t Li, nat_t *kL)
  { int j;
    MatEntry *Lij;
    for (j = 0; j < Li.ne; j++) { Li.e[j] = 0.0; }
    while (((*kL) < Lents->ne) && ((Lij = &(Lents->e[(*kL)]))->row == row))
      { Li.e[Lij->col] = Lij->va;
        affirm(Lij->col <= row, "matrix is not lower triangular");
        (*kL)++;
      }
  }

void SOMatrix_PackRow(double_vec_t Ci, nat_t row, MatEntry_vec_t *Cents, nat_t *nC)
  { int j;
    for (j = 0; j < Ci.ne; j++)
      { if (Ci.e[j] != 0.0)
          { MatEntry_vec_expand(Cents, *nC);
            { MatEntry *Cij = &(Cents->e[*nC]);
              Cij->row = row; Cij->col = j;
              Cij->va = Ci.e[j];
              (*nC)++;
            }
          }
      }
  }

/* Arrays of {MatEntry}: */

vec_typeimpl(MatEntry_vec_t,MatEntry_vec,MatEntry);
