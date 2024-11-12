/* See SPMatrix.h */
/* Last edited on 2023-02-12 07:50:16 by stolfi */

#include <SPMatrix.h>
#include <SPVector.h>

#include <SPBasic.h>

#include <vec.h>
#include <filefmt.h>
#include <fget.h>
#include <nget.h>
#include <affirm.h>
#include <nat.h>
#include <bool.h>

// #include <gsl/gsl_vector_double.h>
// #include <gsl/gsl_matrix_double.h>
// #include <gsl/gsl_linalg.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#include <limits.h>

/* !!! Replace by {spmat_double.c} !!! */
  
/* INTERNAL PROTOTYPES */

void SPMatrix_UnpackRow(MatEntry_vec_t *Aents, nat_t row, SPVector Ai, nat_t *kA);
  /* Unpacks row {i} of the matrix {A} as an ordinary vector of double
    {Ai}. On input, the variable {*kA} must be the index into {Aents ==
    A.ents} of the first element of that row. On output, {*kA} will be
    the index of the first element of row {i+1}. */

void SPMatrix_UnpackHalfRow(MatEntry_vec_t *Lents, nat_t row, SPVector Li, nat_t *kL);
  /* Same as {UnpackRow}, but expects and checks that {L} is lower triangular ---
    meaning that there are no elements {L[i,j]} with {j > i}. */

void SPMatrix_PackRow(SPVector Ci, nat_t row, MatEntry_vec_t *Cents, nat_t *nC);
  /* Appends the non-zero entries of {Ci} as one more row of elements
    of a matrix {C}. Will reallocate the element list {Cents == C.ents} if
    necessary. On input, the variable {*nC} must be the count of
    elements already stored in {Cents^}; this count will be updated on
    exit. Don't forget to trim the element list after adding the last
    row. */

/* IMPLEMENTATIONS */

SPMatrix SPMatrix_Null(nat_t rows, nat_t cols)
  { MatEntry_vec_t ents = MatEntry_vec_new(0);
    return (SPMatrix){rows, rows, ents};
  }

SPMatrix SPMatrix_BinaryOp
  ( SPMatrix A, SPMatrix B,
    int dRow, int dCol,
    double func(int row, int col, double va, double vb),
    bool_t skip00
  )
  { MatEntry *Ael = A.ents.e; int nA = A.ents.ne;
    MatEntry *Bel = B.ents.e; int nB = B.ents.ne;
    MatEntry_vec_t Cents = MatEntry_vec_new(nA + nB);
    nat_t kA = 0, kB = 0, nC = 0;
    nat_t row, col;
    row = -1; col = A.cols-1;
    while (TRUE)
      { /* Get indices {Arow,Acol} next {A} element: */
        int Arow, Acol;
        if (kA >= nA)
          { Arow = A.rows; Acol = 0; }
        else
          { Arow = Ael[kA].row; Acol = Ael[kA].col; }
        
        /* Get indices {Brow,Bcol} of next valid element of shifted {B}: */
        int Brow, Bcol;
        while (TRUE)
          { if (kB >= nB)
              { Brow = A.rows; Bcol = 0; break; }
            else
              { Brow = Bel[kB].row - dRow; 
                Bcol = Bel[kB].col - dCol;
                if (Brow >= A.rows) 
                  { /* Exhausted {B}: */ break; }
                else if ((Brow >= 0) && (Bcol >= 0) && (Bcol < A.cols))
                  { /* Found a valid {B} element: */ break; }
                else
                  { /* Out-of-bounds element, skip: */ kB++; }
              }
          }
        
        /* Get indices {[row,col]} of next element to be computed: */
        if (! skip00)
          { /* Must evaluate even elements which are missing in {A} and {B}: */
            col++; 
            if (col >= A.cols) { row++; col = 0; }
          }
        else
          { /* Pick the nearest element of {A} or {B}: */
            if (Arow < Brow)
              { row = Arow; col = Acol; }
            else if (Brow < Arow) 
              { row = Brow; col = Bcol; }
            else if (Acol < Bcol) 
              { row = Arow; col = Acol; }
            else 
              { row = Brow; col = Bcol; }
          }
          
        /* Are we done already? */
        if (row >= A.rows) { break; }
       
        /* Get element values from {B} and {A}, and skip over them: */
        double va, vb;
        if ((row == Arow) && (col == Acol))
          { va = Ael[kA].va; kA++; }
        else
          { va = 0.0; }
        if ((row == Brow) && (col == Bcol))
          { vb = Bel[kB].va; kB++; }
        else
          { vb = 0.0; }
        
        /* Compute element of result: */
        double vc;
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
    return (SPMatrix){A.rows, A.cols, Cents};
  }

SPMatrix SPMatrix_MixBlock 
  ( double alpha, 
    SPMatrix A, 
    double beta, 
    SPMatrix B, 
    int dRow, 
    int dCol
  )
  { 
    auto double mix(int row, int col, double va, double vb);
    
    double mix(int row, int col, double va, double vb)
      { return alpha*va + beta*vb; }
      
    return SPMatrix_BinaryOp(A, B, dRow, dCol, mix, TRUE);
  }
 
SPMatrix SPMatrix_Mix(double alpha, SPMatrix A, double beta, SPMatrix B)
  { 
    auto double mix(int row, int col, double va, double vb);
    
    double mix(int row, int col, double va, double vb)
      { return alpha*va + beta*vb; }
      
    affirm(A.rows == B.rows, "incompatible row count");
    affirm(A.cols == B.cols, "incompatible col count");
    return SPMatrix_BinaryOp(A, B, 0, 0, mix, TRUE);
  }
 
SPMatrix SPMatrix_RelDiff(SPMatrix A, SPMatrix B)
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
      
    return SPMatrix_BinaryOp(A, B, 0, 0, rdiff, TRUE);
  }

SPMatrix SPMatrix_Scale(double alpha, SPMatrix A)
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
    return (SPMatrix){A.rows, A.cols, Cents};
  }

void SPMatrix_GetDiagonal(SPMatrix A, SPVector d)
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

SPMatrix SPMatrix_Mul(SPMatrix A, SPMatrix B)
  { int nA = A.ents.ne;
    int nB = B.ents.ne;
    nat_t maxnel = (nA > nB ? nA : nB);
    MatEntry_vec_t Cents = MatEntry_vec_new(maxnel);
    SPVector a = double_vec_new(A.rows);
    SPVector c = double_vec_new(B.cols);
    int i; nat_t kA = 0, nC = 0;

    affirm(A.cols == B.rows , "incompatible sizes");
    for (i = 0; i < A.rows; i++)
      { SPMatrix_UnpackRow(&A.ents, i, a, &kA);
        SPMatrix_MulRow(a, B, c);
        SPMatrix_PackRow(c, i, &Cents, &nC);
      }
    MatEntry_vec_trim(&Cents, nC);
    free(a.e); free(c.e);
    return (SPMatrix){A.rows, B.cols, Cents};
  }

void SPMatrix_MulRow(SPVector x, SPMatrix A, SPVector y)
  { MatEntry *Ael = A.ents.e; int nA = A.ents.ne;
    int j, kA;
    affirm(x.ne == A.rows, "vector x has wrong size");
    affirm(y.ne == A.cols, "vector y has wrong size");
    for (j = 0; j < A.cols; j++) { y.e[j] = 0.0; }
    for (kA = 0; kA < nA; kA++)
      { MatEntry *Aij = &(Ael[kA]);
        nat_t i = Aij->row, j = Aij->col;
        y.e[j] += x.e[i] * Aij->va;
      }
  }

void SPMatrix_MulCol(SPMatrix A, SPVector x, SPVector y)
  { MatEntry *Ael = A.ents.e; int nA = A.ents.ne;
    int i, kA;
    affirm(x.ne == A.cols, "vector x has wrong size");
    affirm(y.ne == A.rows, "vector y has wrong size");
    for (i = 0; i < y.ne; i++) { y.e[i] = 0.0; }
    for (kA = 0; kA < nA; kA++)
      { MatEntry *Aij = &(Ael[kA]);
        nat_t i = Aij->row, j = Aij->col;
        y.e[i] += Aij->va * x.e[j];
      }
  }

void SPMatrix_TrMulCol(SPMatrix A, SPVector x, SPVector y)
  { MatEntry *Ael = A.ents.e; int nA = A.ents.ne;
    int i, kA;
    affirm(x.ne == A.rows, "vector x has wrong size");
    affirm(y.ne == A.cols, "vector y has wrong size");
    for (i = 0; i < y.ne; i++) { y.e[i] = 0.0; }
    for (kA = 0; kA < nA; kA++)
      { MatEntry *Aij = &(Ael[kA]);
        nat_t i = Aij->row, j = Aij->col;
        y.e[j] += Aij->va * x.e[i];
      }
  }

void SPMatrix_DivRow(SPVector y, SPMatrix L, SPVector x)
  { MatEntry *Lel = L.ents.e; int nL = L.ents.ne;
    int i, kL;

    affirm(L.rows == L.cols, "matrix is not square");
    affirm(y.ne == L.cols, "vector y has wrong size");
    affirm(x.ne == L.cols, "vector x has wrong size");

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

void SPMatrix_DivCol(SPMatrix L, SPVector y, SPVector x)
  { int i, j, kL;
    SPVector Li = double_vec_new(L.cols);
    
    affirm(L.rows == L.cols, "matrix is not square");
    affirm(y.ne == L.cols, "vector y has wrong size");
    affirm(x.ne == L.cols, "vector x has wrong size");
    
    kL = 0;
    for (i = 0; i < L.rows; i++)
      { SPMatrix_UnpackHalfRow(&L.ents, i, Li, &kL);
        x.e[i] = y.e[i];
        for (j = 0; j < i; j++)
          { x.e[i] -= Li.e[j]*x.e[j]; }
        x.e[i] /= Li.e[i];
      }
    free(Li.e);
  }

SPMatrix SPMatrix_PosDiv(SPMatrix L, SPMatrix A)
  { int i, kA, nX;
    MatEntry_vec_t Xents = MatEntry_vec_new(A.ents.ne);
    SPVector Ai = double_vec_new(A.cols);
    SPVector Xi = double_vec_new(A.cols);
    
    affirm(L.rows == L.cols, "L matrix is not square");
    affirm(A.cols == L.rows, "incompatible matrix sizes");
    
    kA = 0; nX = 0;
    for (i = 0; i < A.rows; i++)
      { SPMatrix_UnpackRow(&(A.ents), i, Ai, &kA);
        SPMatrix_DivRow(Ai, L, Xi);
        SPMatrix_PackRow(Xi, i, &Xents, &nX);
      }
    free(Ai.e); free(Xi.e);
    MatEntry_vec_trim(&Xents, nX);
    return (SPMatrix){A.rows, A.cols, Xents}; 
  }

SPMatrix SPMatrix_PreDiv(SPMatrix L, SPMatrix A)
  { int i, kL, kA, kX, nX;
    MatEntry_vec_t Xents = MatEntry_vec_new(A.ents.ne);
    SPVector Li = double_vec_new(L.cols);
    SPVector Xi = double_vec_new(A.cols);
    
    affirm(L.rows == L.cols, "L matrix is not square");
    affirm(A.rows == L.cols, "A matrix has wrong rows");
    
    kL = 0; kA = 0; nX = 0;
    for (i = 0; i < L.rows; i++)
      { SPMatrix_UnpackRow(&(L.ents), i, Li, &kL);
        SPMatrix_UnpackRow(&(A.ents), i, Xi, &kA);
        kX = 0;
        while(kX < nX)
          { MatEntry *Xjk = &(Xents.e[kX]);
            int j = Xjk->row, k = Xjk->col;
            Xi.e[k] -= Li.e[j] * Xjk->va;
          }
        { int k; for (k = 0; k < A.cols; k++) { Xi.e[k] /= Li.e[i]; } }
        SPMatrix_PackRow(Xi, i, &Xents, &nX);
      }
    free(Li.e); free(Xi.e);
    MatEntry_vec_trim(&Xents, nX);
    return (SPMatrix){A.rows, A.cols, Xents}; 
  }

void SPMatrix_Cholesky(SPMatrix A, double minVal, SPMatrix *L)
  { 
    /* Andre-Louis Cholesky (spelled with a 'y'), born in France in
      1875, was a geodesist in the French military. He developed his
      computational procedure to solve geodetic problems, among which
      was, in the words of his eulogy, the {problem of levelling in
      Morocco.} He died in battle in 1918. Benoit published the
      computational method in {Bulletin geodesique} in 1924.

      The most likely French pronunciation of Cholesky is that heard
      in 'Chopin' -- the English 'sh' sound. 
        -- David Pattison, Washington, DC */
  
    MatEntry_vec_t Lents = MatEntry_vec_new(A.rows);
    int i, j, nL = 0, kL, kA = 0;
    double Lii;
    
    SPVector Ai = double_vec_new(A.rows);
    SPVector Li = double_vec_new(A.rows);
    SPVector dg = double_vec_new(A.rows); /* Diagonal elements */
    
    affirm(A.rows == A.cols , "matrix must be square");
    for (i = 0; i < A.rows; i++)
      { SPMatrix_UnpackRow(&A.ents, i, Ai, &kA);
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
        SPMatrix_PackRow(Li, i, &Lents, &nL);
      }
    free(Ai.e); free(Li.e); free(dg.e);
    MatEntry_vec_trim(&Lents, nL);
    (*L) = (SPMatrix){A.rows, A.cols, Lents};
  }
  
// void SPMatrix_SVD(SPMatrix A, double minVal, SPMatrix *L, SPVector *D, SPMatrix *R)
//   { 
//     /* Unpack matrix {A} as a GSL matrix {GA}: */
//     gsl_matrix *GA = SPMatrix_To_GSLMatrix(&A);
//      
//     /* Allocate the factors {GD,GV} of the SVD, and the work vector {GW}: */
//     gsl_matrix *GV = gsl_matrix_alloc(A.cols, A.cols); /* Yes, {cols} by {cols}. */
//     gsl_vector *GS = gsl_vector_alloc(A.cols);
//     gsl_vector *GW = gsl_vector_alloc(A.cols);
// 
//     /* Compute the SVD decomposition: */
//     int res = gsl_linalg_SV_decomp (GA, GV, GS, GW);
//     affirm(res == 0, "gsl_linalg_SV_decomp failed");
//     gsl_vector_free(GW);
//     
//     /* Eliminate the small singular values from {GD} and find the  rank {p}: */
//     int p = A.cols;
//     if (minVal > 0.0)
//       { 
//         double sMax = gsl_vector_get(GS,0);
//         affirm(sMax >= 0, "singular value is negative");
//         while(p > 0)
//           { double sMin = gsl_vector_get(GS,p-1);
//             if ((sMin > 0.0) && (sMin >= minVal*sMax)) { break; }
//             p--; 
//           }
//       } 
//     
//     /* Copy the first {p} columns of {GA} to {L}: */
//     { gsl_matrix_view GAp = gsl_matrix_submatrix(GA, 0, 0, A.rows, p);
//       (*L) = SPMatrix_From_GSLMatrix(&(GAp.matrix));
//     }
//     gsl_matrix_free(GA);
//     
//     /* Copy the first {p} columns of {GV} to {R}: */
//     { gsl_matrix_view GVp = gsl_matrix_submatrix(GV, 0, 0, A.cols, p);
//       (*R) = SPMatrix_From_GSLMatrix(&(GVp.matrix));
//     }
//     gsl_matrix_free(GV);
//     
//     /* Copy the first {p} elements of vector {GS} to {D}: */
//     { gsl_vector_view GSp = gsl_vector_subvector(GS, 0, p);
//       (*D) = SPVector_From_GSLVector(&(GSp.vector));
//     }
//     gsl_vector_free(GS);
//   }
  
void SPMatrix_GaussSeidel(SPMatrix A, SPVector b, double omega, SPVector x)
  { 
    /* Gauss needs no introduction, right? */
  
    MatEntry *Ael = A.ents.e; int nA = A.ents.ne;
    
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

double SPMatrix_ConjugateGradient(SPMatrix A, SPVector b, SPVector x)
  {
    double rsq;
    SPVector xi = double_vec_new(A.cols);
    SPVector xj = double_vec_new(A.cols);
    SPVector g = double_vec_new(A.cols);
    SPVector h = double_vec_new(A.cols);
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
        SPMatrix_MulCol(A, x, xi);
        for (j = 0; j < A.cols; j++)
          { bsq += b.e[j]*b.e[j];
            xi.e[j] -= b.e[j];
            rp += xi.e[j]*xi.e[j];
          }
        SPMatrix_MulRow(xi, A, g);
        for (j = 0; j < A.cols; j++)
          { g.e[j] = -g.e[j];
            h.e[j] = g.e[j];
          }
        for (iter = 1; iter <= 10*A.cols; iter++)
          { fprintf(stderr, "*");
            SPMatrix_MulCol(A, h, xi);
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
            SPMatrix_MulCol(A, x, xj);
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
            SPMatrix_MulRow(xj, A, xi);
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

nat_t SPMatrix_Find(MatEntry *r, nat_t n, nat_t row, nat_t col)
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

double SPMatrix_GetValue(SPMatrix A, nat_t row, nat_t col)
  {
    nat_t k = SPMatrix_Find(A.ents.e, A.ents.ne, row, col);
    return (k == INT_MAX ? 0.0 : A.ents.e[k].va);
  }

#define SPMatrix_FileFormat "2002-11-20"
    
void SPMatrix_Write(FILE *wr, SPMatrix A)
  { MatEntry *Ael = A.ents.e; int nA = A.ents.ne;
    int kA;
    filefmt_write_header(wr, "SPMatrix", SPMatrix_FileFormat);
    fprintf(wr, "rows = %d\n", A.rows);
    fprintf(wr, "cols = %d\n", A.cols);
    fprintf(wr, "elems = %d\n", nA);
    for (kA = 0; kA < nA; kA++) 
      { MatEntry *Aij = &(Ael[kA]);
        fprintf(wr, "%d %d %24.16e\n", Aij->row, Aij->col, Aij->va);
      }
    filefmt_write_footer(wr, "SPMatrix");
    fflush(wr);
  }

SPMatrix SPMatrix_Read(FILE *rd)
  { SPMatrix A;
    int nA;
    filefmt_read_header(rd, "SPMatrix", SPMatrix_FileFormat);
    A.rows = nget_int32(rd, "rows"); fget_eol(rd); 
    A.cols = nget_int32(rd, "cols"); fget_eol(rd);
    nA = nget_int32(rd, "elems"); fget_eol(rd);
    { MatEntry_vec_t Aents = MatEntry_vec_new(nA);
      int kA;
      for (kA = 0; kA < nA; kA++)
        { Aents.e[kA].row = fget_int32(rd); 
          Aents.e[kA].col = fget_int32(rd); 
          Aents.e[kA].va = fget_double(rd);
          fget_eol(rd);
        }
      A.ents = Aents;
    }  
    filefmt_read_footer(rd, "SPMatrix");
    return A;
  }

/* PACKING AND UNPACKING */

void SPMatrix_UnpackRow(MatEntry_vec_t *Aents, nat_t row, SPVector Ai, nat_t *kA)
  { int j;
    MatEntry *Aij;
    for (j = 0; j < Ai.ne; j++) { Ai.e[j] = 0.0; }
    while (((*kA) < Aents->ne) && ((Aij = &(Aents->e[(*kA)]))->row == row))
      { Ai.e[Aij->col] = Aij->va; (*kA)++; }
  }

void SPMatrix_UnpackHalfRow(MatEntry_vec_t *Lents, nat_t row, SPVector Li, nat_t *kL)
  { int j;
    MatEntry *Lij;
    for (j = 0; j < Li.ne; j++) { Li.e[j] = 0.0; }
    while (((*kL) < Lents->ne) && ((Lij = &(Lents->e[(*kL)]))->row == row))
      { Li.e[Lij->col] = Lij->va;
        affirm(Lij->col <= row, "matrix is not lower triangular");
        (*kL)++;
      }
  }

void SPMatrix_PackRow(SPVector Ci, nat_t row, MatEntry_vec_t *Cents, nat_t *nC)
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

/* HANDLING ARRAYS OF {MatEntry}: */

vec_typeimpl(MatEntry_vec_t,MatEntry_vec,MatEntry);
