/* SPMatrix.h -- sparse matrices */
/* Last edited on 2009-02-09 20:43:50 by stolfi */

#ifndef SPMatrix_H
#define SPMatrix_H

#include <SPVector.h>
#include <vec.h>
#include <bool.h>
#include <nat.h>

/* #include <gsl/gsl_matrix_double.h> */
#include <stdio.h>

typedef struct MatEntry /* Entry (usually nonzero) in a sparse matrix */
  { nat_t row;    /* Row index (from 0). */
    nat_t col;    /* Colum index (from 0). */
    double va;  /* Element value. */
  } MatEntry;
  
vec_typedef(MatEntry_vec_t,MatEntry_vec,MatEntry);
  
typedef struct SPMatrix  
  { nat_t rows;            /* Number of rows. */
    nat_t cols;            /* Number of columns. */
    MatEntry_vec_t ents; /* Non-zero entries, sorted by {row} then {col} */
  } SPMatrix;

SPMatrix SPMatrix_Null(nat_t rows, nat_t cols);  
  /* Returns the null matrix with given dimensions. */
  
void SPMatrix_GetDiagonal(SPMatrix A, SPVector d);
  /* Stores in {d[0..} the diagonal of {A}, where {d.ne = min(A.rows,A.cols)}. */

SPMatrix SPMatrix_Scale(double alpha, SPMatrix A);
  /* Returns the product of {alpha} by {A} */
     
SPMatrix SPMatrix_Mix(double alpha, SPMatrix A, double beta, SPMatrix B);
  /* Returns the matrix  {alpha*A + beta*B}. The two matrices must have
    the same index bounds. */
   
SPMatrix SPMatrix_MixBlock
  ( double alpha, 
    SPMatrix A, 
    double beta, 
    SPMatrix B, 
    int dRow, 
    int dCol
  );
  /* Returns the matrix  {alpha*A + beta*B}, where {B}
    is shifted by {dRow,dCol} and trimmed and/or padded
    with zeros so as to have the same bounds as {A}. */
  
SPMatrix SPMatrix_RelDiff(SPMatrix A, SPMatrix B);  
  /* Returns the matrix {M[i,j] = rdiff(A[i,j],B[i,j])}, where
    {rdiff(x,y) = (x - y)/sqrt(x^2 + y^2)/2}. This number is {0} when
    {x == y}, is {0.5} when only one of them is zero, and has extremal
    values {-1} and {-1} when {x} and {y} are equal and opposite). */
  
SPMatrix SPMatrix_Mul(SPMatrix A, SPMatrix B);
  /* Returns the matrix product {A*B}. */
  
SPMatrix SPMatrix_BinaryOp
  ( SPMatrix A, SPMatrix B,
    int dRow, int dCol,
    double func(int row, int col, double va, double vb),
    bool_t skip00
  );
  /* A generic elementwise binary operation that combines {A} with {B}
    shifted by {dRow,dCol}. Namely, returns the matrix {C}, with the
    same size as {A}, such that {C[i,j] = func(i,j,A[i,j],B[i',j'])},
    where {i'= i-dRow} and {j' = j-dCol)}. Elements {B[i',j']} that
    don't exist are assumed to be zero, and those outside the bounds
    of {A} are ignored. If {skip00 == TRUE}, assumes that {C[i,j]} is
    zero when {A[i,j]} and {B[i',j']} are both zero. */

void SPMatrix_MulRow(SPVector x, SPMatrix A, SPVector y);
  /* Stores in {y} the product of the row vector {x} by the matrix {A}. */
    
void SPMatrix_MulCol(SPMatrix A, SPVector x, SPVector y);
  /* Stores in {y} the product of the matrix {A} by the column vector {x}. */

void SPMatrix_TrMulCol(SPMatrix A, SPVector x, SPVector y);
  /* Stores in {y} the product of the transpose {A^T} of matrix {A} 
    by the column vector {x}. */

/* MATRIX FACTORIZATION */

void SPMatrix_Cholesky(SPMatrix A, double minVal, SPMatrix *L);
 /* Factors the positive definite nxn matrix {A} into {L*L^t} where {L} 
   is lower triangular and {L^t} is the transpose of {L}.
   
   If {minVal} is positive, deletes any off-diagonal elements {L[i,j]}
   which are less than {sqrt(L[i,i]*L[j,j])} in absolute value.
   
   Named after André-Louis Cholesky (1875--1918), French military
   geodesist. Although the name is originally Polish or Russian, he
   probably pronounced it himself in the French fashion, i.e.
   with an {sh} sound, as in {Chopin} and {Chostakovich}; and not
   with a {tsh} sound --- which in French would be spelled {Tch},
   as in {Tchaikovski} or {Tchekov}. */

void SPMatrix_SVD(SPMatrix A, double minVal, SPMatrix *L, SPVector *D, SPMatrix *R);
 /* Factors the {m} by {n} matrix {A} into {L*D*R^T} where {L} and {R} 
   are orthogonal matrices, of size {m} by {p} and {n} by {p},
   respectively; and {D} is a {p} by {p} diagonal matrix (represented as 
   a vector of {p} elements).
   
   The elements of {D} are non-negative and sorted in decreasing order.
   Their number {p} is usually {min(m,n)}. However, if {minVal} is positive, 
   any elements of {D} which are zero, or less than {minVal} times the 
   largest element, are discarded, and {p} is reduced acordingly. */
   
/* OPERATIONS ON LOWER TRIANGULAR MATRICES */

/* The operations in this section require the matrix {L} to 
   be lower triangular, with nonzero entries in the diagonal. */

void SPMatrix_DivRow(SPVector y, SPMatrix L, SPVector x); 
  /* Stores in {x} the row vector {y} times the inverse of {L}.
    In other words, finds the solution {x} to {x*L = y}. */
    
void SPMatrix_DivCol(SPMatrix L, SPVector y, SPVector x); 
  /* Stores in {x} the inverse of {L} times the column vector {y}.
    In other words, finds the solution {x} to {L*x = y}. */
    
SPMatrix SPMatrix_PosDiv(SPMatrix A, SPMatrix L);
  /* Returns the matrix {A} times the inverse of matrix {L}. */

SPMatrix SPMatrix_PreDiv(SPMatrix L, SPMatrix A);
  /* Returns the inverse of matrix {L} times the matrix {A}. */

/* OTHER OPERATIONS ON MATRICES */

void SPMatrix_GaussSeidel(SPMatrix A, SPVector b, double omega, SPVector x);
  /* Performs one iteration of the Gauss-Seidel algorithm for solving
    the linear system {A x == b}. Specifically, recomputes {x[i] =
    x[i] + omega*(b[i] - (SUM){a[i,j]*x[j]})/a[i,i]} for each {i},
    increasing. Note that {x[j]} in the summation is partly old and
    partly new. 
    
    The matrix {A} must be square. Setting {omega == 1.0} gives the
    ordinary method; {omega < 1.0} slows it down, {omega > 1.0} speeds
    it up. */

double SPMatrix_ConjugateGradient(SPMatrix A, SPVector b, SPVector x);
  /* Solves the linear system {A x  == b} using the Conjugate Gradient algorithm 
    [Numerical Recipes section 2.10, subroutine SPARSE].
    Returns the squared norm of the residual {A x - b}. A large
    value means that {A} is singular and therefore {x} is only a 
    least squares solution. */  
  
nat_t SPMatrix_Find(MatEntry *r, nat_t n, nat_t row, nat_t col);
  /* Finds an index {p} in {0..n-1} such that {r[p].row = row} 
    and {r[p].col = col}. Assumes the elements in {r} are stored
    in row-major order. Returns a huge number ({INT_MAX}) if not found. */

double SPMatrix_GetValue(SPMatrix A, nat_t row, nat_t col);
  /* Returns the value of the element {A[row,col]}. In particular,
     returns 0 if there is no {MatEntry} with those indices. */

void SPMatrix_Write(FILE *wr, SPMatrix A);
  /* Write  the  matrix {A} to {wr}, in a format that can be read back. */
  
SPMatrix SPMatrix_Read(FILE *rd);
  /* Reads a matrix from {rd}, assuming it was written by {SPMatrix_Write}. */

#endif
 
