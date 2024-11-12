/* SOMatrix.h -- sparse matrices */
/* Last edited on 2007-01-04 00:20:32 by stolfi */

#ifndef SOMatrix_H
#define SOMatrix_H

#include <stdio.h>
#include <vec.h>

typedef struct MatEntry /* Entry (usually nonzero) in a sparse matrix */
  { nat_t row;  /* Row index (from 0). */ /* Row index (from 0). */
    nat_t col;  /* Colum index (from 0). */
    double va;  /* Element value. */
  } MatEntry;
  
vec_typedef(MatEntry_vec_t,MatEntry_vec,MatEntry);
  
typedef struct SOMatrix  
  { nat_t rows;    /* Number of rows. */
    nat_t cols;    /* Number of columns. */
    MatEntry_vec_t ents;  /* Non-zero entries, sorted by {row} then {col} */
  } SOMatrix;

SOMatrix SOMatrix_Null(nat_t rows, nat_t cols);  
  /* Returns the null matrix with given dimensions. */
  
void SOMatrix_GetDiagonal(SOMatrix A, double_vec_t d);
  /* Stores in {d[0..ne-1]} the diagonal of {A}, where 
    {d.ne = min(A.rows,A.cols)}. */

SOMatrix SOMatrix_Scale(double alpha, SOMatrix A);
  /* Returns the product of {alpha} by {A} */
    
SOMatrix SOMatrix_Mix(double alpha, SOMatrix A, double beta, SOMatrix B);  
  /* Returns the matrix  {alpha*A + beta*B}. */
  
SOMatrix SOMatrix_RelDiff(SOMatrix A, SOMatrix B);  
  /* Returns the matrix {M[i,j] = rdiff(A[i,j],B[i,j])}, where
    {rdiff(x,y) = (x - y)/sqrt(x^2 + y^2)/2}. This number is {0} when
    {x == y}, is {0.5} when only one of them is zero, and has extremal
    values {-1} and {-1} when {x} and {y} are equal and opposite). */
  
SOMatrix SOMatrix_Mul(SOMatrix A, SOMatrix B);
  /* Returns the matrix product {A*B}. */
  
SOMatrix SOMatrix_BinaryOp
  ( SOMatrix A, SOMatrix B,
    double func(int row, int col, double va, double vb),
    bool_t skip00
  );
  /* Computes the matrix {C[i,j] = func(i,j,A[i,j],B[i,j])}. 
    If {skip00 == TRUE}, assumes that {C[i,j]} is zero when
    {A[i,j]} and {B[i,j]} are both zero. */

void SOMatrix_MulRow(double_vec_t x, SOMatrix A, double_vec_t y);
  /* Stores in {y} the product of the row vector {x} by the matrix {A}. */
    
void SOMatrix_MulCol(SOMatrix A, double_vec_t x, double_vec_t y);
  /* Stores in {y} the product of the matrix {A} by the column vector {x}. */

/* MATRIX FACTORIZATION */

SOMatrix SOMatrix_Cholesky(SOMatrix A, double minVal);
 /* Factors the positive definite nxn matrix {A} into {L*L^t} where {L} 
   is lower triangular and {L^t} is the transpose of {L}.
   
   If {minVal} is positive, deletes any elements {L[i,j]}
   which are less than {sqrt(L[i,i]*L[j,j])} in absolute value.
   
   Named after André-Louis Cholesky (1875--1918), French military
   geodesist. Although the name is originally Polish or Russian, he
   probably pronounced it himself in the French fashion, i.e. with an
   `sh' sound, as in `Chopin' and `Chostakovich'; and not with a `tsh'
   sound --- which in French would be spelled `Tch', as in
   `Tchaikovski' or `Tchekov'. */

/* OPERATIONS ON LOWER TRIANGULAR MATRICES */

/* The operations in this section require the matrix {L} to 
   be lower triangular, with nonzero entries in the diagonal. */

void SOMatrix_DivRow(double_vec_t y, SOMatrix L, double_vec_t x); 
  /* Stores in {x} the row vector {y} times the inverse of {L}.
    In other words, finds the solution {x} to {x*L = y}. */
    
void SOMatrix_DivCol(SOMatrix L, double_vec_t y, double_vec_t x); 
  /* Stores in {x} the inverse of {L} times the column vector {y}.
    In other words, finds the solution {x} to {L*x = y}. */
    
SOMatrix SOMatrix_PosDiv(SOMatrix A, SOMatrix L);
  /* Returns the matrix {A} times the inverse of matrix {L}. */

SOMatrix SOMatrix_PreDiv(SOMatrix L, SOMatrix A);
  /* Returns the inverse of matrix {L} times the matrix {A}. */

/* OTHER OPERATIONS ON MATRICES */

void SOMatrix_GaussSeidel(SOMatrix A, double_vec_t b, double omega, double_vec_t x);
  /* Performs one iteration of the Gauss-Seidel algorithm for solving
    the linear system {A x == b}. Specifically, recomputes {x[i] =
    x[i] + omega*(b[i] - (SUM){a[i,j]*x[j]})/a[i,i]} for each {i},
    increasing. Note that {x[j]} in the summation is partly old and
    partly new. 
    
    The matrix {A} must be square. Setting {omega = 1.0} gives the
    ordinary method; {omega < 1.0} slows it down, {omega > 1.0} speeds
    it up. */

double SOMatrix_ConjugateGradient(SOMatrix A, double_vec_t b, double_vec_t x);
  /* Solves the linear system {A x == b} for {x} using the Conjugate
    Gradient algorithm [Numerical Recipes section 2.10, subroutine
    SPARSE]. Returns the squared norm of the residual {A x - b}. A
    large value means that {A} is singular and therefore {x} is only a
    least squares solution. */
  
nat_t SOMatrix_Find(MatEntry *r, nat_t n, nat_t row, nat_t col);
  /* Finds an index {p} in {0..n-1} such that {r[p].row = row} 
    and {r[p].col = col}. Assumes the elements in {r} are stored
    in row-major order. Returns a huge number ({INT_MAX}) if not found. */

double SOMatrix_GetValue(SOMatrix A, nat_t row, nat_t col);
  /* Returns the value of the element {A[row,col]}. In particular,
     returns 0 if there is no {MatEntry} with those indices. */

void SOMatrix_Write(FILE *wr, SOMatrix A);
  /* Write the matrix {A} to {wr}, in a format that can be read back. */
  
SOMatrix SOMatrix_Read(FILE *rd);
  /* Reads a matrix from {rd}, assuming it was written by {SOMatrix_Write}. */

#endif
 
