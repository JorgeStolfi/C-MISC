
#ifndef SparseSymetricMatrix_H
#define SparseSymetricMatrix_H



/*
   Defines a data structure for a sparse symetric square matrix.
   The structure is a vector of vectors, one for each row.
   Each entry of a secondary vector contains the column index (>= 0)
   and value of a non-zero element of the corresponding matrix row.
   The master vector also contains, for each row
   the index in the secondary vector containing the diagonal element
   of that row.
 */
 
#include <Wr.h>
#include <Rd.h>

typedef 
  Vector == ARRAY OF double;
  MatrixElement == struct ?? { nat j, a: double; }
  Row == ARRAY OF MatrixElement;
  T == ARRAY OF  struct ?? { nat d, row: REF Row; }
 
T New(nat n);
  /*
    Allocates a new "n" by "n" matrix, initialized to the identity matrix. */

double GetValue(READONLY T m, i, nat j);
  /*
    Gets element "[i,j]" of the matrix (by linear scanning of row i). */
    
void PutRow(VAR T m, nat i, READONLY Vector row);
  /*
    Replaces row "i" of "m" by a (compressed) copy of "row". */

void Remove(VAR T m, bool_vec sel);
  /*
    Sets to zero the rows and columns of "m" with indices "k"
    such that "sel[k] == TRUE"; except for the diagonal elements,
    which are set to 1. */
    
void SubMul(READONLY T m, READONLY Vector v, VAR Vector x); 
  /* 
    Computes x = x - m * v */
    
PROCEDURE Factor(READONLY T m): REF T;
  /*
    Returns the Cholesky factorization of "m == LDL^T" 
    Valid only if "m" is positive-definite. */
    
void Solve(READONLY T mf, READONLY Vector v, VAR Vector x, i0 = 0);
  /*
    Solves "m*x == v" for "x", given the Choleski factorization "mf" of "m".
    If "i0" > 0, ... ? */

void Print(FILE *wr, READONLY T m);
  /*
    Prints "m" to "wr". */
    
void Write(FILE *wr, READONLY T m, double fp);
  /*
    Writes the matrix "m" to "wr". The "fp" argument should be
    some numerical ``fingerprint'' of the matrix source data. */
    
PROCEDURE Read(FILE *rd, double fp): REF T;
  /*
    Reads a matrix from "rd".  If "fp" is not zero, checks whether
    the fingerprint stored in "rd" matches "fp". */
;
} /* SparseSymetricMatrix */.

#endif
