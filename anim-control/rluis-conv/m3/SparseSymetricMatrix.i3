INTERFACE SparseSymetricMatrix;

(*
   Defines a data structure for a sparse symetric square matrix.
   The structure is a vector of vectors, one for each row.
   Each entry of a secondary vector contains the column index (>= 0)
   and value of a non-zero element of the corresponding matrix row.
   The master vector also contains, for each row
   the index in the secondary vector containing the diagonal element
   of that row.
 *)
 
IMPORT Wr, Rd;

TYPE 
  Vector = ARRAY OF LONGREAL;
  MatrixElement = RECORD j: CARDINAL; a: LONGREAL END;
  Row = ARRAY OF MatrixElement;
  T = ARRAY OF RECORD d: CARDINAL; row: REF Row END;
 
PROCEDURE New(n: CARDINAL): T;
  (*
    Allocates a new "n" by "n" matrix, initialized to the identity matrix. *)

PROCEDURE GetValue(READONLY m: T; i, j: CARDINAL): LONGREAL;
  (*
    Gets element "[i,j]" of the matrix (by linear scanning of row i). *)
    
PROCEDURE PutRow(VAR m: T; i: CARDINAL; READONLY row: Vector);
  (*
    Replaces row "i" of "m" by a (compressed) copy of "row". *)

PROCEDURE Remove(VAR m: T; READONLY sel: ARRAY OF BOOLEAN);
  (*
    Sets to zero the rows and columns of "m" with indices "k"
    such that "sel[k] = TRUE"; except for the diagonal elements,
    which are set to 1. *)
    
PROCEDURE SubMul(READONLY m: T; READONLY v: Vector; VAR x: Vector); 
  (* 
    Computes x := x - m * v *)
    
PROCEDURE Factor(READONLY m: T): REF T;
  (*
    Returns the Cholesky factorization of "m = LDL^T" 
    Valid only if "m" is positive-definite. *)
    
PROCEDURE Solve(READONLY mf: T; READONLY v: Vector; VAR x: Vector; i0 := 0);
  (*
    Solves "m*x = v" for "x", given the Choleski factorization "mf" of "m".
    If "i0" > 0, ... ? *)

PROCEDURE Print(wr: Wr.T; READONLY m: T);
  (*
    Prints "m" to "wr". *)
    
PROCEDURE Write(wr: Wr.T; READONLY m: T; fp: LONGREAL);
  (*
    Writes the matrix "m" to "wr". The "fp" argument should be
    some numerical ``fingerprint'' of the matrix source data. *)
    
PROCEDURE Read(rd: Rd.T; fp: LONGREAL): REF T;
  (*
    Reads a matrix from "rd".  If "fp" is not zero, checks whether
    the fingerprint stored in "rd" matches "fp". *)

END SparseSymetricMatrix.
