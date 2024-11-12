INTERFACE GaussElim;

(* Gauss elimination on non-singular linear systems *)
 
PROCEDURE Solve(
    VAR A: ARRAY OF ARRAY OF LONGREAL; 
    VAR x, b: ARRAY OF LONGREAL;
    n: CARDINAL;
  );
  (*
    Solves a general linear system A x = b in x, by Gaussian elimination 
    with partial pivoting.  Considers only the first n rows/columns/elements
    of A, B, and x. *)

END GaussElim.
