
#ifndef GaussElim_H
#define GaussElim_H



/* Gauss elimination on non-singular linear systems */
 
PROCEDURE Solve(
    ARRAY_vec AOF double; 
    VAR x, b: ARRAY OF double;
    nat n;
  );
  /*
    Solves a general linear system A x == b in x, by Gaussian elimination 
    with partial pivoting.  Considers only the first n rows/columns/elements
    of A, B, and x. */
;
} /* GaussElim */.

#endif
