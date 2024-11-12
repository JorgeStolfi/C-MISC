/* aasolve.h - solves M equations on N unknowns by affine arithmetic */
/* Last edited on 2002-12-01 18:23:39 by stolfi */

/* By Luiz Henrique de Figueiredo and Jorge Stolfi, nov/2002 */

#ifndef aasolve_H
#define aasolve_H

#include <aa.h>

int aa_solve(int N, AAform *X, int M, AAform *Y, int P, int *S);
/* Returns an affine form describing all solutions of a system of 
  M equations on N variables.
  
  The procedure receives: the affine forms X[0..N-1], describing the
  "a priori" joint range of the N variables; the affine forms
  Y[0..M-1] (with M <= N) for M functions of those variables, computed
  over the forms X[0..N-1]; and a list S[0..P-1] of distinct and non-zero
  noise symbol indices. The procedure selects K (K <= M) noise symbols from
  this list; uses the equations Y[i] == 0 to find affine forms for
  those symbols in terms of the remaining symbols; and then uses these
  forms to eliminate those noises from the forms X[r].
  
  The precedure returns the number K of noises that have been actually
  eliminated; which may be less than min(P,M), if the equations Y[i]
  == 0 are not sufficiently independent. In any case, the list
  S[0..P-1] is permuted so that S[0..K-1] are the noises that were
  actually eliminated.
  
  In the process, the forms Y[0..M-1] are modified (by linear
  combination) to produce equivalent, linearly independent equations,
  such that Y[0..K-1] are cardinal on the eliminated noises. */
  
#endif
