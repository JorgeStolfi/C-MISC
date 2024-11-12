#ifndef logsys_sort_H2_H
#define logsys_sort_H2_H

/* Last edited on 2012-12-20 13:49:50 by stolfilocal */
/* General system solving. */

#define _GNU_SOURCE
#include <stdint.h>

#include <bool.h>

#include <logsys.h>
#include <logsys_def.h>

void logsys_sort_variables_for_solver_H2
  ( int heur, 
    logsys_t *S, 
    int nv, 
    logsys_va_t *va[], 
    int vix[], 
    bool_t guess[]
  );
    /* Heuristic H2 ({heur==2}) for sorting the global variables of {S}
      computes an upper bound to the fraction {fr[iv,b]}
      of all assignments to {na[0..nv-1]} that satisfy all equations and
      have the global variable {va[iv]} set to {b} (0 or 1). 
      
      The heuristic then sorts the variables by increasing entropy of the
      distribution {(fr[iv,0],fr[iv,1])}.  For each variable {va[iv]},
      the first guess {guess[iv]} is set to {b} if {fr[iv,b] > fr[iv,1-b]}.
      If {fr[iv,b] == fr[iv,1-b]}, leaves {guess[iv]} unchanged.
      
      The fractions {fr} are estimated iteratively by looking at each
      equation {eq} and computing upper bound to the fraction of all
      solutions that may fit each case accepted by {eq}, then adding
      those fractions for each argument variable {iv} and value {b} to
      get an upper bound on {fr[iv,b]}. 
      
      Note that {fr[iv,0]+fr[iv,1]} is an upper bound to the fraction of
      global assignments that are valid solutions; therefore it may be
      considerably less than 1.
      
      On the other hand, the formulas are such that the estimated
      fractions are always integer multiples of {emin=1/2^nmax} where
      {nmax} is the maximum arity among the equations in the system.
      While this constraint makes it unnecessary to worry about rounding
      errors, it limits the effectiveness of this heuristic, since there
      will be only a small number of distinct entropy values. */

#endif
