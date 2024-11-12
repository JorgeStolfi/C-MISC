#ifndef logsys_sort_H0_H1_H
#define logsys_sort_H0_H1_H

/* Last edited on 2012-12-20 13:49:03 by stolfilocal */
/* General system solving. */

#define _GNU_SOURCE
#include <stdint.h>

#include <bool.h>

#include <logsys.h>
#include <logsys_def.h>

void logsys_sort_variables_for_solver_H0_H1
  ( int heur, 
    logsys_t *S, 
    int nv, 
    logsys_va_t *va[], 
    int vix[], 
    bool_t guess[]
  );
    /* These variable sorting heuristics consider each equation {eq},
      and compute a local score {sc(eq,j)} for each argument slot {j}
      that is not fixed as specified in {lo,hi}. These local scores are
      accumulated for each variable {va[i]} to get its score {score[i]}.
      The permutation {vix} will list the variables in order of
      decreasing score.
      
      Heuristic H0: if {eq} (with arguments fixed as in {lo,hi})
      determines argument {j} functionally, the local score {sc(eq,j)}
      is {(n-f)/n}, otherwise it is {f/n}; where {n} is the number of
      arguments, and {f} the number of functionally determined
      arguments.
      
      Heuristic H1: the local score {sc(eq,j)} is the reduction in the
      entropy of the set of allowed configurations when argument {j} is
      fixed, either to 0 or to 1.
      
      These heuristics do not set {guess[0..nv-1]}. */

#endif
