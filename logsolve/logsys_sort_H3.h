#ifndef logsys_sort_H3_H
#define logsys_sort_H3_H

/* Last edited on 2012-12-20 13:50:22 by stolfilocal */
/* General system solving. */

#define _GNU_SOURCE
#include <stdint.h>

#include <bool.h>

#include <logsys.h>
#include <logsys_def.h>

void logsys_sort_variables_for_solver_H3
  ( int heur, 
    logsys_t *S, 
    int nv, 
    logsys_va_t *va[], 
    int vix[], 
    bool_t guess[]
  );
    /* Heuristic H3 ({heur==3}) computes an estimate {pr0[iv]} of the
      global variable {va[iv]} being zero in a random *valid* assignment
      to {va[0..nv-1]} (an assignment that satisfies all equations).
      The heuristics then sorts the
      variables by increasing entropy of the distribution
      {(pr0[iv],1-pr0[iv])}.  
      
      Each probability {pr0[iv]} is estimated iteratively by looking at
      each equation {eq} that uses the variable {iv}, estimating the
      probability of a valid solution falling into each case of those
      equations (assuming independent a priori probabilities), and then
      combining the probabilities of the cases where variable {iv} is
      zero.
      
      These probabilities can become arbitrarily small so a small bias
      is added to compensate roundoff errors.  For this reason the 
      entropy may never become zero.
      
      The initial guess {guess[iv]} is set to 0 if {pr0[iv] > 1/2},
      set to 1 if {pr0[iv] < 1/2}, and left unchanged if {pr0[iv] == 1/2}. */


#endif
