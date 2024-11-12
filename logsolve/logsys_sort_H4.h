#ifndef logsys_sort_H4_H
#define logsys_sort_H4_H

/* Last edited on 2012-12-20 13:50:37 by stolfilocal */
/* General system solving. */

#define _GNU_SOURCE
#include <stdint.h>

#include <bool.h>

#include <logsys.h>
#include <logsys_def.h>

void logsys_sort_variables_for_solver_H4
  ( int heur, 
    logsys_t *S, 
    int nv, 
    logsys_va_t *va[], 
    int vix[], 
    bool_t guess[]
  );
    /* Heuristic H4 ({heur==4}) tries to sort the variables 
      in a breadth-first order.
      
      The initial guess {guess[iv]} is not changed. */


#endif
