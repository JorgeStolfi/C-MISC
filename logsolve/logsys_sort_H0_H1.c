/* See {logsys_sort_H0_H1.h}. */
/* Last edited on 2012-12-20 13:47:06 by stolfilocal */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <values.h>

#include <bool.h>
#include <affirm.h>
#include <jsmath.h>

#include <logsys.h>
#include <logsys_def.h>
#include <logsys_sort_H0_H1.h>

void logsys_sort_variables_for_solver_H0_H1
  ( int heur, 
    logsys_t *S, 
    int nv, 
    logsys_va_t *va[], 
    int vix[], 
    bool_t guess[]
  )
  {
    assert((heur == 0) || (heur == 1));
    fprintf(stderr, "sorting %d variables by heuristic H%d\n", nv, heur);

    auto void score_eq_arg_H0(logsys_op_t op, int n, double sc[]);
    auto void score_eq_arg_H1(logsys_op_t op, int n, double sc[]);
      /* Compute local score {sc[j]} for argument {j} of {n}-ary operation {op}
        for each {j} in {0..n-1}. */
    
    int i;
    double score[nv];
    for (i = 0; i < nv; i++) { score[i] = 0; }
    if (S->eq != NULL)
      { logsys_eq_t *eq = S->eq;
        do 
          { int n = eq->n;
            assert((n >= 0) && (n <= logsys_n_MAX));
            /* Compute local scores {sc[0..n-1]} of {eq}'s arguments: */
            double sc[n];
            if (heur == 0)
              { /* Compute {sc[0..n-1]} based on functional character of {eq} on args: */
                score_eq_arg_H0(eq->op, n, sc);
              }
            else if (heur == 1)
              { /* Compute {sc[0..n-1]} based on information gain for fixing each arg: */
                score_eq_arg_H1(eq->op, n, sc);
              }
            else
              { assert(FALSE); }
            /* Accumulate local scores {sc[0..n-1]} on global scores {score[0..nv-1]}: */
            int j;
            for (j = 0; j < n; j++)
              { logsys_va_t *vj = eq->arg[j].va;
                assert(vj != NULL);
                assert(vj->sys == S);
                i = (int)(vj->id);
                assert((i >= 0) && (i < nv));
                score[i] += sc[j];
              }
            eq = eq->nxeq;
          }
        while (eq != S->eq);
      }
    
    /* Index-sort variables in {vix} by score decreasing: */
    logsys_sort_variables_by_score(nv, va, score, -1, vix);

    void score_eq_arg_H0(logsys_op_t op, int n, double sc[])
      {
        /* Compute {sc[0..n-1]} based on functional character of {eq} on args: */
        /* Count in {nfun} the args that are functionally defined by {eq}: */
        /* Count in {nfre} the args that are not functionally defined by {eq}: */
        bool_t isfun[n];
        int nfun = 0, nfre = 0;
        int j;
        for (j = 0; j < n; j++) 
          { isfun[j] = logsys_op_is_functional(op, n, j);
            if (isfun[j]) { nfun++; } else { nfre++; }
          }
        /* Compute the total scores to add for functional and non-functional args: */
        double score_fun = -((double)nfre)/((double)nfre+nfun+1.0e-100);
        double score_nuf = +((double)nfun)/((double)nfre+nfun+1.0e-100);
        for (j = 0; j < n; j++) 
          { sc[j] = (isfun[j] ? score_fun : score_nuf); }
      }
      
    void score_eq_arg_H1(logsys_op_t op, int n, double sc[])
      {
        /* Compute {sc[0..n-1]} based on information gain for fixing each arg: */
        int j;
        for (j = 0; j < n; j++) 
          { sc[j] = logsys_op_entropy_reduction(op, n, j); }
      }
  }
  
