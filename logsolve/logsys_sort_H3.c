/* See {logsys_sort_H3.h}. */
/* Last edited on 2012-12-20 13:47:58 by stolfilocal */

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
#include <logsys_sort_H3.h>

double logsys_cook_probability(double p);
  /* Return {p} clipped away from 0 and 1. */
  
double logsys_cook_probability(double p)
  {
    double tiny = 1.0e-14;
    double huge = 1 - tiny;
    assert(huge < 1);
    if (p < tiny) 
      { return tiny; }
    else if (p > huge)
      { return huge; }
    else 
      { return p; }
  }

void logsys_sort_variables_for_solver_H3
  ( int heur, 
    logsys_t *S, 
    int nv, 
    logsys_va_t *va[], 
    int vix[], 
    bool_t guess[]
  )
  {
    bool_t debug = FALSE;
    assert(heur == 3);
    fprintf(stderr, "sorting %d variables by heuristic H%d\n", nv, heur);
    
    /* {pr0[iv]} is exactly 0.0 iff there is absolutely no solution with {va[iv]==0}. */
    /* {pr0[iv]} is exactly 1.0 iff there is absolutely no solution with {va[iv]==1}. */
    double pr0[nv];
    
    auto double recompute_variable_pr0(int iv);
      /* Looks at all equations that include {va[iv]} and
        recomputes the probability {pr0[iv]} of {va[iv]} being zero.
        Returns the RELATIVE change in {pr0[iv]}. */
    
    /* Initially, set all probabilities to 0.5: */
    int iv;
    for (iv = 0; iv < nv; iv++) { pr0[iv] = 0.5; }
    
    /* Now cycle over the variables improving that estimate: */
    int max_passes = 4;
    double tol = 0.02; /* If *relative* change in all probabilities is less than this, give up. */
    int pass = 0;
    double change = +INF;
    while ((pass < max_passes) && (change > tol))
      { change = 0.0;
        for (iv = 0; iv < nv; iv++)
          { double chva = recompute_variable_pr0(iv);
            change = fmax(change, chva);
          }
        fprintf(stderr, "pass %d - max relative prob change %17.15f\n", pass, change);
        pass++;
      }
      
    /* Compute the entropy of each variable, and pick the initial guess: */
    double entropy[nv]; 
    { int iv;
      for (iv = 0; iv < nv; iv++) 
        { double p0 = pr0[iv];
          double p1 = 1 - p0;
          entropy[iv] = logsys_entropy(p0, p1);
          if (guess != NULL)
            { if (p0 > p1)
                { guess[iv] = FALSE; }
              else if (p1 > p0)
                { guess[iv] = TRUE; }
              else
                { /* Leave the default guess. */ } 
            }
        }
    }
    
    /* Index-sort variables in {vix} by entropy increasing: */
    logsys_sort_variables_by_score(nv, va, entropy, +1, vix);
    
    if (debug)
      { /* Show some entries: */
        int k;
        int max_vars = 100; /* Max vars to show. */
        for (k = 0; k < imin(nv,max_vars); k++) 
          { int iv = vix[k];
            logsys_print_variable_id(stderr, "  variable v", iv, 4, NULL);
            if (va[iv] == NULL)
              { fprintf(stderr, " is null\n"); }
            else
              { assert(va[iv]->id == iv);
                double p0 = pr0[iv];
                double p1 = 1 - p0;
                double ent = entropy[iv];
                fprintf(stderr, "  Pr(v=0) = %17.15f", p0); 
                fprintf(stderr, "  Pr(v=1) = %17.15f", p1); 
                fprintf(stderr, "  entropy = %17.15f\n", ent); 
              }
          }
      }
    
    return;
    
    auto double compute_arg_pr0(logsys_eq_t *eq, int j);
      /* Returns the probability of argument {j} of {eq} being zero 
        given that the equation is satisfied, and all variables
        other than arg {j} have independent value distributions defined
        by {pr0[0..nv-1]}.  Uses Bayes's formula with a priori 
        probability {1/2} for arg {j} being zero.
        
        Returns exactly 0 (resp exactly 1) iff there is no assignment
        {X} to the arguments of {eq} where argument {j} has value 0
        (resp.1) and and the values assigned to each arg (excluding {j})
        have nonzero probabilities. */
    
    double recompute_variable_pr0(int iv)
      { 
        /* Enumerate all equations that use {va[iv]}. For every such equation 
          {eq}, assume that it is satisfied, and estimate the probability
          that variable {iv} is 0 assuming that all arguments other
          than {va[iv]} are independently 0 with probabilites 
          given in {pr0}.  Then combine the estimates of all equations
          into a single estimate. */
          
        bool_t debug = FALSE; 
        logsys_va_t *ua = va[iv];
        if (ua == NULL) { /* Nothing to do: */ return 0.0; }
        if (debug) { logsys_print_variable_id(stderr, "  recomputing Pr(v=0) for variable v", ua->id, 4, "\n"); }
        if ((pr0[iv] == 0) || (pr0[iv] == 1)) { /* Variable is totally fixed: */ return 0.0; }
        logsys_eq_t *eq = ua->use;
        double log2_p0 = 0; /* Log2 of product of {Pr(va[iv]=0|OK(eq))} over all {eq}. */
        double log2_p1 = 0; /* Log2 of product of {Pr(va[iv]=1|OK(eq))} over all {eq}. */
        do 
          { /* Find the index {j} of variable {ua} in the args of {eq}: */
            int j = logsys_find_arg(eq, ua);
            assert((j >= 0) && (j < eq->n));
            /* Compute the probability {ep0} of argument {j} being 0, assuming {eq} is satisfied: */
            double ep0 = compute_arg_pr0(eq, j);
            double ep1 = 1 - ep0;
            assert((ep1 < 1) || (ep0 == 0)); /* Roundoff should not cause {ep1==1}. */
            /* Accumulate {ep0} : */
            if (ep0 == 0) { log2_p0 = -INF; } else { log2_p0 += log2(ep0); }
            if (ep1 == 0) { log2_p1 = -INF; } else { log2_p1 += log2(ep1); }
            /* Move to next equation on this same variable: */
            eq = eq->arg[j].nxeq;
          }
        while (eq != ua->use);
        
        /* Now invert the logs and normalize to get the probabilities {p0,p1}: */
        assert(log2_p0 <= 0);
        assert(log2_p1 <= 0);
        demand((log2_p0 != -INF) || (log2_p1 != -INF), "system is not solvable");
        double p0, p1;
        if (log2_p0 == -INF)
          { /* Variable cannot be 0: */
            p0 = 0.0; p1 = 1.0;
          }
        else if (log2_p1 == -INF)
          { /* Variable cannot be 1: */
            p0 = 1.0; p1 = 0.0;
          }
        else
          { /* Convert to relative probs, taking care to avoid underflow in both: */
            double log2_pmax = fmax(log2_p0, log2_p1);
            p0 = exp2(log2_p0 - log2_pmax);
            p1 = exp2(log2_p1 - log2_pmax);
            /* Now at most one of them underflowed. */
            double ptot = p0 + p1;
            p0 = logsys_cook_probability(p0/ptot);
            p1 = 1 - p0;
            assert((p0 > 0) && (p0 < 1));
            assert((p1 > 0) && (p1 < 1));
          }
        /* Update {pr0[iv]} and compute the relative change {relch}: */
        double p0_old = pr0[iv];
        pr0[iv] = p0;
        double relch = fabs(p0_old - p0)/hypot(p0_old, p0)*M_SQRT2;
        if (debug) 
          { logsys_print_variable_id(stderr, "  changing Pr(v=0) for v", ua->id, 4, NULL);
            fprintf(stderr, " from %17.15f to %17.15f (change %17.15f)\n", p0_old, p0, relch);
          }
        return relch;
      }

    double compute_arg_pr0(logsys_eq_t *eq, int j)
      {
        bool_t debug = FALSE; 
        if (debug) { logsys_print_equation_id(stderr, "    considering equation e", eq->id, 4, NULL); }
        
        int n = eq->n;
        assert((n >= 0) && (n <= logsys_n_MAX));
        /* Get max fraction {tfr[2*j+b]} of global assgns that are valid and have arg {j} of {eq} set to {b}: */
        uint32_t N = (1 << n);
        uint32_t X;
        double ep0 = 0; /* Prob. of {eq} satisfied when arg {j} is 0. */
        double ep1 = 0; /* Prob. of {eq} satisfied when arg {j} is 1. */
        for (X = 0; X < N; X++)
          { if ((eq->op & (ONE64 << X)) != 0)
              { /* Local assignment {X} is compatible with {eq}. */
                /* Compute prob {epX} of args being {X} when arg {j} is fixed as in {X}. */
                double log2_epX = 0; /* Log of {epX}. */
                int k;
                for (k = 0; k < n; k++)
                  { if (k != j)
                      { /* Get the index {it} of the variable that is arg {i} of {eq}: */
                        logsys_va_t *ta = eq->arg[j].va;
                        logsys_va_id_t it = ta->id;
                        /* Get the value {b} assigned to that variable in case {X}: */
                        int bt = (X >> k) & 1;
                        /* Get the current probability {pt} of that arg being {bt}: */
                        double log2_ptb;
                        double pt0 = pr0[it];
                        double pt1 = 1 - pt0; 
                        assert((pt1 < 1) || (pt0 == 0)); /* Roundoff should not cause {pt1==0}. */
                        if (bt == 0)
                          { log2_ptb = (pt0 == 0 ? -INF : (pt0 == 1 ? 0 : log2(pt0))); }
                        else if (bt == 1)
                          { log2_ptb = (pt1 == 0 ? -INF : (pt1 == 1 ? 0 : log2(pt1))); }
                        else
                          { assert(FALSE); }
                        log2_epX += log2_ptb;
                      }
                  }
                /* Convert {log2_epX} to actual prob in order to add: */
                double epX = (log2_epX == -INF ? 0 : exp2(log2_epX));
                /* Add to the appropriate bin: */
                int bj = (X >> j) & 1;
                if (bj == 0) { ep0 += epX; } else { ep1 += epX; }
              }
          }
        /* Now {ep0} is zero iff it is impossible to satisfy {eq} */
        /* with arg {j} set to zero; and similarly for {ep1}. */
        /* Normalize {ep0,ep1} to add 1: */
        double ep = ep0 + ep1;
        demand(ep > 0, "equation cannot be satisfied at all");
        ep0 /= ep;
        ep1 /= ep;
        if (debug) 
          { fprintf(stderr, " Pr(x%d=0) = %17.15f", j, ep0);
            fprintf(stderr, " Pr(x%d=1) = %17.15f", j, ep1);
            fprintf(stderr, "\n");
          }
        return ep0;
      }
  }
