/* See {logsys_sort_H2.h}. */
/* Last edited on 2012-12-20 13:47:35 by stolfilocal */

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
#include <logsys_sort_H2.h>


void logsys_sort_variables_for_solver_H2
  ( int heur, 
    logsys_t *S, 
    int nv, 
    logsys_va_t *va[], 
    int vix[], 
    bool_t guess[]
  )
  {
    bool_t debug = TRUE;
    assert(heur == 2);
    fprintf(stderr, "sorting %d variables by heuristic H%d\n", nv, heur);
    
    double fr[2*nv];
    
    auto double consider_equation(logsys_eq_t *eq);
      /* Uses equation {eq} to tighten the bounds {fr[iv,b]}. Returns the maximum 
         absolute change made to any entry of {fr}. */
    
    /* The fraction of all assignments that are solutions and
      have variable {iv} set to {b} are at most 0.5: */
    { int iv;
      for (iv = 0; iv < nv; iv++) 
        { fr[2*iv + 0] = 0.5;
          fr[2*iv + 1] = 0.5;
        }
    }
    if (S->eq != NULL)
      { /* Now cycle over the equations improving that estimate: */
        int max_passes = 20;
        int pass = 0;
        double change = +INF;
        while ((pass < max_passes) && (change > 0.0))
          { logsys_eq_t *fq = S->eq;
            change = 0.0;
            do 
              { double cheq = consider_equation(fq);
                change = fmax(change, cheq);
                fq = fq->nxeq;
              }
            while (fq != S->eq);
            fprintf(stderr, "pass %d - max fraction change %17.15f\n", pass, change);
            pass++;
          }
      }
      
    /* Compute the entropy of each variable and pick the first {guess}: */
    double entropy[nv]; 
    { int iv;
      for (iv = 0; iv < nv; iv++) 
        { double fr0 = fr[2*iv + 0];
          double fr1 = fr[2*iv + 1];
          entropy[iv] = logsys_entropy(fr0, fr1);
          if (guess != NULL)
            { if (fr0 > fr1)
                { guess[iv] = FALSE; }
              else if (fr1 > fr0)
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
                double fr0 = fr[2*iv + 0];
                double fr1 = fr[2*iv + 1];
                double ent = entropy[iv];
                fprintf(stderr, "  Fr(v=0) = %17.15f", fr0); 
                fprintf(stderr, "  Fr(v=1) = %17.15f", fr1); 
                fprintf(stderr, "  entropy = %17.15f\n", ent);
              }
          }
      }
    
    return;
    
    double consider_equation(logsys_eq_t *eq)
      {
        bool_t debug = (eq->n == 6); 
        logsys_print_equation_id(stderr, "    considering equation e", eq->id, 4, "\n");
        
        int n = eq->n;
        assert((n >= 0) && (n <= logsys_n_MAX));
        /* Get max fraction {tfr[2*j+b]} of global assgns that are valid and have arg {j} of {eq} set to {b}: */
        uint32_t N = (1 << n);
        double tfr[2*n];
        int j;
        for (j = 0; j < n; j++) { tfr[2*j + 0] = tfr[2*j + 1] = 0; }
        uint32_t X;
        for (X = 0; X < N; X++)
          { if ((eq->op & (ONE64 << X)) != 0)
              { /* Local assignment {X} is compatible with {eq}: */
                /* Compute {efrX}, max fraction of global assignments that are valid and fall in case {X} of {eq}. */
                /* At most {1/N} of global assignments will give local assignment {X}: */
                double efrX = 1.0/N;
                /* Note that {efrX} is an integer multiple of {1/2^logsys_n_MAX}. */
                /* But the fraction cannot be larger than the bound implied by {fr}: */
                for (j = 0; j < n; j++)
                  { /* Get the index {iv} of the variable that is arg {j} of {eq}: */
                    logsys_va_t *ua = eq->arg[j].va;
                    logsys_va_id_t iu = ua->id;
                    /* Get the value {b} assigned to that variable in case {X}: */
                    int b = (X >> j) & 1;
                    /* Get the max fraction {ufr} of global assgns that are valid and have {ua} set to {b}: */
                    double ufr = fr[2*iu + b]; 
                    /* That is an upper bound to {mfr[X]}: */
                    if (ufr < efrX) { efrX = ufr; }
                  }
                if (debug) { fprintf(stderr, "    case %0*X - max frac %17.15f\n", (int)imax(1,(n+3)/4), X, efrX); }
                /* Now add {efrX} to the relevant facet counts {tfr}: */
                for (j = 0; j < n; j++)
                  { /* Get the value {b} assigned to argument {j} in case {X}: */
                    int b = (X >> j) & 1;
                    tfr[2*j + b] += efrX;
                  }
              }
          }
        /* Now update {fr[iv,b]} for each variable used in {eq}: */
        double ch = 0; /* max absolute change: */
        for (j = 0; j < n; j++)
          { /* Get the index {iv} of the variable that is arg {j} of {eq}: */
            logsys_va_t *ua = eq->arg[j].va;
            logsys_va_id_t iu = ua->id;
            if (debug) 
              { fprintf(stderr, "    argument x%d = ", j); 
                logsys_print_variable_id(stderr, "v", iu, 4, " max fracs");
              }
            int b;
            for (b = 0; b < 2; b++)
              { double tjb = tfr[2*j + b];
                assert (tjb <= 0.5);
                double fub = fr[2*iu + b];
                if (debug) { fprintf(stderr, "  %17.15f", fub); }
                if (fub > tjb) 
                  { if (debug) { fprintf(stderr, " -> %17.15f", tjb); }
                    ch = fmax(ch, fub - tjb);
                    fr[2*iu + b] = tjb;
                  }
              }
            if (debug) { fprintf(stderr, "\n"); }
          }
        return ch;
      }
  }
  
