/* See {logsys_sort_H4.h}. */
/* Last edited on 2012-12-20 17:51:11 by stolfilocal */

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
#include <logsys_sort_H4.h>

void logsys_sort_variables_for_solver_H4
  ( int heur, 
    logsys_t *S, 
    int nv, 
    logsys_va_t *va[], 
    int vix[], 
    bool_t guess[]
  )
  {
    bool_t debug = FALSE;
    assert(heur == 4);
    fprintf(stderr, "sorting %d variables by heuristic H%d\n", nv, heur);
    
    if (S->va == NULL) { /* No variables to sort: */ return; }
    
    /* Fill {vix} with the identity perm: */
    int iv;
    for (iv = 0; iv < nv; iv++) { assert(va[iv]->id == iv); vix[iv] = iv; }
    /* Move all null variables to the end and set {nvok} to their count: */
    int nvok = nv;
    iv = 0;
    while (iv < nvok) 
      { if (va[vix[iv]] == NULL)
          { int tmp = vix[iv]; vix[iv] = vix[nvok-1]; vix[nvok-1] = tmp; nvok--; }
        else
          { iv++; }
      }
    if (nvok == 0) { return; }
    
    /* Gather a list {eq} of all equations, analogous to {va}: */
    logsys_eq_t **eq;
    int ne;
    logsys_get_equations(S, &ne, &eq);
    /* Create an equation index list {eix}, like {vix}, fill it with identity perm: */
    int *eix = notnull(malloc(ne*sizeof(logsys_eq_t *)), "no mem");
    int ie;
    for (ie = 0; ie < ne; ie++) { assert(eq[ie]->id == ie); eix[ie] = ie; }
    /* Move all null equations to the end and set {neok} to their count: */
    int neok = nv;
    ie = 0;
    while (ie < neok) 
      { if (eq[eix[ie]] == NULL)
          { int tmp = eix[ie]; eix[ie] = eix[neok-1]; eix[neok-1] = tmp; neok--; }
        else
          { ie++; }
      }
    if (neok == 0) { return; }
    
    /* Now run a breadth-first Dijkstra-like distance sort on vars and eqs: */
    /* Distances of variables from root node. Indexed by variable ID: */
    double *va_dist = notnull(malloc(nvok*sizeof(double)), "no mem");
    for (iv = 0; iv < nvok; iv++) { va_dist[iv] = +INF; }
    /* Distances of equations from root node. Indexed by equation ID. */
    double *eq_dist = notnull(malloc(neok*sizeof(double)), "no mem");
    for (ie = 0; ie < neok; ie++) { eq_dist[ie] = +INF; }
    
    /* How many nodes have been fixed already: */
    int nvfix = 0; /* Variables already fixed are {va[vix[0..nvfix-1]]}. */
    int nefix = 0; /* Equations already fixed are {eq[eix[0..nefix-1]]}. */

    auto bool_t fixed_by_system(logsys_va_t *ua);
      /* TRUE iff there is some nontrivial unary equation in {S} that fixes {ua}. */ 
    
    auto void scan_variable_neighbors(void);
      /* Assumes all dists {eq_dist[eix[0..nefix-1]]} and {va_dist[vix[0...nvfix-1]]}
        are correct and those are the {nefix+nvfix} nodes closest to the root var.
        Assumes also that {iv = va[vix[nvfix-1]]} is the last node (variable
        or equation) added to that set: so that {va_dist[iv]} is the maximum
        among those distances.  Updates {eq_dist[eix[nefix..neok-1]]} with all the 
        arcs out of {va[iv]}.  */

    auto void scan_equation_neighbors(void);
      /* Symmetrical to {scan_variable_neighbors}. */

    auto void get_next_node(int nfix, int nok, int ix[], double dist[]);
      /* Rearanges {ix[nfix..nok-1]} so that {dist[ix[nfix]]}
        is minimum.  */



    /* Move all variables that are used in unary operations to the beginning and fix them : */
    iv = 0;
    while (iv < nvok) 
      { if ((va[vix[iv]] != NULL) && (fixed_by_system(va[vix[iv]])))
          { int tmp = vix[iv]; vix[iv] = vix[nvfix]; vix[nvfix] = tmp;
            va_dist[vix[nvfix]] = 0.0;
            nvfix++;
            scan_variable_neighbors();
          }
        else
          { iv++; }
      }
    /* Make sure we have at least one fixed variable: */
    if (nvfix == 0)
      { /* If no vars got fixed, fix an arbitrary initial variable: */
        { assert(nvok >= 1); int k = nvok/2; int tmp = vix[0]; vix[0] = vix[k]; vix[k] = tmp; }
        va_dist[vix[0]] = 0;
        nvfix++;
        scan_variable_neighbors();
      }
    /* Make sure that the next pivots are in place: */
    get_next_node(nefix, neok, eix, eq_dist);
    get_next_node(nvfix, nvok, vix, va_dist);
    if (debug) 
      { for (iv = 0; iv < nvfix; iv++) 
          { logsys_print_variable_id(stderr, "  starting with variable v", vix[iv], 4, "\n"); }
      }
    
    while ((nvfix < nvok) || (nefix < neok))
      { /* Decide which is next to be fixed, equation ('e') or variable ('v'): */
        int which;
        double dv = (nvfix < nvok ? va_dist[vix[nvfix]] : +INF);
        double de = (nefix < neok ? eq_dist[eix[nefix]] : +INF);
        if (de < dv)
          { which = 'e'; }
        else if (dv < de)
          { which = 'v'; }
        else 
          { which = "ev"[(nvfix + nefix) % 2]; }
        /* Fix that node and update the distances on the other side accordingly: */
        if (which == 'v')
          { /* Fix {va[vix[nvfix]]}: */
            if (debug) 
              { int iv = vix[nvfix];
                logsys_print_variable_id(stderr, "  fixing v", iv, 4, NULL);
                fprintf(stderr, " dist = %17.15f\n", va_dist[iv]);
              }
            nvfix++;
            /* Scan the equations that use {va[vix[nvfix-1]]}, update their distances: */
            scan_variable_neighbors();
            /* Make sure the nearest of {eq[eix[nefix..neok-1]]} is {eq[eix[nefix]]}: */
            get_next_node(nefix, neok, eix, eq_dist);
          }
        else if (which == 'e')
          { /* Fix {eq[eix[nefix]]}: */
            if (debug) 
              { int ie = eix[nefix];
                logsys_print_equation_id(stderr, "  fixing e", ie, 4, NULL);
                fprintf(stderr, " dist = %17.15f\n", eq_dist[ie]);
              }
            nefix++;
            /* Scan the variables that are used in {eq[eix[nefix-1]]}, update their distances: */
            scan_equation_neighbors();
            /* Make sure the nearest of {va[vix[nvfix..nvok-1]]} is {va[vix[nvfix]]}: */
            get_next_node(nvfix, nvok, vix, va_dist);
          }
        else
          { assert(FALSE); }
      }

    free(eq);
    free(eix);
    free(va_dist);
    free(eq_dist);
    
    return;
    
    auto double equation_to_variable_cost(logsys_eq_t *fq, double df, logsys_va_t *ua);
      /* Defines the length of the edge from a pivot equation {fq} (whose fixed distance is {df})
        and one of its argument variables {ua}. */
      
    auto double variable_to_equation_cost(logsys_va_t *ua, double da, logsys_eq_t *fq);
      /* Defines the length of the edge from a pivot variable {ua} (whose fixed distance is {da})
        and one of the equations {fq} that uses it. */
      
    void scan_variable_neighbors(void)
      {
        assert((nvfix > 0) && (nvfix <= nvok));
        /* Get next pivot variable {ua} from {vix[nvfix-1]} and its distance {du}: */
        int ia = vix[nvfix-1];
        logsys_va_t *ua = va[ia];
        assert(ua->id == ia);
        double da = va_dist[ia];
        /* Enumerate equations that use {ua}: */
        logsys_eq_t *fq = ua->use;
        if (fq != NULL)
          { do
              { int n = fq->n; /* Arg count. */
                /* Get the length {daf} of the step from {ua} to {fq}: */
                double daf = variable_to_equation_cost(ua, da, fq);
                /* Now update the distance of {fq}: */
                int kf = fq->id;
                assert(eq[kf] == fq);
                double dfnew = da + daf;
                if (debug) 
                  { logsys_print_variable_id(stderr, "    arc from v", ia, 4, NULL);
                    logsys_print_equation_id(stderr, " to e", kf, 4, NULL);
                    fprintf(stderr, " daf = %17.15f  dfold = %17.15f dfnew = %17.15f\n", daf, eq_dist[kf], dfnew);
                  }
                if (eq_dist[kf] > dfnew) { eq_dist[kf] = dfnew; }
                /* Go to the next equation that uses {ua}: */
                int j = logsys_find_arg(fq, ua);
                assert((j >= 0) && (j < n));
                fq = fq->arg[j].nxeq;
              }
            while (fq != ua->use);
          }
      }
      
    void scan_equation_neighbors(void)
      { assert((nefix > 0) && (nefix <= nvok));
        /* Get next pivot equation {fq} from {eix[nefix-1]} and its distance {df}: */
        int kf = eix[nefix-1];
        logsys_eq_t *fq = eq[kf];
        assert(fq->id == kf);
        double df = eq_dist[kf];
        /* Enumerate variables used by {fq}: */
        int n = fq->n; /* Arg count. */
        int k;
        for (k = 0; k < n; k++)
          { /* Get an argument variable {ua} and its index {ia}: */
            logsys_va_t *ua = fq->arg[k].va;
            int ia = ua->id; assert(va[ia] == ua);
            /* Get the length {dfa} of the step from {fq} to {ua}: */
            double dfa = equation_to_variable_cost(fq, df, ua);
            /* Now update the distance of {ua}: */
            double danew = df + dfa;
            if (debug) 
              { logsys_print_equation_id(stderr, "    arc from e", kf, 4, NULL);
                logsys_print_variable_id(stderr, " to v", ia, 4, NULL);
                fprintf(stderr, " dfa = %17.15f  daold = %17.15f danew = %17.15f\n", dfa, va_dist[ia], danew);
              }
            if (va_dist[ia] > danew) { va_dist[ia] = danew; } 
          }
      }
                
    double variable_to_equation_cost(logsys_va_t *ua, double da, logsys_eq_t *fq)   
      { /* Let {m} be the count of args of {fq} that appear to be fixed. */
        /* Let {s[0..m-1]} be {(1 + va_dist[ib]/da)/2} for those args. */
        /* The edge length {daf} from {ua} to {fq} is{1/m} of the harmonic mean of {s[0..m-1]}. */
        /* This sort of falsifies Dijkstra's invariant, but still a reasonable heuristic, maybe. */
        int m = 0;
        double sum_invs = 0; /* Sum of {1/s[0..m-1]}. */ 
        int k;
        for (k = 0; k < fq->n; k++)
          { logsys_va_t *ub = fq->arg[k].va;
            int ib = ub->id; assert(va[ib] == ub);
            double darg = va_dist[ib];
            if (darg <= da) 
              { /* Assume fixed: */ 
                sum_invs += (da == 0 ? 1.0 : 2.0/(1.0 + darg/da));
                m++; 
              }
            else
              { assert(fq->arg[k].va != ua); }
          }
        /* Now compute cost {daf} of edge from {ua} to {fq}: */
        assert(m > 0);
        assert(sum_invs > 0); /* May be {+oo}. */
        double daf = 1.0/sum_invs; /* {= (1/m)*havg = (1/m)(1/(sum_invs/m))}. */
        return daf;
      }
            
    double equation_to_variable_cost(logsys_eq_t *fq, double df, logsys_va_t *ua)
      { /* Let {m} be the count of equations that use {ua} and appear to be fixed. */
        /* Let {s[0..m-1]} be {eq_dist[ib]/df} for those args. */
        /* The edge length {dfa} from {fq} to {ua} to is{1/m} of the harmonic mean of {s[0..m-1]}. */
        /* This sort of falsifies Dijkstra's invariant, but still a reasonable heuristic, maybe. */
        int m = 0;
        double sum_invs = 0; /* Sum of {1/s[0..m-1]}. */ 
        logsys_eq_t *gq = ua->use;
        assert(gq != NULL); /* Since {ua} is supposed to be an argument of {fq}. */
        do
          { int ig = gq->id;
            assert(eq[ig] == gq);
            double duse = eq_dist[ig];
            if (duse <= df)
              { /* {gq} appears to be fixed: */ 
                sum_invs += (df == 0 ? 1.0 : 2.0/(1 + duse/df));
                m++;
              }
            else
              { assert(gq != fq); }
            /* Go to next user of {ua}: */
            int j = logsys_find_arg(gq, ua);
            assert((j >= 0) && (j < gq->n));
            assert(gq->arg[j].va == ua);
            gq = gq->arg[j].nxeq;
          }
        while (gq != ua->use);
        /* Now compute cost {dfa} of edge from {fq} to {ua}: */
        assert(m > 0);
        assert(sum_invs > 0); /* May be {+oo}. */
        double dfa = 1.0/sum_invs; /* {= (1/m)*havg = (1/m)(1/(sum_invs/m))}. */
        return dfa;
      }

    bool_t fixed_by_system(logsys_va_t *ua)
      {
        /* Enumerate equations that use {ua}: */
        logsys_eq_t *fq = ua->use;
        if (fq != NULL)
          { do
              { if ((fq->n == 1) && (fq->op != logsys_op1_TRU)) { return TRUE; }
                int j = logsys_find_arg(fq, ua);
                assert((j >= 0) && (j < fq->n));
                fq = fq->arg[j].nxeq;
              }
            while (fq != ua->use);
          }
        return FALSE;
      }
          
    void get_next_node(int nfix, int nok, int ix[], double dist[])
      { 
        /* Bubble sort pass to bring the minimum to {vix[nvfix]}: */
        /* !!! Shoudl use a heap !!! */
        int i = nok - 1;
        while (i > nfix) 
          { int j = i - 1;
            double di = dist[ix[i]];
            double dj = dist[ix[j]];
            if (di < dj) { /* Swap: */ int tmp = ix[i]; ix[i] = ix[j]; ix[j] = tmp; }
            i = j;
          }
      }
  }
