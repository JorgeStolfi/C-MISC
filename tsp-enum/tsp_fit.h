/* Model fitting for TSP/Cycle-cover costs. */
/* Last edited on 2023-03-31 04:18:58 by stolfi */

#ifndef tsp_fit_H
#define tsp_fit_H

/* Must define _GNU_SOURCE in order to get {asprintf} */
#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <bool.h>

#include <tsp_lib.h>

/* 
   ASYMPTOTIC FITTING NEAR OPTIMUM TOURS
   
    The following procedures take the costs of {np} random tours, 
    look at the {ns} smallest ones {Z[0..ns-1]} in increasing order,
    and try to fit power models of the form

      {Z[k] ~ Zref + Coef*T[k]**(1/Expt)}   (1)

    where {T(k)} is the relative rank {(k+0.5)/np}, and
    {Expt,Coef,Zref} are model parameters to be determined so that the
    approximation is good near the smallest end of the list. */

void compute_asymp_model_costs(int32_t ns, double *ZM, int32_t np, double Expt, double Zref, double Coef);
/* Computes costs {Z[0..np-1]} of {np} permutations, according to the 
  asymptotic power model with parameters {Expt,Zref,Coef}. */
    
double compute_asymp_cost_diff_sqr(int32_t ns, double *Z, double *ZM);
/* Computes the total weighted square discrepancy between the the
  actual costs {Z[k]} and the costs predicted by the power model
  {ZM[k]}. Assumes that {Z} is increasing. */
   
void fit_asymp_model
  ( int32_t ns,         /* Number of perms in sub-sample. */   
    double *Z,      /* Their costs, in increasing order. */
    int32_t np,         /* Total sample size. */
    int32_t nv,         /* Number of vertices */
    double useExpt, /* Value of {Expt} to use; if zero, adjusts {Expt} too. */
    double *Expt,   /* (OUT) Exponent of model. */
    double *Zref,   /* (OUT) Super-optimal cost. */
    double *Coef    /* (OUT) Scale coefficient. */
  );
/* Fits the power model (1) to the tour values {Z[0..ns-1]}.
  If {useExpt} is zero, tries several integer and half-integer
  values, and selects the best one. */

/* 
  GLOBAL FITTING
  
    The procedures below try to find parameters {Sdim,Zmid,Zrad} that
    yield the best approximation
    
      {T[k] ~ F(SDim,a[k])} (2)
      
    where {a[k] = (Z[k] - Zmid)/Zrad}, and F(Sdim,x) is the fraction 
    of the volume of {d}-dimensional sphere of unit radius 
    that is contained in the slice between {x=-1} and {x=u}, 
    that is, {F(Sdim,u) = S(Sdim,u)/S(Sdim,+1)}, and 

      {S(Sdim,x) = integral(sqrt(1-u**2)**(Sdim-1), u \in [-1 _ x])}
  
  */

void compute_global_model_ranks(int32_t np, double *Z, double *TM, int32_t Sdim, double Zmid, double Zrad);
/* Computes the ranks {T[0..np-1]} of {np} permutations from their
  costs {Z[0..np-1]}, according to the spherical slice model with
  parameters {Sdim,Zmid,Zrad}.  Fails if any actual {Z} value lies 
  outside the range {[Zmid-Zrad _ Zmid+Zrad]}. */
    
double compute_global_rank_diff_sqr(int32_t np, int32_t *ix, double *TM);
/* Computes the total square discrepancy between the actual ranks {(j+0.5)/np}
  and the ranks {T[ix[j]]} predicted by the global sphere-slice model. The ranks 
  are mapped by {t -> log(t) - log(1-t)} before comparison. */

void fit_global_model
  ( int32_t np,         /* Number of perms in sample. */
    int32_t *ix,        /* Indices of perms in increasing order of cost. */
    double *Z,      /* Their costs, in increasing order. */
    int32_t nv,         /* Number of vertices */
    bool_t tour,    /* TRUE assumes tours, FALSE assumes spins. */
    int32_t useSdim,    /* Value of {Sdim} to use; if 0, adjusts {Sdim} too. */
    int32_t *Sdim,      /* (OUT) Dimension of model ball. */
    double *Zmid,   /* (OUT) Mean cost. */ 
    double *Zrad    /* (OUT) Half-span of costs. */
  );
/* Fits the global sphere-slice model (2) to the tour values
  {Z[0..np-1]}. Assumes that {Z[ix[j]]} is increasing with {j}. 
  If {useSdim} is zero, tries several integer dimensions
  and selects the best one. */


#endif
