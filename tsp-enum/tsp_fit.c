/* See tsp_fit.h */
/* Last edited on 2024-11-23 06:20:14 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <assert.h>
#include <stdlib.h>
#include <values.h>
#include <math.h>

#include <bool.h>
#include <jsmath.h>
#include <ball.h>
#include <sign.h>
#include <affirm.h>
#include <minu_brent.h>

#include <tsp_fit.h>
#include <tsp_lib.h>

/* 
   INTERNAL PROTOTYPES */

void fit_asymp_zref_coef
  ( int32_t ns,         /* Number of perms in sub-sample. */   
    double *Z,      /* Their costs, in increasing order. */
    int32_t np,         /* Total sample size. */               
    double Expt,                 
    double *Zref,
    double *Coef
  );
  /* Adjusts the parameters {*Zref,*Coef} of the power model (1), with
    given exponent {Expt}, for best fit to the values {Z[0..ns-1]}. */

typedef void fit_asymp_proc_t
  ( int32_t ns,         /* Number of perms in sub-sample. */   
    double *Z,      /* Their costs, in increasing order. */
    int32_t np,         /* Total sample size. */               
    double minExpt,
    double maxExpt,
    double *Expt, 
    double *Zref, 
    double *Coef
  );
  /* Type of a procedure that adjusts the parameters
    {*Expt,*Zref,*Coef} of the power model (1) so as to obtain the
    best fit to the values {Z[0..ns-1]}. */
  
void fit_brute_asymp_expt_zref_coef
  ( int32_t ns,         /* Number of perms in sub-sample. */   
    double *Z,      /* Their costs, in increasing order. */
    int32_t np,         /* Total sample size. */               
    double minExpt,
    double maxExpt,
    double *Expt, 
    double *Zref, 
    double *Coef
  );
  /* Adjusts the parameters {*Expt,*Zref,*Coef} of the power model (1)
    so as to obtain the best fit to the values {Z[0..ns-1]}.
    The {Expt} parameter is adjusted in half-integer steps
    between {minExpt} and {maxExpt}. */

void fit_brent_asymp_expt_zref_coef
  ( int32_t ns, 
    double *Z,
    int32_t np, 
    double minExpt,
    double maxExpt,
    double *Expt, 
    double *Zref, 
    double *Coef
  );
  /* Adjusts the parameters {*Expt,*Zref,*Coef} of the power model (1)
    so as to obtain the best fit to the values {Z[0..ns-1]}.
    The optimal {Expt} parameter is computed by Brent's nonlinear
    minimizer. */

void fit_global_zmid_zrad
  ( int32_t np,         /* Number of perms in sub-sample. */   
    int32_t *ix,        /* Indices of perms in increasing order of cost. */
    double *Z,      /* Perm costs. */
    int32_t Sdim,                 
    double *Zmid,
    double *Zrad
  );
  /* Adjusts the parameters {*Zmid,*Zrad} of the global sphere-slice model
    (2), with given sphere dimension {Sdim}, for best fit to the 
    values {Z[0..np-1]}.  Assumes that {Z[ix[j]]} is increasing with {j}. */

void fit_global_sdim_zmid_zrad
  ( int32_t np,         /* Number of perms in sub-sample. */   
    int32_t *ix,        /* Indices of perms in increasing order of cost. */
    double *Z,      /* Perm costs. */
    int32_t minSdim,
    int32_t maxSdim,
    int32_t *Sdim, 
    double *Zmid, 
    double *Zrad
  );
  /* Adjusts the parameters {*Sdim,*Zmid,*Zrad} of the global sphere-slice model (2)
    so as to obtain the best fit to the values {Z[0..np-1]}.
    Assumes that {Z[ix[j]]} is increasing with {j}.
    The {Sdim} parameter is adjusted in integer steps
    between {minSdim} and {maxSdim}. */

#define DEBUG TRUE
/* Define this as TRUE to print various diagnostics. */

void compute_asymp_model_costs(int32_t ns, double *ZM, int32_t np, double Expt, double Zref, double Coef)
  { int32_t k;
    for (k = 0; k < ns; k++)
      { /* Relative rank in range {[0 _ 1]}: */
        double Tk = ((double)k + 0.5)/((double)np);
        /* Predicted cost: */
        ZM[k] = Zref + Coef*exp(log(Tk)/Expt);
      }
  }

double compute_asymp_cost_diff_sqr(int32_t ns, double *Z, double *ZM)
  {
    double S2 = 0;
    int32_t k;
    for (k = 0; k < ns; k++)
      { /* Actual cost: */
        double Zk = Z[k]; 
        /* Predicted cost: */
        double ZMk = ZM[k];
        /* Discrepancy: */
        double dZ = ZMk - Zk; 
        /* Sample weight: */
        double x = ((double)k + 0.5)/((double)ns);
        double W = x*(1.0 - x);
        S2 += W*dZ*dZ;
      }
    return S2;
  }

void fit_asymp_model
  ( int32_t ns,         /* Number of perms in sub-sample. */
    double *Z,      /* Their costs, in increasing order. */
    int32_t np,         /* Total sample size. */
    int32_t nv,         /* Number of vertices */
    double useExpt, /* Value of {Expt} to use (-1 = brute fit, -2 = Brent fit). */
    double *Expt,   /* (OUT) Exponent of model. */
    double *Zref,   /* (OUT) Super-optimal cost. */
    double *Coef    /* (OUT) Scale coefficient. */
  )
  {
    if (ns < 3)
      { /* Not enough data for fitting, print a guess: */
        fprintf(stderr, "not enough data for asymptotic model fitting\n");
        (*Expt) = 1.0; (*Zref) = 0.0; (*Coef) = 1.0;
      }
    else if (useExpt > 0)
      { /* Use given {Expt}, fit {Zref,Coef} to {Z[0..ns-1]}: */
        (*Expt) = useExpt;
        fit_asymp_zref_coef(ns, Z, np, (*Expt), Zref, Coef);
      }
    else
      { /* Fit {Expt,Zref,Coef} to {Z[0..ns-1]}: */
        double minExpt = 0.5;
        double maxExpt = (double)nv;
        fit_asymp_proc_t *fit_proc;
        if (useExpt == -1)
          { fit_proc = &fit_brute_asymp_expt_zref_coef; }
        else if (useExpt == -2)
          { fit_proc = &fit_brent_asymp_expt_zref_coef; }
        else
          { demand(FALSE, "invalid {useExpt}"); }
        fit_proc(ns, Z, np, minExpt, maxExpt, Expt, Zref, Coef);
      }

    if (DEBUG) 
      { 
        fprintf(stderr, "fit_asymp_model:");
        fprintf(stderr, " Expt = %6.1f", *Expt);
        fprintf(stderr, " Zref = %24.16e", *Zref);
        fprintf(stderr, " Coef = %24.16e\n", *Coef);
      }
  }

void fit_asymp_zref_coef
  ( int32_t ns, 
    double *Z, 
    int32_t np,
    double Expt,
    double *Zref,
    double *Coef
  )
  {
    if (DEBUG) { fprintf(stderr, "%s: Expt = %6.1f", __FUNCTION__, Expt); } 
    /* We replace {T[k]} by {R[k] = T[k]**(1/Expt)}, and fit a linear
      model {Z[k] ~ Coef*R[k] + Zref}. */

    /* Set up the normal system for the basis {R,1} and the 
      scalar product { <f|g> = SUM{ W(k)*f(k)g(k) : k \in 0..ns-1 }}. */
    double SRR = 0; double SR1 = 0; double S11 = 0;
    double SZR = 0; double SZ1 = 0;
    int32_t k;
    for (k = 0; k < ns; k++)
      { /* Relative rank in {[0 _ 1]}: */
        double Tk = ((double)k + 0.5)/((double) np);
        /* Compute {Rk = R[k]}: */
        double Rk = exp(log(Tk)/Expt);
        /* Actual cost: */
        double Zk = Z[k]; 
        /* Compute weight: */
        double xk = ((double)k + 0.5)/((double)ns);
        double Wk = xk*(1.0 - xk);
        /* Accumulate scalar products: */
        SRR += Wk*Rk*Rk; SR1 += Wk*Rk; S11 += Wk;
        SZR += Wk*Zk*Rk; SZ1 += Wk*Zk;
      }
    /* Solve normal system {((SRR, SR1),(SR1,S11)) * (CR,C1) = (SZR,SZ1)}: */
    double D = SRR*S11 - SR1*SR1;
    double DR = SZR*S11 - SZ1*SR1;
    double D1 = SRR*SZ1 - SR1*SZR;
    double CR = DR/D;
    double C1 = D1/D;
    (*Coef) = CR;
    (*Zref) = C1;
    /* Consistency: */
    if ((*Coef) <= 0.0) { (*Coef) = 0.001; }
    if (DEBUG)
      { fprintf(stderr, " Zref = %24.16e Coef = %24.16e\n", (*Zref), (*Coef)); }
  }

void fit_brute_asymp_expt_zref_coef
  ( int32_t ns,         /* Number of perms in sub-sample. */   
    double *Z,      /* Their costs, in increasing order. */
    int32_t np,         /* Total sample size. */               
    double minExpt,
    double maxExpt,
    double *Expt, 
    double *Zref, 
    double *Coef
  )
  {
    /* Trial exponents are {EE/2} where {EE} ranges in {EEmin .. EEmax}. */
    int32_t EEmin = (int32_t)floor(2*minExpt + 0.5);
    int32_t EEmax = (int32_t)floor(2*maxExpt + 0.5);
    /* Try all exponents, remember best-fitting model: */
    double bestEx = 0.0, bestZr = 0.0, bestCf = 0.0; 
    double bestS2 = DBL_MAX;
    int32_t EE;
    double ZM[ns];
    for (EE = EEmin; EE <= EEmax; EE++)
      { /* Compute trial exponent: */
        double Ex = ((double)EE)/2;
        /* Find best values {Zr,Cf} for parameters {Zref,Coef}, fixing {Expt=Ex}: */
        double Zr, Cf;
        fit_asymp_zref_coef(ns, Z, np, Ex, &Zr, &Cf);
        /* Check quality of fit: */
        compute_asymp_model_costs(ns, ZM, np, Ex, Zr, Cf);
        double S2 = compute_asymp_cost_diff_sqr(ns, Z, ZM);
        if (S2 < bestS2) { bestEx = Ex; bestZr = Zr; bestCf = Cf; bestS2 = S2; }
      }
    assert(bestS2 < DBL_MAX);
    (*Expt) = bestEx;
    (*Zref) = bestZr;
    (*Coef) = bestCf;
  }

void fit_brent_asymp_expt_zref_coef
  ( int32_t ns, 
    double *Z,
    int32_t np, 
    double minExpt,
    double maxExpt,
    double *Expt, 
    double *Zref, 
    double *Coef
  )
  {
    double xa = minExpt; double xb = maxExpt;
    double tol = 0.001;
    double dist = xb - xa;
    
    auto bool_t d_eval(void *parms, double x, double *fx, double *dfx);
    
    double x = (xa + xb)/2;
    double fx; 
    double dfx;
    
    d_eval(NULL, x, &fx, &dfx);

    minu_brent_minimize
      ( NULL, 
        d_eval, 
        &x, 
        &fx, 
        &dfx, 
        tol, 
        dist, 
        &xa, 
        &xb, 
        NULL, 
        FALSE
      );

    /* Set {Expt} and compute {Zref,Coef}, just to be sure: */
    (*Expt) = x;
    fit_asymp_zref_coef(ns, Z, np, (*Expt), Zref, Coef);
    
    return;

    bool_t d_eval 
      ( void *parms, 
        double x, 
        double *fx, 
        double *dfx
      )
      { double Ex = x;
        double Zr, Cf;
        fit_asymp_zref_coef(ns, Z, np, Ex, &Zr, &Cf);
        compute_asymp_model_costs(ns, ZM, np, Ex, Zr, Cf);
        (*fx) = compute_asymp_cost_diff_sqr(ns, Z, ZM);
        (*dfx) = NAN;
        return FALSE;
      }
    
  }

void compute_global_model_ranks(int32_t np, double *Z, double *TM, int32_t Sdim, double Zmid, double Zrad)
  { int32_t k;
    double eps = 0.1/np; 
    for (k = 0; k < np; k++)
      { /* Actual cost of tour: */
        double Zk = Z[k];
        /* Compute relative cost {Uk} in range {[-1 _ 1]}: */
        double Uk = (Zk - Zmid)/Zrad;
        /* Clip to valid model range: */
        if (Uk < -1) { Uk = -1; }
        if (Uk > +1) { Uk = +1; }
        /* Shrink {Uk} slightly towards center of range: */
        Uk = (1-eps)*Uk;
        /* Compute rank predicted by model: */
        TM[k] = ball_cap_vol_frac_pos(Sdim, Uk);
      }
  }

double compute_global_rank_diff_sqr(int32_t np, int32_t *ix, double *TM)
  {
    double S2 = 0;
    int32_t j; /* Actual rank. */
    for (j = 0; j < np; j++)
      { /* Index of element with rank {j}: */
        int32_t k = ix[j];
        /* Its actual relative rank: */
        double Tk = ((double)j + 0.5)/((double) np);
        /* Its model-predicted relative rank: */
        double TMk = TM[k];
        demand(TMk > 0.0, "invalid predicted rank");
        demand(TMk < 1.0, "invalid predicted rank");
        /* Discrepancy after bi-log map: */
        double tk = log(Tk) - log(1.0 - Tk);
        double tMk = log(TMk) - log(1.0 - TMk);
        double dt = tMk - tk; 
        /* Weight: */
        double W = 1.0;
        /* Accumulate: */
        S2 += W*dt*dt;
      }
    return S2;
  }

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
  )
  {
    if (np < 3)
      { /* Not enough data for fitting, print a guess: */
        fprintf(stderr, "not enough data for global model fitting\n");
        (*Sdim) = 1; (*Zmid) = 1.0; (*Zrad) = 1.0;
      }
    else if (useSdim > 0)
      { /* Use given {Sdim}, fit {Zmid,Zrad} to {Z[0..ns-1]}: */
        (*Sdim) = useSdim;
        fit_global_zmid_zrad(np, ix, Z, (*Sdim), Zmid, Zrad);
      }
    else
      { /* Fit {Zmid,Zrad,Sdim} to {Z[0..ns-1]}: */
        double minSdim = nv/2;
        double maxSdim = 3*nv/2;
        fit_global_sdim_zmid_zrad(np, ix, Z, minSdim, maxSdim, Sdim, Zmid, Zrad);
      }

    if (DEBUG) 
      {
        fprintf(stderr, "fit_global_model:");
        fprintf(stderr, " Sdim = %4d", *Sdim);
        fprintf(stderr, " Zmid = %24.16e", *Zmid);
        fprintf(stderr, " Zrad = %24.16e\n", *Zrad);
      }
  }

void fit_global_zmid_zrad
  ( int32_t np, 
    int32_t *ix,
    double *Z,
    int32_t Sdim, 
    double *Zmid, 
    double *Zrad
  )
  {
    if (DEBUG) { fprintf(stderr, "%s: Sdim = %d", __FUNCTION__, Sdim); } 
    /* We compute {U[k]} such that {ball_cap_vol_frac_pos(Sdim, U[k]) = T[k]},
      then fit a linear model {Z[k] ~ Zmid + Zrad*U[k]}. */

    /* Set up the normal system for the basis {U,1} and the 
      scalar product { <f|g> = SUM{ W(k)*f(k)g(k) : k \in 0..ns-1 }}. */
    double SUU = 0; double SU1 = 0; double S11 = 0;
    double SZU = 0; double SZ1 = 0;
    int32_t j;
    double Uk = -1.0; /* Previous value of {Uk} (see below). */
    double TMk = 0.0; /* {= ball_cap_vol_frac_pos(Sdim,Uk);} */
    double Ustep = 10.0/np;  /* Initial step in {Uk}. */
    double Teps = 1.0e-5; /* Error tolerance for {TMk}. */
    for (j = 0; j < np; j++)
      { /* Index of sample with rank {j}: */
        int32_t k = ix[j];
        /* Its actual cost: */
        double Zk = Z[k];
        /* Its relative rank in {[0 _ 1]}: */
        double Tk = ((double)j + 0.5)/((double) np);
        /* Compute {Uk} such that ball_cap_vol_frac_pos(Sdim,Uk) == Tk}: */
        /* fprintf(stderr, "  Tk = %.8f\n", Tk); */
        { 
          double Uold = Uk;
          /* First, find interval {U0,U1} that brackets {Uk}: */
          double U0 = Uk, U1 = Uk;
          double TM0 = TMk, TM1 = TMk;
          while (TM1 < Tk) 
            { double dU = (U1 == U0 ? Ustep : U1 - U0);
              U0 = U1; TM0 = TM1; 
              U1 = U1 + dU;
              TM1 = ball_cap_vol_frac_pos(Sdim,U1);
              /* fprintf(stderr, "+"); */
              /* fprintf(stderr, "  U1 = %.8f TM1 = %.8f\n", U1, TM1); */
            }
          assert(TM0 <= Tk);
          assert(TM1 >= Tk);
          /* Binary search on {Uk}: */
          while(TRUE)
            { Uk = (U0 + U1)/2;
              TMk = ball_cap_vol_frac_pos(Sdim,Uk);
              /* fprintf(stderr, "  Uk = %.8f TMk = %.8f\n", Uk, TMk); */
              if ((TM1 - TM0) <= Teps) { break; }
              if (Tk < TMk)
                { U1 = Uk; TM1 = TMk; /* fprintf(stderr, "<"); */ }
              else
                { U0 = Uk; TM0 = TMk; /* fprintf(stderr, ">"); */ }
            }
          /* Adjust {Ustep}: */
          Ustep = Uk - Uold;
          if (Ustep < 0.001/np) { Ustep = 0.001/np; }
        }
        /* Compute weight: */
        double Wk = 1.0;
        /* Accumulate scalar products: */
        SUU += Wk*Uk*Uk; SU1 += Wk*Uk; S11 += Wk;
        SZU += Wk*Zk*Uk; SZ1 += Wk*Zk;
      }
    /* Solve normal system {((SUU, SU1),(SU1,S11)) * (CU,C1) = (SZU,SZ1)}: */
    double D = SUU*S11 - SU1*SU1;
    double DU = SZU*S11 - SZ1*SU1;
    double D1 = SUU*SZ1 - SU1*SZU;
    double CU = DU/D;
    double C1 = D1/D;
    (*Zrad) = CU;
    (*Zmid) = C1;
    /* Consistency: */
    if ((*Zrad) <= 0.0) { (*Zrad) = 0.001; }
    if (DEBUG)
      { fprintf(stderr, " Zmid = %24.16e Zrad = %24.16e\n", (*Zmid), (*Zrad)); }
  }

void fit_global_sdim_zmid_zrad
  ( int32_t np, 
    int32_t *ix,        
    double *Z,
    int32_t minSdim,
    int32_t maxSdim,
    int32_t *Sdim, 
    double *Zmid, 
    double *Zrad
  )
  {
    /* Try all dimensions, remember best-fitting model: */
    int32_t bestSd = -1; double bestZm = 0.0, bestZr = 0.0; 
    double bestS2 = DBL_MAX;
    int32_t Sd;
    double TM[np]; /* Model-predicted ranks. */
    for (Sd = minSdim; Sd <= maxSdim; Sd++)
      { /* Compute trial dimension: */
        /* Find best values {Zm,Zr} for parameters {Zmid,Zrad}, fixing {Sdim=Sd}: */
        double Zm, Zr;
        fit_global_zmid_zrad(np, ix, Z, Sd, &Zm, &Zr);
        /* Check quality of fit: */
        compute_global_model_ranks(np, Z, TM, Sd, Zm, Zr);
        double S2 = compute_global_rank_diff_sqr(np, ix, TM);
        if (S2 < bestS2) { bestSd = Sd; bestZm = Zm; bestZr = Zr; bestS2 = S2; }
      }
    assert(bestS2 < DBL_MAX);
    assert(bestSd >= 0);
    (*Sdim) = bestSd;
    (*Zmid) = bestZm;
    (*Zrad) = bestZr;
  }
  
