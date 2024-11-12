/* See SPTimeOnlySolveG3C1Aux.h */
/* Last edited on 2005-10-23 19:43:04 by stolfi */ 

#include <SPTimeOnlySolveG3C1Aux.h>

#include <SPVector.h>
#include <SPFunction.h>
#include <SPTimeOnlyFunction.h> 
#include <SPTimeOnlyFuncMap.h> 
#include <SPTimeOnlyIntegral.h> 
#include <SPTimeBasisG3C1.h>
#include <SPPlot.h>
#include <SPSys.h>

#include <rn.h>
#include <r2x2.h>
#include <vec.h>
#include <bool.h>
#include <affirm.h>

#include <stdio.h>
#include <math.h>

void SPTimeOnlyG3C1_Solve2x2System(r2x2_t *M, SPVector y, SPVector x);
  /* Solves the {2×2} system {M x = y}, stores the solution in {x}. */

/* 
  BASIS PRODUCTS
  
  The time elements implemented so far are {tau[0,u](t)}, 
  {u=0,1]} where 
  {
    tau[0,0](t) = 
      (z + 1)^2(1 - 2z)   for z in [-1 __ 0]
      (z - 1)^2(2z + 1)   for z in [0 __ +1]
      
    tau[0,1](t) = 
      z(z + 1)^2*tStep    for z in [-1 __ 0]
      z(z - 1)^2*tStep    for z in [0 __ +1]
  }
  where {z = t/tStep}. The gauge elements are {pi[0,v]},
  {v=0,1}, where
  {
    pi[0,0](t) = 1       for z in [-1 __ 0]
    pi[0,1](t) = 2z + 1  for z in [-1 __ 0]
  }
*/

void SPTimeOnlyG3C1_ComputeNextFrame
  ( double tj, 
    double tStep,
    SPVector a[],
    r2x2_t N[],
    SPTimeOnlyFuncMap RHS,
    TBasis *tbas, 
    int smpOrder,
    SPSys_GenOptions_t *gso,   /* Parameters for non-linear system solver. */
    bool_t verbose
  )
  { int frameDim = a[1].ne;
    SPVector a0New = double_vec_new(frameDim); /* Work vector. */
    double a0Change, a0Norm; 
    int iter = 0;
    if (verbose) { fprintf(stderr, "guessing initial approximation...\n"); }
    SPSys_GuessSol(a[0]); 
    while (TRUE)
      { if (verbose) 
          { fprintf(stderr, "=== iteration %d", iter); 
            fprintf(stderr, " time = [%8g __ %8g] ===\n", tj-tStep, tj); 
          }
        
        /* Check for termination: */
        a0Norm = rn_norm(frameDim, a[0].e);
        if (SPSys_GenStopCondition(gso, 0.0, a0Norm, iter, 1.0, a0Change, TRUE))
          { break; }
        
        if (verbose) { fprintf(stderr, "refining approximation...\n"); }
        SPTimeOnlyG3C1_RefineFrame
          (tj, tStep, a, N, RHS, tbas, a0New, smpOrder, verbose);
        if (verbose) { fprintf(stderr, "updating solution coefficients...\n"); } 
        a0Change = SPSys_UpdateSolution(a0New, a[0]);

        iter++;
      }
    free(a0New.e);
  }

void SPTimeOnlyG3C1_RefineFrame
  ( double tj, 
    double tStep,
    SPVector a[],
    r2x2_t N[],
    SPTimeOnlyFuncMap RHS,
    TBasis *tbas, 
    SPVector a0New,
    int smpOrder,
    bool_t verbose
  )
  { int frameDim = a[0].ne;
    SPVector d = double_vec_new(frameDim);
    if (verbose) { fprintf(stderr, "computing right-hand side of system...\n"); } 
    SPTimeOnlyG3C1_ComputeSystemRightHandSide
      (tj, tStep, a, N, RHS, tbas, smpOrder, d, verbose);
    if (verbose) { fprintf(stderr, "solving system...\n"); }
    if (verbose)
      { int u, v;
        fprintf(stderr, "  linear system =\n");
        for (v = 0; v < 2; v++) 
          { for (u = 0; u < 2; u++)
              { fprintf(stderr, " %16g", N[0].c[v][u]); }
            fprintf(stderr, "     %16g", d.e[v]); 
            fprintf(stderr, "\n");
          }
      }
    SPTimeOnlyG3C1_Solve2x2System(&(N[0]), d, a0New);
    if (verbose)
      { int i;
        fprintf(stderr, "  aNew = (");
        for (i = 0; i < a0New.ne; i++) 
          { fprintf(stderr, " %16g", a0New.e[i]); }
        fprintf(stderr, " )\n");
      }
    free(d.e);
  }

void SPTimeOnlyG3C1_Solve2x2System(r2x2_t *M, SPVector y, SPVector x)
  {
    affirm(y.ne == 2, "y wrong size");
    affirm(x.ne == 2, "x wrong size");
    double det = M->c[0][0]*M->c[1][1] - M->c[1][0]*M->c[0][1];
    x.e[0] = (+ M->c[1][1]*y.e[0] - M->c[0][1]*y.e[1])/det;
    x.e[1] = (- M->c[1][0]*y.e[0] + M->c[0][0]*y.e[1])/det;
  }

void SPTimeOnlyG3C1_ComputeSystemRightHandSide
  ( double tj,
    double tStep,
    SPVector a[], 
    r2x2_t N[], 
    SPTimeOnlyFuncMap RHS, 
    TBasis *tbas,
    int smpOrder,
    SPVector d,
    bool_t verbose
  )
  { 
    if (verbose)
      { fprintf(stderr, "  tStep = %g", tStep);
        fprintf(stderr, "  smpOrder = %d", smpOrder);
        fprintf(stderr, "  tj = %g\n", tj);
      }
      
    affirm(tbas->tFam.c == 1, "bad continuity");
    affirm(tbas->tFam.g == 3, "bad degree");
    affirm(tbas->tFam.nmp = 2, "bad nmp");
    
    double_vec_t st = double_vec_new(smpOrder);
    double_vec_t wt = double_vec_new(smpOrder);
    int ns = 0;
    SPTimeOnlyIntegral_GaussSampleInterval(tj-tStep, tj, smpOrder, &st, &wt, &ns);  
    if (FALSE) 
      { int i;
        fprintf(stderr, "  sampling = (");
        for (i = 0; i < st.ne; i++)
          { fprintf(stderr, " %5g(%5g)", st.e[i], wt.e[i]); }
        fprintf(stderr, " )\n");
      }
    int u, v;
    for (v = 0; v < 2; v++)
      { 
        /* Compute {d[v] = INTEGRAL(RHS(app(t),t)·pi[j,v](t))}
          for {t} between {tj-tStep} and {tj}. */

        auto double integrand(double t);
        /* Computes {h(t) = pi[0,v](t-tj)·rhs(t)}. */
        
        double integrand(double t)
          { double appt = 
              SPTimeOnlyG3C1_EvalApproxSolution(a, t, tj, tStep, tbas, 0, FALSE);
            double rhst = RHS.map(appt,t);
            double pit = SPTimeBasisG3C1_GaugeEval(v, t, tj, tStep);
            return rhst*pit;
          }
        
        double sum = 0, corr = 0;
        SPTimeOnlyIntegral_BySamples(integrand, st, wt, &sum, &corr);
        d.e[v] = sum;
        if (verbose) { fprintf(stderr, "  v = %d integral = %16g\n", v, sum); }
        
        /* Subtract from {d} the product {N[1] a[0]} */
        
        for (u = 0; u <= 1; u++) { d.e[v] -= N[1].c[v][u]*a[1].e[u]; }
      }
    free(st.e); free(wt.e);
  }
  
void SPTimeOnlyG3C1_PlotSolution
  ( FILE *pltWr,
    SPVector a[], 
    double tj, 
    double tStep, 
    TBasis *tbas,
    double cMass,
    double cFriction,
    double cSpring,
    SPTimeOnlyFuncMap RHS, 
    SPTimeOnlyFunction *sol,
    int plotSteps
  )
  {
    double eps = 1.0e-7*tStep;
    int i;
    for (i = 0; i <= plotSteps; i++)
      { 
        double z = ((double)i)/((double)plotSteps);
        double t = (tj - tStep) + z*tStep;
        /* Adjust {t} to lie strictly inside the time interval: */
        if (i == 0) { t += eps; }
        if (i == plotSteps) { t -= eps; }
        /* Evaluate solutions and other things at time {t}: */
        double appt = 
          SPTimeOnlyG3C1_EvalApproxSolution(a, t, tj, tStep, tbas, 0, FALSE);
        double Dappt = 
          SPTimeOnlyG3C1_EvalApproxSolution(a, t, tj, tStep, tbas, 1, FALSE);
        double DDappt = 
          SPTimeOnlyG3C1_EvalApproxSolution(a, t, tj, tStep, tbas, 2, FALSE);
        double lhst = cMass*DDappt + cFriction*Dappt + cSpring*appt;
        double rhst = RHS.map(appt,t);
        double rest = lhst - rhst;
        double solt = sol->m->eval(sol, t);
        double errt = appt - solt;
        
        /* Write plot file line: */
        fprintf(pltWr, "%16g",  t);
        fprintf(pltWr, "  %16g %16g %16g",  appt, Dappt, DDappt);
        fprintf(pltWr, "  %16g %16g %16g",  lhst, rhst, rest);
        fprintf(pltWr, "  %16g %16g",  solt, errt);
        fprintf(pltWr, "\n");
      }
  }

void SPTimeOnlyG3C1_CheckSolution
  ( SPVector a[], 
    double tj, 
    double tStep, 
    TBasis *tbas,
    double cMass,
    double cFriction,
    double cSpring,
    SPTimeOnlyFuncMap RHS, 
    SPTimeOnlyFunction *sol,
    int checkSteps,
    double *errMax,
    double *resMax
  )
  {
    double eps = 1.0e-7*tStep;
    int i;
    for (i = 0; i <= checkSteps; i++)
      { 
        double z = ((double)i)/((double)checkSteps);
        double t = (tj - tStep) + z*tStep;
        /* Adjust {t} to lie strictly inside the time interval: */
        if (i == 0) { t += eps; }
        if (i == checkSteps) { t -= eps; }
        /* Evaluate solutions and other things at time {t}: */
        double appt = 
          SPTimeOnlyG3C1_EvalApproxSolution(a, t, tj, tStep, tbas, 0, FALSE);
        double Dappt = 
          SPTimeOnlyG3C1_EvalApproxSolution(a, t, tj, tStep, tbas, 1, FALSE);
        double DDappt = 
          SPTimeOnlyG3C1_EvalApproxSolution(a, t, tj, tStep, tbas, 2, FALSE);
        double lhst = cMass*DDappt + cFriction*Dappt + cSpring*appt;
        double rhst = RHS.map(appt,t);
        double rest = lhst - rhst;
        double solt = sol->m->eval(sol, t);
        double errt = appt - solt;
        
        /* Update max errors: */
        double aerrt = fabs(errt);
        if (aerrt > (*errMax)) { (*errMax) = aerrt; }
        double arest = fabs(rest);
        if (arest > (*resMax)) { (*resMax) = arest; }
      }
  }

double SPTimeOnlyG3C1_EvalApproxSolution
  ( SPVector a[], 
    double t, 
    double tj, 
    double tStep, 
    TBasis *tbas,
    int diff,
    bool_t verbose
  )
  {
    double z = (t - tj)/tStep + 1;
    if (verbose) { fprintf(stderr, "  t = %8g  z = %8g\n", t, z); }
    /* Compute {appt = app(t)}: */
    int u, k;
    double appt = 0.0;
    for (u = 0; u <= 1; u++)
      { for (k = 0; k <= 1; k++)
          { if (verbose) { fprintf(stderr, "    u = %d  k = %d", u, k); }
            double aku = a[k].e[u];
            /* Compute {tau[j-k,u](t) = tau[0,u](t-tj-k*tStep))} */
            double taut = SPTimeBasisG3C1_BasisEval(u, t, tj - k*tStep, tStep, diff);
            if (verbose) { fprintf(stderr, "  tau(t) = %16g\n", taut); }
            /* Accumulate in {appt}: */
            appt += aku*taut;
          }
      }
    return appt;
  }
  
void SPTimeOnlyG3C1_GetMatrices
  ( double cMass, 
    double cFriction,
    double cSpring, 
    double tStep,
    TBasis *tbas,
    r2x2_t N[],
    bool_t verbose
  )
  { 
    /* 
      We must build { N[k][v][u] == <D(tau[j-k,u]) | pi[j,v]> } 
      that is 

        { < M·TDD(tau')(t) + R·TD(tau')(t) + K·tau'(t) | pi[j,v]> }

      where {M = cMass}, {R = cFriction}, {K = cSpring}, and {tau'} is
      short for {tau[j-k,u]}. Therefore we can write
      
        { N[k] == M·NTDD[k] + R*NTD[k] + K·NT[k] }
        
      where 
        
        { NTDD[k][v][u] := < TDD(tau[j-k,u]) | pi[j,v] >
          NTD[k][v][u] := < TD(tau[j-k,u]) | pi[j,v] >
          NT[k][v][u] := < tau[j-k,u] | pi[j,v] >
        }
     
     These matrices are small and can be quickly computed.
    */

    int tDim = tbas->tFam.nmp;      /* Number of temporal elements per epoch. */

    /* Computing the matrices: */
    int u, v, k;
    for (k = 0; k < 2; k++)
      { if (verbose) { fprintf(stderr, "  matrix N[%d] =\n", k); }
        for (u = 0; u < tDim; u++)
          { if (verbose) { fprintf(stderr, "    "); }
            for (v = 0; v < tDim; v++)
              { /* Compute the time-domain scalar products: */
                double NTkuv = SPTimeBasisG3C1_BasisGaugeDot(k,u,v, tStep, 0, verbose);
                double NTDkuv = SPTimeBasisG3C1_BasisGaugeDot(k,u,v, tStep, 1, verbose);
                double NTDDkuv = SPTimeBasisG3C1_BasisGaugeDot(k,u,v, tStep, 2, verbose);
                N[k].c[v][u] = cMass*NTDDkuv + cFriction*NTDkuv + cSpring*NTkuv;
                if (verbose) { fprintf(stderr, " %16g", N[k].c[v][u]); }
              }
            if (verbose) {fprintf(stderr, "\n"); }
          }
      }
  }
