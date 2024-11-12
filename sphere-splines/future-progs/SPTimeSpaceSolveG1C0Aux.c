/* See SPTimeSolveG1C0Aux.h */
/* Last edited on 2005-10-29 06:11:26 by stolfi */ 

#include <SPTimeSpaceSolveG1C0Aux.h>

#include <SPTriang.h>
#include <SPMatrix.h>
#include <SPVector.h>
#include <SPSys.h>
#include <SPFunction.h>
#include <SPRange.h>
#include <SPIntegral.h>
#include <SPTimeBasisG1C0.h>
#include <SPTimeSpaceFunction.h> 
#include <SPTimeSpaceFuncMap.h> 
#include <SPTimeSpaceIntegral.h> 
#include <SPApprox.h>

#include <rn.h>

#include <vec.h>
#include <bool.h>
#include <affirm.h>

#include <math.h>
#include <stdio.h>

void SPTimeSpaceG1C0_SummarizeFunction
  ( ScalarField func,
    Triangulation *tri, 
    int smpOrder,
    double *funcMax,
    double *funcSum
  );
/* Evaluates {func(e)} for many points {e} on the sphere, and updates
  {funcMax} with the maximum absolute value seen. Also computes the
  integral of {func(e)^2} over the sphere, and adds that to {funcSum}.
  The client must initialize {funcMax} and {funcSum} before the call.
  
  The sample points for both operations are the union of triangular
  arrays of {smpOrder*(smpOrder+1)/2} points in each triangle of the
  triangulation {tri}. */ 

/* 
  BASIS PRODUCTS
  
  The time basis we use consists of trinagular tent functions.
  There is only one element per epoch (that is, {q == 1}),
  namely
  {
    tau[0,0](t) = 
      (1 + z)   for z in [-1 __ 0]
      (1 - z)   for z in [0 __ +1]
  }
  where {z = t/tStep}. Similarly, there is only one
  gauge element per epoch,
  {
    pi[0,0](t) = 1       for z in [-1 __ 0]
  }
*/

void SPTimeSpaceG1C0_ComputeNextFrame
  ( double tj, 
    double tStep,
    SPVector a[],
    SPMatrix N[],
    SPMatrix N0L,  /* Left factor of {N[0]}. */  
    SPVector N0D,  /* Middle factor of {N[0]}. */
    SPMatrix N0R,  /* Right factor of {N[0]}. */ 
    SPTimeSpaceFuncMap RHS,
    STBasis *stbas, 
    Triangulation *tri,  
    int smpOrderTime,
    SPSys_LinOptions_t *lso,  /* Parameters for linear system solver. */
    SPSys_GenOptions_t *gso,  /* Parameters for non-linear system solver. */
    bool_t verbose
  )
  { int frameDim = a[1].ne;
    SPVector a0New = double_vec_new(frameDim); /* Work vector. */
    double a0Change, a0Norm;
    int iter = 0;
    while (TRUE)
      { fprintf(stderr, "=== iteration %d ===\n", iter);
        if (iter == 0)
          { fprintf(stderr, "guessing initial approximation...\n"); 
            SPSys_GuessSol(a[0]); 
          }
        else
          { fprintf(stderr, "refining approximation...\n"); 
            SPTimeSpaceG1C0_RefineFrame
              ( tj, tStep, a, N, N0L, N0D, N0R, RHS, stbas, tri, a0New, 
                smpOrderTime, lso, verbose
              );
            fprintf(stderr, "updating solution coefficients...\n"); 
            a0Change = SPSys_UpdateSolution(a0New, a[0]);
            if (a0Change <= gso->relTol*a0Norm) { break; }
            if (a0Change <= gso->absTol) { break; }
            if (iter >= gso->maxIter) { break; }
          }
        a0Norm = rn_norm(a[0].ne, a[0].e);
        iter++;
      }
    free(a0New.e);
  }

void SPTimeSpaceG1C0_RefineFrame
  ( double tj, 
    double tStep,
    SPVector a[],
    SPMatrix N[],
    SPMatrix N0L,
    SPVector N0D,
    SPMatrix N0R,
    SPTimeSpaceFuncMap RHS,
    STBasis *stbas, 
    Triangulation *tri, 
    SPVector a0New,
    int smpOrderTime,
    SPSys_LinOptions_t *lso,
    bool_t verbose
  )
  { int frameDim = a[0].ne;
    SPVector d = double_vec_new(frameDim);
    SPVector y = double_vec_new(frameDim); /* Work vector. */
    fprintf(stderr, "computing system's right-hand side...\n"); 
    SPTimeSpaceG1C0_ComputeSystemRightHandSide
      (tj, tStep, a, N, RHS, stbas, tri, smpOrderTime, d, verbose);
    fprintf(stderr, "solving system...\n");
    SPSys_LinSolve(N[0], N0L, N0D, N0R, d, lso, a[0], a0New, TRUE);
    free(d.e); free(y.e);
  }

void SPTimeSpaceG1C0_ComputeSystemRightHandSide
  ( double tj,
    double tStep,
    SPVector a[], 
    SPMatrix N[], 
    SPTimeSpaceFuncMap RHS, 
    STBasis *stbas,
    Triangulation *tri, 
    int smpOrderTime,
    SPVector d,
    bool_t verbose
  )
  { 
    if (verbose)
      { fprintf(stderr, "  tStep = %g", tStep);
        fprintf(stderr, "  smpOrderTime = %d", smpOrderTime);
        fprintf(stderr, "  tj = %g\n", tj);
      }
      
    int sDim = stbas->sBas.ne;             /* Dimension of space basis. */
    
    TimeSpaceFuncMap TSFzero = SPTimeSpaceFuncMap_FromName("TSFzero");

    affirm(stbas->tFam.c == 0, "bad continuity");
    affirm(stbas->tFam.g == 1, "bad degree");
    affirm(stbas->tFam.nmp = 1, "bad nmp");
    
    int nv = 1; /* Number of gauge elements per epoch. */

    /* Compute the product {h = N[1] a[1]} */
    SPVector h = double_vec_new(nv*sDim); /* Work area. */
    SPMatrix_MulCol(N[1], a[1], h);
    
    int r, v;
    
    /* Nest we must compute the dot products
        {d[rv] = INTEGRAL(RHS(app(e,t),e,t)·rho[r](e)·pi[j,v](t))}
      for {t} between {tj-tStep} and {tj}. The expensive part here is
      evaluating {app(e,t)} at many samples points {e} and sample
      times {t}. In order to speed up this evaluation, we build two
      snapshots {app0} and {app1} of the approximate solution, such
      that {app0(e) = app(e,tj)} and {app1(e) = app(e,tj-tStep)}.
      Then {app(e,t)} can be obtained by linear interpolation
      between {app0(e)} and {app1(e)}. */
    
    SPFunction *app0 = SPTimeSpaceG1C0_TimeSliceSingle
      (a[0], tj, tStep, stbas, tri, 0, verbose);
    SPFunction *app1 = SPTimeSpaceG1C0_TimeSliceSingle
      (a[1], tj-tStep, tStep, stbas, tri, 0, verbose);
    
    for (r = 0; r < sDim; r++)
      {
        for (v = 0; v < 1; v++)
          { 
            if (verbose) { fprintf(stderr, "  r = %d v = %d", r, v); }
            
            int rv = sDim*v + r;
            
            /* Compute the dot product
                {rhs_theta_rv = INTEGRAL(RHS(app(e,t),e,t)·rho[r](e)·pi[j,v](t))}
              for {t} between {tj-tStep} and {tj}. */
            double rhs_theta_rv;
              
            if (RHS.map == TSFzero.map) 
              { /* Right-hand side is identically zero: */
                rhs_theta_rv = 0;
              }
            else
              { /* Right-hand side is not identically zero, must integrate: */
              
                auto double integrand(S2Point *e, double t);
                /* Computes {h(t) = pi[0,v](t-tj)·rhs(t)}. */

                double integrand(S2Point *e, double t)
                  { double app_v1 = app1->m->eval(app1, e); /* {app(e,tj-tStep)}. */
                    double app_v0 = app0->m->eval(app0, e); /* {app(e,tj)}. */
                    double s1 = (tj - t)/tStep, s0 = 1 - s1;
                    double app_v = s0*app_v0 + s1*app_v1; /* {app(e,t)}. */
                    double rhs_v = RHS.map(app_v,e,t);
                    double pit_v = SPTimeBasisG1C0_GaugeEval(v, t, tj, tStep);
                    return rhs_v*pit_v;
                  }

                rhs_theta_rv = SPTimeSpaceIntegral_OnSphereInterval
                  (integrand, tri, smpOrderTime, tj-tStep, tj);  
              }              
            if (verbose) { fprintf(stderr, " RHS·theta[rjv] = %16g", rhs_theta_rv); }
            /* Subtract the corresponding element of {N[1] a[1]}: */
            d.e[rv] = rhs_theta_rv - h.e[rv];
            if (verbose) { fprintf(stderr, " d[rv] = %16g\n", d.e[rv]); }
         }
      }
    app0->m->free(app0);
    app1->m->free(app1);
  }

void SPTimeSpaceG1C0_TimeSliceCoeffs
  ( SPVector a[], 
    double t, 
    double tj, 
    double tStep, 
    STBasis *stbas,
    Triangulation *tri, 
    int diff,
    bool_t verbose,
    double c[]
  )
  { 
    /* Get parameters of temporal basis: */
    affirm(stbas->tFam.g == 1, "bad degree");
    affirm(stbas->tFam.c == 0, "bad continuity");
    int tDim = stbas->tFam.nmp;  /* Number {q} of temporal elements per epoch. */
    int sDim = stbas->sBas.ne;  /* Number {m} of spatial elements. */
    
    int i;
    /* Evaluate the time basis elements at time {t}: */
    double tauku[2*tDim];
    { int k;
      for (k = 0; k < 2; k++)
        { int u;
          for (u = 0; u < tDim; u++)
            { int ku = 2*u + k;
              tauku[ku] = SPTimeBasisG1C0_BasisEval(u, t, tj-k*tStep, tStep, diff); 
            }
        }
    }
    /* Compute the space basis coefficients: */
    for (i = 0; i < sDim; i++)
      { double s = 0.0;
        int k;
        for (k = 0; k < 2; k++)
          { int u;
            for (u = 0; u < tDim; u++)
              { int iu = sDim*u + i;
                int ku = 2*u + k;
                s += a[k].e[iu] * tauku[ku]; 
              }
          }
        c[i] = s;
      }
  }

void SPTimeSpaceG1C0_TimeSliceCoeffsSingle
  ( SPVector a, 
    double tj, 
    double tStep, 
    STBasis *stbas,
    Triangulation *tri, 
    int diff,
    bool_t verbose,
    double c[]
  )
  { 
    /* Get parameters of temporal basis: */
    affirm(stbas->tFam.g == 1, "bad degree");
    affirm(stbas->tFam.c == 0, "bad continuity");
    int tDim = stbas->tFam.nmp;  /* Number {q} of temporal elements per epoch. */
    int sDim = stbas->sBas.ne;  /* Number {m} of spatial elements. */
    
    int i;
    /* Evaluate the time basis elements at time {t}: */
    double tau0u[tDim];
    { int u;
      for (u = 0; u < tDim; u++)
        { tau0u[u] = SPTimeBasisG1C0_BasisEval(u, tj, tj, tStep, diff); }
    }
    /* Compute the space basis coefficients: */
    for (i = 0; i < sDim; i++)
      { double s = 0.0;
        int u;
        for (u = 0; u < tDim; u++)
          { int iu = sDim*u + i;
            s += a.e[iu] * tau0u[u]; 
          }
        c[i] = s;
      }
  }

SPFunction *SPTimeSpaceG1C0_TimeSlice
  ( SPVector a[], 
    double t, 
    double tj, 
    double tStep, 
    STBasis *stbas,
    Triangulation *tri, 
    int diff,
    bool_t verbose
  )
  { int sDim = stbas->sBas.ne;  /* Number {m} of spatial elements. */
    double c[sDim];
    SPTimeSpaceG1C0_TimeSliceCoeffs
      (a, t, tj, tStep, stbas, tri, diff, verbose, c);
    return SPFunction_LinComb(c, stbas->sBas);
  }

SPFunction *SPTimeSpaceG1C0_TimeSliceSingle
  ( SPVector a, 
    double tj, 
    double tStep, 
    STBasis *stbas,
    Triangulation *tri, 
    int diff,
    bool_t verbose
  )
  { int sDim = stbas->sBas.ne;  /* Number {m} of spatial elements. */
    double c[sDim];
    SPTimeSpaceG1C0_TimeSliceCoeffsSingle
      (a, tj, tStep, stbas, tri, diff, verbose, c);
    return SPFunction_LinComb(c, stbas->sBas);
  }

double SPTimeSpaceG1C0_EvalApproxSolution
  ( SPVector a[], 
    S2Point *e,
    double t, 
    double tj, 
    double tStep, 
    STBasis *stbas,
    Triangulation *tri, 
    int diff,
    bool_t verbose
  )
  { int sDim = stbas->sBas.ne;  /* Number {m} of spatial elements. */
    double c[sDim];
    SPTimeSpaceG1C0_TimeSliceCoeffs
      (a, t, tj, tStep, stbas, tri, diff, verbose, c);
    int i;
    double s = 0;
    for (i = 0; i < sDim; i++)
      { 
        SPFunction *sigmai = stbas->sBas.e[i];
        s += c[i] * sigmai->m->eval(sigmai, e);
      }
    return s;
  }

void SPTimeSpaceG1C0_GetMatrices
  ( char *matName, 
    double cDiffuse, 
    double cDrift,
    double cDecay, 
    double tStep,
    STBasis *stbas,
    SPMatrix N[],
    bool_t verbose
  )
  { 
    /* 
      We must build { N[k][rv,iu] == <D(phi[i,j-k,u]) | theta[r,j,v]> } 
      that is 

        { < TD(phi) - F*SLap(phi') + v¤SGrd(phi') - L*f | theta[r,j,v]> }

      where {F = cDiffuse}, {V = cDrift}, {L = cDecay}, and {phi'} is
      short for {phi[i,j-k,u] = sigma[i](e)·tau[j-k,u](t)}. Therefore
      we can write
      
        { N[k] == NTD[k] + F*NSL[k] + V·NSG[k] - L*NF[k] }
        
      where 
        
        { NTD[k][rv,iu]
            := < TD(phi[i,j-k,u]) | theta[r,j,v] >
            == < TD(sigma[i]·tau[j-k,u]) | rho[r]·pi[j,v] >
            == < sigma[i]·TD(tau[j-k,u]) | rho[r]·pi[j,v] >
            == < sigma[i] | rho[r] > · < TD(tau[j-k,u]) | pi[j,v] > 
            == SE[r,i] · TDW[k][v,u]
            
          NSL[k][rv,iu] 
            := < SLap(phi[i,j-k,u]) | theta[r,j,v] > 
            == < SLap(sigma[i]·tau[j-k,u]) | rho[r]·pi[j,v] >
            == < SLap(sigma[i])·tau[j-k,u] | rho[r]·pi[j,v] >
            == < SLap(sigma[i]) | rho[r] > · < tau[j-k,u] | pi[j,v] > 
            == SL[r,i] · TW[k][v,u]
            
          NSG[k][rv,iu] 
            := < v¤SGrd(phi[i,j-k,u]) | theta[r,j,v] > 
            == < v¤SGrd(sigma[i]·tau[j-k,u]) | rho[r]·pi[j,v] >
            == < v¤SGrd(sigma[i])·tau[j-k,u] | rho[r]·pi[j,v] >
            == < (v/V)¤SGrd(sigma[i]) | rho[r] > · < tau[j-k,u] | pi[j,v] > 
            == SV[r,i] · TW[k][v,u]
            
          NF[k][rv,iu] 
            := < phi[i,j-k,u] | theta[r,j,v] > 
            == < sigma[i]·tau[j-k,u] | rho[r]·pi[j,v] >
            == < sigma[i] | rho[r] > · < tau[j-k,u] | pi[j,v] > 
            == SE[r,i] · TW[k][v,u]
        }
        
     and
     
        { SE[r,i] := < sigma[i] | rho[r] >
          SL[r,i] := < SLap(sigma[i]) | rho[r] >
          SV[r,i] := < V¤SGrd(sigma[i]) | rho[r] >
          
          TW[k][v,u] := < tau[-k,u] | pi[0,v] >
          TDW[k][v,u] := < TD(tau[-k,u]) | pi[0,v] >
        }
     
     If the velocity field {v} is known beforehand, the matrices
     {SE,SL,SV} can be precomputed. The matrices {TW[k]} and {TDW[k]}
     are small and can be quickly computed.
    */

    Basis sBas = stbas->sBas;        /* Spatial basis. */
    int sDim = sBas.ne;             /* Dimension of space basis. */
    affirm(stbas->tFam.g == 1, "bad time basis degree");
    affirm(stbas->tFam.c == 0, "bad time basis continuity");
    int tDim = stbas->tFam.nmp;      /* Number of temporal elements per epoch. */
    int colsN = sDim*tDim, rowsN = sDim*tDim;
    N[0] = SPMatrix_Null(colsN, rowsN); 
    N[1] = SPMatrix_Null(colsN, rowsN);

    SPMatrix SE = SPApprox_ReadMatrix(txtcat(matName, "-ev.mat"));
    SPMatrix SL = SPApprox_ReadMatrix(txtcat(matName, "-sl.mat"));
    SPMatrix SV = SPApprox_ReadMatrix(txtcat(matName, "-vg.mat"));
    
    /* Computing the matrices: */
    int u, v, k;
    for (u = 0; u < tDim; u++)
      { for (v = 0; v < tDim; v++)
          { for (k = 0; k < 2; k++)
              { /* Compute the time-domain scalar prods {TW[k,u,v],TDW[k,u,v]}: */
                double TWkvu = SPTimeBasisG1C0_BasisGaugeDot(k,u,v, tStep, 0, FALSE);
                double TDWkvu = SPTimeBasisG1C0_BasisGaugeDot(k,u,v, tStep, 1, FALSE);
                /* Compute the submatrix {N[k][rv,iu]} for current {u,v}: */
                double SEcoef = TDWkvu - cDecay*TWkvu;
                double SLcoef = cDiffuse*TWkvu;
                double SVcoef = cDrift*TWkvu;
                SPMatrix Wkvu = SPMatrix_Mix(SEcoef, SE, SLcoef, SL);
                SPMatrix Zkvu = SPMatrix_Mix(1.0, Wkvu, SVcoef, SV);
                /* Insert {Zkvu} in {N[k]}: */
                SPMatrix Mk = SPMatrix_MixBlock(1.0, N[k], 1.0, Zkvu, u*sDim, v*sDim);
                free(N[k].ents.e); 
                N[k] = Mk;
                free(Wkvu.ents.e); 
                free(Zkvu.ents.e);
              }
          }
      }
  }
  
void SPTimeSpaceG1C0_CheckSolution
  ( SPVector a[], 
    double tj, 
    double tStep, 
    STBasis *stbas,
    Triangulation *tri, 
    double cDiffuse,
    double cDrift,
    double cDecay,
    SPTimeSpaceFuncMap RHS, 
    SPTimeSpaceFunction *sol,
    int smpOrderTime,
    int smpOrderSpace,
    double *errMax, double *errAvg,
    double *resMax, double *resAvg
  )
  {
    double eps = 1.0e-7*tStep;
    int i;
    /* Initialize error accumulators: */
    (*errMax) = 0; (*resMax) = 0;
    double errSum = 0, resSum = 0; /* Integrals of squared error and residue. */
    for (i = 0; i <= smpOrderTime; i++)
      { 
        double z = ((double)i)/((double)smpOrderTime);
        double t = (tj - tStep) + z*tStep;
        /* Adjust {t} to lie strictly inside the time interval: */
        if (i == 0) { t += eps; }
        if (i == smpOrderTime) { t -= eps; }
        /* Evaluate solutions and other things at time {t}: */
        SPFunction *appt = 
          SPTimeSpaceG1C0_TimeSlice(a, t, tj, tStep, stbas, tri, 0, FALSE);

        SPFunction *dift = 
          SPTimeSpaceG1C0_TimeSlice(a, t, tj, tStep, stbas, tri, 1, FALSE);
        
        auto double eval_err(S2Point *e);
        auto double eval_res(S2Point *e);
          /* Return the error {app(e,t) - sol(e,t)} and the 
            residue {(LHS(app)-RHS(app))(e,t)} at the given
            point {e} and the current time {t}. */

        double eval_err(S2Point *e) 
          { double app_p = appt->m->eval(appt, e);
            double sol_p = sol->m->eval(sol, e, t);
            double err_p = app_p - sol_p;
            return err_p;
          }
        
        double eval_res(S2Point *e) 
          { 
            return 0.0;
            
            /* Compute time and space derivatives of approximate sol: */
            double app_p = appt->m->eval(appt, e);   /* {app(e,t)}. */
            double dif_p = dift->m->eval(dift, e);   /* {TD(app,1)(e,t)}. */
            r3_t   grd_p = SPFunction_SGrd(appt, e); /* {SGrd(app)(e,t)}. */
            double lap_p = SPFunction_SLap(appt, e); /* {SLap(app)(e,t)}. */
            
            /* Compute velocity term: */
            r3_t velocity = (r3_t){{ -e->c[1], e->c[0], 0 }}; /* {V(e)}. */
            double gdv_p = r3_dot(&grd_p, &velocity); /* {V(e)¤SGrd(f)(e,t)}. */
            
            /* Compute LHS, RHS, and residue: */
            double lhs_p = dif_p - cDiffuse*lap_p + cDrift*gdv_p - cDecay*app_p;
            double rhs_p = RHS.map(app_p,e,t);
            double res_p = lhs_p - rhs_p;
            
            return res_p;
          }
        SPTimeSpaceG1C0_SummarizeFunction
          (eval_err, tri, smpOrderSpace, errMax, &errSum);
        SPTimeSpaceG1C0_SummarizeFunction
          (eval_res, tri, smpOrderSpace, resMax, &resSum);
        
      }
      
    /* Compute RMS average over all slices. */
    (*errAvg) = sqrt(errSum/(smpOrderTime+1)/FOURPI);
    (*resAvg) = sqrt(resSum/(smpOrderTime+1)/FOURPI);
  }

void SPTimeSpaceG1C0_SummarizeFunction
  ( ScalarField func,
    Triangulation *tri, 
    int smpOrder,
    double *funcMax,
    double *funcSum
  )
  {
    auto double func_sqr(S2Point *e);
      /* The square of {func}. */

    double func_mx = SPRange_OnTriangulation(func, tri, smpOrder);
    if (func_mx > (*funcMax)) { (*funcMax) = func_mx; }

    (*funcSum) += SPIntegral_OnSphere(func_sqr, tri);

    double func_sqr(S2Point *e)
      { double fe = func(e); return fe*fe; }
  }
