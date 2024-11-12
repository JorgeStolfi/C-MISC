/* Tools for ODE integration. */
/* Last edited on 2005-10-23 19:37:36 by stolfi */

#ifndef SPTimeOnlySolveG3C1Aux_H
#define SPTimeOnlySolveG3C1Aux_H

#include <SPVector.h>
#include <SPTimeOnlyFunction.h> 
#include <SPTimeOnlyFuncMap.h> 
#include <SP1DSpline.h>
#include <SPSys.h>

#include <r2x2.h>
#include <bool.h>
#include <stdio.h>

/* MATRICES OF THE LINEAR SYSTEM */

void SPTimeOnlyG3C1_GetMatrices
  ( double cMass, 
    double cFriction,
    double cSpring, 
    double tStep,
    TBasis *tbas,
    r2x2_t N[],
    bool_t verbose
  );
  /* Obtains the matrices {N[0],N[1]} of the linear system (6) for the
    diffusion-decay equation (1,2). The arguments {cSpring},
    {cFriction}, and {cMass} are the constants {K}, {R} and {M} of
    equation (2), respectively. The matrices are allocated as needed. */

/* INTEGRATION LOOP */

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
  );
  /* Computes the next frame {a[0]} of the approximate solution, assumed
    to be centered at the given time {tj}, given the previous frame {a[1]}
    and the right-hand side operator {RHS}. */

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
  );
  /* Given the frame {a[1]} for the previous epoch, and a tentative
    guess {a[0]} for the current frame, computes an improved guess
    {a0New}, by solving the system of equations (6). The right-hand
    side {d} is computed as per equation (4), where the function {h} is
    defined by {a[0..1]} through the operator {RHS}. */
    
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
  );
  /* Given frames {a[0]} and {a[1]} for the epochs {j} and {j-1}, respectively,
    computes the right-hand side {d} of equation system (6). Uses equation (4),
    where the function {h} is defined by {a[0]} and {a[1]} through the
    operator {RHS}. Assumes that the gauge functions are confined to the 
    interval between the two epochs. */

/* EVALUATION OF APPROX SOLUTION */

double SPTimeOnlyG3C1_EvalApproxSolution
  ( SPVector a[], 
    double t, 
    double tj, 
    double tStep, 
    TBasis *tbas,
    int diff,
    bool_t verbose
  );
  /* Evaluates the derivative of order {diff}
    of the approximate solution {app} at time {t}.
    
    Assumes that {t} is in the interval {[tj-tStep __ tj]},
    and that {a[k][u]} is the coefficient 
    of {tau[j-k][u]}, the basis function of index {u} at
    time {tj-k*tStep}, for {k=0,1} and {u = 0,1}. */

/* PLOTTING AND CHECKING THE APPROX SOLUTION */

void SPTimeOnlyG3C1_PlotSolution
  ( FILE *plotWr,
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
  );
  /* Given the two frames {a[0],a[1]} bracketing the time interval
    {tj-tStep} to {tj}, evaluates the approximate solution {app} and
    the true solution {sol} at {plotSteps+1} equally spaced times {t}
    in that interval, and writes that data to {pltWr}. Writes also the
    derivatives of {app}, the left-hand and the right-hand sides of
    the differential equation, the residue {lhs-rhs}, and the error
    {app(t)-sol(t)}. The parameters {a,tj,tStep,tbas} are interpreted
    as in {EvalApproxSolution}. */

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
  );
  /* Given the two frames {a[0],a[1]} bracketing the time interval
    {tj-tStep} to {tj}, evaluates the approximate solution {app} and
    the true solution {sol} at {checkSteps+1} equally spaced times {t}
    in that interval, and computes the error {err(t) = app(t)-sol(t)}. 
    Computes also the residue {res(t)}, that is, the left-hand 
    side of the differential equation, minus the right-hand side. 
    Updates the variables {errMax} and {resMax} to show
    the maximum of {abs(err(t))} and {abs(res(t))}, respectively;
    these variables must be initialized by the client. 
    The parameters {a,tj,tStep,tbas} are interpreted
    as in {EvalApproxSolution}. */

#endif
