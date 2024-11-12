/* Tools for time-space PDE integration. */
/* Last edited on 2005-10-29 06:01:08 by stolfi */

#ifndef SPTimeSpaceSolveG1C0Aux_H
#define SPTimeSpaceSolveG1C0Aux_H

#include <SPBasic.h>
#include <SPTriang.h>
#include <SPMatrix.h>
#include <SPVector.h>
#include <SPSys.h>
#include <SPFunction.h>
#include <SPTimeSpaceFunction.h> 
#include <SPTimeSpaceFuncMap.h> 
#include <SP1DSpline.h>

#include <bool.h>
#include <stdio.h>
  
/* MATRICES OF LINEAR SYSTEM */

void SPTimeSpaceG1C0_GetMatrices
  ( char *matName, 
    double cDiffuse, 
    double cDrift,
    double cDecay, 
    double tStep,
    STBasis *stbas,
    SPMatrix N[],
    bool_t verbose
  );
  /* Obtains the matrices {N[0],N[1]} of the linear system (6) for the
    diffusion-decay equation (1,2). The arguments {cDiffuse},
    {cDrift}, and {cDecay} are the constants {F}, {R} and {L} of
    equation (2), respectively. The matrices are allocated as needed.
    
    When computing scalar products of sphere-time elements, the
    space-domain factor is obtained from the scalar product matrices
    {<sigma[j]|sigma[i]>}, {<V¤SGrd(sigma[j])|sigma[i]>}, and
    {<SLap(sigma[j])|sigma[i]>}, which are read from files
    "{matName}-ev.mat", "{matName}-vg.mat", "{matName}-sl.mat",
    respectively. */
  
/* PROGRESSIVE INTEGRATION */

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
  );
  /* Computes the next frame {a[0]} of the approximate solution, assumed
    to be centered at the given time {tj}, given the previous frame {a[1]}
    and the right-hand side operator {RHS}. */

void SPTimeSpaceG1C0_RefineFrame
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
    SPVector a0New,
    int smpOrderTime,
    SPSys_LinOptions_t *lso,
    bool_t verbose
  );
  /* Given the frame {a[1]} for the previous epoch, and a tentative
    guess {a[0]} for the current frame, computes an improved guess
    {a0New}, by solving the system of equations (6). The right-hand
    side {d} is computed as per equation (4), where the function {h} is
    defined by {a[0..1]} through the operator {RHS}. */
    
/* EVALUATING THE APPROXIMATE SOLUTION */

/* The procedures below are given frames {a[0]} for epoch {j} (time
  {tj}) and {a[1]} for epoch {j-1} (time {tj-tStep}), for a space-time
  function {f} described in terms of basis {stbas}.
  The time {t} must lie in the interval
  {[tj-tStep__tj]}. Frame {a[0]} is not used if {t == tj-tStep}, and
  frame {a[1]} is not used if {t==tj}. */

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
  );
  /* Computes the right-hand side {d} of equation system (6). Uses equation (4),
    where the function {h} is defined by {a[0]} and {a[1]} through the
    operator {RHS}. Assumes that the gauge functions are confined to the 
    interval between the two epochs. */

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
  );
  /* Evaluates {TD(f,diff)(e,t)}. */
    
SPFunction *SPTimeSpaceG1C0_TimeSlice
  ( SPVector a[], 
    double t, 
    double tj, 
    double tStep, 
    STBasis *stbas,
    Triangulation *tri, 
    int diff,
    bool_t verbose
  );
  /* Returns a spherical function {gt} (independent of time) such that
    {gt(e) = TD(f,diff)(e,t)} for all points {e} and the given time {t}. */

SPFunction *SPTimeSpaceG1C0_TimeSliceSingle
  ( SPVector a, 
    double t, 
    double tStep, 
    STBasis *stbas,
    Triangulation *tri, 
    int diff,
    bool_t verbose
  );
  /* Same as {SPTimeSpaceG1C0_TimeSlice}, but uses a single frame {a}
    for epoch {j}, and returns the snapshot at time {t==tj}. */

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
  );
  /* Returns the coefficients {c[0..m-1]}, relative to the spatial
    basis {stbas->sBas}, of the spherical function {gt}
    returned by {SPTimeSpaceG1C0_TimeSlice}. */

void SPTimeSpaceG1C0_TimeSliceCoeffsSingle
  ( SPVector a, 
    double tj, 
    double tStep, 
    STBasis *stbas,
    Triangulation *tri, 
    int diff,
    bool_t verbose,
    double c[]
  );
  /* Same as {SPTimeSpaceG1C0_TimeSliceCoeffs}, but uses a single frame {a}
    for epoch {j}, and returns the snapshot at time {t==tj}. */

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
  );
  /* Evaluates the approximate solution {app} and the true solution
    {sol} at {smpOrderTime+1} equally spaced times {t} in the interval
    {[tj-tStep _ tj]}, and computes the error {err(e,t) =
    app(e,t)-sol(e,t)}. Computes also the residue {res(t)}, that is,
    the left-hand side of the differential equation, minus the
    right-hand side. Updates the variables {errMax} and {resMax} to
    show the maximum of {abs(err(t))} and {abs(res(t))}, respectively;
    these variables must be initialized by the client. The parameters
    {a,tj,tStep,stbas,tri} are interpreted as in 
    {SPTimeSpaceG1C0_EvalApproxSolution}. */

#endif
