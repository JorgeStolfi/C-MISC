/* SPTimeSpaceFunction.h -- General functions of the sphere × time */
/* Last edited on 2009-02-09 21:45:23 by stolfi */

#ifndef SPTimeSpaceFunction_H
#define SPTimeSpaceFunction_H

#include <SPTimeSpaceFuncMap.h>
#include <SPFunction.h>
#include <SPBasic.h>
#include <SPTriang.h>
#include <SPMatrix.h>

#include <udg_pulse.h>
#include <r3.h>
#include <SPH3.h>
#include <vec.h>
#include <bool.h>
#include <nat.h>
#include <stdio.h>
 
/* An {SPTimeSpaceFunction} object represents a /sphere-time function/, a real
  function whose domain is {S^2 × R}. The second parameter could be,
  for example, time in dynamic shallow-water models, or radius above
  ground in static geophisical modeling. */
  
#define OBJ void /* Actually {SPTimeSpaceFunction} or a subclass */

typedef double SPTimeSpaceFunction_EvalMth(OBJ *f, R3Point *p, double t);     
typedef void SPTimeSpaceFunction_WriteMth(OBJ *f, FILE *wr);        
typedef void SPTimeSpaceFunction_FreeMth(OBJ *f);        

typedef struct SPTimeSpaceFunction_Methods    /* Methods for any SPTimeSpaceFunction: */
  { SPTimeSpaceFunction_EvalMth *eval;    /* Evaluates the function at {p}. */
    SPTimeSpaceFunction_WriteMth *write;  /* Writes the function to {wr}, for reading back. */
    SPTimeSpaceFunction_FreeMth *free;    /* Reclaims the function's private storage. */
  } SPTimeSpaceFunction_Methods;
  /* Some of these methods may be missing for certain subclasses. */
  
typedef struct SPTimeSpaceFunction_Data    /* Data fields for any SPTimeSpaceFunction: */
  { } SPTimeSpaceFunction_Data;            
  
typedef struct SPTimeSpaceFunction /* Abstract class: a function defined on the sphere. */
  { char *type;              /* Type identifier */
    SPTimeSpaceFunction_Methods *m;  /* TimeSpaceFunction methods. */
    SPTimeSpaceFunction_Data *d;     /* TimeSpaceFunction parameters. */
  } SPTimeSpaceFunction;
  
#define SPTimeSpaceFunction_TypeId "STF."

SPTimeSpaceFunction *SPTimeSpaceFunction_Cast(OBJ *f);
  /* If {f} (according to {f->type}) is a subclass of {SPTimeSpaceFunction},
    returns {f} cast to that type; otherwise returns NULL. */

void SPTimeSpaceFunction_M_Write(SPTimeSpaceFunction *f, FILE *wr); 
  /* The default {write} method for a generic sphere-time function.
    Writes a standard {SPTimeSpaceFunction} header line and the specific
    subclass ID, then calls the subclass-specific {write} method, and
    closes off with the generic {SPTimeSpaceFunction} footer. */

SPTimeSpaceFunction *SPTimeSpaceFunction_Read(FILE *rd);  
  /* Reads a function's description from {rd}. */

/* BASES FOR SPHERE-TIME FUNCTION SPACES */

/* For sphere-time functions, we only consider tensor-product type bases 
  of the form {phi[i,j](p,t) == sigma[i](p)*tau[j](t)} where {sigma} is a 
  sherical function basis (independent of time) and {tau} is a 
  time function basis (idependent of space).  For the time being,
  the time basis is determined by a single parameter {timeCont},
  the continuity class; its elements are the splines 
  {S1DSpline_Eval(cont,r,z-k)} for {r} in {0..c} and integer {k}. */

typedef struct SPTimeSpaceFunction_Basis /* A tensor-product basis for sphere-time fns. */
  { SPFunction_Basis sBas;    /* The spatial basis. */
    udg_pulse_family_t tFam;  /* The time basis. */
  } SPTimeSpaceFunction_Basis;

typedef SPTimeSpaceFunction_Basis STBasis;

void SPTimeSpaceFunction_WriteBasis(FILE *wr, STBasis stbas);
  /* Writes a basis {stbas} to {wr}, in a format that can be 
    read back by {Read}.  */

STBasis SPTimeSpaceFunction_ReadBasis(FILE *rd);
  /* Reads from {rd} a basis for a space of sphere-time functions.
    Assumes the contents of {rd} was created by {write} above. */

/* INTEGRALS OF SPHERE-TIME FUNCTIONS */

double SPTimeSpaceFunction_Integral
  ( SPTimeSpaceFunction *f,
    TimeSpaceFuncMap FMap, 
    SPTimeSpaceFunction *w,
    Triangulation *tri,
    int smpOrderTime,
    double tMin, 
    double tMax, 
    bool_t verbose
  );
  /* Computes the integral of {w(p,t)*FMap(f(p,t),p,t)} for {p}
    ranging over the whole sphere and {t} ranging over the 
    interval {[tMin _ tMax]}.
    
    ??? If the function {f} happens to be an instance of {SPSpline},
    and {tri} is either NULL or {f}'s underlying triangulation, then
    the procedure reduces to {SPSpline_IntegralPW} (q.v.).
    Otherwise, the integral is computed with {SPTimeSpaceIntegral_OnSphere},
    with the same {tri} argument. */

double SPTimeSpaceFunction_RawIntegral
  ( SPTimeSpaceFunction *f,
    TimeSpaceFuncMap FMap, 
    SPTimeSpaceFunction *w,
    Triangulation *tri,
    int smpOrderTime,
    double tMin, 
    double tMax, 
    bool_t verbose
  );
  /* Like {SPTimeSpaceFunction_Integral}, but always computes the
    integral with {SPTimeSpaceIntegral_OnSphereInterval}, with the
    given {tri} argument, ignoring any underlying triangulation of
    {f}. */

double SPTimeSpaceFunction_CustomIntegral
  ( SPTimeSpaceFunction *f, 
    TimeSpaceFuncMap FMap, 
    SPTimeSpaceFunction *w,
    S2Point_vec_t sp,
    double_vec_t wp,
    double_vec_t st,
    double_vec_t wt
  );
  /* Computes the sample-based integral of {WF(p,t) =
    w(p,t)*FMap(f(p,t),p,t)}, defined as
    {SUM(wp[k]*wt[j]*WF(sp[k],st[j]), k=0..M-1, j=0..N-1)}, where {sp}
    is a client-specified set of sample points on the sphere, {wp} are
    the corresponding weights (areas of surface elements), {st} is a
    set of sample times, and {wt} are the corresponding weights
    (interval widths). */
    
/* INNER PRODUCT OF SPHERE-TIME FUNCTIONS */

/* These procedures compute weighted dot product of spherical functions,
  defined as integrals of {w(p,t)*H(f(p,t),g(p,t))} for {p} ranging
  over the sphere and {t} over some interval {[tMin _ tMax]},
  where {w} is a weight function, and {H()} is some pointwise operator.
  The conventions are similar to those of the corresponding integration
  procedures above.  In particular, all these procedures 
  assume {w(p,t) == 1.0} for all {p}, when {w == NULL}.*/

double SPTimeSpaceFunction_Dot
  ( SPTimeSpaceFunction *f,
    TimeSpaceFuncMap FMap, 
    SPTimeSpaceFunction *g,
    TimeSpaceFuncMap GMap,
    SPTimeSpaceFunction *w,
    Triangulation *tri,
    int smpOrderTime,
    double tMin, 
    double tMax, 
    bool_t verbose
  );
  /* Computes the weighted dot product of {F(p,t)=FMap(f(p,t),p,t)}
    and {G(p,t)=GMap(g(p,t),p,t)}, that is the integral of
    {w(p,t)*F(p,t)*G(p,t)} over the whole sphere and the interval
    {[tMin _ tMax]}. Uses Gaussian quadrature of order {smpOrderTime} in
    the time direction. (If {FMap == NULL}, then {F(p,t)} reduces to
    {f(p,t)}, and similarly for {GMap}.)
    
    ??? If both functions {f,g} happen to be instances of {SPSpline}
    with the same triangulation, and {tri} is either NULL or that same
    triangulation, then the procedure reduces to
    {SPSpline_DotBothPW} (q.v.).
    
    ??? Else, if either {f} or {g} happen to be an instance of {SPSpline}
    and {tri} is either NULL or that element's underlying triangulation,
    then the procedure reduces to {SPSpline_DotSinglePW} (q.v.).
    
    ??? In all other cases, the integral is computed with {SPTimeSpaceIntegral_OnSphere},
    with same {tri} argument. */

double SPTimeSpaceFunction_RawDot
  ( SPTimeSpaceFunction *f,
    TimeSpaceFuncMap FMap, 
    SPTimeSpaceFunction *g,
    TimeSpaceFuncMap GMap,
    SPTimeSpaceFunction *w,
    Triangulation *tri,
    int smpOrderTime,
    double tMin, 
    double tMax, 
    bool_t verbose
  );
  /* Like {SPTimeSpaceFunction_Dot}, but always computes the integral with
    {SPTimeSpaceIntegral_OnSphereInterval}, with the given {tri} argument,
    ignoring any underlying triangulation of {f} and {g}. */
    
double SPTimeSpaceFunction_CustomDot
  ( SPTimeSpaceFunction *f, 
    TimeSpaceFuncMap FMap, 
    SPTimeSpaceFunction *g,
    TimeSpaceFuncMap GMap,
    SPTimeSpaceFunction *w,
    S2Point_vec_t sp,
    double_vec_t wp,
    double_vec_t st,
    double_vec_t wt
  );
  /* Analogous to {SPTimeSpaceFunction_Dot}, but uses {SPTimeSpaceIntegral_Custom} to compute
    the integral, with the given weights and sample points.
    
    If {FMap == NULL}, then {F(p,t) == f(p,t)}, and similarly for {GMap}. */

#endif
