/* SPTimeOnlyFunction.h -- General functions of time */
/* Last edited on 2009-02-10 09:00:05 by stolfi */

#ifndef SPTimeOnlyFunction_H
#define SPTimeOnlyFunction_H

#include <SPTimeOnlyFuncMap.h>
#include <SPBasic.h>
#include <SPMatrix.h>

#include <udg_pulse.h>
#include <vec.h>
#include <bool.h>
#include <nat.h>
#include <stdio.h>
 
/* An {SPTimeOnlyFunction} object represents a a real
  function of a real variable (called /time/ here). */
  
#define OBJ void /* Actually {SPTimeOnlyFunction} or a subclass */

typedef double SPTimeOnlyFunction_EvalMth(OBJ *f, double t);     
typedef void SPTimeOnlyFunction_WriteMth(OBJ *f, FILE *wr);        
typedef void SPTimeOnlyFunction_FreeMth(OBJ *f);        

typedef struct SPTimeOnlyFunction_Methods    /* Methods for any SPTimeOnlyFunction: */
  { SPTimeOnlyFunction_EvalMth *eval;    /* Evaluates the function at {p}. */
    SPTimeOnlyFunction_WriteMth *write;  /* Writes the function to {wr}, for reading back. */
    SPTimeOnlyFunction_FreeMth *free;    /* Reclaims the function's private storage. */
  } SPTimeOnlyFunction_Methods;
  /* Some of these methods may be missing for certain subclasses. */
  
typedef struct SPTimeOnlyFunction_Data    /* Data fields for any SPTimeOnlyFunction: */
  { } SPTimeOnlyFunction_Data;            
  
typedef struct SPTimeOnlyFunction /* Abstract class: a function defined on the sphere. */
  { char *type;              /* Type identifier */
    SPTimeOnlyFunction_Methods *m;  /* TimeOnlyFunction methods. */
    SPTimeOnlyFunction_Data *d;     /* TimeOnlyFunction parameters. */
  } SPTimeOnlyFunction;
  
#define SPTimeOnlyFunction_TypeId "TOF."

SPTimeOnlyFunction *SPTimeOnlyFunction_Cast(OBJ *f);
  /* If {f} (according to {f->type}) is a subclass of {SPTimeOnlyFunction},
    returns {f} cast to that type; otherwise returns NULL. */

void SPTimeOnlyFunction_M_Write(SPTimeOnlyFunction *f, FILE *wr); 
  /* The default {write} method for a generic sphere-time function.
    Writes a standard {SPTimeOnlyFunction} header line and the specific
    subclass ID, then calls the subclass-specific {write} method, and
    closes off with the generic {SPTimeOnlyFunction} footer. */

SPTimeOnlyFunction *SPTimeOnlyFunction_Read(FILE *rd);  
  /* Reads a function's description from {rd}. */

/* BASES FOR TIME FUNCTION SPACES */

/* For time functions, we only consider finite-element bases 
  of the form {tau[j,u](t)} where the index {j} identifies an /epoch/
  (distinguished time) and {u} is an additional /shape index/. */

typedef struct SPTimeOnlyFunction_Basis /* A tensor-product basis for sphere-time fns. */
  { udg_pulse_family_t tFam;  /* The time basis. */
  } SPTimeOnlyFunction_Basis;

typedef SPTimeOnlyFunction_Basis TBasis;

void SPTimeOnlyFunction_WriteBasis(FILE *wr, TBasis tbas);
  /* Writes a basis {tbas} to {wr}, in a format that can be 
    read back by {Read}.  */

TBasis SPTimeOnlyFunction_ReadBasis(FILE *rd);
  /* Reads from {rd} a basis for a space of sphere-time functions.
    Assumes the contents of {rd} was created by {write} above. */

/* INTEGRALS OF SPHERE-TIME FUNCTIONS */

double SPTimeOnlyFunction_Integral
  ( SPTimeOnlyFunction *f,
    TimeOnlyFuncMap FMap, 
    SPTimeOnlyFunction *w,
    int smpOrder,
    double tMin, 
    double tMax, 
    bool_t verbose
  );
  /* Computes the integral of {w(t)*FMap(f(t),t)} for {t}
    ranging over the interval {[tMin _ tMax]}, using Gaussian
    quadrature of order {smpOrder}. */

double SPTimeOnlyFunction_CustomIntegral
  ( SPTimeOnlyFunction *f, 
    TimeOnlyFuncMap FMap, 
    SPTimeOnlyFunction *w,
    double_vec_t st,
    double_vec_t wt
  );
  /* Computes the sample-based integral of {WF}, namely {WF(t) =
    w(t)*FMap(f(t),t)}, defined as {SUM(wt[j]*WF(st[j]), j=0..N-1)},
    where {st} is a set of sample times, and {wt} are the
    corresponding weights (interval widths). */
    
/* INNER PRODUCT OF SPHERE-TIME FUNCTIONS */

/* These procedures compute weighted dot product of spherical
  functions, defined as integrals of {w(t)*H(f(t),g(t))} for {t}
  ranging over some interval {[tMin _ tMax]}, where {w} is a weight
  function, and {H()} is some pointwise operator. The conventions are
  similar to those of the corresponding integration procedures above.
  In particular, all these procedures assume {w(t) == 1.0}
  for all {t}, when {w == NULL}.*/

double SPTimeOnlyFunction_Dot
  ( SPTimeOnlyFunction *f,
    TimeOnlyFuncMap FMap, 
    SPTimeOnlyFunction *g,
    TimeOnlyFuncMap GMap,
    SPTimeOnlyFunction *w,
    int smpOrder,
    double tMin, 
    double tMax, 
    bool_t verbose
  );
  /* Computes the weighted dot product of {F(t)=FMap(f(t),t)} and
    {G(t)=GMap(g(t),t)}, that is the integral of {w(t)*F(t)*G(t)} 
    over the interval {[tMin _ tMax]}, using Gaussian quadrature
    of order {smpOrder}. (If {FMap == NULL}, then {F(t)} reduces to
    {f(t)}, and similarly for {GMap}.) */
    
double SPTimeOnlyFunction_CustomDot
  ( SPTimeOnlyFunction *f, 
    TimeOnlyFuncMap FMap, 
    SPTimeOnlyFunction *g,
    TimeOnlyFuncMap GMap,
    SPTimeOnlyFunction *w,
    double_vec_t st,
    double_vec_t wt
  );
  /* Analogous to {SPTimeOnlyFunction_Dot}, but uses {SPTimeOnlyIntegral_Custom} 
    to compute the integral, with the given weights and sample times.
    If {FMap == NULL}, then {F(t) == f(t)}, and similarly for {GMap}. */

#endif
