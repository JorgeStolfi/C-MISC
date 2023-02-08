/* SOFunction.h -- General functions on an adaptive dyadic grid */
/* Last edited on 2007-01-04 00:22:13 by stolfi */

#ifndef SOFunction_H
#define SOFunction_H

#include <SOGrid.h>
#include <SOBasic.h>

#include <dg_tree.h>
#include <dg_grid.h>

#include <vec.h>
#include <stdio.h>
 
/* An {SOFunction} object represents some function from some {R^m} to some {R^n}. */

#define OBJ void /* Actually {SOFunction} or a subclass. */

typedef void EvalMth(OBJ *f, double *p, double *fp);  /* Function value at {p}. */
typedef void GradMth(OBJ *f, double *p, double *dp);  /* Gradient at {p}. */
typedef void HessMth(OBJ *f, double *p, double *ddp); /* Hessian at {p}. */
typedef void WriteMth(OBJ *f, FILE *wr);       /* Writes function def to {wr}. */
typedef OBJ *CopyMth(OBJ *f);                  /* Clones {f} (new data). */

typedef struct SOFunction_Methods  /* Methods for any SOFunction: */
  { EvalMth *eval;     /* Stores in {fp[0..n-1]} the function value at {p}. */
    GradMth *grad;     /* Stores in {dp[0..m*n-1]} the gradient at {p}. */
    HessMth *hess;     /* Stores in {ddp[0..m*(m+1)/2*n-1]} the Hessian at {p}. */
    WriteMth *write;   /* Writes the function description to {wr}. */
    CopyMth *copy;     /* Makes a clone of {f} (same meths, new data). */
  } SOFunction_Methods;
  /* Some of these methods may be missing for certain subclasses. */
  
typedef struct SOFunction_Data    /* Data fields for any SOFunction: */
  { nat_t pDim;    /* Dimension of the argument points (domain). */
    nat_t fDim;    /* Dimension of the function values (range). */
  } SOFunction_Data;            
  
typedef struct SOFunction /* Abstract class: a function from {R^m} to {R^n}. */
  { char *type;             /* Type identifier. */
    SOFunction_Methods *m;  /* Function methods (shared through class). */
    SOFunction_Data *d;     /* Function attributes (exclusive of this fn). */
  } SOFunction;
  
typedef SOFunction* SOFunctionRef;

#define SOFunction_TypeId "SOF."

SOFunction *SOFunction_Cast(OBJ *f);
  /* If {f} (according to {f->type}) is a subclass of {SOFunction},
    returns {f} cast to that type; otherwise returns NULL. */

void SOFunction_WriteMth(SOFunction *f, FILE *wr); 
  /* The default {write} method for a generic function. Writes a
    standard {SOFunction} header line and the specific subclass ID,
    then calls the subclass-specific {write} method, and closes off
    with the generic {SOFunction} footer. */

SOFunction *SOFunction_Read(FILE *rd);  
  /* Reads a function's description from {rd}, which must have been
    created by {SOFunction_WriteMth}. */

/* BASES FOR FUNCTION SPACES */

vec_typedef(SOFunctionRef_vec_t,SOFunctionRef_vec,SOFunctionRef);

#define Basis SOFunctionRef_vec_t 
  /* A {Basis} is a list of {SOFunction}s, hopefully linearly independent. */

void SOFunction_WriteBasis(FILE *wr, Basis F);
  /* Writes a basis {F} to {wr}, in a format that can be 
    read back by {Read}.  */

Basis SOFunction_ReadBasis(FILE *rd);
  /* Reads from {rd} a basis for a space of {SOFunction}s. 
    Assumes that the contents of {rd} was created by {write} above. */

Basis SOFunction_ReadBasisCached(char *fileName);

/* INTEGRALS OF FUNCTIONS */

/* These procedures compute weighted and modified integrals of a
  function {f}, defined as integrals of {w(p)*FMap(f(p),p)}
  over some box {b}, where {w} is a weight function, and {FMap()} is some
  local operator (a {FuncMap}). 
  
  The default box {b} is the root (canonical) cell of the dyadic grid.
  All procedures assume that {w(p) = 1.0} for all {p}, when {w = NULL}.
  
  In all procedures, the integrand {f} may have arbitrary
  domain and range dimensions.  The weight function {w} must
  have the same domain dimension as {f}, and it must return
  a scalar result (i.e. {w->d.fDim = 1}). */

typedef void FuncMapProc(double *u, double *p, double *v); 
  /* A procedure that applies a local transformation to 
    the value {u = f(p)} of some function {f}.  The transforma
    tion may depend on the position {p}.  The transformed result
    is returned in {v}. */

typedef struct FuncMap    /* Describes a local mapping of function values. */
  { FuncMapProc *map;  /* The map proper (NULL means the identity map). */
    dg_dim_t uDim;   /* Dimension of the function's value {u = f(p)}. */
    dg_dim_t pDim;   /* Dimension of the argument points {p}. */
    dg_dim_t vDim;   /* Dimension of the {FuncMap}'s result {v = FMap(f(p),p)}. */
    bool_t zeroPres;     /* TRUE implies {map(0,p) = 0} for any {p}. */
  } FuncMap;
  
FuncMap IdentityMap (dg_dim_t uDim, dg_dim_t pDim);
  /* Returns an identity {FuncMap} {F}, such that {F(u,p) = u}
   for any {u} in {R^uDim} and any {p} in {R^pDim}. */

typedef void SOFunction_fgDot(SOFunction *f, SOFunction *g, double *iv);
  /* A generic procedure that returns the dot product of a pair of 
     SOFunctions. iv = <f | g>, iv with the same dimension as f and g. */

void SOFunction_Integral
  ( SOFunction *f,
    FuncMap *FMap, 
    SOFunction *w,
    SOGrid_Tree *tree,
    bool_t verbose,
    double *rv
  );
  /* Computes the integral of {w(p)*FMap(f(p),p)} over the 
    root cell.  Returns the result in {rv}.
    
    If the function {f} happens to be an instance of {SOSpline},
    and {tree} is either NULL or {f}'s underlying grid, then
    the procedure reduces to {SOSpline_IntegralPW} (q.v.).
    Otherwise, the integral is computed with {SOIntegral_OnRootCell}
    on the root cell, with the same {tree} argument. */

void SOFunction_RawIntegral
  ( SOFunction *f,
    FuncMap *FMap, 
    SOFunction *w,
    SOGrid_Tree *tree,
    bool_t verbose,
    double *rv
  );
  /* Like {SOFunction_Integral}, but always computes the integral with
    {SOIntegral_OnRootCell}, with the given {tree} argument, ignoring any
    underlying grid of {f}. */
    
/* INNER PRODUCT OF FUNCTIONS */

/* These procedures compute weighted dot product of functions, defined
  as integrals of {w(p)*H(f(p),g(p))} over the root cell, where {w} is
  a weight function (a scalar field), and {H()} is some local
  operator. The conventions are similar to those of the corresponding
  integration procedures above. In particular, all of these procedures
  assume {w(p) == 1.0} for all {p}, when {w == NULL}.
  
  The mapped function results {FMap(f(p),p)} and {GMap(g(p),p)}
  must have the same dimension. */

void SOFunction_Dot
  ( SOFunction *f,
    FuncMap *FMap, 
    SOFunction *g,
    FuncMap *GMap,
    SOFunction *w,
    SOGrid_Tree *tree,
    double *iv,
    bool_t verbose
  );
  /* Computes the weighted dot product of {F(p)=FMap(f(p),p)} and
    {G(p)=GMap(g(p),p)}, that is the integral of {w(p)*vdot(F(p),G(p))} over
    the whole root cell, where {vdot} is the ordinary dot product.
    
    If both functions {f,g} happen to be instances of {SOSpline}
    with the same grid, and {tree} is either NULL or that same
    grid, then the procedure reduces to
    {SOSpline_DotBothPW} (q.v.).
    
    Else, if either {f} or {g} happen to be an instance of {SOSpline}
    and {tree} is either NULL or that element's underlying grid,
    then the procedure reduces to {SOSpline_DotSinglePW} (q.v.).
    
    In all other cases, the integral is computed with {SOIntegral_OnRootCell},
    with same {tree} argument. */

void SOFunction_GradDot
  ( SOFunction *f,
    SOFunction *g, 
    SOFunction *w,
    SOGrid_Tree *tree,
    double *iv,
    bool_t verbose
  );
  /* Computes the dot product of the {Grad(f)} and {Grad(g)}, that is
    the integral of {w(p)*vdot(Grad(f)(p),Grad(g)(p))} over the root
    cell, where {Grad} is the gradient operator.
    
    The integral is approximated as described for {SOFunction_Dot}. */

double SOFunction_RawDot
  ( SOFunction *f,
    FuncMap *FMap, 
    SOFunction *g,
    FuncMap *GMap,
    SOFunction *w,
    SOGrid_Tree *tree,
    bool_t verbose
  );
  /* Like {SOFunction_Dot}, but always computes the integral with
    {SOIntegral_OnRootCell}, with the given {tree} argument,
    ignoring any underlying grid of {f} and {g}. */
    
double SOFunction_RawGradDot
  ( SOFunction *f,
    SOFunction *g, 
    SOFunction *w,
    SOGrid_Tree *tree,
    bool_t verbose
  );
  /* Like {SOFunction_GradDot}, but always computes the integral with
    {SOIntegral_OnRootCell}, with the given {tree} argument,
    ignoring any underlying grid of {f} and {g}. */
  
/* {Basis} is a nicer name for an array of {SOFunction*}: */

#define Basis_new SOFunctionRef_vec_new 
#define Basis_trim SOFunctionRef_vec_trim 
#define Basis_expand SOFunctionRef_vec_expand 

#endif
