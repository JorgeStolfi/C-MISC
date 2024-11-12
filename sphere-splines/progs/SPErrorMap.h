/* SPErrorMap.h -- Approximation error map as an {SPFunction} */
/* Last edited on 2003-08-17 18:01:05 by stolfi */
 
#ifndef SPErrorMap_H
#define SPErrorMap_H

#include <SPBasic.h> 
#include <SPFunction.h> 
#include <stdio.h>
#include <vec.h>
#include <js.h>

/* SPHERICAL ERROR MAP FUNCTION OBJECT */

/* An {SPErrorMap} gives the worst-case error at a given point {p}
  when a unit-norm function from a function space {\CF} (the `gauge'
  or `target space') is optimally approximated by some other function
  space {\CA} (the `approximation space').
  
  The space {\CF} must be finite-dimensional, and is given by a basis
  {F[0..N-1]}, which is assumed to be orthonormal in the metric of
  {\CF} (and therefore defines that metric).
  
  The worst-case approximation error at {p} turns out to be the
  root-mean-square of the errors {e[i](p) = F[i](p) - G[i](p)}, where
  {G[i]} is the best approximation (orthogonal projection) of {F[i]}
  in the space {\CA}. See
  
    A. Gomide and J. Stolfi, "Approximation Error Maps", in Proc. of
    the 2001 International Symposium Algorithms for Approximation IV,
    446--453 (University of Huddersfield, UK, 2002).
  
  */

#define OBJ void /* Actualy {SPErrorMap} or a subclass */

typedef struct SPErrorMap_Methods  /* Methods for {SPErrorMap}: */
  { SPFunction_Methods fn;  
      /* Methods inherited from superclass, with some overrides. 
        The {eval} method returns the root-mean-square approximation
        error at the given point. The {grad}, {hess}, and {scale}
        methods are defined. The {copy} method works and makes a
        copy of all functions {f->d->F[i]} and {f->d->G[i]}. 
        The {add} and {maple} methods are undefined. */
    
    SPFunction_WriteMth *write;
      /* Writes the function's def to {wr}, as an {SPErrorMap}. */
        
  } SPErrorMap_Methods;

typedef struct SPErrorMap_Data  /* Data fields for {SPErrorMap}: */
  { SPFunction_Data fn; /* Data fields inherited from superclass. */
    Basis F;            /* An orthonormal basis for the gauge space {\CF}. */
    Basis G;            /* Best approximation of each {F[i]} in {\CA}. */
  } SPErrorMap_Data;

typedef struct SPErrorMap  /* Logically a subclass of SPFunction. */
  { char *type;             /* Type identifier. */
    SPErrorMap_Methods *m;  /* Function methods. */
    SPErrorMap_Data *d;     /* Function parameters. */
  } SPErrorMap;
  /* The data record {*d}, the basis vctors {*(d->F.e),*(d->G.e)},
    and the functions in those bases are private a single {SPErrorMap}
    object, and are all reclaimed by the {free} method --- that also 
    calls the {free} method of each element {*(d->F.e[i]), *(d->G.e[i])} */

#define SPErrorMap_TypeId "SF.emap."

SPErrorMap *SPErrorMap_Cast(OBJ *f);
  /* If {f} (according to {f->type}) is a subclass of {SPErrorMap},
    returns {f} cast to that type; otherwise returns NULL. */

/* ====================================================================== */

SPErrorMap *SPErrorMap_FromBases(Basis F, Basis G);
  /* Creates an {SPErrorMap} given an orthonormal basis {F}
    for the gauge space {\CF}, and its projection {G} onto the
    approximation space {\CA}. Note: the basis vectors {F.e} and
    {G.e}, as well as the functions in them, will become 
    private to the new {ErrorMap} object, and should not
    be used elsewhere. */

SPErrorMap *SPErrorMap_Read(FILE *rd);
  /* Reads an {SPErrorMap} from {rd}. It must have been
    created by calling the {write} method of an {SPErrorMap}. */

#endif
