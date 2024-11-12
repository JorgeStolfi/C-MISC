/* SPHarmonic.h -- Spherical harmonic functions. */
/* Last edited on 2009-01-09 22:04:11 by stolfi */

/* SPHERICAL ERROR MAP FUNCTION OBJECT */

/* A {SPHarmonic} is a spherical harmonic function of given 
  order and degree, packaged as a {SPFunction} object. */

#ifndef SPHarmonic_H
#define SPHarmonic_H

#include <SPBasic.h> 
#include <SPFunction.h> 

#include <vec.h>
#include <nat.h>
#include <js.h>

#include <stdio.h>

#define OBJ void /* Actualy {SPHarmonic} or a subclass */

typedef struct SPHarmonic_Methods  /* Methods for {SPHarmonic}: */
  { SPFunction_Methods fn;  
      /* Methods inherited from superclass, with some overrides. The
        {eval}, {grad}, {hess}, {scale}, and {copy} methods are
        defined. The {add} and {maple} methods are undefined. */
    
    SPFunction_WriteMth *write;
      /* Writes the function's def to {wr}, as an {SPHarmonic}. */
        
  } SPHarmonic_Methods;

typedef struct HarmonicTerm 
  { nat_t degree;  /* Degree of harmonic. */
    int order;     /* In the range [-degree..+degree]. */;
    double coeff;  /* A scale factor. */
  } HarmonicTerm;
  /* A spherical harmonic {f} with given order and degree, scaled by
    {coeff}.  The term is always real-valued, with a quarter-wave
    shift (i.e. the the value at {lat = lon = 0} is {1/sqrt(2)} times
    the maximum amplitude).  The function is normalized so that
    {SPFunction_Dot(f,f) = coeff^ 2}. */

vec_typedef(HarmonicTerm_vec_t,HarmonicTerm_vec,HarmonicTerm);

typedef struct SPHarmonic_Data  /* Data fields for {SPHarmonic}: */
  { SPFunction_Data fn;         /* Data fields inherited from superclass. */
    HarmonicTerm_vec_t t;       /* Harmonic terms, in lex order. */
  } SPHarmonic_Data;

typedef struct SPHarmonic  /* Logically a subclass of SPFunction. */
  { char *type;             /* Type identifier. */
    SPHarmonic_Methods *m;  /* Function methods. */
    SPHarmonic_Data *d;     /* Function parameters. */
  } SPHarmonic;
  /* The data record {*d} and the sequence of terms {d->t} are 
    private to a single {SPHarmonic} object, and and are reclaimed 
    by the {free} method.  The {add} and {scale} methods may 
    change the size of {d->t}. */

#define SPHarmonic_TypeId "SF.harm."

SPHarmonic *SPHarmonic_Cast(OBJ *f);
  /* If {f} (according to {f->type}) is a subclass of {SPHarmonic},
    returns {f}, cast to that type; otherwise returns NULL. */

/* ====================================================================== */

/* INPUT/OUTPUT */

SPHarmonic *SPHarmonic_Read(FILE *rd);
  /* Reads an {SPHarmonic} from {rd}. It must have been
    created by calling the {write} method of an {SPHarmonic}. */

SPHarmonic *SPHarmonic_FromTerms(HarmonicTerm_vec_t t);
  /* Creates an {SPHarmonic} function object {f} from a sequence of terms.
    The vector {t} becomes private to {f} and should not be reused. */

HarmonicTerm SPHarmonic_MakeTerm(nat_t degree, nat_t order);
  /* Creates an {SPHarmonic} term {f} of the the specified degree and
    order, with unit {coeff}. */

SPHarmonic *SPHarmonic_FromTerm(nat_t degree, nat_t order);
  /* Creates an {SPHarmonic} function object {f} with a single term
    of the given {degree} and {order}, and unit {coeff}.  */

Basis SPHarmonic_MakeBasis(nat_t degree);
  /* Creates a basis consisting of all spherical harmonics 
    with degrees in {0 .. degree}. */

/* HARMONIC TERM EVALUATION */

void SPHarmonic_EvalTerms
  ( HarmonicTerm_vec_t t, 
    R3Point *p, 
    int diff,
    double *v,
    r3_t *dv,
    r6_t *ddv
  );
/* Computes the sum of the terms {t[0..t.ne-1]} evaluated 
  at the point {p}. If {diff>=0}, returns the value in {v};
  if {diff>=1}, returns the gradient in {dv};
  if {diff>=2}, return the hessian in {ddv}. */ 

void SPHarmonic_EvalTerm
  ( HarmonicTerm *t, 
    double clon, 
    double slon, 
    double slat, 
    double R, 
    int diff,
    double *v,
    r3_t *dv,
    r6_t *ddv
  );
  /* Computes the harmonic term {t} at a point {p=(x,y,z)},
    given {clon = x/r}, {slon = y/r}, {slat = z/R}, and {R},
    where {r = sqrt(x^2+y^2)} and {R = sqrt(x^2+y^y+z^2)}.
    
    The terms are defined as real-valued and rotated 
    by {PI/4} around Z. */
    
/* TOOLS FOR SPHERICAL HARMONICS */

void SPHarmonic_Trig(int m, double *x, double *y);
  /* Computes the complex power {(x + I*y)^m}. 
    Returns result in {x, y}. */

double SPHarmonic_Legendre(int d, int m, double z);
  /* Computes the Legendre polynomial {P_d^m(z)}. 
    Requires {abs(z) <= 1}. */
    
double SPHarmonic_LegendreOld(int d, int m, double z);
  /* The old implementation of {SPHarmonic_Legendre},
    from Numerical Recipes. */

#endif
