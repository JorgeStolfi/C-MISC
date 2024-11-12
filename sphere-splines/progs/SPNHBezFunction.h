/* SPNHBezFunction.h -- Non-homogeneous spherical polynomial in Bezier form */
/* Last edited on 2004-08-18 00:06:35 by stolfi */

#ifndef SPNHBezFunction_H
#define SPNHBezFunction_H

#include <SPDeCasteljau.h>
#include <SPHBezFunction.h>
#include <SPBasic.h>
#include <vec.h>
#include <js.h>

/* NON-HOMOGNEOUS SPHERICAL POLYNOMIALS IN BEZIER FORM */

/* An {SPNHBezFunction} object describes a real-valued function {f}
  defined on {S^3} as a non-homogeneous polynomial of total 
  degree {<= d} in the {x,y,z} coordinates. 
  
  As shown by Gomide and Stolfi, such polynomial is equivalent on the
  sphere to the direct sum of homogeneous polynomials of degree {d}
  and {d-1}. The representation of {f} consists by the Bezier
  coefficients of these two components. */

typedef struct SPNHBezFunction_Methods  /* Methods for SPNHBezFunction: */
  { SPFunction_Methods fn;  
      /* Methods inherited from superclass, with some overrides. 
        The {add} method only works with other functions of the same
        type. The {scale} and {copy} methods work. */
    
    SPFunction_WriteMth *write;
      /* Writes the function's def to {wr}, as an {SPNHBezFunction}. */
        
  } SPNHBezFunction_Methods;

typedef struct SPNHBezFunction_Data  /* Data fields for SPNHBezFunction: */
  { SPFunction_Data fn;     /* Data fields inherited from superclass. */
    int deg;                /* Degree of polynomial. */
    BezCoeff_vec_t c0;      /* Bezier coeficients of the {deg} component. */
    BezCoeff_vec_t c1;      /* Bezier coeficients of the {deg-1} component. */
  } SPNHBezFunction_Data;

typedef struct SPNHBezFunction  /* Logically a subclass of SPFunction. */
  { char *type;                   /* Type identifier. */
    SPNHBezFunction_Methods *m;   /* Function methods. */
    SPNHBezFunction_Data *d;      /* Function parameters. */
  } SPNHBezFunction;
  /* The data record {*d} and the coefficient vectors {*(d->c0),*(d->c1)}
    are private to a single {SPNHBezFunction} object, and are reclaimed
    by the {free} method. */
  
#define SPNHBezFunction_TypeId "SF.NHBez."

SPNHBezFunction *SPNHBezFunction_Cast(OBJ *f);
  /* If {f} (according to {f->type}) is a subclass of {SPNHBezFunction},
    returns {f} cast to that type; otherwise returns NULL. */

/* ====================================================================== */

SPNHBezFunction *SPNHBezFunction_FromCoeffs
  ( int deg, 
    BezCoeff_vec_t c0,
    BezCoeff_vec_t c1
  );
  /* Creates an {SPNHBezFunction} with given degree and coefficient
    vectors. The latter will become private to the new
    {SPNHBezFunction}, and must not be used for other purposes. */

SPNHBezFunction *SPNHBezFunction_Read(FILE *rd);
  /* Reads an {SPNHBezFunction} from {rd} (mainly, the Bezier coeffs of
    its {deg} and {deg-1} homogeneous components). It must have been
    created by calling the {write} method of an {SPNHBezFunction}. */

SPNHBezFunction *SPNHBezFunction_FromHomo(SPHBezFunction *f, int deg);
  /* Converts a homogeneous spherical polynomial {f},
    in Bezier form, to a non-homogenous spherical polynomial,
    also in Bezier form.  The degree of {f} must be {deg} or {deg-1}.
    The coefficient vector of {f} is copied, not shared. */

#endif
