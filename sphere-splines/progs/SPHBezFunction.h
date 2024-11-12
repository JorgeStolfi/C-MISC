/* SPHBezFunction.h -- Homogeneous spherical polynomial in Bezier form */
/* Last edited on 2005-10-06 21:47:12 by stolfi */

#ifndef SPHBezFunction_H
#define SPHBezFunction_H

#include <SPFunction.h>
#include <SPDeCasteljau.h>
#include <SPBasic.h>

#include <vec.h>
#include <js.h>
#include <nat.h>

/* HOMOGNEOUS SPHERICAL POLYNOMIALS IN BEZIER FORM */

/* An {SPHBezFunction} object describes a real-valued function
  defined on {S^3} as a homogeneous polynomial in the {x,y,z} 
  coordinates given by its Bezier coefficients. */

typedef struct SPHBezFunction_Methods  /* Methods for SPHBezFunction: */
  { SPFunction_Methods fn;  
      /* Methods inherited from superclass, with some overrides. 
        The {add} method only works with other functions of the same
        type. The {scale} and {copy} methods work.*/
    
    SPFunction_WriteMth *write;
      /* Writes the function's def to {wr}, as an {SPHBezFunction}. */
        
  } SPHBezFunction_Methods;

typedef struct SPHBezFunction_Data  /* Data fields for SPHBezFunction: */
  { SPFunction_Data fn;   /* Data fields inherited from superclass. */
    int deg;              /* Degree of polynomial. */
    BezCoeff_vec_t c;     /* Bezier coeficients. */
  } SPHBezFunction_Data;

typedef struct SPHBezFunction  /* Logically a subclass of SPFunction. */
  { char *type;                  /* Type identifier. */
    SPHBezFunction_Methods *m;   /* Function methods. */
    SPHBezFunction_Data *d;      /* Function parameters. */
  } SPHBezFunction;
  /* The data record {*d} and the coefficient vector {*(d->c)}
    are private to a single {SPHBezFunction} object, and are reclaimed
    by the {free} method. */
  
#define SPHBezFunction_TypeId "SF.HBez."

SPHBezFunction *SPHBezFunction_Cast(OBJ *f);
  /* If {f} (according to {f->type}) is a subclass of {SPHBezFunction},
    returns {f} cast to that type; otherwise returns NULL. */

/* ====================================================================== */

nat_t SPHBezFunction_NumCoeffs(int deg);
  /* Number of Bezier coefficients in a homogeneous spherical
    polynomial of degree {deg}. */

SPHBezFunction *SPHBezFunction_FromCoeffs(int deg, BezCoeff_vec_t c);
  /* Creates an {SPHBezFunction} object with given degree and coefficient
    vector. The latter will become private to the new {SPHBezFunction}, 
    and must not be used for other purposes. */

SPHBezFunction *SPHBezFunction_Read(FILE *rd);
  /* Reads an {SPHBezFunction} from {rd} (mainly, its Bezier coeffs).
    It must have been created by calling the {write} method
    of an {SPHBezFunction}. */

#endif
 
