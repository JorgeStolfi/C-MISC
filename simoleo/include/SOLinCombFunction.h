/* SOLinCombFunction.h -- Linear combinations of other functions. */
/* Last edited on 2007-01-04 00:20:29 by stolfi */

#ifndef SOLinCombFunction_H
#define SOLinCombFunction_H

#include <SOBasic.h>
#include <SOFunction.h>

#include <dg_grid.h>

#include <vec.h>

/* An {SOLinCombFunction} object describes a linear combination of
  functions from {R^m} to {R}, with coefficients in {R^n}, resulting
  in a single {SOFunction} object with from {R^m} to {R^n}.  */

typedef struct SOLinCombFunction_Methods  /* Methods for SOLinCombFunction: */
  { SOFunction_Methods fn;  
      /* Methods inherited from superclass, with some overrides. 
        The {add} method only works with other functions of the same
        type. The {scale} and {copy} methods work. */
    
    WriteMth *write;
      /* Writes the function's def to {wr}, as an {SOLinCombFunction}. */
        
  } SOLinCombFunction_Methods;

typedef struct SOLinCombFunction_Data  /* Data fields for SOLinCombFunction: */
  { SOFunction_Data fn;   /* Data fields inherited from superclass. */
    double_vec_t coef;    /* Coefficients of the linear combination. */
    Basis bas;            /* The basis functions. */
    char *basFile;        /* Name of basis file (with extension). */
  } SOLinCombFunction_Data;
    /* In an {SOLinCombFunction} {f}, all elements of the basis {bas}
      must be functions from the same domain {R^d} to {R}, where {d =
      f->pDim}. The {coef} vector should contain {n*bas.ne} elements,
      where {n = f->fDim}. Coefficient {coef[i*n + j]} is the
      contribution of basis element {bas[i]} to component {f[j]} of
      the function. */

typedef struct SOLinCombFunction  /* Logically a subclass of SOFunction. */
  { char *type;                 /* Type identifier. */
    SOLinCombFunction_Methods *m;   /* Function methods. */
    SOLinCombFunction_Data *d;      /* Function parameters. */
  } SOLinCombFunction;
  
#define SOLinCombFunction_TypeId "SOF.LC."

SOLinCombFunction *SOLinCombFunction_Cast(OBJ *f);
  /* If {f} (according to {f->type}) is a subclass of {SOLinCombFunction},
    returns {f} cast to that type; otherwise returns NULL. */

/* ====================================================================== */

SOLinCombFunction *SOLinCombFunction_Read(FILE *rd, dg_dim_t pDim);
  /* Reads an {SOLinCombFunction} from {rd}.
    It must have been created by calling the {write} method
    of an {SOLinCombFunction}. */

SOLinCombFunction *SOLinCombFunction_Make
  ( dg_dim_t pDim, 
    dg_dim_t fDim,
    char *basFile, 
    Basis bas,
    double_vec_t a
  );
  /* Returns the linear combination {SUM{a[i]*bas[i],
    i=0..bas.ne-1}}, as a new {SOLinCombFunction}. All the functions
    {bas[i]} must have domain dimension {pDim} and range dimension
    {fDim}. The {basFile} must be the name (with extension) 
    of a file where {bas} was read from. */

#endif
 
