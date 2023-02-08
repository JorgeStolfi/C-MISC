/* SOShiftedFunction.h -- a scaled and shifted copy of another function. */
/* Last edited on 2007-01-04 00:23:00 by stolfi */

#ifndef SOShiftedFunction_H
#define SOShiftedFunction_H

#include <SOBasic.h>
#include <SOFunction.h>

#include <vec.h>

/* PROCEDURALLY DEFINED FUNCTIONS */

/* A {SOShiftedFunction} object describes a function {f} from {R^m} to {R^n}, 
  such that {f(p) = g(q)} where {g} is some other {SOFunction} and {q} 
  is the coordinate vector of point {p} relative to a specified box {B}. */

typedef struct SOShiftedFunction_Methods  /* Methods for SOShiftedFunction: */
  { SOFunction_Methods fn;  
      /* Methods inherited from superclass, with some overrides. */
      
    WriteMth *write;       
      /* Writes the function's def to {wr}, as an {SOShiftedFunction}. */
        
  } SOShiftedFunction_Methods;

typedef struct SOShiftedFunction_Data  /* Data fields for SOShiftedFunction: */
  { SOFunction_Data fn;    /* Data fields inherited from superclass. */
    char *descr;              /* Function's description = "Sh(B,g)". */
    SOFunction *g;            /* The shifted function. */
    interval_t B[MAX_PDIM]; /* Reference box */
  } SOShiftedFunction_Data;

typedef struct SOShiftedFunction  /* Logically a subclass of SOFunction. */
  { char *type;                     /* Type identifier. */
    SOShiftedFunction_Methods *m;   /* Function methods. */
    SOShiftedFunction_Data *d;      /* Function parameters. */
  } SOShiftedFunction;
  
#define SOShiftedFunction_TypeId "SOF.Shf."

SOShiftedFunction *SOShiftedFunction_Cast(OBJ *f);
  /* If {f} (according to {f->type}) is a subclass of {SOShiftedFunction},
    returns {f} cast to that type; otherwise returns NULL. */

/* ====================================================================== */

SOShiftedFunction *SOShiftedFunction_Make(SOFunction *g, interval_t B[]);
  /* Returns ashifted and scaled version of {g}, such that 
    {f(p) = g(q)} for every point {p}, where {q} is the 
    coordinate vector of {p} relative to box {B}
    The function will have domain {R^g->d->pDim} and range {R^g->d->fDim}. */

SOShiftedFunction *SOShiftedFunction_Read(FILE *rd);
  /* Reads a procedural function from {rd}.
    It must have been created by calling the {write} method
    of an {SOShiftedFunction}. */

#endif
