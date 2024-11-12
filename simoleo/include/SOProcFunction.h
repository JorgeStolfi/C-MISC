/* SOProcFunction.h -- Procedurally defined   functions. */
/* Last edited on 2007-01-04 00:20:54 by stolfi */

#ifndef SOProcFunction_H
#define SOProcFunction_H

#include <SOBasic.h>
#include <SOFunction.h>

#include <vec.h>

/* PROCEDURALLY DEFINED FUNCTIONS */

/* An {SOProcFunction} object describes a real-valued function
  from {R^m} to {R^n}, defined by a C procedure. */

typedef struct SOProcFunction_Methods  /* Methods for SOProcFunction: */
  { SOFunction_Methods fn;  
      /* Methods inherited from superclass, with some overrides. */
      
    WriteMth *write;       
      /* Writes the function's def to {wr}, as an {SOProcFunction}. */
        
  } SOProcFunction_Methods;

typedef struct SOProcFunction_Data  /* Data fields for SOProcFunction: */
  { SOFunction_Data fn;   /* Data fields inherited from superclass. */
    char *descr;          /* Function's description (e.g. formula). */
  } SOProcFunction_Data;

typedef struct SOProcFunction  /* Logically a subclass of SOFunction. */
  { char *type;                  /* Type identifier. */
    SOProcFunction_Methods *m;   /* Function methods. */
    SOProcFunction_Data *d;      /* Function parameters. */
  } SOProcFunction;
  
#define SOProcFunction_TypeId "SOF.Proc."

SOProcFunction *SOProcFunction_Cast(OBJ *f);
  /* If {f} (according to {f->type}) is a subclass of {SOProcFunction},
    returns {f} cast to that type; otherwise returns NULL. */

/* ====================================================================== */

SOProcFunction *SOProcFunction_FromName(char *name, nat_t pDim, nat_t fDim);
  /* Returns a procedural function with given {name},
    which must match one of the predefined names.
    The function will have domain {R^pDim} and range {R^fDim}. */

SOProcFunction *SOProcFunction_Read(FILE *rd, nat_t pDim, nat_t fDim);
  /* Reads a procedural function from {rd}.
    It must have been created by calling the {write} method
    of an {SOProcFunction}.  The domain dimension {pDim} and the 
    range dimension {fDim} must be provided by the caller. */

#endif
