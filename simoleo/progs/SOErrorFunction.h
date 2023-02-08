/* SOErrorFunction.h -- Special function for Error Map determination. */
/* Last edited on 2007-01-04 00:21:54 by stolfi */

#ifndef SOErrorFunction_H
#define SOErrorFunction_H

#include <SOBasic.h>
#include <SOFunction.h>

#include <dg_grid.h>

#include <vec.h>

/* An {SOErrorFunction} object describes a ???*** {R^m} to {R}???.  */

typedef struct SOErrorFunction_Methods  /* Methods for SOErrorFunction: */
  { SOFunction_Methods fn;  
      /* Methods inherited from superclass, with some overrides. 
        The {add} method only works with other functions of the same
        type. */
    
    WriteMth *write;
      /* Writes the function's def to {wr}, as an {SOErrorFunction}. */
        
  } SOErrorFunction_Methods;

typedef struct SOErrorFunction_Data  /* Data fields for SOErrorFunction: */
  { SOFunction_Data fn;       /* Data fields inherited from superclass. */

    double_vec_t weight;      /* Weights for error determination. */

    Basis testbas;            /* The test basis functions. */
    Basis appbas;             /* The approximation basis functions. */

    char *testbasFile;        /* Name of test basis file (with extension). */
    char *appbasFile;         /* Name of approximation basis file (with extension). */
  } SOErrorFunction_Data;

typedef struct SOErrorFunction  /* Logically a subclass of SOFunction. */
  { char *type;                 /* Type identifier. */
    SOErrorFunction_Methods *m;   /* Function methods. */
    SOErrorFunction_Data *d;      /* Function parameters. */
  } SOErrorFunction;
  
#define SOErrorFunction_TypeId "SOF.Error."

SOErrorFunction *SOErrorFunction_Cast(OBJ *f);
  /* If {f} (according to {f->type}) is a subclass of {SOErrorFunction},
    returns {f} cast to that type; otherwise returns NULL. */

/* ====================================================================== */

SOErrorFunction *SOErrorFunction_Read(FILE *rd, dg_dim_t pDim);
  /* Reads an {SOErrorFunction} from {rd}.
    It must have been created by calling the {write} method
    of an {SOErrorFunction}. */


SOErrorFunction *SOErrorFunction_Make
  ( dg_dim_t pDim, 
    dg_dim_t fDim,
    char *testbasFile,   /* Name of test basis file (with extension). */
    char *appbasFile,    /* Name of approximation basis file (with extension). */
    Basis testbas,
    Basis appbas,
    double_vec_t w
  );

#endif
