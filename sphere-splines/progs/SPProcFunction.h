/* SPProcFunction.h -- Procedurally defined spherical functions. */
/* Last edited on 2005-08-14 16:43:22 by stolfi */

#ifndef SPProcFunction_H
#define SPProcFunction_H

#include <SPFunction.h>
#include <r3.h>
#include <js.h>

/* PROCEDURALLY DEFINED SPHERICAL FUNCTIONS */

/* An {SPProcFunction} object describes a real-valued function
  defined on {S^2} (actually, {R^3}) by a black-box procedure. */

typedef struct SPProcFunction_Methods  /* Methods for SPProcFunction: */
  { SPFunction_Methods fn;  
      /* Methods inherited from superclass, with some overrides. 
        The {add} method is generally undefined for these functions,
        except between two functions with the same final type. 
        The {scale} and {copy} methods generally work. */
      
    SPFunction_WriteMth *write;       
      /* Writes the function's def to {wr}. */
        
  } SPProcFunction_Methods;

typedef struct SPProcFunction_Data  /* Data fields for SPProcFunction: */
  { SPFunction_Data fn;   /* Data fields inherited from superclass. */
    char *descr;          /* Function's description (e.g. formula). */
    double scale;         /* Scale factor for function values. */
  } SPProcFunction_Data;

typedef struct SPProcFunction  /* Logically a subclass of SPFunction. */
  { char *type;                  /* Type identifier. */
    SPProcFunction_Methods *m;   /* Function methods. */
    SPProcFunction_Data *d;      /* Function parameters. */
  } SPProcFunction;
  /* The data record {*d} is private to a single {SPProcFunction} object,
    and is reclaimed by the {free} method. */
  
#define SPProcFunction_TypeId "SF.Proc."

SPProcFunction *SPProcFunction_Cast(OBJ *f);
  /* If {f} (according to {f->type}) is a subclass of {SPProcFunction},
    returns {f} cast to that type; otherwise returns NULL. */

/* ====================================================================== */

SPProcFunction *SPProcFunction_FromName(char *name);
  /* Returns a procedural spherical function with given {name},
    which must match one of the predefined names. */

SPProcFunction *SPProcFunction_Read(FILE *rd);
  /* Reads an procedural spherical function from {rd}.
    It must have been created by calling the {write} method
    of an {SPProcFunction}. */

#endif
