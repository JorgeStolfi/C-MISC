/* SPTimeOnlyProcFunction.h -- Procedurally defined spherical functions. */
/* Last edited on 2005-08-18 18:14:32 by stolfi */

#ifndef SPTimeOnlyProcFunction_H
#define SPTimeOnlyProcFunction_H

#include <SPTimeOnlyFunction.h>
#include <js.h>

/* PROCEDURALLY DEFINED FUNCTIONS OF TIME */

/* An {SPTimeOnlyProcFunction} object describes a real-valued function
  of a real variable, defined by a black-box procedure. */

typedef struct SPTimeOnlyProcFunction_Methods  /* Methods for SPTimeOnlyProcFunction: */
  { SPTimeOnlyFunction_Methods fn;  
      /* Methods inherited from superclass, with some overrides. 
        The {add} method is generally undefined for these functions,
        except between two functions with the same final type. 
        The {scale} and {copy} methods generally work. */
      
    SPTimeOnlyFunction_WriteMth *write;       
      /* Writes the function's def to {wr}. */
        
  } SPTimeOnlyProcFunction_Methods;

typedef struct SPTimeOnlyProcFunction_Data  /* Data fields for SPTimeOnlyProcFunction: */
  { SPTimeOnlyFunction_Data fn;   /* Data fields inherited from superclass. */
    char *descr;           /* TimeOnlyFunction's description (e.g. formula). */
    double scale;          /* Scale factor for function values. */
  } SPTimeOnlyProcFunction_Data;

typedef struct SPTimeOnlyProcFunction  /* Logically a subclass of SPTimeOnlyFunction. */
  { char *type;                  /* Type identifier. */
    SPTimeOnlyProcFunction_Methods *m;   /* TimeOnlyFunction methods. */
    SPTimeOnlyProcFunction_Data *d;      /* TimeOnlyFunction parameters. */
  } SPTimeOnlyProcFunction;
  /* The data record {*d} is private to a single {SPTimeOnlyProcFunction} object,
    and is reclaimed by the {free} method. */
  
#define SPTimeOnlyProcFunction_TypeId SPTimeOnlyFunction_TypeId "Proc."

SPTimeOnlyProcFunction *SPTimeOnlyProcFunction_Cast(OBJ *f);
  /* If {f} (according to {f->type}) is a subclass of {SPTimeOnlyProcFunction},
    returns {f} cast to that type; otherwise returns NULL. */

/* ====================================================================== */

SPTimeOnlyProcFunction *SPTimeOnlyProcFunction_FromName(char *name);
  /* Returns a procedural spherical function with given {name},
    which must match one of the predefined names.
    The data record will be brand-new, with the {scale} field 
    set to 1. */

SPTimeOnlyProcFunction *SPTimeOnlyProcFunction_Read(FILE *rd);
  /* Reads a procedural spherical function from {rd}. It must have
    been created by calling the {write} method of an
    {SPTimeOnlyProcFunction}. */

#endif
