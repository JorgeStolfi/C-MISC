/* SPTimeSpaceProcFunction.h -- Procedurally defined spherical functions. */
/* Last edited on 2005-08-14 16:43:10 by stolfi */

#ifndef SPTimeSpaceProcFunction_H
#define SPTimeSpaceProcFunction_H

#include <SPTimeSpaceFunction.h>
#include <r3.h>
#include <js.h>

/* PROCEDURALLY DEFINED SPHERICAL FUNCTIONS */

/* An {SPTimeSpaceProcFunction} object describes a real-valued function
  defined on {S^2} (actually, {R^3}) by a black-box procedure. */

typedef struct SPTimeSpaceProcFunction_Methods  /* Methods for SPTimeSpaceProcFunction: */
  { SPTimeSpaceFunction_Methods fn;  
      /* Methods inherited from superclass, with some overrides. 
        The {add} method is generally undefined for these functions,
        except between two functions with the same final type. 
        The {scale} and {copy} methods generally work. */
      
    SPTimeSpaceFunction_WriteMth *write;       
      /* Writes the function's def to {wr}. */
        
  } SPTimeSpaceProcFunction_Methods;

typedef struct SPTimeSpaceProcFunction_Data  /* Data fields for SPTimeSpaceProcFunction: */
  { SPTimeSpaceFunction_Data fn;   /* Data fields inherited from superclass. */
    char *descr;           /* TimeSpaceFunction's description (e.g. formula). */
    double scale;          /* Scale factor for function values. */
  } SPTimeSpaceProcFunction_Data;

typedef struct SPTimeSpaceProcFunction  /* Logically a subclass of SPTimeSpaceFunction. */
  { char *type;                  /* Type identifier. */
    SPTimeSpaceProcFunction_Methods *m;   /* TimeSpaceFunction methods. */
    SPTimeSpaceProcFunction_Data *d;      /* TimeSpaceFunction parameters. */
  } SPTimeSpaceProcFunction;
  /* The data record {*d} is private to a single {SPTimeSpaceProcFunction} object,
    and is reclaimed by the {free} method. */
  
#define SPTimeSpaceProcFunction_TypeId "STF.Proc."

SPTimeSpaceProcFunction *SPTimeSpaceProcFunction_Cast(OBJ *f);
  /* If {f} (according to {f->type}) is a subclass of {SPTimeSpaceProcFunction},
    returns {f} cast to that type; otherwise returns NULL. */

/* ====================================================================== */

SPTimeSpaceProcFunction *SPTimeSpaceProcFunction_FromName(char *name);
  /* Returns a procedural spherical function with given {name},
    which must match one of the predefined names.
    The data record will be brand-new, with the {scale} field 
    set to 1. */

SPTimeSpaceProcFunction *SPTimeSpaceProcFunction_Read(FILE *rd);
  /* Reads a procedural spherical function from {rd}. It must have
    been created by calling the {write} method of an
    {SPTimeSpaceProcFunction}. */

#endif
