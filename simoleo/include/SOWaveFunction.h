/* SOWaveFunction.h -- Functions defined in terms of waves */
/* Last edited on 2007-01-04 00:21:06 by stolfi */

#ifndef SOWaveFunction_H
#define SOWaveFunction_H

#include <SOBasic.h>
#include <SOFunction.h>
#include <vec.h>

/* FUNCTIONS DEFINED IN TERMS OF WAVES */

/* An {SOWaveFunction} object describes a Hartley wave with 
  integer frequency vector {k[0..m-1]} in the unit cube, i.e.
  a real-valued function {w_k} from {R^m} to {R} such that
    
    {w_k(x) = sin( PI/4 + 2 * PI * SUM { x[i]*k[i] : i=0..m } ) }
  
  */

typedef struct SOWaveFunction_Methods  /* Methods for SOWaveFunction: */
  { SOFunction_Methods fn;  
      /* Methods inherited from superclass, with some overrides.
        The {copy} method will also copy the frequency vector
        {f->d->freq}. */
      
    WriteMth *write;       
      /* Writes the function's def to {wr}, as an {SOWaveFunction}. */
        
  } SOWaveFunction_Methods;

typedef struct SOWaveFunction_Data  /* Data fields for SOWaveFunction: */
  { SOFunction_Data fn;   /* Data fields inherited from superclass. */
    int freq[MAX_PDIM];   /* Frequency vector. */
  } SOWaveFunction_Data;

typedef struct SOWaveFunction  /* Logically a subclass of SOFunction. */
  { char *type;                  /* Type identifier. */
    SOWaveFunction_Methods *m;   /* Function methods. */
    SOWaveFunction_Data *d;      /* Function parameters. */
  } SOWaveFunction;
  
#define SOWaveFunction_SubTypeId "Wave."
#define SOWaveFunction_TypeId (SOFunction_TypeId SOWaveFunction_SubTypeId)

SOWaveFunction *SOWaveFunction_Cast(OBJ *f);
  /* If {f} (according to {f->type}) is a subclass of {SOWaveFunction},
    returns {f} cast to that type; otherwise returns NULL. */

/* ====================================================================== */

SOWaveFunction *SOWaveFunction_FromFreq(nat_t pDim, int *freq);
  /* Returns a wave function with frequency vector 
    {k[0..pDim-1] = freq[0..pDim-1]}. The function will have
    domain {R^pDim} and range {R} (i.e. {f->d->fDim = 1}). */

SOWaveFunction *SOWaveFunction_Read(FILE *rd, nat_t pDim);
  /* Reads a {SOWaveFunction} from {rd}. It must have been created by
    calling the {write} method of an {SOWaveFunction}. The domain
    dimension {pDim} must be provided by the caller, and must 
    be that of the original function that was written to the file. */

#endif
