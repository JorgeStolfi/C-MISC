/* SOTentFunction.h -- Spline functions in a diadic grid domain */
/* Last edited on 2007-01-04 00:21:02 by stolfi */

#ifndef SOTentFunction_H
#define SOTentFunction_H

#include <SOTent.h>
#include <SOBasic.h>
#include <SOFunction.h>

#include <dg_grid.h>

#include <vec.h>

/* An {SOTentFunction} object describes a tent-like multiaffine {C0}
  spline function from {R^m} to {R}, whose support consists of {2^d}
  congruent cells. The function varies from 1 at the ``tent's pole''
  (center of the support) to 0 on the periphery of the support. Within
  each cell, the tent is a polynomial of degree 1 in each
  coordinate. */

typedef struct SOTentFunction_Methods  /* Methods for SOTentFunction: */
  { SOFunction_Methods fn;  
      /* Methods inherited from superclass, with some overrides. */
    
    WriteMth *write;
      /* Writes the function's def to {wr}, as an {SOTentFunction}. */
        
  } SOTentFunction_Methods;

typedef struct SOTentFunction_Data  /* Data fields for SOTentFunction: */
  { SOFunction_Data fn;   /* Data fields inherited from superclass. */
    dg_cell_index_t index;   /* Index of lower-numbered brick in tent's domain. */
    dg_rank_t rank;         /* Rank of bricks. */
  } SOTentFunction_Data;

typedef struct SOTentFunction  /* Logically a subclass of SOFunction. */
  { char *type;                 /* Type identifier. */
    SOTentFunction_Methods *m;   /* Function methods. */
    SOTentFunction_Data *d;      /* Function parameters. */
  } SOTentFunction;
  
#define SOTentFunction_TypeId "SOF.Tent."

SOTentFunction *SOTentFunction_Cast(OBJ *f);
  /* If {f} (according to {f->type}) is a subclass of {SOTentFunction},
    returns {f} cast to that type; otherwise returns NULL. */

/* ====================================================================== */

SOTentFunction *SOTentFunction_Read(FILE *rd, dg_dim_t pDim);
  /* Reads an {SOTentFunction} from {rd}.
    It must have been created by calling the {write} method
    of an {SOTentFunction}. */

SOTentFunction *SOTentFunction_FromHighBrick(dg_dim_t pDim, dg_cell_index_t k);
  /* Creates an {SPTentFunction} for a tent of dimension {d} and rank {r}, 
   whose high brick is cell {k}. */

/* INNER PRODUCTS */

double SOTentFunction_DotBoth(SOTentFunction *ftf, SOTentFunction *gtf);
  /* Returns the dot product of 2 Tent Functions (ftf, gtf), using 
    {SOTent_tt_dot}. */
void SOTentFunction_DotSingle(SOTentFunction *tf, SOFunction *f, FuncMap *FMap, double *sum);
  /* Returns the dot product of a Tent Function (tf) against a general
    (f) function using {SOTent_tf_dot}. */
double SOTentFunction_GradDotBoth(SOTentFunction *ftf, SOTentFunction *gtf);
/* Returns the dot product of the gradient of 2 Tent Functions (ftf, gtf). */

#endif
 
