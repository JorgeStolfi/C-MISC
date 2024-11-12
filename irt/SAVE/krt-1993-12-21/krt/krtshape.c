#include "krtshape.h"
#include <r3.h>

#define FUDGE_FACTOR 0.99999

double intsolid ( r3_t pos, r3_t ray, t_solid solid )
{
  return ( (*(solid->shape->intp)) ( pos, ray, solid.parms ) );
}

void nrmsolid ( r3_t hit, t_solid solid, r3_t nrm )
{
  (*(solid->shape->nrmp)) ( hit, solid.parms, nrm );
}

void prnsolid ( FILE *f, t_solid solid )
{
  (*(solid->shape->prnp)) (f, solid.parms);
}


