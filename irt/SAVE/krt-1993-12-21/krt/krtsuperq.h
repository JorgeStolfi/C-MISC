#ifndef KRTSUPERQ_H
#define KRTSUPERQ_H

#include "krtshape.h"
#include <r3.h>

t_solid maksup(
    double x, double y, double z,
    double xs, double ys, double zs,
    double power
);

double intsup( 
    r3_t *pos,          /* origin of ray */
    r3_t *ray,          /* ray vector */
    t_parms *parms      /* superquadric description */
);

void nrmsup(
    r3_t *pos,       /* point of intersection */
    t_parms *parms,  /* superquadric description */
    r3_t *nrm        /* return surface normal */
);

#endif
