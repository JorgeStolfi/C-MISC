#ifndef KRTINTERSECT_H
#define KRTINTERSECT_H

#include <krtshape.h>
#include <r3.h>

double krt_intersect(
    r3_t *pos, r3_t *ray, /* Origin and direction of ray */
    t_solid *solid[],     /* Objects to try */
    int nsolid,           /* Number of objects */
    r3_t *hit, 
    r3_t *nrm
);

#endif
