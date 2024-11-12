#include <math.h>
#include "krtintersect.h"
#include "krtshape.h"
#include <r3.h>

#define FAR_AWAY 99.99E+20
#define FUDGE_FACTOR 0.000001

double krt_intersect(
    r3_t *pos, r3_t *ray, /* Origin and direction of ray */
    t_solid *solid[],     /* Objects to try */
    int nsolid,           /* Number of objects */
    r3_t *hit, 
    r3_t *nrm
)
{
    int isldhit, isld;
    double s, ss, eps;

    isldhit = -1;
    ss = FAR_AWAY;
    /* check for intersection of ray with all objects */
    for( isld=0; isld<nsolid; isld++ ){
	
	/* special check used for reflections */
	s = intsolid ( pos, ray, solid[isld] );
	/* keep track of closest intersection */
	if( ( s > 0.0 ) && ( s <= ss ) ){
	    isldhit = isld;
	    ss = s;
	}
    }

    if( isldhit < 0 ) return(0); /* ray hit no objects */

    /* find point of intersection */
    hit->x = pos->x + ss * ray->x;
    hit->y = pos->y + ss * ray->y;
    hit->z = pos->z + ss * ray->z;

    krt_nrmsolid ( hit, solid[isldhit], nrm );
    
    /* Fudge intersection: */
    eps = FUDGE_FACTOR;
    if ( dotp (nrm, ray) > 0.0 ) eps = -eps;
    hit->x += eps * nrm->x;
    hit->y += eps * nrm->y;
    hit->z += eps * nrm->z;

    return(ss);
}



    
