#include <math.h>
#include <malloc.h>
#include "krtsuperq.h"
#include <r3.h>

#define ERRCONST 1.05

/* superquadratic class: */

static t_shape supshape = (t_shape) { intsup, nrmsup, prnsup };

/* superquadratic parameters: */

typedef struct
  {
    t_shape *shape;      /* Shape (must be &supshape) */
    double  xs, ys, zs;             /* size of sides        */
    double  x, y, z;                /* center position      */
    double  pow;                    /* n in formula         */
    double  a, b, c, r;             /* coefficients         */
    double  err;                    /* error measure        */
  }
  o_sup;

/* create superquadratic object */

t_solid maksup(
    double x, double y, double z,
    double xs, double ys, double zs,
    double power
)
{
    double max;
    o_sup *sup = (o_sup *) malloc(sizeof(o_sup));

    sup->shape = &supshape;
    sup->x = x;
    sup->y = y;
    sup->z = z;
    sup->xs = xs;
    sup->ys = ys;
    sup->zs = zs;
    sup->pow = power;
    max = xs;
    if( ys > max ) max = ys;
    if( zs > max ) max = zs;
    sup->a = xs / max;
    sup->b = ys / max;
    sup->c = zs / max;
    sup->r = pow( max, power );
    sup->err = pow( ( ERRCONST * max ), power ) - sup->r;
    return ( (t_solid) { &supshape, (t_parms *) sup } );
}

/* intersection calculation for ray and superquadratic */

double intsup( 
    t_3d *pos,          /* origin of ray */
    t_3d *ray,          /* ray vector */
    t_parms *parms      /* superquadratic parameters */
)
{
    double xsiz, usiz, zsiz;
    double s, s1, t;
    double xadj, yadj, zadj;
    double old, result, p;

    o_sup *sup = (o_sup *) parms;

    assert (sup->shape == &supshape, "intsup: wrong parms!");
    
    /* find box intersection */
    s = intbbox( pos, ray, sup->x, sup->y, sup->z, sup->xs, sup->ys, sup->zs );
    if( s == 0 ) return( 0.0 );
    xadj = pos->x - sup->x;
    yadj = pos->y - sup->y;
    zadj = pos->z - sup->z;

    /* special case - ray origin within box */
    if( ( fabs( xadj ) < sup->xs ) &&
        ( fabs( yadj ) < sup->ys ) &&
        ( fabs( zadj ) < sup->zs ) ) s = 0;

    /* initial solution */
    p = sup->pow;
    result = pow( fabs( ( xadj + ray->x * s ) / sup->a ), p ) +
             pow( fabs( ( yadj + ray->y * s ) / sup->b ), p ) +
	     pow( fabs( ( zadj + ray->z * s ) / sup->c ), p ) -
	     sup->r;

    if( result < sup->err ) return( s );

    s1 = s;
    s = s + 0.001;
    /* interactive refinament */
    while( result > sup->err ){

	old = result;
	result = pow( fabs( ( xadj + ray->x * s ) / sup->a ), p ) +
                 pow( fabs( ( yadj + ray->y * s ) / sup->b ), p ) +
	         pow( fabs( ( zadj + ray->z * s ) / sup->c ), p ) -
	         sup->r;
	if( result >= old ) return( 0.0 );
	t = ( result * ( s - s1 ) ) / ( result - old );
	s1 = s;
	s -=t;
    }
    return( s );
}

/* normal calculation for superquadratic */

void nrmsup(
    t_3d *pos,       /* point of intersection */
    t_parms *parms,  /* superquadratic description */
    t_3d *nrm        /* return surface normal */
)
{
    double k;
    o_sup *sup = (o_sup *) parms;

    assert (sup->shape == &supshape, "nrmsup: wrong parms!");

    nrm->x = ( pos->x - sup->x ) / sup->a;
    nrm->y = ( pos->y - sup->y ) / sup->b;
    nrm->z = ( pos->z - sup->z ) / sup->c;
    k = sup->pow - 1;
    if( nrm->x > 0 ) nrm->x = pow( nrm->x, k );
    else nrm->x = -pow( -nrm->x, k );
    if( nrm->y > 0 ) nrm->y = pow( nrm->y, k );
    else nrm->y = -pow( -nrm->y, k );
    if( nrm->z > 0 ) nrm->z = pow( nrm->z, k );
    else nrm->z = -pow( -nrm->z, k );
    normalize( nrm );
}

void prnsup ( FILE *f, t_parms * parms )
{
    o_sup  *sup = (o_sup *) parms;

    assert (sup->shape == &supshape, "prnsup: wrong parms!");

    fprintf (f, "(superquadratic (%lf %lf %lf) (%lf %lf %lf) %lf)", 
      sup->x, sup->y, sup->z
      sup->xs, sup->ys, sup->zs,
      sup->pow
    );
}
    
double intbbox( 
    t_3d *pos,       /* origin of ray */
    t_3d *ray,       /* ray vector */
    double x, double y, double z,
    double xs, double ys, double zs
)
{
    double s, ss, xhit, yhit, zhit;
    double xadj, yadj, zadj;

    ss = FAR_AWAY;
    /* translate ray origin to object's space */
    xadj = pos->x - x;
    yadj = pos->y - y;
    zadj = pos->z - z;

    if( ray->x != 0 ){   /* check x faces */

	s = ( xs - xadj ) / ray->x;
	if( ( s > 0 ) && ( s < ss ) ){
	    yhit = fabs( yadj + s * ray->y );
	    zhit = fabs( zadj + s * ray->z );
	    if( ( yhit < ys ) && ( zhit < zs ) ){
		ss = s;
	    }
	}
	s = ( -xs - xadj ) / ray->x;
	if( ( s > 0 ) && ( s < ss ) ){
	    yhit = fabs( yadj + s * ray->y );
	    zhit = fabs( zadj + s * ray->z );
	    if( ( yhit < ys ) && ( zhit < zs ) ){
		ss = s;
	    }
	}
    }
    
    if( ray->y != 0 ){   /* check y faces */

	s = ( ys - yadj ) / ray->y;
	if( ( s > 0 ) && ( s < ss ) ){
	    xhit = fabs( xadj + s * ray->x );
	    zhit = fabs( zadj + s * ray->z );
	    if( ( xhit < xs ) && ( zhit < zs ) ){
		ss = s;
	    }
	}
	s = ( -ys - yadj ) / ray->y;
	if( ( s > 0 ) && ( s < ss ) ){
	    xhit = fabs( xadj + s * ray->x );
	    zhit = fabs( zadj + s * ray->z );
	    if( ( xhit < xs ) && ( zhit < zs ) ){
		ss = s;
	    }
	}
    }
       
    if( ray->z != 0 ){   /* check z faces */

	s = ( zs - zadj ) / ray->z;
	if( ( s > 0 ) && ( s < ss ) ){
	    xhit = fabs( xadj + s * ray->x );
	    yhit = fabs( yadj + s * ray->y );
	    if( ( xhit < xs ) && ( yhit < ys ) ){
		ss = s;
	    }
	}
	s = ( -zs - zadj ) / ray->z;
	if( ( s > 0 ) && ( s < ss ) ){
	    xhit = fabs( xadj + s * ray->x );
	    yhit = fabs( yadj + s * ray->y );
	    if( ( xhit < xs ) && ( yhit < ys ) ){
		ss = s;
	    }
	}
    }

    if( ss == FAR_AWAY ) return( 0.0 );
    return( ss );
}

