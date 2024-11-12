#include <malloc.h>
#include "krttriangle.h"
#include <r3.h>

/* triangle class: */

static t_shape trishape = (t_shape) { inttri, nrmtri, prntri };

/* triangle parameters: */

typedef struct
  {
    t_shape *shape;      /* Shape (must be &trishape) */
    t_3d    nrm;         /* triangle normal      */
    double  d;           /* plane constant       */
    t_3d    e1, e2, e3;  /* edge vectors         */
    double  d1, d2, d3;  /* plane constants      */
    t_3d    v1, v2, v3;  /* vertices */
  } 
  o_tri;

/* create triangle object */

t_solid maktri(
    t_3d *p1, t_3d *p2, t_3d *p3
)
{
    t_3d vc1, vc2, vc3;
    o_tri *tri = (o_tri *) malloc(sizeof(o_tri));

    tri->shape = &trishape;

    /* vertices: */
    tri->v1 = p1;
    tri->v2 = p2;
    tri->v3 = p3;
    
    /* edge vectors */
    vc1.x = p2->x - p1->x;
    vc1.y = p2->y - p1->y;
    vc1.z = p2->z - p1->z;
    vc2.x = p3->x - p2->x;
    vc2.y = p3->y - p2->y;
    vc2.z = p3->z - p2->z;
    vc3.x = p1->x - p3->x;
    vc3.y = p1->y - p3->y;
    vc3.z = p1->z - p3->z;

    /* plane of triangle */
    crossp( &tri->nrm, &vc1, &vc2 );
    normalize( &tri->nrm );
    tri->d = dotp( &tri->nrm, p1 );

    /* edge planes */
    crossp( &tri->e1, &tri->nrm, &vc1 );
    normalize( &tri->e1 );
    tri->d1 = dotp( &tri->e1, p1 );

    crossp( &tri->e2, &tri->nrm, &vc2 );
    normalize( &tri->e2 );
    tri->d2 = dotp( &tri->e2, p2 );

    crossp( &tri->e3, &tri->nrm, &vc3 );
    normalize( &tri->e3 );
    tri->d3 = dotp( &tri->e3, p3 );
    
    return ( (t_solid) { &trishape, (t_parms *) tri } );
}

/* intersection calculation for ray and triangle */

double inttri( 
    t_3d *pos,          /* origin of ray */
    t_3d *ray,          /* ray vector */
    t_parms *parms      /* triangle parameters */
)
{
    double s, k;
    t_3d point;
    o_tri *tri = (o_tri *) parms;

    assert (tri->shape == &trishape, "inttri: wrong parms!");
    
    /* plane intersection */
    k = dotp( &tri->nrm, ray );
    if( k == 0 ) return( 0.0 );
    s = ( tri->d - dotp( &tri->nrm, pos ) ) / k;
    if( s <= 0 ) return( 0.0 );

    point.x = pos->x + ray->x *s;
    point.y = pos->y + ray->y *s;
    point.z = pos->z + ray->z *s;

    /* edge checks */
    k = dotp( &tri->e1, &point ) - tri->d1;
    if( k < 0 ) return( 0.0 );
    k = dotp( &tri->e2, &point ) - tri->d2;
    if( k < 0 ) return( 0.0 );
    k = dotp( &tri->e3, &point ) - tri->d3;
    if( k < 0 ) return( 0.0 );

    return( s );
}

/* normal calculation for triangle */

void nrmtri(
    t_3d *pos,       /* point of intersection */
    t_parms *parms,  /* triangle parameters */
    t_3d *nrm        /* return surface normal */
)
{
    o_tri *tri = (o_tri *) parms;

    assert (tri->shape == &trishape, "nrmtri: wrong parms!");

    nrm->x = tri->nrm.x;
    nrm->y = tri->nrm.y;
    nrm->z = tri->nrm.z;
    return;
}

void prntri ( FILE *f, t_parms * parms )
{
    o_tri  *tri = (o_tri *) parms;

    assert (tri->shape == &trishape, "prntri: wrong parms!");

    fprintf (f, "(triangle (%lf %lf %lf) (%lf %lf %lf) (%lf %lf %lf))", 
      tri->v1.x, tri->v1.y, tri->v1.z,
      tri->v2.x, tri->v2.y, tri->v2.z,
      tri->v3.x, tri->v3.y, tri->v3.z,
    );
}
