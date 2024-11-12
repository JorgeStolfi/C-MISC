#include "krtsphere.h"
#include <math.h>
#include <malloc.h>
#include <r3.h>
#include <krtshape.h>

/* sphere class: */

static t_shape sphshape = (t_shape) {krt_int_sph, krt_nrm_sph, krt_prn_sph};

/* sphere parameters: */

typedef struct
  {
    t_shape *shape;      /* Shape (must be &sphshape) */
    double  r;           /* radius       */
    double  x;           /* x - position */
    double  y;           /* y - position */
    double  z;           /* z - position */
  } 
  o_sph;

/* create sphere object */

t_solid krt_make_sph(
    double r, 
    double x, double y, double z
)
{
    o_sph *sph = (o_sph *) malloc(sizeof(o_sph));

    sph->shape = &sphshape;
    sph->r = r;
    sph->center.c[0] = x;
    sph->center.c[1] = y;
    sph->center.c[2] = z;

    return ((t_solid) { &sphshape, (t_parms *) sph });
}

double krt_int_sph(
    t_3d *org,          /* origin of ray */
    t_3d *dir,          /* ray vector */
    t_parms *parms,     /* sphere parameters */
    t_3d *hit           /* (out) hit point */
)
{
    double b, t, s1, s2, s;
    double xadj, yadj, zadj;
    o_sph *sph = (o_sph *) parms;

    assert (sph->shape == &sphshape, "krt_int_sph: wrong parms!");
    
    /* translate ray origin to object's space */
    xadj = org->c[0] - sph->center.c[0];
    yadj = org->c[1] - sph->center.c[1];
    zadj = org->c[2] - sph->center.c[2];

    /* solve quadratic equation */
    b = xadj * dir->c[0] + yadj * dir->c[1] + zadj * dir->c[2];
    t = b * b - xadj * xadj - yadj * yadj - zadj * zadj + sph->r * sph-> r;
    if(t < 0) return(0.0);
    s = -b - sqrt(t);   /* try smaller solution */
    if(s > 0.0) return(s);
    s = -b + sqrt(t);   /* try larger solution */
    if(s > 0.0) return(s);
    return(0.0); /* both solutions are negative */
}

void krt_nrm_sph(
    t_3d *hit,       /* point of intersection */
    t_parms *parms,  /* sphere description */
    t_3d *nrm        /* return surface normal */
)
{
    o_sph *sph = (o_sph *) parms;

    assert (sph->shape == &sphshape, "krt_nrm_sph: wrong parms!");

    nrm->c[0] = (hit->c[0] - sph->center.c[0]) / sph->r;
    nrm->c[1] = (hit->c[1] - sph->center.c[1]) / sph->r;
    nrm->c[2] = (hit->c[2] - sph->center.c[2]) / sph->r;
}

void krt_prn_sph (FILE *f, t_parms * parms)
{
    o_sph  *sph = (o_sph *) parms;

    assert (sph->shape == &sphshape, "krt_prn_sph: wrong parms!");

    fprintf (f, "(sphere %lf (%lf %lf %lf))", 
      sph->r, sph->center.c[0], sph->center.c[1], sph->center.c[2]
    );
}
    
