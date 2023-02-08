/* See {intersect_triangle.h}. */
/* Last edited on 2005-09-15 19:07:10 by stolfi */

#include <basic.h>
#include <intersect_triangle.h>

void intersect_triangle 
  ( r3_t *u, r3_t *v, 
    r3_t *a, r3_t *b, r3_t *c, 
    /*out:*/ 
    r3_t *q, r3_t *nq, 
    double *alpha
  )
  {

    /* Convert all points to homogeneous coordinates: */
    r4_t U = (r4_t){{ 1, u->c[0], u->c[1], u->c[2] }};
    r4_t V = (r4_t){{ 1, v->c[0], v->c[1], v->c[2] }};
    r4_t A = (r4_t){{ 1, a->c[0], a->c[1], a->c[2] }};
    r4_t B = (r4_t){{ 1, b->c[0], b->c[1], b->c[2] }};
    r4_t C = (r4_t){{ 1, c->c[0], c->c[1], c->c[2] }};

    r4_t P;
    r4_cross(&A, &B, &C, &P); /* {P} are the coeffs of the plane through {A,B,C}. */

    double d1 = r4_dot(&P,&U);
    double d2 = r4_dot(&P,&V);
    if ((d1*d2) > 0.0) { *alpha = INF; return; }

    bool_t flipuv = (d1 < 0);
    if (flipuv) { r4_t T = U; U = V; V = T; d1 = -d1; d2 = -d2; }

    if (r4_det(&A,&B,&V,&U) < 0) { *alpha = INF; return; }
    if (r4_det(&B,&C,&V,&U) < 0) { *alpha = INF; return; }
    if (r4_det(&C,&A,&V,&U) < 0) { *alpha = INF; return; }

    /* Compute the fraction {r}: */
    if ((d1 - d2) == 0) 
      { /* Points are too close to plane: */
        *alpha = 0.5; return;
      }
    double r = d1/(d1 - d2);
    if (r < 0.0) { r = 0.0; }
    if (r > 1.0) { r = 1.0; }
    if (flipuv) { r = 1.0 - r; }
    (*alpha) = r;

    /* Compute point of intersection: */
    r4_t Q;
    r4_mix(-d2, &U, d1, &V, &Q);

    /* Convert back to Cartesian coordinates: */
    (*q) = (r3_t){{ Q.c[1]/Q.c[0], Q.c[2]/Q.c[0], Q.c[3]/Q.c[0] }};

    /* Get plane's normal: */
    (*nq) = (r3_t){{ P.c[1], P.c[2], P.c[3] }};
    r3_dir(nq, nq);
  }

