/* See {interact.h}. */
/* Last edited on 2005-09-15 19:11:05 by stolfi */

#include <interact.h>

#include <basic.h>
#include <wavefront.h>
#include <geomodel.h>
#include <triangulate.h>
#include <intersect_triangle.h>

bool_t inside_bbox ( r3_t *p, interval_t bb[] )
  {
    int k;
    for (k = 0; k < 3; k++)
      { double pk = p->c[k];
        if ((pk > bb[k].end[1]) || (pk < bb[k].end[0])) 
          { return FALSE; }
      }
   return TRUE;
  }

sign_t reflector_side ( r3_t *p, reflector_t *rf )
  {
    /* Compute the index of XY projected cells containing {p}'s projection: */
    int ipos[2];
    double step[2];
    int k;
    for (k = 0; k < 2; k++)
      { 
        /* Reflector's span, number of nodes, and step along axis {k}: */
        double rinf = rf->bb[k].end[0];
        double rsup = rf->bb[k].end[1];
        int nk = rf->nv[k];
        step[k] = (rsup - rinf)/(nk-1);

        /* Get cell index: */
        ipos[k] = (int)floor((p->c[k] - rinf)/step[k]);
        if (ipos[k] < 0) { ipos[k] = 0; }
        if (ipos[k] >= nk) { ipos[k] = nk - 1; }
      }

    /* Test grid cell for location: */
    int ix = ipos[0], iy = ipos[1];
    double x = rf->bb[0].end[0] + ix*step[0];
    double y = rf->bb[1].end[0] + iy*step[1];

    /* Mesh cell vertices in CCW order: */
    int nx = rf->nv[0];
    r3_t v1 = (r3_t){{ x + 0*step[0], y + 0*step[1], rf->z[(ix+0) + nx*(iy+0)] }};
    r3_t v2 = (r3_t){{ x + 1*step[0], y + 0*step[1], rf->z[(ix+1) + nx*(iy+0)] }};
    r3_t v3 = (r3_t){{ x + 1*step[0], y + 1*step[1], rf->z[(ix+1) + nx*(iy+1)] }};
    r3_t v4 = (r3_t){{ x + 0*step[0], y + 1*step[1], rf->z[(ix+0) + nx*(iy+1)] }};

    if (orient_xy(&v1, p, &v3) > 0)
      { return plane_side(&v1, &v2, &v3, p); }
    else 
      { return plane_side(&v3, &v4, &v1, p); }
  }
  
int containing_medium ( r3_t *p, geomodel_t *geo )
  {
    /* Check bounding box: */
    if (! inside_bbox(p, geo->bb)) { return -1; } 
    
    /* Binary search on reflectors: */
    int i = 0, j = geo->nmd - 1;
    while (i < j)
      { 
        /* At this point, {p} lies in media {i..j} or intervening reflectors. */
        /* Get an intermediate reflector: */
        int k = (i + j)/2;
        /* Check against it: */
        sign_t s = reflector_side(p, &(geo->rf[k]));
        if (s < 0)
          { /* Below {rf[k]}, i.e. in layers {k+1..j}: */ i = k+1; }
        else if (s > 0)
          { /* Above {rf[k]}, i.e. in layers {i..k}: */ j = k; }
        else
          { /* On {rf[k]}, pretend it is in layer {k}: */  i = k; j = k; }
          
        /* Paranoia: */
        assert(i <= j);
      }
    return i;
  } 
  
void first_intersection_1 
  ( r3_t *p, 
    r3_t *c, 
    double time,
    reflector_t *rf, 
    /*Out*/ 
    r3_t *q, 
    r3_t *nq, 
    double *alpha
  )
  {
    int k;
    interval_t sbb[3]; /* Bounding box of segment. */
    for (k = 0; k < 3; k++)
      { double ck = c->c[k], pk = p->c[k];
        sbb[k] = (interval_t){{ (ck < pk ? ck : pk), (ck > pk ? ck : pk) }};
        if (sbb[k].end[0] > rf->bb[k].end[1]) { (*alpha) = INF; return; }
        if (sbb[k].end[1] < rf->bb[k].end[0]) { (*alpha) = INF; return; }
      }

    /* Compute the index range of XY projected cells intersected by the XY segment: */
    int inf[2], sup[2];
    double step[2];
    for (k = 0; k < 2; k++)
      { interval_t sbbk = sbb[k];
        double ainf = sbbk.end[0];
        double asup = sbbk.end[1];

        double rinf = rf->bb[k].end[0];
        double rsup = rf->bb[k].end[1];

        int nk = rf->nv[k];
        step[k] = (rsup - rinf)/(nk-1);

        inf[k] = (int)floor((ainf - rinf)/step[k]);
        if (inf[k] < 0) { inf[k] = 0; }
        sup[k] = (int)ceil((asup - rinf)/step[k]);
        if (sup[k] >= nk) { sup[k] = nk - 1; }
      }

    /* Test all grid cells for intersection: */
    (*alpha) = INF;
    int ix, iy;
    for (ix = inf[0]; ix < sup[0]; ix++)
      {
        double x = rf->bb[0].end[0] + ix*step[0];
        for (iy = inf[1]; iy < sup[1]; iy++)
          {
            double y = rf->bb[1].end[0] + iy*step[1];

            /* Mesh cell vertices in CCW order: */
            r3_t v1 = (r3_t){{ x + 0*step[0], y + 0*step[1], rf->z[(ix+0) + rf->nv[0]*(iy+0)] }};
            r3_t v2 = (r3_t){{ x + 1*step[0], y + 0*step[1], rf->z[(ix+1) + rf->nv[0]*(iy+0)] }};
            r3_t v3 = (r3_t){{ x + 1*step[0], y + 1*step[1], rf->z[(ix+1) + rf->nv[0]*(iy+1)] }};
            r3_t v4 = (r3_t){{ x + 0*step[0], y + 1*step[1], rf->z[(ix+0) + rf->nv[0]*(iy+1)] }};

            double alphaa, alphab; /* Fractions up to each intersection. */
            r3_t qa, qb;           /* Positions of intersections. */
            r3_t nqa, nqb;         /* Normal directions at intersections. */

            intersect_triangle(p, c, &v1, &v2, &v3, /*out*/ &qa, &nqa, &alphaa);
            intersect_triangle(p, c, &v3, &v4, &v1, /*out*/ &qb, &nqb, &alphab);

            if (alphaa <= alphab)
              { if (alphaa < (*alpha)) { (*q) = qa; (*nq) = nqa; (*alpha) = alphaa; } }
            else
              { if (alphab < (*alpha)) { (*q) = qb; (*nq) = nqb; (*alpha) = alphab; } }
          }
      }
  }

void first_intersection 
  ( r3_t *p, 
    r3_t *c, 
    double time,
    geomodel_t *geo, 
    /*Out*/ 
    r3_t *q, 
    r3_t *nq, 
    double *alpha, 
    int *irf 
  )
  {
    int i;
    (*alpha) = INF;
    for (i = 0; i < geo->nrf; i++)
      { 
        double alphai; /* Time of intersection with reflector {i}. */
        r3_t qi, nqi;  /* Position of and normal direction at intersection. */

        first_intersection_1(p, c, time, &(geo->rf[i]), &qi, &nqi, &alphai);
        if (alphai < (*alpha)) 
          { (*q) = qi; (*nq) = nqi; (*alpha) = alphai; (*irf) = i; }
      }
  }

void interact
  ( sample_t *u,
    r3_t *ndir,
    int irf,
    reflector_t *rf, 
    double zeps
  )
  {      
    bool_t live = TRUE;
    assert(u->sgn[0] != '\000');
    if (u->sgn[1] == '\000')
      { /* End of signature, ray dies: */
        live = FALSE;
      }
    else
      { /* Get signature of interaction: */
        assert(u->sgn[2] != '\000');
        assert(u->sgn[3] != '\000');
        char imode = u->sgn[0];  /* Input wave mode: 'P' = pressure, 'S' = shear. */
        char xirf = u->sgn[1];   /* Reflector number in ASCII decimal. */
        char side = u->sgn[2];   /* Which side: 'r' = reflected, 't' = transmitted. */
        char omode = u->sgn[3];  /* Output wave mode: 'P' = pressure, 'S' = shear. */

        if ((xirf < '0') || (xirf > '9'))
          { affirm(FALSE, "invalid signature"); }
        else if (xirf != '0' + irf)
          { /* There is no continuation ray. */
            live = FALSE;
          }
        else
          { /* Perform reflection or refraction: */
            int dimd = snell_law(&(u->vel), imode, ndir, rf, side, omode);

            /* Update current medium index: */
            u->imd += dimd;

            /* Perturb position slightly to avoid bogus re-intersection: */
            double dotvn = r3_dot(&(u->vel), ndir);
            u->pos.c[2] += (dotvn > 0 ? +zeps : -zeps);
          }
        /* Update signature of remaining ray: */
        u->sgn += 3;
      }
      
    if (! live) { kill_sample(u); }
  }    

int snell_law
  ( r3_t *v, 
    char imode, 
    r3_t *ndir, 
    reflector_t *rf, 
    char side, 
    char omode
  )
  {
    /* Which side of the reflector are we? */
    double dotvn = r3_dot(v, ndir);

    /* Get the input medium, and relative co-sine sign: */
    medium_t *medium_in = (dotvn > 0 ? rf->m_below : rf->m_above );
    medium_t *medium_ot;
    double sign;
    if (side == 'r')
        { 
          /* Input and output ray are in the same medium: */
          medium_ot = medium_in;
          sign = -1;
        }
      else if (side == 't')
        { /* Follow transmitted (refracted) ray into other medium: */
          medium_ot = (dotvn > 0 ? rf->m_above : rf->m_below );
          sign = +1;
        }
      else
        { /* error - invalid signature */
          assert(FALSE);
        }

    if (medium_ot == NULL)
      { 
        /* Ray went out of model: */
        r3_zero(v);
      }
    else
      { /* Get input and output scalar speeds: */
        double sp_in = (imode == 'P' ?  medium_in->v_P : medium_in->v_S );
        double sp_ot = (omode == 'P' ?  medium_ot->v_P : medium_ot->v_S );

        /* Snell's law: */
        double sin_in = r3_sin(v, ndir);
        double sin_ot = sin_in*sp_ot/sp_in;
        if (sin_ot > 0.999) 
          { /* There is practically no outgoing ray of mode {omode}: */
            r3_zero(v);
          }
        else
          { 
            /* Co-sine of output angle w.r.t {ndir}: */
            double cos_ot = sqrt(1-sin_ot*sin_ot);
            if (sign*dotvn < 0) { cos_ot = -cos_ot; }

            r3_t tdir; /* Direction of tangential component of input {v}. */
            r3_mix(1.0, v, -dotvn, ndir, &tdir);
            r3_dir(&tdir, &tdir); 
            r3_t v_ot;
            r3_mix(sin_ot, &tdir, cos_ot, ndir, &v_ot);
            r3_scale(sp_ot, &v_ot, &v_ot);
            (*v) = v_ot;
          }
      }
    /* Position of output medium relative to input one: */
    if (medium_ot == NULL) 
      { return 0; }
    else
      { return (side == 'r' ? 0 : 1) * (dotvn > 0 ? -1 : +1); }
  }

