/* See SPH3.h */
/* Last edited on 2003-01-17 01:57:37 by anamaria */

/* Created 93-05-18 by Marcos C. Carrard.
  Based on H3.pas by J. Stolfi. */

#include <SPH3.h>
#include <r3.h>
#include <r4.h>
#include <r4x4.h>
#include <r6.h>
#include <math.h>
#include <affirm.h>
#include <nat.h>

SPH3_Point SPH3_FromCartesian(r3_t *c)
  { return (SPH3_Point){(r4_t){{1.0, c->c[0], c->c[1], c->c[2]}}}; }
    
r3_t SPH3_ToCartesian(SPH3_Point *p)
  { double w = p->h.c[0];
    affirm(w != 0.0, "null weight");
    return (r3_t){{p->h.c[1]/w, p->h.c[2]/w, p->h.c[3]/w}};
  }

Sign SPH3_Side(SPH3_Point *p, SPH3_Plane *Q)
  { double s = r4_dot(&(p->h), &(Q->f));
    if (s > 0.0)
      { return 1; }
    else if (s < 0.0)
      { return -1; }
    else
      { return 0; }
  }     

Sign SPH3_Orient(SPH3_Point *p, SPH3_Point* q, SPH3_Point *r, SPH3_Point *s)
  { double dt = r4_det(&(p->h), &(q->h), &(r->h), &(s->h));
    if (dt > 0.0)
      { return 1; }
    else if (dt < 0.0)
      { return -1; }
    else
      { return 0; }
  }

SPH3_Line SPH3_LineFromTwoPoints(SPH3_Point *p, SPH3_Point *q)
  { double a0 = p->h.c[0];
    double a1 = p->h.c[1];
    double a2 = p->h.c[2];
    double a3 = p->h.c[3];

    double b0 = q->h.c[0];
    double b1 = q->h.c[1];
    double b2 = q->h.c[2];
    double b3 = q->h.c[3];

    return (SPH3_Line){
      (r6_t){{
          a0*b1 - a1*b0,
          a0*b2 - a2*b0,
          a1*b2 - a2*b1,
          a0*b3 - a3*b0,
          a1*b3 - a3*b1,
          a2*b3 - a3*b2
      }}
    };
  }

SPH3_Plane SPH3_PlaneFromThreePoints(SPH3_Point *p, SPH3_Point *q, SPH3_Point *r)
  { SPH3_Plane P;
    r4_cross(&(p->h), &(q->h), &(r->h), &(P.f));
    return P;
  }   

SPH3_Point SPH3_PointFromThreePlanes(SPH3_Plane *P, SPH3_Plane *Q, SPH3_Plane *R)
  { SPH3_Point p;
    r4_cross(&(Q->f), &(P->f), &(R->f), &(p.h));
    return p;
  } 

SPH3_Plane SPH3_PlaneFromLineAndPoint(SPH3_Line *n, SPH3_Point *r)
  { double a01 = n->k.c[0];
    double a02 = n->k.c[1];
    double a12 = n->k.c[2];
    double a03 = n->k.c[3];
    double a13 = n->k.c[4];
    double a23 = n->k.c[5];

    double b0 = r->h.c[0];
    double b1 = r->h.c[1];
    double b2 = r->h.c[2];
    double b3 = r->h.c[3];

    SPH3_Plane P;

    P.f = (r4_t){{
      - a12*b3 + a13*b2 - a23*b1,
        a02*b3 - a03*b2 + a23*b0,
      - a01*b3 + a03*b1 - a13*b0,
        a01*b2 - a02*b1 + a12*b0
    }};
    return P;
  }

SPH3_Line SPH3_LineFromTwoPlanes(SPH3_Plane *P, SPH3_Plane *Q)
  { double a0 = P->f.c[0];
    double a1 = P->f.c[1];
    double a2 = P->f.c[2];
    double a3 = P->f.c[3];

    double b0 = Q->f.c[0];
    double b1 = Q->f.c[1];
    double b2 = Q->f.c[2];
    double b3 = Q->f.c[3];

    SPH3_Line m;

    m.k = (r6_t){{
        a2*b3 - a3*b2,
      - a1*b3 + a3*b1,
        a0*b3 - a3*b0,
        a1*b2 - a2*b1,
      - a0*b2 + a2*b0,
        a0*b1 - a1*b0
    }};
    return m;
  }

SPH3_Point SPH3_PointFromLineAndPlane(SPH3_Line *n, SPH3_Plane *R)
  { double a01 = n->k.c[0];
    double a02 = n->k.c[1];
    double a12 = n->k.c[2];
    double a03 = n->k.c[3];
    double a13 = n->k.c[4];
    double a23 = n->k.c[5];

    double b123 = R->f.c[0];
    double b023 = R->f.c[1];
    double b013 = R->f.c[2];
    double b012 = R->f.c[3];

    SPH3_Point p;

    p.h = (r4_t){{
        a01*b023 + a02*b013 + a03*b012,
      - a01*b123 + a12*b013 + a13*b012,
      - a02*b123 - a12*b023 + a23*b012,
      - a03*b123 - a13*b023 - a23*b013
    }};
    return p;
  }

SPH3_Point SPH3_MapPoint(SPH3_Point *p, SPH3_PMap *m)
  { SPH3_Point q;
    r4x4_map_row(&(p->h), &(m->dir), &(q.h));
    return q;
  }

SPH3_Point SPH3_InvMapPoint(SPH3_Point *p, SPH3_PMap *m)
  { SPH3_Point q;
    r4x4_map_row(&(p->h), &(m->inv), &(q.h));
    return q;
  } 

SPH3_Plane SPH3_MapPlane(SPH3_Plane *P, SPH3_PMap *m)
  { SPH3_Plane Q;
    r4x4_map_col(&(m->inv), &(P->f), &(Q.f));
    return Q;
  }

SPH3_Plane SPH3_InvMapPlane(SPH3_Plane *P, SPH3_PMap *m)
  { SPH3_Plane Q;
    r4x4_map_col(&(m->dir), &(P->f), &(Q.f));
    return Q;
  }
    
SPH3_PMap SPH3_CompMap(SPH3_PMap *m, SPH3_PMap *n)
  { SPH3_PMap mn;
    r4x4_mul(&(m->dir), &(n->dir), &(mn.dir));
    r4x4_mul(&(n->inv), &(m->inv), &(mn.inv));
    return mn;
  }   
  
SPH3_PMap SPH3_InvMap(SPH3_PMap *m)
  { SPH3_PMap n;
    n.dir = m->inv; n.inv = m->dir;
    return n;
  }

r3_t SPH3_Normal(SPH3_Plane *P)
  { double nx = P->f.c[1];
    double ny = P->f.c[2];
    double nz = P->f.c[3];
    double length = hypot(hypot(nx, ny), nz);
    return (r3_t){{nx/length, ny/length, nz/length}};
  }

r3_t SPH3_Dir(SPH3_Point *frm, SPH3_Point *tto)
  { double fw = frm->h.c[0];
    double tw = tto->h.c[0];

    double fx = frm->h.c[1];
    double tx = tto->h.c[1];
    double dx = fw * tx - tw * fx;

    double fy = frm->h.c[2];
    double ty = tto->h.c[2];
    double dy = fw * ty - tw * fy;

    double fz = frm->h.c[3];
    double tz = tto->h.c[3];
    double dz = fw * tz - tw * fz;

    double length = hypot(hypot(dx, dy), dz);

    return (r3_t){{dx/length, dy/length, dz/length}};
  }

double SPH3_Dist(SPH3_Point *a, SPH3_Point *b)
  { double aw = 1.0/a->h.c[0];
    double bw = 1.0/b->h.c[0];

    double ax = a->h.c[1];
    double bx = b->h.c[1];
    double dx = ax*aw - bx*bw;

    double ay = a->h.c[2];
    double by = b->h.c[2];
    double dy = ay*aw - by*bw;

    double az = a->h.c[3];
    double bz = b->h.c[3];
    double dz = az*aw - bz*bw;

    return hypot(hypot(dx, dy), dz);
  }

double SPH3_DistSqr(SPH3_Point *a, SPH3_Point *b)
  { double aw = 1.0/a->h.c[0];
    double bw = 1.0/b->h.c[0];

    double ax = a->h.c[1];
    double bx = b->h.c[1];
    double dx = ax*aw - bx*bw;

    double ay = a->h.c[2];
    double by = b->h.c[2];
    double dy = ay*aw - by*bw;

    double az = a->h.c[3];
    double bz = b->h.c[3];
    double dz = az*aw - bz*bw;

    return dx * dx + dy * dy + dz * dz;
  }

SPH3_PMap SPH3_PerspMap(SPH3_Point *obs, SPH3_Point *foc, SPH3_Point *upp)
  { SPH3_PMap m;
    r3_t r, s, t, u, v;
    r4x4_t mt, mr, mc, mrt;
    int i, j;
    affirm(foc->h.c[0] > 0.0, "bad foc");  /* Focus must be finite and hither */
    affirm(obs->h.c[0] >= 0.0, "bad obs"); /* Observer must be hither or infinite */
    affirm(upp->h.c[0] >= 0.0, "bad upp"); /* Zenith must be hither or infinite */

    /* Compute translation matrix */
    for (i = 0; i < 4; i++){
      for (j = 0; j < 4; j++){
        if (i == j)
          { mt.c[i][j] = foc->h.c[0]; }
        else if (i == 0)
          { mt.c[i][j] = -foc->h.c[j]; }
        else
          { mt.c[i][j] = 0.0; }
      }
    }
    
    /* Compute rotation matrix */
    t = SPH3_Dir(foc, obs);
    u = SPH3_Dir(foc, upp);
    r3_decomp(&u, &t, &v, &s);
    /* Zenith reference point {upp} must not be on imagesys Z axis: */
    affirm(r3_norm_sqr(&s) >= 1.0e-12, "bad zenith"); 
    r3_dir(&s, &s);
    r3_cross(&s, &t, &r);
    r3_dir(&r, &r);

    mr.c[0][0] = 1.0;
    for (i = 1; i < 4; i++)
      { mr.c[0][i] = 0.0;
        mr.c[i][0] = 0.0;
        mr.c[i][1] = r.c[i-1];
        mr.c[i][2] = s.c[i-1];
        mr.c[i][3] = t.c[i-1];
      }
    
    /* Compose the two matrices: */
    r4x4_mul(&mt, &mr, &mrt);

    if (obs->h.c[0] == 0.0)
      { /* Observer is at infinity; cilindrical projection. */
        m.dir = mrt;
      }
    else
      { /* Observer is finite; add conical projection step. */
        double d = SPH3_Dist(obs, foc);
        double uno = -1.0;
        if (fabs(d) > 1.0)
          { uno = -1.0/d; d = 1.0; }
        for (i = 0; i < 4; i++)
          { for (j = 0; j < 4; j++)
              { if (i == j)
                  { mc.c[i][j] = d; }
                else
                  { mc.c[i][j] = 0.0; }
              }
          }
        mc.c[3][0] = uno;
        r4x4_mul(&mrt, &mc, &(m.dir));
      }

    /* Compute inverse matrix: */
    r4x4_inv(&m.dir, &m.inv);
    return m;
  }
