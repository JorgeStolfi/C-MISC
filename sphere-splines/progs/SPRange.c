/* See SPRange.h. */
/* Last edited on 2008-05-24 12:29:54 by stolfi */

#include <SPRange.h>
#include <SPTriang.h>
#include <SPBasic.h>
#include <r3.h>
#include <js.h>
#include <values.h>
#include <stdio.h>
#include <math.h>

double SPRange_OnTriangle
  ( S2Point *p, S2Point *q, S2Point *r, 
    ScalarField func, 
    int smpOrder
  )
  { double fMax = -INFINITY;
    int mm  = smpOrder;
    double fmm = ((double)mm);
    int j,k;
    
    for (k = 0; k <= mm; k++)
      { for (j = 0; j <= mm-k; j++)
          { int i = mm - (k+j);
            double fi = ((double)i);
            double fj = ((double)j);
            double fk = ((double)k);
            r3_t t = (r3_t){{
              (p->c[0]*fi + q->c[0]*fj + r->c[0]*fk)/fmm,
              (p->c[1]*fi + q->c[1]*fj + r->c[1]*fk)/fmm,
              (p->c[2]*fi + q->c[2]*fj + r->c[2]*fk)/fmm
            }};
            double ft, absft;
            r3_dir(&t, &t);
            ft = func(&t);
            absft = fabs(ft);
            if (absft > fMax) { fMax = absft; }
          }
      }
    return fMax;
  }

double SPRange_OnTriangulation(ScalarField func, Triangulation *tri, int smpOrder)
  { int k;
    int NT = tri->side.ne;
    double fMax = -INFINITY;
    for (k = 0; k < NT; k++)
      { Arc e = tri->side.e[k];
        r3_t *p = &(Org(e)->pos);
        r3_t *q = &(Org(Lnext(e))->pos);
        r3_t *r = &(Org(Lprev(e))->pos);
        double fMaxTri = SPRange_OnTriangle(p, q, r, func, smpOrder);
        if (fMaxTri > fMax) { fMax = fMaxTri; }
      }
    return fMax;
  }

double SPRange_OnSphere(ScalarField func, int smpOrder)
  { double fMax = -INFINITY;
    int x, y, z;
    for (x = -1; x <= 1; x += 2)
      { for (y = -1; y <= 1; y += 2)
          { for (z = -1; z <= 1; z += 2)
              { r3_t p = (r3_t){{((double)x), 0.0, 0.0}};
                r3_t q = (r3_t){{0.0, ((double)y), 0.0}};
                r3_t r = (r3_t){{0.0, 0.0, ((double)z)}};
                double fMaxTri = SPRange_OnTriangle(&p, &q, &r, func, smpOrder);
                if (fMaxTri > fMax) { fMax = fMaxTri; }
              }
          }
      }
    return fMax;
  }
