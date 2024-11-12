/* See SPIntegral.h */
/* Last edited on 2005-06-06 11:36:15 by stolfi */

#include <SPIntegral.h>
#include <SPTriang.h>
#include <SPBasic.h>
#include <vec.h>
#include <r3.h>
#include <r3x3.h>
#include <affirm.h>
#include <nat.h>
#include <math.h>
#include <stdlib.h>

typedef struct TmpTriangle { S2Point v[3]; } TmpTriangle;
typedef r3_t *TmpPlaneRef;

static double GWeight[13] = 
  {
   -0.149570044467670, 
    0.175615257433204, 
    0.175615257433204,
    0.175615257433204, 
    0.053347235608839,
    0.053347235608839,
    0.053347235608839, 
    0.077113760890257, 
    0.077113760890257,
    0.077113760890257,
    0.077113760890257,
    0.077113760890257,
    0.077113760890257
  };

static r3_t GPos[13] =
  { (r3_t){{
      0.3333333333333333, 
      0.3333333333333333, 
      0.3333333333333333
    }},
    (r3_t){{
      0.479308067841923,
      0.260345966079038,
      0.260345966079038
    }},
    (r3_t){{
      0.260345966079038,
      0.479308067841923,
      0.260345966079038
    }},
    (r3_t){{
      0.260345966079038,
      0.260345966079038,
      0.479308067841923
    }},
    (r3_t){{
      0.869739794195568,
      0.065130102902216,
      0.065130102902216
    }},
    (r3_t){{
      0.065130102902216,
      0.869739794195568,
      0.065130102902216
    }},
    (r3_t){{
      0.065130102902216,
      0.065130102902216,
      0.869739794195568
    }},
    (r3_t){{
      0.638444188569809,
      0.312865496004875,
      0.048690315425316
    }},
    (r3_t){{
      0.638444188569809,
      0.048690315425316,
      0.312865496004875
    }},
    (r3_t){{
      0.312865496004875,
      0.638444188569809,
      0.048690315425316
    }},
    (r3_t){{
      0.312865496004875,
      0.048690315425316,
      0.638444188569809
    }},
    (r3_t){{
      0.048690315425316,
      0.312865496004875,
      0.638444188569809
    }},
    (r3_t){{
      0.048690315425316,
      0.638444188569809,
      0.312865496004875
    }}
  }; 

void SPIntegral_OnSphericalTriangle
  ( S2Point *u, S2Point *v, S2Point *w, 
    double func(S2Point *p),
    double *sum, 
    double *corr
  )
  { r3x3_t b2c;
    r3_t uv, uw, n;
    double A = TriangleArea(u,v,w);
    int i, j;

    for (j = 0; j < 3; j++)
      { b2c.c[j][0] = u->c[j];
        b2c.c[j][1] = v->c[j];
        b2c.c[j][2] = w->c[j];
      }
    r3_sub(v, u, &uv);
    r3_sub(w, u, &uw);
    r3_cross(&uv, &uw, &n);
    r3_dir(&n, &n);
    for (i = 0;  i < 13;  i++)
      { R3Point q;
        S2Point p;
        double sqra, dp, na, term;
        r3_t *GPi = &(GPos[i]);
        r3x3_map_col(&b2c, GPi, &q);
        r3_dir(&q, &p);
        sqra = r3_norm_sqr(&q);
        na = r3_dot(&n, &p);
        dp = A*fabs(na/sqra);
        /* fprintf(stderr, "  i = %02d dp = %22.15e\n", i, dp); */

        term = GWeight[i]*func(&p)*dp;
        { /* Kahan's summation formula: */
          double tcorr = term - *corr;
          double newSum = *sum + tcorr;
          *corr = (newSum - *sum) - tcorr;
          *sum = newSum;
        }
      }
  }

void SPIntegral_OnFlatTriangle
  ( R3Point *u, R3Point *v, R3Point *w,
    double func(R3Point *p),
    double *sum, 
    double *corr
  )
  { r3x3_t b2c;
    double A = TriangleArea(u,v,w);
    int i, j;

    for (j = 0; j < 3; j++)
      { b2c.c[j][0] = u->c[j];
        b2c.c[j][1] = v->c[j];
        b2c.c[j][2] = w->c[j];
      }
    for (i = 0;  i < 13;  i++)
      { R3Point p;
        r3_t *GPi = &(GPos[i]);
        double term;
        r3x3_map_col(&b2c, GPi, &p);
        term = A*GWeight[i]*func(&p);

        { /* Kahan's summation formula: */
          double tcorr = term - *corr;
          double newSum = *sum + tcorr;
          *corr = (newSum - *sum) - tcorr;
          *sum = newSum;
        }
      }
  }

void SPIntegral_BySamples
  ( double func(S2Point *p),
    S2Point_vec_t sp,
    double_vec_t wp,
    double *sum, 
    double *corr
  )
  { int i;
    for (i = 0;  i < sp.ne;  i++)
      { S2Point *p = &(sp.e[i]);
        double w = wp.e[i];
        double term = w*func(p);
        /* Kahan's summation formula: */
        double tcorr = term - *corr;
        double newSum = *sum + tcorr;
        *corr = (newSum - *sum) - tcorr;
        *sum = newSum;
      }
  }
 
double SPIntegral_OnSphere
  ( double func(S2Point *p),
    Triangulation *tri
  )
  { int i;
    double sum = 0.0, corr = 0.0;
    if (tri == NULL) { tri = SPIntegral_GetDefaultTriangulation(); }
    for (i = 0;  i < tri->side.ne;  i++)
      { Arc e = tri->side.e[i];
        Face *f = Left(e);
        SPIntegral_BySamples(func, f->sp, f->wp, &sum, &corr);
      }
    return sum;
  }


/* GENERATING SAMPLE POINTS AND WEIGHTS */

void SPIntegral_GaussSampleTriangle
  ( S2Point *u, S2Point *v, S2Point *w, 
    S2Point_vec_t *sp, 
    double_vec_t *wp, 
    int *ns
  )
  { r3x3_t b2c;
    r3_t uv, uw, n;
    int i, j;
    int m = *ns;

    for (j = 0; j < 3; j++)
      { b2c.c[j][0] = u->c[j];
        b2c.c[j][1] = v->c[j];
        b2c.c[j][2] = w->c[j];
      }
    r3_sub(v, u, &uv);
    r3_sub(w, u, &uw);
    r3_cross(&uv, &uw, &n); /* Normal to triangle, scaled by 2*area. */
    S2Point_vec_expand(sp, m+12);
    double_vec_expand(wp, m+12);
    for (i = 0;  i < 13;  i++)
      { R3Point t; /* Point on flat triangle. */
        S2Point s; /* Its central projection on the sphere. */
        r3_t *GPi = &(GPos[i]);
        r3x3_map_col(&b2c, GPi, &t);
        r3_dir(&t, &s);
        { double t2 = r3_norm_sqr(&t);
          double dsdb = 0.5 * fabs(r3_dot(&n, &s))/t2;
          double ds = GWeight[i]*dsdb;
          /* fprintf(stderr, "  i = %02d ds = %22.15e\n", i, ds); */
          sp->e[m] = s;
          wp->e[m] = ds;
          m++;
        }
      }
    (*ns) = m;
  }

void SPIntegral_SuperSampleTriangle
  ( S2Point *u, S2Point *v, S2Point *w,
    int smpOrder, 
    S2Point_vec_t *sp, 
    double_vec_t *wp, 
    int *ns
  )
  { 
    int N = smpOrder;
    double fN = (double)N;
    S2Point_vec_t s = S2Point_vec_new(N + 1); /* Saved mesh corners. */
    int i;
    affirm(N >= 1, "invalid sampling order");
    for (i = 0; i <= N; i++)
      { int j;
        for (j = 0; j <= N-i; j++) 
          { int k = N - i - j;
            r3_t t = (r3_t)
              {{(u->c[0]*i + v->c[0]*j + w->c[0]*k)/fN,
                (u->c[1]*i + v->c[1]*j + w->c[1]*k)/fN,
                (u->c[2]*i + v->c[2]*j + w->c[2]*k)/fN
              }};
            r3_dir(&t, &t); 
            if (i > 0)
              { /* Sample triangles using points from previous row: */
                r3_t *p = &(s.e[j]);
                r3_t *q = &(s.e[j+1]);
                if (j > 0)
                  { r3_t *r = &(s.e[j-1]);
                    SPIntegral_GaussSampleTriangle(p, r, &t, sp, wp, ns);
                  }
                SPIntegral_GaussSampleTriangle(q, p, &t, sp, wp, ns);
              }
            /* Save point for next row: */
            s.e[j] = t;
          }
      }
    free(s.e);
  }

void SPIntegral_RecursiveSampleTriangle
  ( S2Point *u, S2Point *v, S2Point *w, 
    int nMin,
    S2Point_vec_t *sp, 
    double_vec_t *wp, 
    int *ns
  )
  { if (nMin <= 1)
      { double A = TriangleArea(u, v, w);
        int m = *ns;
        R3Point ctr;
        r3_add(u, v, &ctr);
        r3_add(w, &ctr, &ctr);
        S2Point_vec_expand(sp, m); 
        double_vec_expand(wp, m); 
        r3_dir(&ctr, &(sp->e[m]));
        wp->e[m] = A;
        (*ns) = m + 1;
      }
    else if (nMin <= 3)
      { double A = TriangleArea(u, v, w);
        R3Point ctr;
        int m = *ns;
        r3_add(u, v, &ctr);
        r3_add(w, &ctr, &ctr);
        r3_dir(&ctr, &ctr);
        S2Point_vec_expand(sp, m+2); 
        double_vec_expand(wp, m+2); 
        r3_mix(2.0, u, 1.0, &ctr, &(sp->e[m+0])); wp->e[m+0] = A/3.0;
        r3_mix(2.0, v, 1.0, &ctr, &(sp->e[m+1])); wp->e[m+1] = A/3.0;
        r3_mix(2.0, w, 1.0, &ctr, &(sp->e[m+2])); wp->e[m+2] = A/3.0;
        (*ns) = m + 3;
      }
    else
      { int NQ = (nMin + 3) / 4;
        S2Point uv, vw, wu;
        r3_add(u, v, &uv); r3_dir(&uv, &uv);
        r3_add(v, w, &vw); r3_dir(&vw, &vw);
        r3_add(w, u, &wu); r3_dir(&wu, &wu);
        SPIntegral_RecursiveSampleTriangle(  u, &uv, &wu, NQ, sp, wp, ns);
        SPIntegral_RecursiveSampleTriangle(  v, &vw, &uv, NQ, sp, wp, ns);
        SPIntegral_RecursiveSampleTriangle(  w, &wu, &vw, NQ, sp, wp, ns);
        SPIntegral_RecursiveSampleTriangle(&vw, &wu, &uv, NQ, sp, wp, ns);
      }
  }

void SPIntegral_SampleOctants
  ( int smpOrder, 
    S2Point_vec_t *sp, 
    double_vec_t *wp, 
    int *ns
  )
  { int x, y, z;
    for (x = -1;  x <= +1;  x+= 2)
      { for (y = -1;  y <= +1;  y += 2)
          { for (z = -1;  z <= +1;  z += 2)
              { S2Point p = (S2Point){{(double)x, 0.0, 0.0}};
                S2Point q = (S2Point){{0.0, (double)y, 0.0}};
                S2Point r = (S2Point){{0.0, 0.0, (double)z}};
                SPIntegral_SuperSampleTriangle(&p, &q, &r, smpOrder, sp, wp, ns);
              }
          }
      }
  }
  
void SPIntegral_SampleTriangulation
  ( Triangulation *tri, 
    int smpOrder,
    S2Point_vec_t *sp, 
    double_vec_t *wp, 
    int *ns
  )
  { int k;
    int NT = tri->side.ne;
    for (k = 0; k < NT; k++)
      { Arc e = tri->side.e[k];
        r3_t *p = &(Org(e)->pos);
        r3_t *q = &(Org(Lnext(e))->pos);
        r3_t *r = &(Org(Lprev(e))->pos);
        SPIntegral_SuperSampleTriangle(p, q, r, smpOrder, sp, wp, ns);
      }
  }

/* DEFAULT SAMPLING ORDER AND DEFAULT TRIANGULATION */

static int smpOrderDefault = 0;

static Triangulation *triDefault = NULL;

void SPIntegral_SetDefaultSamplingOrder(int smpOrder)
  { int old = smpOrderDefault;
    affirm(smpOrder > 0, "bad smpOrder");
    affirm((old == 0) || (smpOrder == old), "cannot change smpOrder");
    smpOrderDefault = smpOrder;
  }

int SPIntegral_GetDefaultSamplingOrder(void)
  { affirm(smpOrderDefault > 0, "sampling order was never set");
    return smpOrderDefault;
  }

Triangulation *SPIntegral_GetDefaultTriangulation(void)
  { if (triDefault == NULL)
      { affirm(smpOrderDefault > 0, "sampling order was never set");
        triDefault = SPTriang_RegularOctahedron(smpOrderDefault);
      }
    return triDefault;
  }

/* OBSOLETE */

double SPIntegral_OnTriangleXXX
  ( S2Point *u, S2Point *v, S2Point *w,
    double func(S2Point *p),
    int depth
  )
  { 
    double sum = 0.0, corr = 0.0;
    
    auto void process_triangle(S2Point *u, S2Point *v, S2Point *w, int depth);
      /* Subdivides {u,v,w} into {4^depth} spherical triangles
        and calls {SPIntegral_OnSphericalTriangle} on each of them. */
    
    void process_triangle(S2Point *u, S2Point *v, S2Point *w, int depth)
      { /* PrintTriangle(depth, u, v, w); */
        if (depth == 0)
          { SPIntegral_OnSphericalTriangle(u, v, w, func, &sum, &corr); }
        else
          { r3_t uv, vw, wu;
            r3_add(u, v, &uv); r3_dir(&uv, &uv);
            r3_add(v, w, &vw); r3_dir(&vw, &vw);
            r3_add(w, u, &wu); r3_dir(&wu, &wu);
            process_triangle(&uv,   v, &vw, depth-1);
            process_triangle(&wu,   u, &uv, depth-1);
            process_triangle(&vw,   w, &wu, depth-1);
            process_triangle(&uv, &vw, &wu, depth-1);
          }
      }
      
    process_triangle(u, v, w, depth);
    return sum;
  }

void CutOneTriangle
  ( S2Point *p, S2Point *q, S2Point *r,
    TmpPlaneRef pl,
    TmpTriangle *Parts,
    int *np
  );  
  /* Cuts triangle {p q r} by plane {*pl}. Breaks the part on the
    positive side of {*pl} into {k} (zero or more) triangles. Stores
    those triangles in {Parts[*np..*np+k-1]}, and increments {*np} by {k}. */

double SPIntegral_OnTwoTrianglesXXX
  ( S2Point *pa, S2Point *qa, S2Point *ra,
    S2Point *pb, S2Point *qb, S2Point *rb,
    double func(S2Point *p),
    int depth
  )
  { TmpTriangle R[9], Parts[9];
    int m, n;
    double sum = 0.0, corr = 0.0;
    r3_t u, v, w;
    TmpPlaneRef pl[3] = {&u, &v, &w};
    int i, j;

    r3_cross(qa, ra, &u);
    r3_cross(ra, pa, &v);
    r3_cross(pa, qa, &w);
    
    /* Decompose the intersection into triangles: */
    R[0] = (TmpTriangle){{*pb,*qb,*rb}}; m = 1;
    for (i = 0;  i < 3;  i++)
      { n = 0;
        for (j = 0; j < m; j++)
          { S2Point *v = &(R[j].v[0]);
            CutOneTriangle(&(v[0]), &(v[1]), &(v[2]), pl[i], Parts, &n);
          }
        for (j = 0; j < n; j++)
          { R[j] = Parts[j]; }
        m = n;
      }
      
    /* Integrate over all parts: */
    for (i = 0; i < m; i++)
      { S2Point *v = &(R[i].v[0]);
        { double term = SPIntegral_OnTriangleXXX(&(v[0]), &(v[1]), &(v[2]), func, depth);
          /* Kahan's summation formula: */
          double tcorr = term - corr;
          double newSum = sum + tcorr;
          corr = (newSum - sum) - tcorr;
          sum = newSum;
        }
      }
    return sum;
  }

void CutOneTriangle
  ( S2Point *p, S2Point *q, S2Point *r,
    r3_t *pl,
    TmpTriangle *Parts,
    int *np
  )
  { double sp, sq, sr;
    sp = r3_dot(p, pl);
    sq = r3_dot(q, pl);
    sr = r3_dot(r, pl);

    /* Sort {p, q, r} so that {sp >= sq >= sr}. */
    if (sp < sq )
      { { double aux = sq; sq = sp; sp = aux; }
        { r3_t *aux = q; q =  p; p = aux; }
      }
    if (sq < sr )
      { { double aux = sr; sr = sq; sq = aux; }
        { r3_t *aux = r; r = q; q = aux; }
      }
    if (sp < sq )
      { { double aux = sq; sq = sp; sp = aux; }
        { r3_t *aux = q; q = p; p = aux; }
      }

    /* Cuts {p q r} by plane: */
    if (sp <= 0.0 )
      { /* All negative or zero */
        return;
      }
    else if (sr >= 0.0 )
      { /* All positive */
        Parts[(*np)] = (TmpTriangle){{*p, *q, *r}};
        (*np)++;
      }
    else if (sq <= 0.0 )
      { r3_t u, v;
        r3_mix(fabs(sq), p, fabs(sp), q, &u); r3_dir(&u, &u);
        r3_mix(fabs(sr), p, fabs(sp), r, &v); r3_dir(&v, &v);
        Parts[(*np)] = (TmpTriangle){{*p, u, v}};
        (*np)++;
      }
    else
      { r3_t u, v;
        r3_mix(fabs(sr), q, fabs(sq), r, &u); r3_dir(&u, &u);
        r3_mix(fabs(sr), p, fabs(sp), r, &v); r3_dir(&v, &v);
        Parts[(*np)] = (TmpTriangle){{*p, *q, v}};
        (*np)++;
        Parts[(*np)] = (TmpTriangle){{*q, v, u}};
        (*np)++;
      }
  }

/*
void PrintTriangle(int depth, Point *p, Point *q, Point *r)
  {
    fprintf(stderr, "depth = %d\n",  depth);
    #define PC SPTriang_PrintPoint
    fprintf(stderr, "p = "); PC(stderr, p); fprintf(stderr, "\n");
    fprintf(stderr, "q = "); PC(stderr, q); fprintf(stderr, "\n");
    fprintf(stderr, "r = "); PC(stderr, r); fprintf(stderr, "\n");
    fprintf(stderr, "%s",  "\n");
  }
*/
