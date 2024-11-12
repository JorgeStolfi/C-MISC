/* See SPHarmonic.h */
/* Last edited on 2023-02-12 07:50:23 by stolfi */

#include <SPHarmonic.h>
#include <SPFunction.h>
#include <SPBasic.h>

#include <fget.h>
#include <nget.h>
#include <filefmt.h>
#include <affirm.h>
#include <nat.h>

#include <math.h>
// #include <gsl/gsl_sf_legendre.h>
#include <stdio.h>
#include <stdlib.h>
#include <values.h>

#define T SPHarmonic

double SPHarmonic_M_Eval(T *f, R3Point *p);
R3Gradient SPHarmonic_M_Grad(T *f, R3Point *p);
R3Hessian SPHarmonic_M_Hess(T *f, R3Point *p);
void SPHarmonic_M_Write(T *f, FILE *wr);
void SPHarmonic_M_Add(T *f, double a, T *h);
void SPHarmonic_M_Scale(T *f, double a);
void SPHarmonic_M_Maple(T *f, FILE *wr);
SPHarmonic *SPHarmonic_M_Copy(T *f);
void SPHarmonic_M_Free(T *f);

SPHarmonic_Methods *SPHarmonic_Methods_New(void);
SPHarmonic_Data *SPHarmonic_Data_New(void);
SPHarmonic *SPHarmonic_New(void);

SPHarmonic *SPHarmonic_FullNew(void);
  /* Allocates a new {SPHarmonic} object {f}, including its data
    record, with the proper {type} fields and
    methods record. The data record fields are left undefined;
    the coefficient vector {f->d->t} is not allocated. */

SPHarmonic *SPHarmonic_Cast(OBJ *f)
  { SPFunction *ff = (SPFunction *)f;
    if ((f != NULL) && isprefix(SPHarmonic_TypeId, ff->type))
      { return (SPHarmonic *)f; }
    else
      { return NULL; }
  }

#define SPHarmonic_FileFormat "2003-01-24"

/* CACHE OF HARMONIC FACTORS */ 

/* The following fields are used to speed up the computation
  of harmonic terms for a given argument point. 
  They are defined only for positive {m}
  and up to indices {curd} and {curm}, and are valid for 
  arguments {curx, cury, curz}. */

/* OVERRIDES FOR PARENT CLASS METHODS */
  
double SPHarmonic_M_Eval(T *f, R3Point *p) 
  { double v;
    SPHarmonic_EvalTerms(f->d->t, p, 0, &v, NULL, NULL);
    return v;
  } 

R3Gradient SPHarmonic_M_Grad(T *f, R3Point *p) 
  { double v;
    R3Gradient dv;
    SPHarmonic_EvalTerms(f->d->t, p, 1, &v, &dv, NULL);
    return dv;
  } 

R3Hessian SPHarmonic_M_Hess(T *f, R3Point *p) 
  { double v;
    R3Gradient dv;
    R3Hessian ddv;
    SPHarmonic_EvalTerms(f->d->t, p, 2, &v, &dv, &ddv);
    return ddv;
  } 

void SPHarmonic_M_Maple(T *f, FILE *wr)
  { affirm(FALSE , "maple method not implemented yet"); }

#define TermSeq(t) (((t)->degree)*((t)->degree) + (t)->degree + (t)->order)

void SPHarmonic_M_Add(T *f, double a, T *h)
  { affirm(isprefix(SPHarmonic_TypeId, f->type), "type/method bug");
    affirm(isprefix(SPHarmonic_TypeId, h->type), "incompatible argument");
    { /* Merge the two term lists: */
      int fn = f->d->t.ne, hn = h->d->t.ne;
      int fi = 0, hi = 0;
      int rn = 0;
      HarmonicTerm_vec_t t = HarmonicTerm_vec_new(f->d->t.ne); 
      /* Merge/add the term lists: */
      while ((fi < fn) || (hi < hn))
        { HarmonicTerm *ft = (fi < fn ? &(f->d->t.e[fi]) : NULL);
          HarmonicTerm *ht = (hi < hn ? &(h->d->t.e[hi]) : NULL);
          int fs = (ft != NULL ? TermSeq(ft) : INT_MAX);
          int hs = (ht != NULL ? TermSeq(ht) : INT_MAX);
          double rv;
          HarmonicTerm r;
          if (fs < hs)
            { r.degree = ft->degree; r.order = ft->order;
              r.coeff = ft->coeff;
              fi++;
            }
          else if (hs < fs)                                        
            { r.degree = ht->degree; r.order = ht->order;  
              r.coeff = a * ht->coeff;
              hi++;
            }
          else
            { r.degree = ht->degree; r.order = ht->order;
              r.coeff = ft->coeff + a * ht->coeff; 
              fi++; hi++;
            }
          if (rv != 0.0)
            { /* Store new term in {t[rn]}: */
              HarmonicTerm_vec_expand(&t, rn);
              t.e[rn] = r;
              rn++;
            }
        }
      /* Install new terms on {f}: */
      HarmonicTerm_vec_trim(&t, rn);
      free(f->d->t.e); f->d->t = t;
    }
  }

void SPHarmonic_M_Scale(T *f, double a)
  { if (a == 0.0)
      { HarmonicTerm_vec_trim(&(f->d->t), 0); }
    else
      { int NT = f->d->t.ne, i;
        HarmonicTerm *t = &(f->d->t.e[0]);
        for (i = 0; i < NT; i++, t++)
          { t->coeff *= a; }
      }
  }

T *SPHarmonic_M_Copy(T *f)
  { affirm(isprefix(SPHarmonic_TypeId, f->type), "type/method bug");
    { int NT = f->d->t.ne;
      SPHarmonic *g = SPHarmonic_FullNew();
      g->d->t = HarmonicTerm_vec_new(NT);
      { int j; for (j = 0; j < NT; j++) { g->d->t.e[j] = f->d->t.e[j]; } }
      return g;
    }
  }
 
void SPHarmonic_M_Free(T *f)
  { affirm(isprefix(SPHarmonic_TypeId, f->type), "type/method bug");
    free(f->d->t.e);
    free(f->d);
    free(f);
  }
  
/* CLASS-SPECIFIC METHODS */
  
void SPHarmonic_M_Write(T *f, FILE *wr)
  { HarmonicTerm_vec_t t = f->d->t;
    int i;
    filefmt_write_header(wr, "SPHarmonic", SPHarmonic_FileFormat);
    fprintf(wr, "terms = %d\n", t.ne);
    for (i = 0; i < t.ne; i++)
      { HarmonicTerm *ti = &(t.e[i]);
        fprintf(wr, "%d %d %22.16e\n", ti->degree, ti->order, ti->coeff);
      }
    filefmt_write_footer(wr, "SPHarmonic");
    fflush(wr);
  }
    
/* OTHER PROCS */
  
SPHarmonic_Methods *SPHarmonic_Methods_New(void)
  { void *v = malloc(sizeof(SPHarmonic_Methods));
    return (SPHarmonic_Methods *)notnull(v, "no mem for SPHarmonic_Methods");
  }

SPHarmonic_Data *SPHarmonic_Data_New(void)
  { void *v = malloc(sizeof(SPHarmonic_Data));
    return (SPHarmonic_Data *)notnull(v, "no mem for SPHarmonic_Data");
  }

SPHarmonic *SPHarmonic_New(void)
  { void *v = malloc(sizeof(SPHarmonic));
    return (SPHarmonic *)notnull(v, "no mem for SPHarmonic");
  }

static SPHarmonic_Methods *HarmonicMths = NULL;

SPHarmonic *SPHarmonic_FullNew(void)
  { SPHarmonic *f = SPHarmonic_New();
    f->type = SPHarmonic_TypeId;
    f->d = SPHarmonic_Data_New();
    f->d->t = (HarmonicTerm_vec_t){0, NULL};
    if (HarmonicMths == NULL)
      { HarmonicMths = SPHarmonic_Methods_New();
        HarmonicMths->fn.eval = (SPFunction_EvalMth *)&SPHarmonic_M_Eval;
        HarmonicMths->fn.grad = (SPFunction_GradMth *)&SPHarmonic_M_Grad;
        HarmonicMths->fn.hess = (SPFunction_HessMth *)&SPHarmonic_M_Hess;
        HarmonicMths->fn.maple = (SPFunction_MapleMth *)&SPHarmonic_M_Maple;
        /* Note: the {fn.write} method is inherited from {SPFunction}! */
        HarmonicMths->fn.write = (SPFunction_WriteMth *)&SPFunction_M_Write;
        HarmonicMths->fn.scale = (SPFunction_ScaleMth *)&SPHarmonic_M_Scale;
        HarmonicMths->fn.add = (SPFunction_AddMth *)&SPHarmonic_M_Add;
        HarmonicMths->fn.copy = (SPFunction_CopyMth *)&SPHarmonic_M_Copy;
        HarmonicMths->fn.free = (SPFunction_FreeMth *)&SPHarmonic_M_Free;
        /* Class-specific methods */
        HarmonicMths->write = (SPFunction_WriteMth *)&SPHarmonic_M_Write;
      }
    f->m = HarmonicMths;
    return f;
  }

/* HARMONIC TERMS */

void SPHarmonic_EvalTerms
  ( HarmonicTerm_vec_t t, 
    R3Point *p, 
    int diff,
    double *v,
    r3_t *dv,
    r6_t *ddv
  )
  { int NT = t.ne;
    
    /* Convert Cartesian cordinates to {cos(lon), sin(lon), sin(lat)}: */
    double clon = p->c[0], slon = p->c[1];
    double r = hypot(clon, slon);
    if (r > 0) { clon /= r; slon /= r; }

    double slat = p->c[2];
    double R = hypot(r, slat);
    if (R > 0) { slat /= R; }

    double term, sum, corr; 
    r3_t dterm, dsum, dcorr; 
    r6_t ddterm, ddsum, ddcorr;
    
    /* Loop over terms: */
    int i;
    sum = corr = 0.0;
    dsum = dcorr = (r3_t) {{ 0.0, 0.0, 0.0 }};
    ddsum = ddcorr = (r6_t) {{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }};
    HarmonicTerm *ti = &(t.e[0]);
    for (i = 0; i < NT; i++, ti++)
      { SPHarmonic_EvalTerm(ti, clon, slon, slat, R, diff, &term, &dterm, &ddterm);
        if (diff >= 0)
          { /* Kahan's summation formula: */
            double tcorr = term - corr;
            double newSum = sum + tcorr;
            corr = (newSum - sum) - tcorr;
            sum = newSum;
          }
        if (diff >= 1)
          {
            int k;
            for (k = 0; k < 3; k++)
              { /* Kahan's summation formula: */
                double tcorr = dterm.c[k] - dcorr.c[k];
                double newSum = dsum.c[k] + tcorr;
                dcorr.c[k] = (newSum - dsum.c[k]) - tcorr;
                dsum.c[k] = newSum;
              }
          }
        if (diff >= 2)
          {
            int k;
            for (k = 0; k < 6; k++)
              { /* Kahan's summation formula: */
                double tcorr = ddterm.c[k] - ddcorr.c[k];
                double newSum = ddsum.c[k] + tcorr;
                ddcorr.c[k] = (newSum - ddsum.c[k]) - tcorr;
                ddsum.c[k] = newSum;
              }
          }
      }
    if (diff >= 0) { (*v) = sum; }
    if (diff >= 1) { (*dv) = dsum; }
    if (diff >= 2) { (*ddv) = ddsum; }
  }

#define INV_SQRT_4_PI (0.28209479177387814)

void SPHarmonic_EvalTerm
  ( HarmonicTerm *t, 
    double clon, 
    double slon, 
    double slat, 
    double R, 
    int diff,
    double *v,
    r3_t *dv,
    r6_t *ddv
  )
  {
    int d = t->degree;
    if (d == 0)
      { if (diff >= 0) { (*v) = t->coeff*INV_SQRT_4_PI; }
        if (diff >= 1) { (*dv) = (r3_t){{ 0.0, 0.0, 0.0 }}; }
        if (diff >= 2) { (*ddv) = (r6_t){{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }}; }
      }
    else
      { int m = t->order;
        m = (m + d) % (2*d+1) - d;
        double xx = clon, yy = slon;
        SPHarmonic_Trig(m, &xx, &yy);
        double leg = SPHarmonic_Legendre(d, m, slat);
        double hpos = t->coeff * leg * (xx + yy);
        double hneg = t->coeff * leg * (xx - yy);
        if (diff >= 0) 
          { 
            (*v) = hpos * pow(R, d);
          }
        if (diff >= 1)
          {
            affirm(FALSE , "gradient not implemented yet");
            double alpha = m * pow(R, d-1);
            (*dv) = (r3_t)
              {{alpha * (clon * hpos - slon * hneg),
                alpha * (slon * hpos + clon * hneg),
                alpha * slat * hpos
              }};
          }
        if (diff >= 2)
          { 
            if (d == 1)
              { (*ddv) = (r6_t){{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }}; }
            else
              { affirm(FALSE , "hessian not implemented yet");
                (*ddv) = (r6_t){{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }};
              }
          }
      }
  }

SPHarmonic *SPHarmonic_FromTerms(HarmonicTerm_vec_t t)
  { SPHarmonic *f = SPHarmonic_FullNew();
    f->d->t = t;
    return f;
  }

HarmonicTerm SPHarmonic_MakeTerm(nat_t degree, nat_t order)
  { HarmonicTerm t;
    t.degree = degree;
    t.order = order;
    t.coeff = 1.0;
    return t;
  }
    
SPHarmonic *SPHarmonic_FromTerm(nat_t degree, nat_t order)
  { SPHarmonic *f = SPHarmonic_FullNew();
    HarmonicTerm_vec_t t = HarmonicTerm_vec_new(1);
    t.e[0] = SPHarmonic_MakeTerm(degree, order);
    f->d->t = t;
    return f;
  }

SPHarmonic *SPHarmonic_Read(FILE *rd)
  { SPHarmonic *f = SPHarmonic_FullNew();
    int NT, i;
    filefmt_read_header(rd, "SPHarmonic", SPHarmonic_FileFormat);
    fget_skip_formatting_chars(rd);
    NT = nget_int32(rd, "terms"); fget_eol(rd);
    f->d->t = HarmonicTerm_vec_new(NT);
    for (i = 0; i < NT; i++)
      { HarmonicTerm *ti = &(f->d->t.e[i]);
        ti->degree = fget_int32(rd);
        ti->order = fget_int32(rd);
        ti->coeff = fget_double(rd);
        fget_eol(rd);
      }
    filefmt_read_footer(rd, "SPHarmonic");
    return f;
  }
  
Basis SPHarmonic_MakeBasis(nat_t degree)
  { int dim = (degree + 1)*(degree + 1);
    int n = 0, m;
    int d;
    Basis H = Basis_new(dim);
    for (d = 0; d <= degree;  d++) 
      { for (m = -d; m <= +d; m++)
          { HarmonicTerm_vec_t t = HarmonicTerm_vec_new(1);
            t.e[0] = SPHarmonic_MakeTerm(d, m);
            H.e[n] = (SPFunction *)SPHarmonic_FromTerms(t); n++;
          }
      }
    affirm (n == dim, "dimension error");
    return H;
  }

/* HARMONIC FACTORS */
  
double trigx = 0.0, trigy = 0.0;
static int trigm = -1;
static double_vec_t trigxm =  (double_vec_t){ 0, NULL };
static double_vec_t trigym =  (double_vec_t){ 0, NULL };
  /* Cache for {SPHarmonic_Trig}. The values {trigxm[m],trigym[m]}
    contain the real and imaginary parts of {(trigx + I*trigy)^m}.
    They are defined up to index {trigm}. */

void SPHarmonic_Trig(int m, double *x, double *y)
  { 
    if (m == 0)
      { (*x) = 1.0; (*y) = 0.0; }
    else
      { double cx, cy; 
        int am = (m < 0 ? -m : m);
        /* Flush cache {trigxm,trigym} when arguments {x,y} change: */
        if (((*x) != trigx) || ((*y) != trigy))
          { trigx = (*x); trigy = (*y); trigm = -1; }
        /* Get element {(x + I*y)^am} from the cache: */
        if (am > trigm) 
          { int i;
            /* Not yet cached, must compute it: */
            /* Make sure that there is space for it: */
            double_vec_expand(&trigxm, am);
            double_vec_expand(&trigym, am);
            /* Make sure that there is at least one entry in cache: */
            if (trigm < 0) { trigxm.e[0] = 1.0; trigym.e[0] = 0.0; trigm = 0; }
            /* Get last entry {cx,cy} in cache: */
            cx = trigxm.e[trigm]; cy = trigym.e[trigm];
            /* Compute missing cache entries and save them: */
            for (i = trigm+1; i <= am; i++)
              { double tx = cx*trigx - cy*trigy;
                double ty = cx*trigy + cy*trigx;
                trigxm.e[i] = cx = tx;
                trigym.e[i] = cy = ty;
              }
            trigm = am;
          }
        else
          { cx = trigxm.e[am]; cy = trigym.e[am]; }
        /* Return element with proper sign: */
        (*x) = cx; (*y) = (m < 0 ? -cy : cy );
      }
  }

//double SPHarmonic_Legendre(int d, int m, double z)
//  {
//    if (m < 0) { m = -m; }
//    return gsl_sf_legendre_sphPlm(d, m, z);
//  }

double SPHarmonic_Legendre(int d, int m, double z)
  { /* From /Numerical Recipes in Fortran/, p.180. */
    double p, p1, p2;
    double sgn = 1.0;
    double fm;
    int i;
    affirm(d >= 0 , "invalid degree");
    if (m < 0) { m = -m; if ((m % 2) == 1) { sgn = -1.0; } }
    fm = (double)m;
    affirm(m <= d , "invalid order");
    affirm(fabs(z) <= 1.0, "invalid argument");
    /* Compute {P_m^m(z)}: */
    p = 1.0;
    if (m > 0)
      { double h = sqrt((1.0 - z)*(1.0 + z));
        p = -p*h;
        if (m > 1) 
          { double fact = 3.0;
            for (i = 2; i <= m; i++) 
              { p = -p*fact*h;
                fact = fact + 2.0;
              }
          }
      }
    if (m == d) { return sgn*p; }
    /* Compute {P_{m+1}^m(z)}: */
    p1 = z*(2.0 * fm + 1.0)*p;
    if (m == d-1) { return sgn*p1; }
    /* Compute {P_i^m(z)} for {i = m+2,..d}: */
    for (i = m+2; i <= d; i++) 
      { double fi = (double)i;
        p2 = (z*(2.0*fi - 1.0)*p1 - (fi + fm - 1.0)*p)/(fi - fm);
        p = p1;
        p1 = p2;
      }
    return sgn*p2;
  }

/* Arrays of {HarmonicTerm}: */

vec_typeimpl(HarmonicTerm_vec_t,HarmonicTerm_vec,HarmonicTerm);
