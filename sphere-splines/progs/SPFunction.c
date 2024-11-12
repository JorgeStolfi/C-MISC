/* See SPFunction.h */
/* Last edited on 2023-02-12 07:50:34 by stolfi */

#include <SPFunction.h>
#include <SPProcFunction.h>
#include <SPSpline.h>
#include <SPHBezFunction.h>
#include <SPNHBezFunction.h>
#include <SPHarmonic.h>
#include <SPErrorMap.h>
#include <SPTriang.h>
#include <SPIntegral.h>
#include <SPMatrix.h>
#include <SPSys.h>
#include <SPApprox.h>
#include <SPBasisMatrices.h>
#include <SPQuad.h>

#include <filefmt.h>
#include <nget.h>
#include <fget.h>
#include <math.h> 
#include <r3.h> 
#include <rmxn.h> 
#include <sym_eigen.h> 
#include <affirm.h> 
#include <nat.h>
#include <bool.h> 

#include <stdio.h> 
#include <stdlib.h> 

#define OnTriangle SPIntegral_OnTriangle

#define T SPFunction

void SPFunction_WriteAsSubclass(SPFunction *f, FILE *wr);
  /* Casts {f} to its widest proper effective subclass, and 
    calls the corresponding {write} method (which must shadow,
    but must not override, the {write} method of {f}. */
    
SPMatrix SPFunction_GetNormalMatrix
  ( Basis g, 
    Metric dot, 
    bool_t verbose,
    SPMatrix *B
  );
  /* If {B} is not NULL, returns {*B}. Otherwise computes the matrix
     of dot products {BB[i,j] = dot(g[j],g[i])}, and returns it. */

SPMatrix SPFunction_GetCholeskyMatrix
  ( Basis g, 
    Metric dot, 
    bool_t verbose,
    SPMatrix *B,
    SPMatrix *L
  );
  /* If {L} is not NULL, returns {*L}. Otherwise, obtains the 
    matrix of dot products {BB[i,j] = dot(g[j],g[i])}
    by calling {SPFunction_GetNormalMatrix(g,dot,verbose,B)},
    and returns its lower Cholesky factor. */

void debugmat(int m, int n, char *name, double *A);
  /* Prints the {m} by {n} matrix {A} on {stderr}. */

/* IMPLEMENTATIONS */

void debugmat(int m, int n, char *name, double *A)
  { if (n != 5) { return; }
    fprintf(stderr, "\n%s:\n", name);
    rmxn_gen_print(stderr, m, n, A, "%16.9e", "", "", "\n", "  ", " ", "\n");
  }

SPFunction *SPFunction_Cast(OBJ *f)
  { SPFunction *ff = (SPFunction *)f;
    if ((f != NULL) && isprefix(SPFunction_TypeId, ff->type))
      { return (SPFunction *)f; }
    else
      { return NULL; }
  }

S2Gradient SPFunction_SGrd(SPFunction *f, S2Point *p)
  { R3Gradient fDp = f->m->grad(f, p);
    double fDr = r3_dot(&fDp, p);
    /* Remove radial component: */
    r3_mix(1.0, &fDp, -fDr, p, &fDp);
    return fDp;
  }

double SPFunction_SLap(SPFunction *f, S2Point *p)
  { 
    double x = p->c[0];
    double y = p->c[1];
    double z = p->c[2];
    
    r6_t H = f->m->hess(f, p); /* Cartesian Hessian. */
    double fxx = H.c[0];
    double fxy = H.c[1];
    double fyy = H.c[2];
    double fxz = H.c[3];
    double fyz = H.c[4];
    double fzz = H.c[5];
    
    r3_t g = f->m->grad(f, p); /* Cartesian gradient. */
    double fx = g.c[0];
    double fy = g.c[1];
    double fz = g.c[2];
    
    return 
      + fxx*(1-x*x) + fyy*(1-y*y) + fzz*(1-z*z) 
      - 2*(fxy*x*y + fxz*x*z + fyz*y*z + fx*x + fy*y + fz*z);
  }

SPH3_Plane SPFunction_GetSupportingPlane(SPFunction *f)
  { SPSpline *fpw;
    if ((fpw = SPSpline_Cast(f)) != NULL)
      { return fpw->d->supp; }
    else
      { return Omega; }
  }

#define SPFunction_FileFormat "2002-11-12"

void SPFunction_M_Write(T *f, FILE *wr)
  { char *t = f->type;
    affirm(isprefix(SPFunction_TypeId, t), "type mismatch");
    filefmt_write_header(wr, "SPFunction", SPFunction_FileFormat);
    fprintf(wr, "type = %s\n", t);
    SPFunction_WriteAsSubclass(f, wr);
    filefmt_write_footer(wr, "SPFunction");
    fflush(wr);
  }
  
void SPFunction_WriteAsSubclass(SPFunction *f, FILE *wr)
  {
    { SPProcFunction *fp = SPProcFunction_Cast(f);
      if (fp != NULL) { fp->m->write(fp, wr); return; } 
    }
    { SPSpline *fp = SPSpline_Cast(f);
      if (fp != NULL) { fp->m->write(fp, wr); return; } 
    }
    { SPHBezFunction *fp = SPHBezFunction_Cast(f);
      if (fp != NULL) { fp->m->write(fp, wr); return; } 
    }
    { SPNHBezFunction *fp = SPNHBezFunction_Cast(f);
      if (fp != NULL) { fp->m->write(fp, wr); return; } 
    }
    { SPHarmonic *fp = SPHarmonic_Cast(f);
      if (fp != NULL) { fp->m->write(fp, wr); return; } 
    }
    { SPErrorMap *fp = SPErrorMap_Cast(f);
      if (fp != NULL) { fp->m->write(fp, wr); return; } 
    }
    fprintf(stderr, "unknown SPFunction subclass \"%s\"\n", f->type);
    affirm(FALSE, "aborted");
  }

SPFunction *SPFunction_Read(FILE *rd)
  { char *t;
    SPFunction *f;
    filefmt_read_header(rd, "SPFunction", SPFunction_FileFormat);
    t = nget_string(rd, "type"); fget_eol(rd);
    affirm(isprefix(SPFunction_TypeId, t), "incorrect SPFunction type");
    if (isprefix(SPProcFunction_TypeId, t))
      { f = (SPFunction *)SPProcFunction_Read(rd); }
    else if (isprefix(SPSpline_TypeId, t))
      { f = (SPFunction *)SPSpline_Read(rd); }
    else if (isprefix(SPHBezFunction_TypeId, t))
      { f = (SPFunction *)SPHBezFunction_Read(rd); }
    else if (isprefix(SPNHBezFunction_TypeId, t))
      { f = (SPFunction *)SPNHBezFunction_Read(rd); }
    else if (isprefix(SPHarmonic_TypeId, t))
      { f = (SPFunction *)SPHarmonic_Read(rd); }
    else if (isprefix(SPErrorMap_TypeId, t))
      { f = (SPFunction *)SPErrorMap_Read(rd); }
    else 
      { fprintf(stderr, "unknown SPFunction subclass \"%s\"\n", t);
        affirm(FALSE, "aborted");
      }
    fget_skip_formatting_chars(rd);
    filefmt_read_footer(rd, "SPFunction");
    free(t);
    return f;
  }

#define SPFunction_Basis_FileFormat "2002-11-19"
  
void SPFunction_WriteBasis(FILE *wr, Basis F)
  { int i;
    filefmt_write_header(wr, "SPFunction_Basis", SPFunction_Basis_FileFormat);
    fprintf(wr, "dimension = %d\n", F.ne);
    for (i = 0; i < F.ne; i++)
      { SPFunction *f = F.e[i]; 
        fprintf(wr, "elemindex = %d\n", i);
        f->m->write(f, wr); fputc('\n', wr);
      }
    filefmt_write_footer(wr, "SPFunction_Basis");
    fflush(wr);
  }

Basis SPFunction_ReadBasis(FILE *rd)
  { int dim;
    Basis F;
    int i;
    filefmt_read_header(rd, "SPFunction_Basis", SPFunction_Basis_FileFormat);
    dim = nget_int32(rd, "dimension"); fget_eol(rd);
    affirm(dim >= 0, "bad dimension");
    F = SPFunctionRef_vec_new(dim);
    for (i = 0; i < dim; i++) 
      { int ii = nget_int32(rd, "elemindex"); 
        affirm(ii == i, "wrong element index");
        fget_eol(rd);
        F.e[i] = SPFunction_Read(rd); fget_skip_formatting_chars(rd); }
    filefmt_read_footer(rd, "SPFunction_Basis");
    return F;
  }

Basis SPFunction_GetSubBasis(Basis F, int_vec_t ix)
  { Basis G = Basis_new(ix.ne);
    int i;
    for (i = 0; i < ix.ne; i++) 
      { int ixi = ix.e[i];
        demand((ixi >= 0) && (ixi < F.ne), " invalid elem index" );
        G.e[i] = F.e[ixi];
      }
    return G;
  }

Basis SPFunction_CopyBasis(Basis F)
  { Basis G = Basis_new(F.ne);
    int i;
    for (i = 0; i < F.ne; i++) 
      { SPFunction *fi = F.e[i]; G.e[i] = fi->m->copy(fi); }
    return G;
  }

SPFunction *SPFunction_LinComb(double *a, Basis F)
  { if (F.ne == 0) 
      { return NULL; }
    else
      { SPFunction *r = F.e[0]->m->copy(F.e[0]);
        int i;
        r->m->scale(r, a[0]);
        for (i = 1; i < F.ne; i++)
          { double ai = a[i];
            if (ai != 0.0) {  r->m->add(r, ai, F.e[i]); }
          }
        return r;
      }
  }

double_vec_t SPFunction_Sample(SPFunction *f, S2Point_vec_t s)
  { int N = s.ne;
    double_vec_t v = double_vec_new(N);
    int i;
    for (i = 0; i < s.ne;  i++)
      { v.e[i] = f->m->eval(f, &(s.e[i])); }
    return v;
  }
  
/* INTEGRALS */

double SPFunction_Integral
  ( SPFunction *f,
    FuncMap FMap, 
    SPFunction *w,
    Triangulation *tri,
    bool_t verbose
  )
  {
    if (SPSpline_IsCompatible(f, tri) && FMap.zeroPres)
      { /* A piecewise function (on {tri}, if given) versus a general function: */
        SPSpline *fpw = SPSpline_Cast(f);
        return SPSpline_IntegralPW(fpw, FMap, w, verbose);
      }
    else
      { /* Use a general integrator (on {tri}, if given): */
        return SPFunction_RawIntegral(f, FMap, w, tri, verbose);
      }
  }

double SPFunction_RawIntegral
  ( SPFunction *f,
    FuncMap FMap, 
    SPFunction *w,
    Triangulation *tri,
    bool_t verbose
  )
  {    
    auto double func(S2Point *p);
      /* The integrand of the dot product. */

    double func(S2Point *p) 
      { double wval = (w == NULL ? 1.0 : w->m->eval(w, p));
        double fval = (f == NULL ? 0.0 : f->m->eval(f, p));
        if (FMap.map != NULL) { fval = FMap.map(fval, p); }
        return wval * fval;
      }

    return SPIntegral_OnSphere(func, tri);
  }

double SPFunction_CustomIntegral
  ( SPFunction *f, 
    FuncMap FMap, 
    SPFunction *w,
    S2Point_vec_t sp,
    double_vec_t wp
  ) 
  { 
    auto double func(S2Point *p);
    
    double func(S2Point *p)
      { double wval = (w == NULL ? 1.0 : w->m->eval(w, p));
        double fval = f->m->eval(f, p);
        if (FMap.map != NULL) { fval = FMap.map(fval, p); }
        return wval * fval;
      }

    double sum = 0.0, corr = 0.0;
    SPIntegral_BySamples(func, sp, wp, &sum, &corr);
    return sum;
  }

R3Point SPFunction_Centroid
  ( SPFunction *f,
    FuncMap FMap, 
    SPFunction *w,
    Triangulation *tri
  )
  { 
    if (SPSpline_IsCompatible(f, tri) && FMap.zeroPres)
      { /* A piecewise function (on {tri}, if given) versus a general function: */
        SPSpline *fpw = SPSpline_Cast(f);
        return SPSpline_CentroidPW(fpw, FMap, w, FALSE);
      }
    else
      { /* Use a general integrator (on {tri}, if given): */
        return SPFunction_RawCentroid(f, FMap, w, tri, FALSE);
      }
  }

R3Point SPFunction_RawCentroid
  ( SPFunction *f,
    FuncMap FMap, 
    SPFunction *w,
    Triangulation *tri,
    bool_t verbose
  )
  {    
    r3_t b;
    int i;
    for (i = 0; i < 3; i++) 
      { 
        auto double CFunc(S2Point *p);
          /* Integrand to compute the weighted mean of coordinate {i}. */ 
        
        double CFunc(S2Point *p)
          { double wval = (w == NULL ? 1.0 : w->m->eval(w, p));
            double fval = (f == NULL ? 0.0 : f->m->eval(f, p));
            if (FMap.map != NULL) { fval = FMap.map(fval, p); }
            return wval * fval * p->c[i]; }
          
        b.c[i] = SPIntegral_OnSphere(CFunc, tri);
      }
    return b;
  }

/* INNER PRODUCTS */

double SPFunction_Dot
  ( SPFunction *f,
    FuncMap FMap, 
    SPFunction *g, 
    FuncMap GMap,
    SPFunction *w,
    Triangulation *tri,
    bool_t verbose
  )
  {
    if (SPSpline_AreCompatible(f, g, tri) && FMap.zeroPres && GMap.zeroPres)
      { /* Two piecewise functions with same triangulation ({tri}, if given): */
        SPSpline *fpw = SPSpline_Cast(f);
        SPSpline *gpw = SPSpline_Cast(g);
        return SPSpline_DotBothPW(fpw, FMap, gpw, GMap, w, verbose);
      }
    else if (SPSpline_IsCompatible(f, tri) && FMap.zeroPres)
      { /* A piecewise function (on {tri}, if given) versus a general function: */
        SPSpline *fpw = SPSpline_Cast(f);
        return SPSpline_DotSinglePW(fpw, FMap, g, GMap, w, verbose);
      }
    else if (SPSpline_IsCompatible(g, tri) && GMap.zeroPres)
      { /* A piecewise function (on {tri}, if given) versus a general function: */
        SPSpline *gpw = SPSpline_Cast(g);
        return SPSpline_DotSinglePW(gpw, GMap, f, FMap, w, verbose);
      }
    else
      { /* Use a general integrator (on {tri}, if given): */
        return SPFunction_RawDot(f, FMap, g, GMap, w, tri, verbose);
      }
  }

double SPFunction_RawDot
  ( SPFunction *f,
    FuncMap FMap, 
    SPFunction *g, 
    FuncMap GMap,
    SPFunction *w,
    Triangulation *tri,
    bool_t verbose
  )
  {    
    auto double func(S2Point *p);
      /* The integrand of the dot product. */

    double func(S2Point *p) 
      { double wval = (w == NULL ? 1.0 : w->m->eval(w, p));
        double fval = (f == NULL ? 0.0 : f->m->eval(f, p));
        double gval = (g == NULL ? 0.0 : g->m->eval(g, p));
        if (FMap.map != NULL) { fval = FMap.map(fval, p); }
        if (GMap.map != NULL) { gval = GMap.map(gval, p); }
        return wval * fval * gval;
      }

    return SPIntegral_OnSphere(func, tri);
  }

double SPFunction_CustomDot
  ( SPFunction *f, 
    FuncMap FMap, 
    SPFunction *g,
    FuncMap GMap,
    SPFunction *w,
    S2Point_vec_t sp,
    double_vec_t wp
  ) 
  { 
    auto double func(S2Point *p);
    
    double func(S2Point *p)
      { double wval = (w == NULL ? 1.0 : w->m->eval(w, p));
        double fval = f->m->eval(f, p);
        double gval = g->m->eval(g, p);
        if (FMap.map != NULL) { fval = FMap.map(fval, p); }
        if (GMap.map != NULL) { gval = GMap.map(gval, p); }
        return wval * fval * gval;
      }

    double sum = 0.0, corr = 0.0;
    SPIntegral_BySamples(func, sp, wp, &sum, &corr);
    return sum;
  }

/* VELOCITY-GRADIENT INNER PRODUCTS */

double SPFunction_VelSGrdDot
  ( SPFunction *f,
    SPFunction *g, 
    SPFunction *w,
    Triangulation *tri,
    bool_t verbose
  )
  { 
    if (SPSpline_AreCompatible(f, g, tri))
      { /* Two piecewise functions with same triangulation ({tri}, if given): */
        SPSpline *fpw = SPSpline_Cast(f);
        SPSpline *gpw = SPSpline_Cast(g);
        return SPSpline_VelSGrdDotBothPW(fpw, gpw, w, verbose);
      }
    else if (SPSpline_IsCompatible(f, tri))
      { /* A piecewise function (on {tri}, if given) versus a general function: */
        SPSpline *fpw = SPSpline_Cast(f);
        return SPSpline_VelSGrdDotSinglePW(fpw, g, w, verbose);
      }
    else if (SPSpline_IsCompatible(g, tri))
      { /* A piecewise function (on {tri}, if given) versus a general function: */
        SPSpline *gpw = SPSpline_Cast(g);
        return SPSpline_VelSGrdDotSinglePW(gpw, f, w, verbose);
      }
    else
      { /* Use a general integrator (on {tri}, if given): */
        return SPFunction_RawVelSGrdDot(f, g, w, tri, verbose);
      }
  }

double SPFunction_RawVelSGrdDot
  ( SPFunction *f,
    SPFunction *g, 
    SPFunction *w,
    Triangulation *tri,
    bool_t verbose
  )
  { 
    auto double func(S2Point *p);
      /* The integrand of the dot product. */

    double func(S2Point *p) 
      { double wval = (w == NULL ? 1.0 : w->m->eval(w, p));
        r3_t ftan, tmp;
        R3Gradient fDp = f->m->grad(f, p);
        r3_decomp(&fDp, p, &tmp, &ftan);
        double dval = p->c[0]*ftan.c[1] - p->c[1]*ftan.c[0];
        double gval = g->m->eval(g, p);
        return wval * dval * gval;
      }

    return SPIntegral_OnSphere(func, tri);
  }

double SPFunction_CustomVelSGrdDot
  ( SPFunction *f, 
    SPFunction *g,
    SPFunction *w,
    S2Point_vec_t sp,
    double_vec_t wp
  ) 
  { 
    auto double func(S2Point *p);
    
    double func(S2Point *p)
      { double wval = (w == NULL ? 1.0 : w->m->eval(w, p));
        r3_t ftan, tmp;
        R3Gradient fDp = f->m->grad(f, p);
        r3_decomp(&fDp, p, &tmp, &ftan);
        double dval = p->c[0]*ftan.c[1] - p->c[1]*ftan.c[0];
        double gval = g->m->eval(g, p);
        return wval * dval * gval;
      }

    double sum = 0.0, corr = 0.0;
    SPIntegral_BySamples(func, sp, wp, &sum, &corr);
    return sum;
  }

/* LAPLACIAN INNER PRODUCTS */

double SPFunction_SLapDot
  ( SPFunction *f,
    SPFunction *g, 
    SPFunction *w,
    Triangulation *tri,
    bool_t verbose
  )
  { 
    if (SPSpline_AreCompatible(f, g, tri))
      { /* Two piecewise functions with same triangulation ({tri}, if given): */
        SPSpline *fpw = SPSpline_Cast(f);
        SPSpline *gpw = SPSpline_Cast(g);
        return SPSpline_SLapDotBothPW(fpw, gpw, w, verbose);
      }
    else if (SPSpline_IsCompatible(f, tri))
      { /* A piecewise function (on {tri}, if given) versus a general function: */
        SPSpline *fpw = SPSpline_Cast(f);
        return SPSpline_SLapDotSinglePW(fpw, g, w, verbose);
      }
    else if (SPSpline_IsCompatible(g, tri))
      { /* A piecewise function (on {tri}, if given) versus a general function: */
        SPSpline *gpw = SPSpline_Cast(g);
        return SPSpline_SLapDotSinglePW(gpw, f, w, verbose);
      }
    else
      { /* Use a general integrator (on {tri}, if given): */
        return SPFunction_RawSLapDot(f, g, w, tri, verbose);
      }
  }

double SPFunction_RawSLapDot
  ( SPFunction *f,
    SPFunction *g, 
    SPFunction *w,
    Triangulation *tri,
    bool_t verbose
  )
  { 
    auto double func(S2Point *p);
      /* The integrand of the dot product. */

    double func(S2Point *p) 
      { double wval = (w == NULL ? 1.0 : w->m->eval(w, p));
        r3_t ftan, gtan, tmp;
        R3Gradient fDp = f->m->grad(f, p);
        R3Gradient gDp = g->m->grad(g, p);
        r3_decomp(&fDp, p, &tmp, &ftan);
        r3_decomp(&gDp, p, &tmp, &gtan);
        return wval * r3_dot(&ftan, &gtan);
      }

    return SPIntegral_OnSphere(func, tri);
  }

double SPFunction_CustomSLapDot
  ( SPFunction *f, 
    SPFunction *g,
    SPFunction *w,
    S2Point_vec_t sp,
    double_vec_t wp
  ) 
  { 
    auto double func(S2Point *p);
    
    double func(S2Point *p)
      { double wval = (w == NULL ? 1.0 : w->m->eval(w, p));
        r3_t ftan, gtan, tmp;
        R3Gradient fDp = f->m->grad(f, p);
        R3Gradient gDp = g->m->grad(g, p);
        r3_decomp(&fDp, p, &tmp, &ftan);
        r3_decomp(&gDp, p, &tmp, &gtan);
        return wval * r3_dot(&ftan, &gtan);
      }

    double sum = 0.0, corr = 0.0;
    SPIntegral_BySamples(func, sp, wp, &sum, &corr);
    return sum;
  }

/* ORTHONORMALIZATION OF SPHERICAL FUNCTIONS - GENERAL METRICS  */

static SPFunction *unit = NULL;

void SPFunction_GenMakePositive(SPFunction *f, Metric dot)
  { double ff;
    if (unit == NULL) 
      { unit = (SPFunction *)SPProcFunction_FromName("unit"); }
    ff = dot(f, unit);
    if (ff < 0.0) { f->m->scale(f, -1.0); }
  }

double SPFunction_GenNormalize(SPFunction *f, Metric dot)
  { double ff = dot(f,f);
    if (ff == 0.0) 
      { return 0.0; }
    else
      { f->m->scale(f, 1.0/sqrt(fabs(ff))); 
        return (ff < 0.0 ? -1.0 : +1.0);
      }
  }

void SPFunction_GenOrthize
  ( SPFunction *f, 
    SPFunction *g, 
    Metric dot,
    double gg
  )
  { int maxIter = 3, iter;
    double fMag;
    for (iter = 1; iter <= maxIter; iter++)
      { double fg = dot(f,g);
        if (fg == 0.0) { break; }
        affirm(gg != 0.0, "orthizing against a null vector");
        f->m->add(f, -fg/gg, g);
        /* Check for convergence: */
        { double df = fabs(fg)/sqrt(fabs(gg));
          if (iter == 1) { fMag = sqrt(dot(f,f)); }
          if (df <= 1.0e-13*fMag) { break; } 
          if (iter >= maxIter) { fprintf(stderr, "!o(%8.2e:%8.2e) ", df, fMag); }
        }
      }
  }

SPFunction *SPFunction_GenProject
  ( SPFunction *f, 
    SPFunction *g, 
    Metric dot,
    double gg
  )
  { int maxIter = 3, iter;
    double fMag, hh = gg;
    SPFunction *h = g->m->copy(g);
    for (iter = 1; iter <= maxIter; iter++)
      { double fh = dot(f,h);
        if (fh == hh) { break; }
        affirm(hh != 0.0, "orthizing against a null vector");
        h->m->scale(h, fh/hh);
        /* Check for convergence: */
        { double df = (fabs(fh - hh))/sqrt(fabs(hh));
          if (df <= 1.0e-13) { break; } 
          if (iter == 1) { fMag = sqrt(dot(f,f)); }
          if (df <= 1.0e-13*fMag) { break; } 
          if (iter >= maxIter) { fprintf(stderr, "!p(%8.2e:%8.2e) ", df, fMag); }
        }
      }
    return h;
  }

void SPFunction_GenOrthizeBasisAgainstBasis
  ( Basis f, 
    Basis g, 
    Metric dot, 
    bool_t verbose,
    SPMatrix *B,
    SPMatrix *L
  )
  { if ((g.ne == 0) || (f.ne == 0))
      { /* Nothing to do. */ }
    else if (g.ne == 1)
      { /* Only one element to orthize against. */
        SPFunction *g0 = g.e[0];
        double g0g0 = dot(g0,g0);
        int i;
        for (i = 0; i < f.ne; i++)
          { SPFunction *fi = f.e[i];
            SPFunction_GenOrthize(fi, g0, dot, g0g0);
          }
      }
    else
      { /* Since the elements of {g} may not be orthogonal among
          themselves, we must solve for each {i} the linear system 
          {B x = y} where {y[s] = dot(f[i],g[s])} and {B[r,s] = dot(g[r],g[s])}.
          The system is solved through Cholesky factorization of {B}.
          Then we make {f[i] = f[i] - SUM{ x[j]*g[j] : j }}.

          We expect the matrix {B} to be large but sparse, so we use the 
          SPMatrix package rather than plain arrays.
        */
        double_vec_t y = double_vec_new(g.ne); /* Right-hand side of system. */
        double_vec_t x = double_vec_new(g.ne); /* Coefficients of projection. */ 
        double_vec_t t = double_vec_new(g.ne); /* Work area for {CholeskySolve}. */ 
        double_vec_t gNorm = double_vec_new(g.ne); /* Norm of {g} vectors. */
        SPMatrix LL = SPFunction_GetCholeskyMatrix(g, dot, verbose, B, L);
        int i;
        
        for (i = 0; i < g.ne; i++)
          { gNorm.e[i] = sqrt(fabs(dot(g.e[i],g.e[i]))); }
        if (verbose) { fprintf(stderr, "s"); }
        for (i = 0; i < f.ne; i++)
          { SPFunction *ci = f.e[i];
            double ciMag;
            int maxIter = 3, iter;
            for (iter = 1; iter <= maxIter; iter++)
              { double dciMax = 0.0;
                int j;
                /* Find coeffs {x[j]} of projection of {ci} on space of {g}: */
                SPBasisMatrices_RecomputeVectorGen(ci, g, dot, y, FALSE);
                SPSys_CholeskySolve(LL, y, t, x, FALSE);
                /* Subtract {x[j]*g[j]} from {f[i]}: */
                for (j = 0; j < g.ne; j++)
                  { double xj = x.e[j];
                    if (xj != 0.0)
                      { ci->m->add(ci, -xj, g.e[j]);
                        { double dci = fabs(xj)*gNorm.e[j]; 
                          if (dci > dciMax) { dciMax = dci; }
                        }
                      }
                  }
                /* if (verbose) { fprintf(stderr, "(%8.2e)", xMax); } */
                if (dciMax == 0.0) { break; }
                if (iter == 1) { ciMag = sqrt(dot(ci,ci)); }
                if (dciMax <= 1.0e-13*ciMag) { break; }
                if (iter >= maxIter)
                  { fprintf(stderr, "!");
                    if (verbose)
                      { fprintf(stderr, "(%8.2e:%8.2e)", dciMax, ciMag); }
                  }
              }
          }
        if (verbose) { fprintf(stderr, " "); }
        if (L == NULL) { free(LL.ents.e); }
        free(gNorm.e);
        free(x.e); free(t.e); free(y.e);
      }
  }

Basis SPFunction_GenProjectBasisOnBasis
  ( Basis f, 
    Basis g, 
    Metric dot, 
    bool_t verbose,
    SPMatrix *B,
    SPMatrix *L
  )
  { Basis h = Basis_new(f.ne);
    if (f.ne == 0)
      { /* Nothing to do. */ }
    else if (g.ne == 0)
      { /* The space generated by {g} is trivial, the projection is null: */
        int i; 
        for (i = 0; i < h.ne; i++) { h.e[i] = NULL; };
      }
    else if (g.ne == 1)
      { /* There is only one {g} element, reduces to simple projection. */
        SPFunction *g0 = g.e[0];
        double g0g0 = dot(g0,g0);
        int i;
        for (i = 0; i < f.ne; i++)
          { h.e[i] = SPFunction_GenProject(f.e[i], g0, dot, g0g0); }
      }
    else
      { /* Since the elements of {g} may not be orthogonal among
          themselves, we must solve for each {i} the linear system 
          {B x = y} where {y[s] = dot(f[i],g[s])} and {B[r,s] = dot(g[r],g[s])}.
          The system is solved through Cholesky factorization of {B}.
          Then we make {h[i] = SUM{ x[j]*g[j] : j }}.

          We expect the matrix {B} to be large but sparse, so we use the 
          SPMatrix package rather than plain arrays.
        */
        double_vec_t y  = double_vec_new(g.ne); /* Right-hand side of system. */
        double_vec_t dy = double_vec_new(g.ne); /* Residual of {y}. */
        double_vec_t x  = double_vec_new(g.ne); /* Coefficients of projection. */ 
        double_vec_t dx = double_vec_new(g.ne); /* Correction to {x}. */ 
        double_vec_t t  = double_vec_new(g.ne); /* Workarea for {CholeskySolve}. */
        double_vec_t gNorm = double_vec_new(g.ne); /* Norm of {g} vectors. */
        SPMatrix LL = SPFunction_GetCholeskyMatrix(g, dot, verbose, B, L);
        int i;
        
        for (i = 0; i < g.ne; i++)
          { gNorm.e[i] = sqrt(fabs(dot(g.e[i],g.e[i]))); }

        if (verbose) { fprintf(stderr, "s"); }
        for (i = 0; i < f.ne; i++)
          { SPFunction *fi = f.e[i];
            SPFunction *hi = NULL;
            double hiMag;
            int maxIter = 3, iter;
            /* Compute coeffs {x[i]} of projection, rel basis {g}: */
            if (verbose) { fprintf(stderr, "."); }
            SPBasisMatrices_RecomputeVectorGen(fi, g, dot, y, FALSE);
            SPSys_CholeskySolve(LL, y, t, x, FALSE);
            hi = SPFunction_LinComb(x.e, g);
            /* Now do some residual correction: */
            for (iter = 1; iter <= maxIter; iter++)
              { double dhiMax = 0.0;
                int j;
                bool_t converged = TRUE;
                /* Find correction {dx[j]} for {x[j]}: */
                SPBasisMatrices_RecomputeVectorGen(hi, g, dot, t, FALSE);
                converged = FALSE;
                for (j = 0; j < g.ne; j++) 
                  { double yj = y.e[j], tj = t.e[j];
                    double dyj = yj - tj;
                    if (fabs(dyj) > 1.0e-13*(fabs(yj) + fabs(tj)))
                      { converged = FALSE; } 
                    dy.e[j] = dyj;
                  }
                if (converged) { break; } 
                SPSys_CholeskySolve(LL, dy, t, dx, FALSE);
                /* Apply {dx} to {x}, recompute {h}, and check for convergence: */
                hi->m->scale(hi, 0.0);
                for (j = 0; j < g.ne; j++)
                  { double dxj = dx.e[j], xj = x.e[j];
                    if (dxj != 0.0)
                      { double dhi = fabs(dxj)*gNorm.e[j]; 
                        if (dhi > dhiMax) { dhiMax = dhi; }
                        xj += dxj; x.e[j] = xj;
                      }
                    if (xj != 0.0) { hi->m->add(hi, xj, g.e[j]); }
                  }
                if (dhiMax == 0.0) { break; }
                /* if (verbose) { fprintf(stderr, "(%8.2e)", dhiMax); } */
                if (iter == 1) { hiMag = sqrt(dot(hi,hi)); }
                if (dhiMax <= 1.0e-13*hiMag) { break; }
                if (iter >= maxIter)
                  { fprintf(stderr, "!");
                    if (verbose)
                      { fprintf(stderr, "(%8.2e:%8.2e)", dhiMax, hiMag); }
                  }
              }
            /* Save projection: */
            h.e[i] = hi;
          }
        if (verbose) { fprintf(stderr, " "); }
        if (L == NULL) { free(LL.ents.e); }
        free(gNorm.e);
        free(x.e); free(t.e); free(y.e);
        free(dx.e); free(dy.e);
      }
    return h;
  }

void SPFunction_GenOrthonizeBasis(Basis bas, Metric dot, bool_t verbose)
  { int i, j;
    for (i = 0; i < bas.ne; i++)
      { if (verbose) { fprintf(stderr, "-"); }
        SPFunction *bi = bas.e[i];
        double dii = SPFunction_GenNormalize(bi, dot);
        for (j = i+1; j < bas.ne; j++)
          { SPFunction *bj = bas.e[j];
            SPFunction_GenOrthize(bj, bi, dot, dii);
          }
      }
  }

void SPFunction_GenEigenFuncBasis
  ( Basis bas, 
    double dot(SPFunction *f, SPFunction *g),
    double mdot(SPFunction *f, SPFunction *g),
    double_vec_t ev,
    bool_t verbose
  )
  { int nZer = 0, i;
    affirm(ev.ne == bas.ne, "wrong vector length"); 
    /* Normalize basis elements; should be good for numerics. */ 
    for (i = 0; i < bas.ne; i++)
      { SPFunction_GenNormalize(bas.e[i], dot); }
    /* Check whether there are any zero-norm vectors, get them out of the way */
    for (i = 0; i < bas.ne; i++)
      { /* Elements {bas[0..nZer-1]} are the pariahs with zero {mdot} norm. */
        SPFunction *bi = bas.e[i];
        double dii = dot(bi, bi);
        double mdii = mdot(bi, bi);
        if (fabs(mdii) < 1.0e-8*fabs(dii)) 
          { /* Move {bas[i]} to the front of the basis: */
            bas.e[i] = bas.e[nZer];
            bas.e[nZer] = bi; ev.e[nZer] = 0.0;
            nZer++;
          }
        else
          { affirm(fabs(dii) > 0.5, "element has zero \"dot\" norm"); }
      }
    /* Orthogonalize the pariahs with plain Gram-Schmidt: */
    { Basis zer = (Basis) { /* nel */ nZer, /* el */ bas.e};
      SPFunction_GenOrthonizeBasis(zer, dot, verbose);
    }
    /* Now compute an eigenfunction basis for the good vectors: */
    { /* Fake a {Basis} descriptor for the good subset {sub} of {bas}: */
      int n = bas.ne - nZer;
      Basis sub = (Basis) { /* nel */ n, /* el */ &(bas.e[nZer])};
      double_vec_t B = double_vec_new(n*n); /* {B[i,j] = dot(bas[j],bas[i])}. */
      double_vec_t L = double_vec_new(n*n); /* Cholesky factor of {B}. */
      double_vec_t G = double_vec_new(n*n); /* {G[i,j] = mdot(bas[j],bas[i])}. */
      double_vec_t J = double_vec_new(n*n); /* Cholesky factor of {G}. */
      double_vec_t K = double_vec_new(n*n);
      double_vec_t M = double_vec_new(n*n); /* Relative metric matrix. */
      double_vec_t R = double_vec_new(n*n); /* Eigenvectors of {M}. */
      double_vec_t e = double_vec_new(n);
      double_vec_t d = double_vec_new(n);
      int neig;
      Basis nbas = SPFunctionRef_vec_new(n);

      /* Compute dot product matrices: */
      if (verbose) { fprintf(stderr, " d"); }
      for (i = 0; i < n; i++)
        { SPFunction *bi = sub.e[i];
          int j;
          for (j = 0; j <= i; j++)
            { SPFunction *bj = sub.e[j];
              double Bij = dot(bj, bi);
              double Gij = mdot(bj, bi);
              B.e[n*i + j] = Bij;
              G.e[n*i + j] = Gij;
              if (j < i) { B.e[n*j + i] = Bij; G.e[n*j + i] = Gij; }
            }
        }
      /* Compute Cholesky factorization {L*L^t} of {B}: */
      if (verbose) { fprintf(stderr, "c"); }
      rmxn_cholesky(n, B.e, L.e);
      rmxn_cholesky(n, G.e, J.e);
      /* Compute {M = (L^{-1})*G*(L^{-1})^t}: */
      if (verbose) { fprintf(stderr, "m"); }
      rmxn_LT_pre_div(n, n, L.e, J.e, K.e);
      rmxn_mul_tr(n, n, n, K.e, K.e, M.e);
      /* Find the eigenvectors of {M}: */
      if (n == 1)
        { R.e[0] = 1.0; d.e[0] = M.e[0]; neig = 1; }
      else
        { if (verbose) { fprintf(stderr, "e"); }
          syei_tridiagonalize(n, M.e, d.e, e.e, R.e);
          syei_trid_eigen(n, d.e, e.e, R.e, &neig, 1);
          affirm(neig == n, "eigenvectors of M could not be determined");
        }
      /* Convert them to basis coefficents: */
      if (verbose) { fprintf(stderr, "b"); }
      rmxn_LT_pos_div(n, n, R.e, L.e, B.e);
      /* Now form the functions: */
      if (verbose) { fprintf(stderr, "x"); }
      for (i = 0; i < n; i++)
        { double *a = &(B.e[n*i]);
          /*
          if (verbose)
            { int k; 
              fprintf(stderr, "    coefs = (");
              for (k = 0; k < n; k++) { fprintf(stderr, " %g", a[k]); }
              fprintf(stderr, " )\n");
            }
          */
          nbas.e[i] = SPFunction_LinComb(a, sub);
        }
      /* Store back: */
      if (verbose) { fprintf(stderr, " "); }
      for (i = 0; i < n; i++)
        { /* Should reclaim storage used by old {bas.e[i]}; */
          /* however we don't know whether it is safe to do so. */ 
          bas.e[nZer + i] = nbas.e[i];
          ev.e[nZer + i] = d.e[i];
        }
      /* Release temporary storage: */
      free(B.e); free(G.e);
      free(J.e); free(K.e); free(L.e);
      free(M.e); free(R.e); 
      free(d.e); free(e.e);
      free(nbas.e);
    }
  }

int SPFunction_GenCheckOrthonization
  ( Basis bas, 
    int ini, int lim, 
    Metric dot, 
    int maxTests, 
    double maxErr,
    int maxBad,
    bool_t verbose
  )
  { 
    int n = lim - ini; /* Size of submatrix to check. */
    if (verbose) { fprintf(stderr, "checking orthonormality...\n"); }
    int i, j;
    double diiWorst = 1.0; /* Worst diagonal element found do far. */
    double dijWorst = 0.0; /* Worst off-diagonal element found do far. */
    double fracTests = maxTests/(0.5*n*(n+1)); /* Obligation to test each element. */
    if (fracTests > 1.0) { fracTests = 1.0; }
    double dueTests = 0.0; /* Tests due minus tests done. */
    int nBad = 0; /* Number of tests that failed so far. */
    
    for (i = lim-1; i >= ini; i--)
      { 
        SPFunction *bi = bas.e[i];
        for (j = i; j < lim; j++)
          { 
            SPFunction *bj = bas.e[j];
            
            /* Decide whether this element is worth testing: */
            dueTests += fracTests;
            if (dueTests >= 0.999999)
              {
                if (verbose) { fprintf(stderr, "[%04d,%04d] ", i, j); }
                double dij = dot(bj, bi), err;
                if (i == j)
                  { /* Check normalization: */
                    err = fabs(dij - 1.0);
                    if (fabs(dij - 1.0) > fabs(diiWorst - 1.0)) { diiWorst = dij; }
                  }
                else
                  { /* Check orthogonality: */
                    err = fabs(dij);
                    if (fabs(dij) > fabs(dijWorst)) { dijWorst = dij; }
                  }
                if (err > maxErr)
                  { fprintf
                      ( stderr, 
                        " ** dot error <%04d|%04d> = %14.8e\n", 
                        i, j, dij
                      );
                    nBad++;
                    if (nBad > maxBad) { return nBad; }
                  }
                /* Account for one more test done, fix roundoff errors: */
                dueTests -= 1.0;
                if ((dueTests < 0.0) && (dueTests >= -0.000001)) { dueTests = 0.0; }
                if (verbose) { fprintf(stderr, "\n"); }
              }
          }
      }
    fprintf(stderr, "worst <b[i]|b[i]> = %14.8e\n", diiWorst);
    fprintf(stderr, "worst <b[j]|b[i]> = %14.8e\n", dijWorst);
    return nBad;
  }

/* NORMALIZATION */

void SPFunction_MakePositive
  ( SPFunction *f, 
    SPFunction *w, 
    Triangulation *tri
  )
  { double ff = SPFunction_Integral(f, NoFMap, w, tri, FALSE);
    if (ff < 0.0) { f->m->scale(f, -1.0); }
  }

double SPFunction_StdNormalize
  ( SPFunction *f, 
    SPFunction *w, 
    Triangulation *tri
  )
  { double ff = SPFunction_Dot(f, NoFMap, f, NoFMap, w, tri, FALSE);
    if (ff != 0.0) 
      { f->m->scale(f, 1.0/sqrt(fabs(ff))); return 1.0; }
    else
      { return 0.0; }
  }

void SPFunction_StdOrthize
  ( SPFunction *f, 
    SPFunction *g, 
    SPFunction *w, 
    Triangulation *tri,
    bool_t gNormal
  )
  { int iter;
    double fgMax = 0.0;
    for (iter = 1; iter <= 3; iter++)
      { double fg = SPFunction_Dot(f, NoFMap, g, NoFMap, w, tri, FALSE);
        if (fg == 0.0) { break; } 
        if (! gNormal)
          { double gg = SPFunction_Dot(g, NoFMap, g, NoFMap, w, tri, FALSE);
            fg = fg/gg; 
          }
        f->m->add(f, -fg, g);
        if ((iter > 1) && (fg < 1.0e-14*fgMax)) { break; } 
        if (fabs(fg) > fgMax) { fgMax = fabs(fg); }
      }
  }

void SPFunction_StdOrthonizeBasis
  ( Basis bas, 
    SPFunction *w, 
    Triangulation *tri
  )
  { int i, j;
    for (i = 0; i < bas.ne; i++)
      { SPFunction *bi = bas.e[i];
        for (j = 0; j < i; j++)
          { SPFunction *bj = bas.e[j];
            SPFunction_StdOrthize(bi, bj, w, tri, TRUE);
          }
        SPFunction_StdNormalize(bi, w, tri);
      }
  }

void SPFunction_SLapDotEigenFuncBasis
  ( Basis bas, 
    double_vec_t ev,
    SPFunction *w, 
    Triangulation *tri,
    bool_t verbose
  )
  { auto double EvalDot(SPFunction *f, SPFunction *g);
      /* Standard metric {<f|g>} */

    auto double SLapDot(SPFunction *f, SPFunction *g);
      /* Metric for eigen-analysis: {<SGrd(f)|SGrd(g)>} */

    double EvalDot(SPFunction *f, SPFunction *g)
      { return SPFunction_Dot(f, NoFMap, g, NoFMap, w, tri, FALSE); }

    double SLapDot(SPFunction *f, SPFunction *g)
      { return SPFunction_SLapDot(f, g, w, tri, FALSE); }

    SPFunction_GenEigenFuncBasis(bas, EvalDot, SLapDot, ev, verbose);
  }
  
int SPFunction_StdCheckOrthonization
  ( Basis bas, 
    SPFunction *w, 
    Triangulation *tri, 
    int maxTests, 
    double maxErr,
    int maxBad,
    bool_t verbose
  )
  { auto double EvalDot(SPFunction *f, SPFunction *g);
      /* Standard metric {<f|g>} */

    double EvalDot(SPFunction *f, SPFunction *g)
      { return SPFunction_Dot(f, NoFMap, g, NoFMap, w, tri, FALSE); }

    return SPFunction_GenCheckOrthonization
      ( bas, 0, bas.ne, EvalDot, maxTests, maxErr, maxBad, verbose );
  }

SPMatrix SPFunction_GetNormalMatrix
  ( Basis g, 
    Metric dot, 
    bool_t verbose,
    SPMatrix *B
  )
  { if (B != NULL)
      { return *B; }
    else
      { SPMatrix BB; 
        if (verbose) { fprintf(stderr, " d"); }
        BB = SPBasisMatrices_BuildMatrixGen(g, g, dot, NULL, FALSE, FALSE, TRUE, 0.0, FALSE);
        if (verbose) { fprintf(stderr, "(%d)", BB.ents.ne); }
        return BB;
      }
  }

SPMatrix SPFunction_GetCholeskyMatrix
  ( Basis g, 
    Metric dot, 
    bool_t verbose,
    SPMatrix *B,
    SPMatrix *L
  )
  { if (L != NULL)
      { return *L; }
    else
      { SPMatrix BB = SPFunction_GetNormalMatrix(g, dot, verbose, B); 
        SPMatrix LL;
        if (verbose) { fprintf(stderr, "c"); }
        SPMatrix_Cholesky(BB, 0.0, &LL);
        if (verbose) { fprintf(stderr, "(%d)", LL.ents.ne); }
        if (B == NULL) { free(BB.ents.e); }
        return LL;
      }
  }

/* Arrays of {SPFunction*}: */

vec_typeimpl(SPFunctionRef_vec_t,SPFunctionRef_vec,SPFunction*);
