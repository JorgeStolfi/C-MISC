/* See SOFunction.h */
/* Last edited on 2023-02-12 07:51:19 by stolfi */

#include <SOGrid.h>
#include <SOFunction.h>
#include <SOTentFunction.h>
#include <SOProcFunction.h>
#include <SOWaveFunction.h>
#include <SOShiftedFunction.h>
#include <SOErrorFunction.h>
#include <SOLinCombFunction.h>

//#include <SOIntegral.h>
//#include <SOMatrix.h>
//#include <SOApprox.h>
//#include <SOBasisMatrices.h>

#include <dg_grid.h>

#include <filefmt.h>
#include <nget.h>
#include <fget.h>
#include <math.h> 
#include <stdlib.h> 
#include <string.h> 
#include <rn.h>
#include <rmxn.h> 
//#include <sym_eigen.h> 
#include <affirm.h> 
#include <nat.h>

#define OnBox SOIntegral_OnBox

#define FourPi 12.56637061435917

#define T SOFunction

void SOFunction_WriteAsSubclass(SOFunction *f, FILE *wr);
  /* Casts {f} to its widest proper effective subclass, and 
    calls the corresponding {write} method (which must shadow,
    but must not override, the {write} method of {f}. */
    
void debugmat(int m, int n, char *name, double *A);
  /* Prints the {m} by {n} matrix {A}, linearized by rows, on {stderr}. */

/* IMPLEMENTATIONS */

void debugmat(int m, int n, char *name, double *A)
  { if (n != 5) { return; }
    fprintf(stderr, "\n%s:\n", name);
    rmxn_gen_print(stderr, m, n, A, "%16.9e", "", "", "\n", "  ", " ", "\n");
  }

SOFunction *SOFunction_Cast(OBJ *f)
  { SOFunction *ff = (SOFunction *)f;
    if ((f != NULL) && isprefix(SOFunction_TypeId, ff->type))
      { return (SOFunction *)f; }
    else
      { return NULL; }
  }


#define SOFunction_FileFormat "2003-01-16"

void SOFunction_WriteMth(T *f, FILE *wr)
  { char *t = f->type;
    affirm(isprefix(SOFunction_TypeId, t), "type mismatch");
    filefmt_write_header(wr, "SOFunction", SOFunction_FileFormat);
    fprintf(wr, "type = %s\n", t);
    fprintf(wr, "domain_dim = %d\n", f->d->pDim);
    fprintf(wr, "range_dim = %d\n", f->d->fDim);
    SOFunction_WriteAsSubclass(f, wr);
    filefmt_write_footer(wr, "SOFunction");
  }
  
void SOFunction_WriteAsSubclass(SOFunction *f, FILE *wr)
  {
    { SOTentFunction *fp = SOTentFunction_Cast(f);
      if (fp != NULL) { fp->m->write(fp, wr); return; } 
    }
    { SOProcFunction *fp = SOProcFunction_Cast(f);
      if (fp != NULL) { fp->m->write(fp, wr); return; } 
    }
    { SOWaveFunction *fp = SOWaveFunction_Cast(f);
      if (fp != NULL) { fp->m->write(fp, wr); return; } 
    }
    { SOLinCombFunction *fp = SOLinCombFunction_Cast(f);
      if (fp != NULL) { fp->m->write(fp, wr); return; } 
    }
    { SOErrorFunction *fp = SOErrorFunction_Cast(f);
      if (fp != NULL) { fp->m->write(fp, wr); return; } 
    }
    { SOShiftedFunction *fp = SOShiftedFunction_Cast(f);
      if (fp != NULL) { fp->m->write(fp, wr); return; } 
    }
    /* 
    { SOSpline *fp = SOSpline_Cast(f);
      if (fp != NULL) { fp->m->write(fp, wr); return; } 
    }
    { SOBezFunction *fp = SOBezFunction_Cast(f);
      if (fp != NULL) { fp->m->write(fp, wr); return; } 
    }
    */
    fprintf(stderr, "unknown SOFunction subclass \"%s\"\n", f->type);
    affirm(FALSE, "aborted");
  }

SOFunction *SOFunction_Read(FILE *rd)
  { char *t;
    SOFunction *f;
    int pDim, fDim;
    filefmt_read_header(rd, "SOFunction", SOFunction_FileFormat);

    t = nget_string(rd, "type"); fget_eol(rd);
    affirm(isprefix(SOFunction_TypeId, t), "incorrect SOFunction type");
    pDim = nget_int32(rd, "domain_dim"); fget_eol(rd);
    affirm((pDim >= 1) && (pDim <= MAX_PDIM), "invalid domain dimension");
    fDim = nget_int32(rd, "range_dim"); fget_eol(rd);
    affirm((pDim >= 1) && (pDim <= MAX_FDIM), "invalid range dimension");
    /* Find the concrete subclass and call the corresponding {Read} function: */
    if (isprefix(SOTentFunction_TypeId, t))
      { f = (SOFunction *)SOTentFunction_Read(rd, pDim); }
    else if (isprefix(SOProcFunction_TypeId, t))
      { f = (SOFunction *)SOProcFunction_Read(rd, pDim, fDim); } 
    else if (isprefix(SOWaveFunction_TypeId, t))
      { f = (SOFunction *)SOWaveFunction_Read(rd, pDim); } 
    else if (isprefix(SOErrorFunction_TypeId, t))
      { f = (SOFunction *)SOErrorFunction_Read(rd, pDim); } 
    else if (isprefix(SOLinCombFunction_TypeId, t))
      { f = (SOFunction *)SOLinCombFunction_Read(rd, pDim);}
    else if (isprefix(SOShiftedFunction_TypeId, t))
      { f = (SOFunction *)SOShiftedFunction_Read(rd);}
    else
      { fprintf(stderr, "unknown SOFunction subclass \"%s\"\n", t);
        affirm(FALSE, "aborted");
      }

    fget_skip_formatting_chars(rd);

    filefmt_read_footer(rd, "SOFunction");
    free(t);
    return f;
  }

#define SOFunction_Basis_FileFormat "2003-01-16"
  
void SOFunction_WriteBasis(FILE *wr, Basis F)
  { int i;
    filefmt_write_header(wr, "SOFunction_Basis", SOFunction_Basis_FileFormat);
    fprintf(wr, "basis_dim = %d\n", F.nel);
    for (i = 0; i < F.nel; i++)
      { SOFunction *f = F.el[i]; 
        fprintf(wr, "elem_index = %d\n", i);
        f->m->write(f, wr); fputc('\n', wr);
      }
    filefmt_write_footer(wr, "SOFunction_Basis");
    fflush(wr);
  }

Basis SOFunction_ReadBasis(FILE *rd)
  { int dim;
    Basis F;
    int i;
    filefmt_read_header(rd, "SOFunction_Basis", SOFunction_Basis_FileFormat);
    dim = nget_int32(rd, "basis_dim"); fget_eol(rd);
    affirm(dim >= 0, "bad dimension");
    F = SOFunctionRef_vec_new(dim);
    fget_skip_formatting_chars(rd);
    for (i = 0; i < dim; i++) 
      { int ii = nget_int32(rd, "elem_index"); 
        affirm(ii == i, "wrong element index");
        fget_eol(rd);
        F.el[i] = SOFunction_Read(rd); fget_skip_formatting_chars(rd);
      }
    filefmt_read_footer(rd, "SOFunction_Basis");
    return F;
  }

static char **nameCache = NULL;
static Basis *basCache = NULL;
static nat_t szCache = 0; /* Allocated size of {nameCache}, {basCache}. */
static nat_t nCached = 0; /* Occupied entries in {nameCache}, {basCache}. */

Basis SOFunction_ReadBasisCached(char *fileName)
  { int i;
    Basis bas;
    /* Look up in cache */
    for (i = 0; i < nCached; i++)
      { if (strcmp(fileName, nameCache[i]) == 0) { return basCache[i]; }
      }
    /* Not found, read from file: */
    FILE *rd = open_read(fileName, TRUE);
    bas = SOFunction_ReadBasis(rd);
    if (rd != stdin) { fclose(rd); }
    /* Ensure there is space in cache: */
    if (nCached >= szCache)
      { /* Expand cache: */
        nat_t sztmp = 2*szCache + 10;
        char **tmpc = (char **)notnull(malloc(sztmp*sizeof(char*)), "no mem for tmpc");         
        Basis *tmpb = (Basis *)notnull(malloc(sztmp*sizeof(Basis)), "no mem for tmpb");
        for (i = 0; i < nCached; i++) 
          { tmpc[i] = nameCache[i]; tmpb[i] = basCache[i]; }
        if (nameCache != NULL) { free(nameCache); }
        if (basCache != NULL) { free(basCache); }
        nameCache = tmpc;
        basCache = tmpb;
        szCache = sztmp;
      }
    /* Store in cache: */
    nameCache[nCached] = fileName;
    basCache[nCached] = bas;
    nCached++;
    return bas;
  }

/* INTEGRALS */

FuncMap IdentityMap (dg_dim_t uDim, dg_dim_t pDim)
  { FuncMap NoMap;
    NoMap.map = (FuncMapProc *)NULL;
    NoMap.uDim = uDim;
    NoMap.pDim = pDim;
    NoMap.vDim = uDim; /* Same as {uDim}. */
    return NoMap;
  }

void SOFunction_Integral
  ( SOFunction *f,
    FuncMap *FMap, 
    SOFunction *w,
    SOGrid_Tree *tree,
    bool_t verbose,
    double *rv
  )
  {
    static double fv[MAX_FDIM];
    //if (SOSpline_IsCompatible(f, tree))
    //{ /* A piecewise function (on {tree}, if given) versus a general function: */
    //  SOSpline *fpw = SOSpline_Cast(f);
    //  return SOSpline_IntegralPW(fpw, FMap, w, verbose);
    //}
    //else
    //{ /* Use a general integrator (on {tree}, if given): */
    SOFunction_RawIntegral(f, FMap, w, tree, verbose, fv);
    //}
  }

void SOFunction_RawIntegral
  ( SOFunction *f,
    FuncMap *FMap, 
    SOFunction *w,
    SOGrid_Tree *tree,
    bool_t verbose,
    double *rv
  )
  { dg_dim_t pDim = f->d->pDim; /* Dimension of domain. */
    dg_dim_t iDim = FMap->vDim; /* Dimension of integrand. */
    static double corr[MAX_FDIM]; /* Correction term for Kahan summation. */
    
    auto void Func(double *p, double *iv);
      /* Computes the integrand {I(p) = w(p)*FMap(f(p),p)}; 
        returns the result in {iv}. */

    void Func(double *p, double *iv) 
      { 
        static double fv[MAX_FDIM];
        double wv;
        /* Compute {iv = FMap(f(p),p)}: */
        if (FMap->map == NULL) 
          { f->m->eval(f, p, iv); }
        else
          { f->m->eval(f, p, fv);
            FMap->map(fv, p, iv);
          }
        /* Apply weight {w(p)}: */
        if(w != NULL) 
          { w->m->eval(w, p, &wv);
            rn_scale(iDim, wv, iv, iv);
          }
      }

    affirm(f != NULL, "null function"); 
    affirm(FMap != NULL, "null FuncMap"); 
    affirm(w->d->fDim == 1, "weight must be a scalar field");
    affirm(w->d->pDim == f->d->pDim, "weight/function domain mismatch");
    affirm(f->d->fDim == FMap->uDim, "function/FuncMap range mismatch"); 
    rn_zero(iDim, rv);
    rn_zero(iDim, corr);
    SOIntegral_OnRootCell(Func, pDim, iDim, tree, rv, corr);
  }

/* INNER PRODUCTS */

void SOFunction_Dot
  ( SOFunction *f,
    FuncMap *FMap, 
    SOFunction *g, 
    FuncMap *GMap,
    SOFunction *w,
    SOGrid_Tree *tree,
    double *iv,
    bool_t verbose
  )
  {
    if (isprefix(SOTentFunction_TypeId, f->type) && isprefix(SOTentFunction_TypeId, g->type))
      { /* Two piecewise functions with same grid ({tree}, if given): */
        SOTentFunction *ftf = SOTentFunction_Cast(f);
        SOTentFunction *gtf = SOTentFunction_Cast(g);
        *iv = SOTentFunction_DotBoth(ftf, gtf);
        return;
      }
    else 
    if (isprefix(SOTentFunction_TypeId, f->type))
      { /* A piecewise function (on {tree}, if given) versus a general function: */
        SOTentFunction *ftf = SOTentFunction_Cast(f);
        SOTentFunction_DotSingle(ftf, g, GMap, iv);
        return;
      }
    else 
    if (isprefix(SOTentFunction_TypeId, g->type))
      { /* A piecewise function (on {tree}, if given) versus a general function: */
        SOTentFunction *gtf = SOTentFunction_Cast(g);
        SOTentFunction_DotSingle(gtf, f, FMap, iv);
        return;
      }

    //    else
    //    { /* Use a general integrator (on {tree}, if given): */
    //      *iv = SOFunction_RawDot(f, FMap, g, GMap, w, tree, verbose);
    //      return;
    //    }
  }


void SOFunction_GradDot
  ( SOFunction *f,
    SOFunction *g, 
    SOFunction *w,
    SOGrid_Tree *tree,
    double *iv,
    bool_t verbose
  )
  { 
    if (isprefix(SOTentFunction_TypeId, f->type) && isprefix(SOTentFunction_TypeId, g->type))
      { /* Two piecewise functions with same grid ({tree}, if given): */
        SOTentFunction *ftf = SOTentFunction_Cast(f);
        SOTentFunction *gtf = SOTentFunction_Cast(g);
        *iv = SOTentFunction_GradDotBoth(ftf, gtf);
        return;
      }
/*     else  */
/*     if (isprefix(SOTentFunction_TypeId, f->type)) */
/*       { /\* A piecewise function (on {tree}, if given) versus a general function: *\/ */
/*         SOTentFunction *ftf = SOTentFunction_Cast(f); */
/*         SOTentFunction_GradDotSingle(ftf, g, iv); */
/*         return; */
/*       } */
/*     else  */
/*     if (isprefix(SOTentFunction_TypeId, g->type)) */
/*       { /\* A piecewise function (on {tree}, if given) versus a general function: *\/ */
/*         SOTentFunction *gtf = SOTentFunction_Cast(g); */
/*         SOTentFunction_GradDotSingle(gtf, f, iv); */
/*         return; */
/*       } */
/*     else */
/*     *iv = SOFunction_RawGradDot(f, g, w, tree, verbose); */
  }


double SOFunction_RawDot
  ( SOFunction *f,
    FuncMap *FMap, 
    SOFunction *g, 
    FuncMap *GMap,
    SOFunction *w,
    SOGrid_Tree *tree,
    bool_t verbose
  )
  { dg_dim_t pDim = f->d->pDim; /* Dim of domain. */
    dg_dim_t mDim = FMap->vDim; /* Dim of {FMap(f(p),p)} and {GMap(g(p),p)}. */
    dg_dim_t iDim = 1;          /* Dimension of integrand. */
    double fv[f->d->fDim], mfv[FMap->vDim];
    double gv[g->d->fDim], mgv[GMap->vDim];
    double rv, corr;
    
    auto void Func(double *p, double *iv);
      /* The integrand of the dot product, namely
        {I(p) = w(p)*vdot(FMap(f(p),p), GMap(g(p),p))} */

    void Func(double *p, double *iv) 
      { 
        if (FMap->map != NULL) 
          { f->m->eval(f, p, fv);
            FMap->map(fv, p, mfv);
          }
        else
          { f->m->eval(f, p, mfv); }
        
        if (GMap->map != NULL) 
          { g->m->eval(g, p, gv);
            GMap->map(gv, p, mgv);
          }
        else
          { g->m->eval(g, p, mgv); }

        iv[0] = rn_dot(mDim, mfv, mgv);
        if(w != NULL) 
          { double wv;
            w->m->eval(w, p, &wv); 
            iv[0] *= wv;
          }
      }

    affirm(f != NULL, "null function f"); 
    affirm(FMap != NULL, "null FuncMap FMap");
    affirm(FMap->uDim == f->d->fDim, "function/FuncMap range mismatch");
    affirm(FMap->pDim == FMap->pDim, "function/FuncMap domain mismatch");
    
    affirm(g != NULL, "null function g"); 
    affirm(GMap != NULL, "null FuncMap GMap");
    affirm(GMap->uDim == g->d->fDim, "function/FuncMap range mismatch"); 
    affirm(GMap->pDim == GMap->pDim, "function/FuncMap domain mismatch");

    rv = 0.0;
    corr = 0.0;
    SOIntegral_OnRootCell(Func, pDim, iDim, tree, &rv, &corr);
    return rv;
  }

double SOFunction_RawGradDot
  ( SOFunction *f,
    SOFunction *g, 
    SOFunction *w,
    SOGrid_Tree *tree,
    bool_t verbose
  )
  { 
    auto void Func(double *p, double *iv);
      /* The integrand of the dot product. */

    void Func(double *p, double *iv) 
      { //double wp = (w == NULL ? 1.0 : w->m->eval(w, p));
        double wv;
        if(w == NULL) wv = 1.0; else w->m->eval(w, p, &wv);
        /* r3_t ftan, gtan, tmp; */
        /* R3Gradient fDp = f->m->grad(f, p); */
        /* R3Gradient gDp = g->m->grad(g, p); */
        /* return wp * r3_dot(&ftan, &gtan);  */
        iv[0] = 0;
      }

    // return SOIntegral_OnSphere(Func, tree);
    return (double)FALSE;// Not converted!
  }
 
/* Arrays of {SOFunction*}: */

vec_typeimpl(SOFunctionRef_vec_t,SOFunctionRef_vec,SOFunctionRef);
