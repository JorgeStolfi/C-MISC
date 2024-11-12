/* See SPTimeSpaceFunction.h */
/* Last edited on 2023-02-12 07:55:40 by stolfi */

#include <SPTimeSpaceFunction.h>
#include <SPTimeSpaceProcFunction.h>
#include <SPTriang.h>
#include <SPTimeSpaceIntegral.h>
#include <SPQuad.h>
#include <filefmt.h>
#include <nget.h>
#include <fget.h>
#include <math.h> 
#include <stdlib.h> 
#include <r3.h> 
#include <rmxn.h> 
#include <affirm.h> 
#include <nat.h>
#include <bool.h> 

#define OnTriangle SPTimeSpaceIntegral_OnTriangle

#define T SPTimeSpaceFunction

void SPTimeSpaceFunction_WriteAsSubclass(SPTimeSpaceFunction *f, FILE *wr);
  /* Casts {f} to its widest proper effective subclass, and 
    calls the corresponding {write} method (which must shadow,
    but must not override, the {write} method of {f}. */
    
void debugmat(int m, int n, char *name, double *A);
  /* Prints the {m} by {n} matrix {A} on {stderr}. */

/* IMPLEMENTATIONS */

SPTimeSpaceFunction *SPTimeSpaceFunction_Cast(OBJ *f)
  { SPTimeSpaceFunction *ff = (SPTimeSpaceFunction *)f;
    if ((f != NULL) && isprefix(SPTimeSpaceFunction_TypeId, ff->type))
      { return (SPTimeSpaceFunction *)f; }
    else
      { return NULL; }
  }

#define SPTimeSpaceFunction_FileFormat "2003-08-17"

void SPTimeSpaceFunction_M_Write(T *f, FILE *wr)
  { char *t = f->type;
    affirm(isprefix(SPTimeSpaceFunction_TypeId, t), "type mismatch");
    filefmt_write_header(wr, "SPTimeSpaceFunction", SPTimeSpaceFunction_FileFormat);
    fprintf(wr, "type = %s\n", t);
    SPTimeSpaceFunction_WriteAsSubclass(f, wr);
    filefmt_write_footer(wr, "SPTimeSpaceFunction");
    fflush(wr);
  }
  
void SPTimeSpaceFunction_WriteAsSubclass(SPTimeSpaceFunction *f, FILE *wr)
  {
    { SPTimeSpaceProcFunction *fp = SPTimeSpaceProcFunction_Cast(f);
      if (fp != NULL) { fp->m->write(fp, wr); return; } 
    }
    fprintf(stderr, "unknown SPTimeSpaceFunction subclass \"%s\"\n", f->type);
    affirm(FALSE, "aborted");
  }

SPTimeSpaceFunction *SPTimeSpaceFunction_Read(FILE *rd)
  { char *t;
    SPTimeSpaceFunction *f;
    filefmt_read_header(rd, "SPTimeSpaceFunction", SPTimeSpaceFunction_FileFormat);
    t = nget_string(rd, "type"); fget_eol(rd);
    affirm(isprefix(SPTimeSpaceFunction_TypeId, t), "incorrect SPTimeSpaceFunction type");
    if (isprefix(SPTimeSpaceProcFunction_TypeId, t))
      { f = (SPTimeSpaceFunction *)SPTimeSpaceProcFunction_Read(rd); }
    else 
      { fprintf(stderr, "unknown SPTimeSpaceFunction subclass \"%s\"\n", t);
        affirm(FALSE, "aborted");
      }
    fget_skip_formatting_chars(rd);
    filefmt_read_footer(rd, "SPTimeSpaceFunction");
    free(t);
    return f;
  }

#define SPTimeSpaceFunction_Basis_FileFormat "2003-08-17"
  
void SPTimeSpaceFunction_WriteBasis(FILE *wr, STBasis stbas)
  { filefmt_write_header(wr, "SPTimeSpaceFunction_Basis", SPTimeSpaceFunction_Basis_FileFormat);
    fprintf(wr, "timePulseKind = %c\n", udg_pulse_kind_to_char(stbas.tFam.pkind));
    fprintf(wr, "timeDegree = %d\n", stbas.tFam.g);
    fprintf(wr, "timeCont = %d\n", stbas.tFam.c);
    SPFunction_WriteBasis(wr, stbas.sBas);
    filefmt_write_footer(wr, "SPTimeSpaceFunction_Basis");
    fflush(wr);
  }

STBasis SPTimeSpaceFunction_ReadBasis(FILE *rd)
  { STBasis stbas;
    filefmt_read_header(rd, "SPTimeSpaceFunction_Basis", SPTimeSpaceFunction_Basis_FileFormat);
    char pkc = nget_char(rd, "timePulseKind"); fget_eol(rd);
    if (pkc == 'B')
      { stbas.tFam.pkind = udg_PK_B; }
    else if (pkc == 'H')
      { stbas.tFam.pkind = udg_PK_H; }
    else if (pkc == 'N')
      { stbas.tFam.pkind = udg_PK_N; }
    else 
      { affirm(FALSE, "bad pulse kind"); }
    stbas.tFam.g = nget_int32(rd, "timeDegree"); fget_eol(rd);
    stbas.tFam.c = nget_int32(rd, "timeCont"); fget_eol(rd);
    stbas.sBas = SPFunction_ReadBasis(rd);
    filefmt_read_footer(rd, "SPTimeSpaceFunction_Basis");
    return stbas;
  }

/* INTEGRALS */

double SPTimeSpaceFunction_Integral
  ( SPTimeSpaceFunction *f,
    TimeSpaceFuncMap FMap, 
    SPTimeSpaceFunction *w,
    Triangulation *tri,
    int smpOrderTime,
    double tMin, 
    double tMax, 
    bool_t verbose
  )
  {
    /* ??? SHOULD INTEGRATE OVER THE SPACE SUPPORT OF {f} */
    /* if (SPSpline_IsCompatible(f, tri) && FMap.zeroPres) {... } */ 
    /* Use a general integrator (on {tri}, if given): */
    return SPTimeSpaceFunction_RawIntegral(f, FMap, w, tri, smpOrderTime, tMin, tMax, verbose);
  }

double SPTimeSpaceFunction_RawIntegral
  ( SPTimeSpaceFunction *f,
    TimeSpaceFuncMap FMap, 
    SPTimeSpaceFunction *w,
    Triangulation *tri,
    int smpOrderTime,
    double tMin, 
    double tMax, 
    bool_t verbose
  )
  {    
    auto double func(S2Point *p, double t);
      /* The integrand of the dot product. */

    double func(S2Point *p, double t) 
      { double wv = (w == NULL ? 1.0 : w->m->eval(w, p, t));
        double fv = (f == NULL ? 0.0 : f->m->eval(f, p, t));
        if (FMap.map != NULL) { fv = FMap.map(fv, p, t); }
        return wv * fv;
      }

    return SPTimeSpaceIntegral_OnSphereInterval(func, tri, smpOrderTime, tMin, tMax);
  }

double SPTimeSpaceFunction_CustomIntegral
  ( SPTimeSpaceFunction *f, 
    TimeSpaceFuncMap FMap, 
    SPTimeSpaceFunction *w,
    S2Point_vec_t sp,
    double_vec_t wp,
    double_vec_t st,
    double_vec_t wt
  ) 
  { 
    auto double func(S2Point *p, double t);
    
    double func(S2Point *p, double t)
      { double wv = (w == NULL ? 1.0 : w->m->eval(w, p, t));
        double fv = f->m->eval(f, p, t);
        if (FMap.map != NULL) { fv = FMap.map(fv, p, t); }
        return wv * fv;
      }

    double sum = 0.0, corr = 0.0;
    SPTimeSpaceIntegral_BySamples(func, sp, wp, st, wt, &sum, &corr);
    return sum;
  }

/* INNER PRODUCTS */

double SPTimeSpaceFunction_Dot
  ( SPTimeSpaceFunction *f,
    TimeSpaceFuncMap FMap, 
    SPTimeSpaceFunction *g, 
    TimeSpaceFuncMap GMap,
    SPTimeSpaceFunction *w,
    Triangulation *tri,
    int smpOrderTime,
    double tMin, 
    double tMax, 
    bool_t verbose
  )
  {
    /* ??? SHOULD INTEGRATE ONLY OVER SUPPORT. */
    /* Use a general integrator (on {tri}, if given): */
    return SPTimeSpaceFunction_RawDot(f, FMap, g, GMap, w, tri, smpOrderTime, tMin, tMax, verbose);
  }

double SPTimeSpaceFunction_RawDot
  ( SPTimeSpaceFunction *f,
    TimeSpaceFuncMap FMap, 
    SPTimeSpaceFunction *g, 
    TimeSpaceFuncMap GMap,
    SPTimeSpaceFunction *w,
    Triangulation *tri,
    int smpOrderTime,
    double tMin, 
    double tMax, 
    bool_t verbose
  )
  {    
    auto double func(S2Point *p, double t);
      /* The integrand of the dot product. */

    double func(S2Point *p, double t) 
      { double wv = (w == NULL ? 1.0 : w->m->eval(w, p, t));
        double fv = (f == NULL ? 0.0 : f->m->eval(f, p, t));
        double gv = (g == NULL ? 0.0 : g->m->eval(g, p, t));
        if (FMap.map != NULL) { fv = FMap.map(fv, p, t); }
        if (GMap.map != NULL) { gv = GMap.map(gv, p, t); }
        return wv * fv * gv;
      }

    return SPTimeSpaceIntegral_OnSphereInterval(func, tri, smpOrderTime, tMin, tMax);
  }

double SPTimeSpaceFunction_CustomDot
  ( SPTimeSpaceFunction *f, 
    TimeSpaceFuncMap FMap, 
    SPTimeSpaceFunction *g,
    TimeSpaceFuncMap GMap,
    SPTimeSpaceFunction *w,
    S2Point_vec_t sp,
    double_vec_t wp,
    double_vec_t st,
    double_vec_t wt
  ) 
  { 
    auto double func(S2Point *p, double t);
    
    double func(S2Point *p, double t)
      { double wv = (w == NULL ? 1.0 : w->m->eval(w, p, t));
        double fv = f->m->eval(f, p, t);
        double gv = g->m->eval(g, p, t);
        if (FMap.map != NULL) { fv = FMap.map(fv, p, t); }
        if (GMap.map != NULL) { gv = GMap.map(gv, p, t); }
        return wv * fv * gv;
      }

    double sum = 0.0, corr = 0.0;
    SPTimeSpaceIntegral_BySamples(func, sp, wp, st, wt, &sum, &corr);
    return sum;
  }
