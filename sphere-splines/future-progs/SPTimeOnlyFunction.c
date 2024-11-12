/* See SPTimeOnlyFunction.h */
/* Last edited on 2023-02-12 07:55:56 by stolfi */

#include <SPTimeOnlyFunction.h>
#include <SPTimeOnlyProcFunction.h>
#include <SPTimeOnlyIntegral.h>
#include <filefmt.h>
#include <nget.h>
#include <fget.h>
#include <math.h> 
#include <stdlib.h> 
#include <affirm.h> 
#include <nat.h>
#include <bool.h> 

#define T SPTimeOnlyFunction

void SPTimeOnlyFunction_WriteAsSubclass(SPTimeOnlyFunction *f, FILE *wr);
  /* Casts {f} to its widest proper effective subclass, and 
    calls the corresponding {write} method (which must shadow,
    but must not override, the {write} method of {f}. */
    
void debugmat(int m, int n, char *name, double *A);
  /* Prints the {m} by {n} matrix {A} on {stderr}. */

/* IMPLEMENTATIONS */

SPTimeOnlyFunction *SPTimeOnlyFunction_Cast(OBJ *f)
  { SPTimeOnlyFunction *ff = (SPTimeOnlyFunction *)f;
    if ((f != NULL) && isprefix(SPTimeOnlyFunction_TypeId, ff->type))
      { return (SPTimeOnlyFunction *)f; }
    else
      { return NULL; }
  }

#define SPTimeOnlyFunction_FileFormat "2005-08-17"

void SPTimeOnlyFunction_M_Write(T *f, FILE *wr)
  { char *t = f->type;
    affirm(isprefix(SPTimeOnlyFunction_TypeId, t), "type mismatch");
    filefmt_write_header(wr, "SPTimeOnlyFunction", SPTimeOnlyFunction_FileFormat);
    fprintf(wr, "type = %s\n", t);
    SPTimeOnlyFunction_WriteAsSubclass(f, wr);
    filefmt_write_footer(wr, "SPTimeOnlyFunction");
    fflush(wr);
  }
  
void SPTimeOnlyFunction_WriteAsSubclass(SPTimeOnlyFunction *f, FILE *wr)
  {
    { SPTimeOnlyProcFunction *fp = SPTimeOnlyProcFunction_Cast(f);
      if (fp != NULL) { fp->m->write(fp, wr); return; } 
    }
    fprintf(stderr, "unknown SPTimeOnlyFunction subclass \"%s\"\n", f->type);
    affirm(FALSE, "aborted");
  }

SPTimeOnlyFunction *SPTimeOnlyFunction_Read(FILE *rd)
  { char *t;
    SPTimeOnlyFunction *f;
    filefmt_read_header(rd, "SPTimeOnlyFunction", SPTimeOnlyFunction_FileFormat);
    t = nget_string(rd, "type"); fget_eol(rd);
    affirm(isprefix(SPTimeOnlyFunction_TypeId, t), "incorrect SPTimeOnlyFunction type");
    if (isprefix(SPTimeOnlyProcFunction_TypeId, t))
      { f = (SPTimeOnlyFunction *)SPTimeOnlyProcFunction_Read(rd); }
    else 
      { fprintf(stderr, "unknown SPTimeOnlyFunction subclass \"%s\"\n", t);
        affirm(FALSE, "aborted");
      }
    fget_skip_formatting_chars(rd);
    filefmt_read_footer(rd, "SPTimeOnlyFunction");
    free(t);
    return f;
  }

#define SPTimeOnlyFunction_Basis_FileFormat "2003-08-17"
  
void SPTimeOnlyFunction_WriteBasis(FILE *wr, TBasis tbas)
  { filefmt_write_header(wr, "SPTimeOnlyFunction_Basis", SPTimeOnlyFunction_Basis_FileFormat);
    fprintf(wr, "timePulseKind = %c\n", udg_pulse_kind_to_char(tbas.tFam.pkind));
    fprintf(wr, "timeDegree = %d\n", tbas.tFam.g);
    fprintf(wr, "timeCont = %d\n", tbas.tFam.c);
    filefmt_write_footer(wr, "SPTimeOnlyFunction_Basis");
    fflush(wr);
  }

TBasis SPTimeOnlyFunction_ReadBasis(FILE *rd)
  { TBasis tbas;
    filefmt_read_header(rd, "SPTimeOnlyFunction_Basis", SPTimeOnlyFunction_Basis_FileFormat);
    char pkc = nget_char(rd, "timePulseKind"); fget_eol(rd);
    if (pkc == 'B')
      { tbas.tFam.pkind = udg_PK_B; }
    else if (pkc == 'H')
      { tbas.tFam.pkind = udg_PK_H; }
    else if (pkc == 'N')
      { tbas.tFam.pkind = udg_PK_N; }
    else 
      { affirm(FALSE, "bad pulse kind"); }
    tbas.tFam.g = nget_int32(rd, "timeDegree"); fget_eol(rd);
    tbas.tFam.c = nget_int32(rd, "timeCont"); fget_eol(rd);
    filefmt_read_footer(rd, "SPTimeOnlyFunction_Basis");
    return tbas;
  }

/* INTEGRALS */

double SPTimeOnlyFunction_Integral
  ( SPTimeOnlyFunction *f,
    TimeOnlyFuncMap FMap, 
    SPTimeOnlyFunction *w,
    int smpOrder,
    double tMin, 
    double tMax, 
    bool_t verbose
  )
  {    
    auto double func(double t);
      /* The integrand of the dot product. */

    double func(double t) 
      { double wv = (w == NULL ? 1.0 : w->m->eval(w, t));
        double fv = (f == NULL ? 0.0 : f->m->eval(f, t));
        if (FMap.map != NULL) { fv = FMap.map(fv, t); }
        return wv * fv;
      }

    return SPTimeOnlyIntegral_OnInterval(func, smpOrder, tMin, tMax);
  }

double SPTimeOnlyFunction_CustomIntegral
  ( SPTimeOnlyFunction *f, 
    TimeOnlyFuncMap FMap, 
    SPTimeOnlyFunction *w,
    double_vec_t st,
    double_vec_t wt
  ) 
  { 
    auto double func(double t);
    
    double func(double t)
      { double wv = (w == NULL ? 1.0 : w->m->eval(w, t));
        double fv = f->m->eval(f, t);
        if (FMap.map != NULL) { fv = FMap.map(fv, t); }
        return wv * fv;
      }

    double sum = 0.0, corr = 0.0;
    SPTimeOnlyIntegral_BySamples(func, st, wt, &sum, &corr);
    return sum;
  }

/* INNER PRODUCTS */

double SPTimeOnlyFunction_Dot
  ( SPTimeOnlyFunction *f,
    TimeOnlyFuncMap FMap, 
    SPTimeOnlyFunction *g, 
    TimeOnlyFuncMap GMap,
    SPTimeOnlyFunction *w,
    int smpOrder,
    double tMin, 
    double tMax, 
    bool_t verbose
  )
  {    
    auto double func(double t);
      /* The integrand of the dot product. */

    double func(double t) 
      { double wv = (w == NULL ? 1.0 : w->m->eval(w, t));
        double fv = (f == NULL ? 0.0 : f->m->eval(f, t));
        double gv = (g == NULL ? 0.0 : g->m->eval(g, t));
        if (FMap.map != NULL) { fv = FMap.map(fv, t); }
        if (GMap.map != NULL) { gv = GMap.map(gv, t); }
        return wv * fv * gv;
      }

    return SPTimeOnlyIntegral_OnInterval(func, smpOrder, tMin, tMax);
  }

double SPTimeOnlyFunction_CustomDot
  ( SPTimeOnlyFunction *f, 
    TimeOnlyFuncMap FMap, 
    SPTimeOnlyFunction *g,
    TimeOnlyFuncMap GMap,
    SPTimeOnlyFunction *w,
    double_vec_t st,
    double_vec_t wt
  ) 
  { 
    auto double func(double t);
    
    double func(double t)
      { double wv = (w == NULL ? 1.0 : w->m->eval(w, t));
        double fv = f->m->eval(f, t);
        double gv = g->m->eval(g, t);
        if (FMap.map != NULL) { fv = FMap.map(fv, t); }
        if (GMap.map != NULL) { gv = GMap.map(gv, t); }
        return wv * fv * gv;
      }

    double sum = 0.0, corr = 0.0;
    SPTimeOnlyIntegral_BySamples(func, st, wt, &sum, &corr);
    return sum;
  }
