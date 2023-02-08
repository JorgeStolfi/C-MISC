/* See SOWaveFunction.h */
/* Last edited on 2005-06-05 21:05:17 by stolfi */

#include <SOWaveFunction.h>
#include <SOFunction.h>

#include <filefmt.h>
#include <fget.h>
#include <nget.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <affirm.h>
#include <nat.h>

#define T SOWaveFunction 
    
void SOWaveFunction_WriteMth(T *f, FILE *wr);
SOWaveFunction *SOWaveFunction_CopyMth(T *f);

SOWaveFunction_Methods *SOWaveFunction_Methods_New(void);
SOWaveFunction_Data *SOWaveFunction_Data_New(void);
SOWaveFunction *SOWaveFunction_New(void);

SOWaveFunction *SOWaveFunction_Cast(OBJ *f)
  { SOFunction *ff = (SOFunction *)f;
    if ((f != NULL) && isprefix(SOWaveFunction_TypeId, ff->type))
      { return (SOWaveFunction *)f; }
    else
      { return NULL; }
  }

static SOWaveFunction_Methods *SOWaveFunction_Mths = NULL;

SOWaveFunction *SOWaveFunction_Make(void);
  /* Creates a new {SOWaveFunction} object (and its data record). 

    If {*mp != NULL}, uses {*mp} as the methods record, ignoring the
    other arguments. All instances of {SOWaveFunction}
    share the same methods vector {*SOWaveFunction_Mths},
    which is created if necessary. */

void SOWaveFunction_EvalMth(T *f, double *p, double *fp);
void SOWaveFunction_GradMth(T *f, double *p, double *dfp);
void SOWaveFunction_HessMth(T *f, double *p, double *ddfp);
  /* Methods of a {SOWaveFunction}. */

void SOWaveFunction_EvalMth(T *f, double *p, double *fp)
  { 
    int *k = &(f->d->freq[0]); /* Frequency vector. */
    double u = 0.0, s = 1.0;
    int i;
    
    affirm(f->d->fn.fDim == 1, "wrong fDim");
    for (i = 0; i < f->d->fn.pDim; i++)
      { u += ((double)k[i])*p[i]; }
    (*fp) = s * sin(TWOPI*u + PIBYFOUR);
  }
  
void SOWaveFunction_GradMth(T *f, double *p, double *dfp)
  {  
    int *k = &(f->d->freq[0]); /* Frequency vector. */
    double u = 0.0, s = 1.0;
    int i;
    
    affirm(f->d->fn.fDim == 1, "wrong fDim");
    for (i = 0; i < f->d->fn.pDim; i++)
      { double fri = (double)k[i]; u += fri*p[i]; dfp[i] = fri; }
    double ds = s * TWOPI * cos(TWOPI*u + PIBYFOUR);
    for (i = 0; i < f->d->fn.pDim; i++)
      { dfp[i] *= ds; }
  }

void SOWaveFunction_HessMth(T *f, double *p, double *ddfp)
  { affirm(FALSE, "not implemented yet"); }

#define SOWaveFunction_FileFormat "2003-07-03"

/* OVERRIDES FOR PARENT CLASS METHODS */
  
SOWaveFunction *SOWaveFunction_CopyMth(T *f)
  { SOWaveFunction *fnew = SOWaveFunction_New();
    fnew->type = f->type; 
    fnew->d = SOWaveFunction_Data_New();
    *(fnew->d) = *(f->d); /* Will copy {pDim,fDim,freq}. */
    fnew->m = f->m;       /* Methods vector is shared. */
    return fnew;
  }

/* CLASS-SPECIFIC METHODS */
  
void SOWaveFunction_WriteMth(T *f, FILE *wr)
  { int *k = &(f->d->freq[0]); /* Frequency vector. */
  int i;
  
    filefmt_write_header(wr, "SOWaveFunction", SOWaveFunction_FileFormat);

    fprintf(wr, "freq ="); 
    for (i = 0; i < f->d->fn.pDim; i++) { fprintf(wr, " %d", k[i]); }
    fputc('\n', wr);

    filefmt_write_footer(wr, "SOWaveFunction");
    fflush(wr);
  }

/* OTHER PROCS */
  
SOWaveFunction_Methods *SOWaveFunction_Methods_New(void)
  { void *v = malloc(sizeof(SOWaveFunction_Methods));
    return (SOWaveFunction_Methods *)notnull(v, "out of mem for SOWaveFunction_Methods");
  }

SOWaveFunction_Data *SOWaveFunction_Data_New(void)
  { void *v = malloc(sizeof(SOWaveFunction_Data));
    return (SOWaveFunction_Data *)notnull(v, "out of mem for SOWaveFunction_Data");
  }

SOWaveFunction *SOWaveFunction_New(void)
  { void *v = malloc(sizeof(SOWaveFunction));
    return (SOWaveFunction *)notnull(v, "no mem for SOWaveFunction");
  }

SOWaveFunction *SOWaveFunction_Make()
  {
    SOWaveFunction *f = SOWaveFunction_New();
    SOWaveFunction_Methods *m = SOWaveFunction_Mths;
    int i;

    f->type = SOWaveFunction_TypeId;
    if (m == NULL)
      { /* Creates the methods vector: */
        m = SOWaveFunction_Methods_New();
        /* Superclass methods: */
        m->fn.eval = (EvalMth *)&SOWaveFunction_EvalMth;
        m->fn.grad = (GradMth *)&SOWaveFunction_GradMth;
        m->fn.hess = (HessMth *)&SOWaveFunction_HessMth;
        m->fn.write = (WriteMth *)&SOFunction_WriteMth; /* NOT SOWave! */
        m->fn.copy = (CopyMth *)&SOWaveFunction_CopyMth;
        /* Class-specific methods */
        m->write = (WriteMth *)&SOWaveFunction_WriteMth;
        SOWaveFunction_Mths = m;
      }
    f->m = m;
    f->d = SOWaveFunction_Data_New();
    f->d->fn.fDim = 1;
    f->d->fn.pDim = 1; /* For now. */
    for (i = 0; i < MAX_PDIM; i++) { f->d->freq[i] = 0; } /* For now. */
    return f;
  }

SOWaveFunction *SOWaveFunction_FromFreq(nat_t pDim, int *freq)
  { SOWaveFunction *f = SOWaveFunction_Make();
    int i;
    f->d->fn.pDim = pDim;
    f->d->fn.fDim = 1;
    for (i = 0; i < pDim; i++) { f->d->freq[i] = freq[i]; }
    return f;
  }

SOWaveFunction *SOWaveFunction_Read(FILE *rd, nat_t pDim)
  { SOWaveFunction *f;
    int k[MAX_PDIM];
    int i;
    
    filefmt_read_header(rd, "SOWaveFunction", SOWaveFunction_FileFormat);
    
    nget_name_eq(rd, "freq"); 
    for (i = 0; i < pDim; i++) { k[i] = fget_int(rd); }
    fget_eol(rd);

    f = SOWaveFunction_FromFreq(pDim, k);
    filefmt_read_footer(rd, "SOWaveFunction");
    return f;
  }
