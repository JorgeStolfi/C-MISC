/* See SOErrorFunction.h */
/* Last edited on 2023-02-12 07:51:24 by stolfi */

#include <SOErrorFunction.h>
#include <SOFunction.h>

#include <dg_grid.h>

#include <filefmt.h>
#include <fget.h>
#include <nget.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
    
#define T SOErrorFunction

void SOErrorFunction_EvalMth(T *f, double *p, double *fp);
void SOErrorFunction_GradMth(T *f, double *p, double *dp);
void SOErrorFunction_HessMth(T *f, double *p, double *ddp);
void SOErrorFunction_WriteMth(T *f, FILE *wr);
SOErrorFunction *SOErrorFunction_CopyMth(T *f);

SOErrorFunction_Methods *SOErrorFunction_Methods_New(void);
SOErrorFunction_Data *SOErrorFunction_Data_New(void);
SOErrorFunction *SOErrorFunction_New(void);

SOErrorFunction *SOErrorFunction_FullNew(void);
  /* Allocates a new {SOErrorFunction} object {f}
    (and its data record), with the proper {type}
    fields and methods record.*/

SOErrorFunction *SOErrorFunction_Cast(OBJ *f)
  { SOFunction *ff = (SOFunction *)f;
    if ((f != NULL) && isprefix(SOErrorFunction_TypeId, ff->type))
      { return (SOErrorFunction *)f; }
    else
      { return NULL; }
  }

#define SOErrorFunction_FileFormat "2003-06-02"

/* OVERRIDES FOR PARENT CLASS METHODS */
  
void SOErrorFunction_EvalMth(T *f, double *p, double *fp) 
  { int i, j, n = f->d->fn.fDim;
    int testdim, appdim;  
    double sum[n], testv[n], appv[n], diff[n];
    double *w = &(f->d->weight.el[0]);
  
    testdim = f->d->testbas.nel;
    appdim = f->d->appbas.nel;

    affirm((testdim == appdim), "test and approximation dimensions mismatch");

    for(j = 0; j < n; j++) sum[j] = 0.0;

    for(i = 0; i < testdim; i++)
      {
        f->d->testbas.el[i]->m->eval(f->d->testbas.el[i], p, testv);
        f->d->appbas.el[i]->m->eval(f->d->appbas.el[i], p, appv);
        for(j = 0; j < n; j++) 
          {
            diff[j] = testv[j] - appv[j];
            sum[j] = sum[j] + (double) w[i] * diff[j] * diff[j];
          }
      }
    for(j = 0; j < n; j++) fp[j] = sqrt(sum[j]);
  }

void SOErrorFunction_GradMth(T *f, double *p, double *dp) 
  { affirm(isprefix(SOErrorFunction_TypeId, f->type), "type mismatch");
    //dg_Error_grad(t, p, fp); 
    affirm(FALSE, "Grad not implemented");
  }
  
void SOErrorFunction_HessMth(T *f, double *p, double *ddp) 
  { affirm(FALSE, "Hess not implemented"); }
  
T *SOErrorFunction_CopyMth(T *f)
  { affirm(isprefix(SOErrorFunction_TypeId, f->type), "type/method bug");
    SOErrorFunction *g = SOErrorFunction_FullNew();
    g->d->fn.pDim = f->d->fn.pDim;
    g->d->fn.fDim = f->d->fn.fDim;

    g->d->testbas = f->d->testbas;
    g->d->appbas = f->d->appbas;

    g->d->testbasFile = f->d->testbasFile;
    g->d->appbasFile = f->d->appbasFile;

    { int N = f->d->weight.nel, i;
      double *gw = &(g->d->weight.el[0]), *fw = &(f->d->weight.el[0]);
      g->d->weight = double_vec_new(N);
      for (i = 0; i < N; i++) { gw[i] = fw[i]; }
    }
    return g;
  }

/* CLASS-SPECIFIC METHODS */
  
void SOErrorFunction_WriteMth(T *f, FILE *wr)
  { int n = f->d->fn.fDim;
    affirm(isprefix(SOErrorFunction_TypeId, f->type), "type mismatch");
    
    filefmt_write_header(wr, "SOErrorFunction", SOErrorFunction_FileFormat);
    fprintf(wr, "domain_dim = %d", f->d->fn.pDim); fputc('\n', wr);
    fprintf(wr, "range_dim = %d", n); fputc('\n', wr);

    fprintf(wr, "basis_size = %d", f->d->testbas.nel); fputc('\n', wr);

    fprintf(wr, "testbasis = %s", f->d->testbasFile); fputc('\n', wr);
    fprintf(wr, "appbasis = %s", f->d->appbasFile); fputc('\n', wr);

    { int i;
      double *w = &(f->d->weight.el[0]);

      for (i = 0; i < f->d->testbas.nel; i++)
        { fprintf(wr, "%.16g", w[i]);
          fputc('\n', wr);
        }
    }
    filefmt_write_footer(wr, "SOErrorFunction");
  }
    
/* OTHER PROCS */
  
SOErrorFunction_Methods *SOErrorFunction_Methods_New(void)
  { void *v = malloc(sizeof(SOErrorFunction_Methods));
    return (SOErrorFunction_Methods *)notnull(v, "no mem for SOErrorFunction_Methods");
  }

SOErrorFunction_Data *SOErrorFunction_Data_New(void)
  { void *v = malloc(sizeof(SOErrorFunction_Data));
    return (SOErrorFunction_Data *)notnull(v, "no mem for SOErrorFunction_Data");
  }

SOErrorFunction *SOErrorFunction_New(void)
  { void *v = malloc(sizeof(SOErrorFunction));
    return (SOErrorFunction *)notnull(v, "no mem for SOErrorFunction");
  }

static SOErrorFunction_Methods *ErrorMths = NULL;

SOErrorFunction *SOErrorFunction_FullNew(void)
  { SOErrorFunction *f = SOErrorFunction_New();
    f->type = SOErrorFunction_TypeId;
    f->d = SOErrorFunction_Data_New();
    if (ErrorMths == NULL)
      { ErrorMths = SOErrorFunction_Methods_New();
        ErrorMths->fn.eval = (EvalMth *)&SOErrorFunction_EvalMth;
        ErrorMths->fn.grad = (GradMth *)&SOErrorFunction_GradMth;
        ErrorMths->fn.hess = (HessMth *)&SOErrorFunction_HessMth;
	ErrorMths->fn.write = (WriteMth *)&SOFunction_WriteMth; /* NOT SOError! */
        ErrorMths->fn.copy = (CopyMth *)&SOErrorFunction_CopyMth;
        /* Class-specific methods */
        ErrorMths->write = (WriteMth *)&SOErrorFunction_WriteMth;
      }
    f->m = ErrorMths;
    return f;
  }

SOErrorFunction *SOErrorFunction_Read(FILE *rd, dg_dim_t dDim)
  { int bSize;
    SOErrorFunction *f = SOErrorFunction_FullNew();
    
    filefmt_read_header(rd, "SOErrorFunction", SOErrorFunction_FileFormat);
    
    f->d->fn.pDim = nget_int32(rd, "domain_dim"); fget_eol(rd);
    f->d->fn.fDim = nget_int32(rd, "range_dim"); fget_eol(rd);
    bSize = nget_int32(rd, "basis_size"); fget_eol(rd);
    
    f->d->testbasFile = nget_string(rd, "testbasis"); fget_eol(rd);
    f->d->appbasFile = nget_string(rd, "appbasis"); fget_eol(rd);

    f->d->testbas = SOFunction_ReadBasisCached(f->d->testbasFile);
    f->d->appbas = SOFunction_ReadBasisCached(f->d->appbasFile);

    affirm(f->d->testbas.nel == bSize, "inconsistent test basis size");
    affirm(f->d->appbas.nel == bSize, "inconsistent approximation basis size");
    
    f->d->weight = double_vec_new(bSize);
    
    // Read weights:
    { int i;
      double *w = &(f->d->weight.el[0]);
      for (i = 0; i < bSize; i++)
        { w[i] = fget_double(rd);
          fget_skip_formatting_chars(rd);
        }
    }

    filefmt_read_footer(rd, "SOErrorFunction");
    return f;
  }

SOErrorFunction *SOErrorFunction_Make
  ( dg_dim_t pDim, 
    dg_dim_t fDim,
    char *testbasFile,   /* Name of test basis file (with extension). */
    char *appbasFile,    /* Name of approximation basis file (with extension). */
    Basis testbas,
    Basis appbas,
    double_vec_t w
  )
  { SOErrorFunction *g = SOErrorFunction_FullNew();
    
    g->d->fn.pDim = pDim;
    g->d->fn.fDim = fDim;

    g->d->testbasFile = testbasFile;
    g->d->appbasFile = appbasFile;

    g->d->testbas = testbas;
    g->d->appbas = appbas;
    
    g->d->weight = w;

    return g;
  }
