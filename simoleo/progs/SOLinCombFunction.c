/* See SOLinCombFunction.h */
/* Last edited on 2023-02-12 07:51:06 by stolfi */

#include <SOLinCombFunction.h>
#include <SOFunction.h>

#include <dg_grid.h>

#include <filefmt.h>
#include <fget.h>
#include <nget.h>
#include <stdio.h>
#include <stdlib.h>
    
#define T SOLinCombFunction

void SOLinCombFunction_EvalMth(T *f, double *p, double *fp);
void SOLinCombFunction_GradMth(T *f, double *p, double *dp);
void SOLinCombFunction_HessMth(T *f, double *p, double *ddp);
void SOLinCombFunction_WriteMth(T *f, FILE *wr);
SOLinCombFunction *SOLinCombFunction_CopyMth(T *f);

SOLinCombFunction_Methods *SOLinCombFunction_Methods_New(void);
SOLinCombFunction_Data *SOLinCombFunction_Data_New(void);
SOLinCombFunction *SOLinCombFunction_New(void);

SOLinCombFunction *SOLinCombFunction_FullNew(void);
  /* Allocates a new {SOLinCombFunction} object {f}
    (and its data record), with the proper {type}
    fields and methods record.*/

SOLinCombFunction *SOLinCombFunction_Cast(OBJ *f)
  { SOFunction *ff = (SOFunction *)f;
    if ((f != NULL) && isprefix(SOLinCombFunction_TypeId, ff->type))
      { return (SOLinCombFunction *)f; }
    else
      { return NULL; }
  }

#define SOLinCombFunction_FileFormat "2003-05-06"

/* OVERRIDES FOR PARENT CLASS METHODS */
  
void SOLinCombFunction_EvalMth(T *f, double *p, double *fp) 
{ int i, j, k;
  int n = f->d->fn.fDim;
  Basis bas = f->d->bas;
  double *c = &(f->d->coef.el[0]);
  affirm(isprefix(SOLinCombFunction_TypeId, f->type), "type mismatch");
  for (j = 0; j < n; j++) { fp[j] = 0.0; }
  k = 0;
  for (i = 0; i < bas.nel; i++)
    { double bp;
      affirm(bas.el[i]->d->fDim == 1, "basis must be scalar");
      bas.el[i]->m->eval(bas.el[i], p, &bp);
      for (j = 0; j < n; j++,k++) { fp[j] += c[k] * bp; }
    }
  }

void SOLinCombFunction_GradMth(T *f, double *p, double *dp) 
  { affirm(isprefix(SOLinCombFunction_TypeId, f->type), "type mismatch");
    //dg_LinComb_grad(t, p, fp); // problema com coordenadas relativas!
    affirm(FALSE, "Grad not implemented");
  }
  
void SOLinCombFunction_HessMth(T *f, double *p, double *ddp) 
  { affirm(FALSE, "Hess not implemented"); }
  
T *SOLinCombFunction_CopyMth(T *f)
  { affirm(isprefix(SOLinCombFunction_TypeId, f->type), "type/method bug");
    SOLinCombFunction *g = SOLinCombFunction_FullNew();
    g->d->fn.pDim = f->d->fn.pDim;
    g->d->fn.fDim = f->d->fn.fDim;
    // ??? Deveria copiar a base ???
    g->d->bas = f->d->bas;
    g->d->basFile = f->d->basFile;
    { int N = f->d->coef.nel, i;
      double *gc = &(g->d->coef.el[0]), *fc = &(f->d->coef.el[0]);
      g->d->coef = double_vec_new(N);
      for (i = 0; i < N; i++) { gc[i] = fc[i]; }
    }
    return g;
  }

/* CLASS-SPECIFIC METHODS */
  
void SOLinCombFunction_WriteMth(T *f, FILE *wr)
  { int n = f->d->fn.fDim;
    affirm(isprefix(SOLinCombFunction_TypeId, f->type), "type mismatch");
    
    filefmt_write_header(wr, "SOLinCombFunction", SOLinCombFunction_FileFormat);
    fprintf(wr, "domain_dim = %d", f->d->fn.pDim); fputc('\n', wr);
    fprintf(wr, "range_dim = %d", n); fputc('\n', wr);
    fprintf(wr, "basis_size = %d", f->d->bas.nel); fputc('\n', wr);
    fprintf(wr, "basis = %s", f->d->basFile); fputc('\n', wr);
    { int i, j, k = 0;
      double *c = &(f->d->coef.el[0]);

      for (i = 0; i < f->d->bas.nel; i++)
        { for (j = 0; j < n; j++,k++)
            { if (j > 0) { fputc(' ', wr); }
	      fprintf(wr, "%.16g", c[k]);
            }
          fputc('\n', wr);
        }
    }
    filefmt_write_footer(wr, "SOLinCombFunction");
  }
    
/* OTHER PROCS */
  
SOLinCombFunction_Methods *SOLinCombFunction_Methods_New(void)
  { void *v = malloc(sizeof(SOLinCombFunction_Methods));
    return (SOLinCombFunction_Methods *)notnull(v, "no mem for SOLinCombFunction_Methods");
  }

SOLinCombFunction_Data *SOLinCombFunction_Data_New(void)
  { void *v = malloc(sizeof(SOLinCombFunction_Data));
    return (SOLinCombFunction_Data *)notnull(v, "no mem for SOLinCombFunction_Data");
  }

SOLinCombFunction *SOLinCombFunction_New(void)
  { void *v = malloc(sizeof(SOLinCombFunction));
    return (SOLinCombFunction *)notnull(v, "no mem for SOLinCombFunction");
  }

static SOLinCombFunction_Methods *LinCombMths = NULL;

SOLinCombFunction *SOLinCombFunction_FullNew(void)
  { SOLinCombFunction *f = SOLinCombFunction_New();
    f->type = SOLinCombFunction_TypeId;
    f->d = SOLinCombFunction_Data_New();
    if (LinCombMths == NULL)
      { LinCombMths = SOLinCombFunction_Methods_New();
        LinCombMths->fn.eval = (EvalMth *)&SOLinCombFunction_EvalMth;
        LinCombMths->fn.grad = (GradMth *)&SOLinCombFunction_GradMth;
        LinCombMths->fn.hess = (HessMth *)&SOLinCombFunction_HessMth;
	LinCombMths->fn.write = (WriteMth *)&SOFunction_WriteMth; /* NOT SOLinComb! */
        LinCombMths->fn.copy = (CopyMth *)&SOLinCombFunction_CopyMth;
        /* Class-specific methods */
        LinCombMths->write = (WriteMth *)&SOLinCombFunction_WriteMth;
      }
    f->m = LinCombMths;
    return f;
  }

SOLinCombFunction *SOLinCombFunction_Read(FILE *rd, dg_dim_t dDim)
  { int bSize;
    SOLinCombFunction *f = SOLinCombFunction_FullNew();
    
    filefmt_read_header(rd, "SOLinCombFunction", SOLinCombFunction_FileFormat);
    
    f->d->fn.pDim = nget_int32(rd, "domain_dim"); fget_eol(rd);
    f->d->fn.fDim = nget_int32(rd, "range_dim"); fget_eol(rd);
    bSize = nget_int32(rd, "basis_size"); fget_eol(rd);
    
    f->d->basFile = nget_string(rd, "basis"); fget_eol(rd);
    f->d->bas = SOFunction_ReadBasisCached(f->d->basFile);
    affirm(f->d->bas.nel == bSize, "inconsistent basis size");
    
    f->d->coef = double_vec_new(f->d->fn.fDim * bSize);
    
    // Read coefficients:
    { int i, j, k = 0;
      int n = f->d->fn.fDim;
      double *c = &(f->d->coef.el[0]);
      for (i = 0; i < bSize; i++)
        { for (j = 0; j < n; j++,k++) 
            { fget_skip_spaces(rd);
              c[k] = fget_double(rd);
            }
          fget_skip_formatting_chars(rd);
        }
    }

    filefmt_read_footer(rd, "SOLinCombFunction");
    return f;
  }

SOLinCombFunction *SOLinCombFunction_Make
  ( dg_dim_t pDim, 
    dg_dim_t fDim,
    char *basFile, 
    Basis bas,
    double_vec_t a
  )
  { SOLinCombFunction *g = SOLinCombFunction_FullNew();
    
    g->d->fn.pDim = pDim;
    g->d->fn.fDim = fDim;

    g->d->basFile = basFile;
    g->d->bas = bas;
    g->d->coef = a;
    
    return g;
  }

// double SOLinCombFunction_DotBoth(SOLinCombFunction *ftf, SOLinCombFunction *gtf)
//   { SOLinComb_Pair tp;   
//     dg_dim_t d;     /* Dimension of domain. */
//     double sum;       /* Accumulator for dot product. */
//     double corr;      /* Low-order bits of {*sum}. */
// 
//     tp.cx[0] = ftf->d->index; 
//     tp.cx[1] = gtf->d->index;  
// 
//     tp.r[0] = ftf->d->rank;
//     tp.r[1] = gtf->d->rank;
// 
//     affirm(ftf->d->fn.pDim == gtf->d->fn.pDim, "LinComb Functions domain dimension must be equal");
// 
//     d = ftf->d->fn.pDim;
// 
//     SOLinComb_tt_dot(&tp, d, &sum, &corr);
// 
//     return(sum);
//   }
// 
// 
// void SOLinCombFunction_DotSingle(SOLinCombFunction *tf, SOFunction *f, double *sum)
//   { SOLinComb t;   
//     dg_dim_t pd, fd;          /* Dimension of domain and (f) general function range */
//     double corr[f->d->fDim];     /* Low-order bits of {*sum}. */
//     int i;
// 
//     t.cx = tf->d->index;   
//     t.r = tf->d->rank;
//     pd = tf->d->fn.pDim;
//     fd = f->d->fDim;
//     for(i = 0; i < fd; i++){ sum[i] = 0.0; corr[i] = 0.0; }
// 
//     affirm(pd == f->d->pDim, "LinComb and General Functions domain dimension must be equal");
// 
//     auto void General_f(double *p, double *iv);
//       /* Encapsulates the (eval) method as a (SOIntegral_Func).*/
// 
//     void General_f(double *p, double *iv){f->m->eval((SOFunction *)f, p, iv);} 
//     
//     SOLinComb_tf_dot(&t, pd, General_f, fd, SOIntegral_GaussOrder, sum, corr);
// 
//   }
