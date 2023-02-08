/* See SOTentFunction.h */
/* Last edited on 2004-06-20 10:34:56 by stolfi */

#include <SOTentFunction.h>
#include <SOFunction.h>

#include <dg_grid.h>

#include <filefmt.h>
#include <fget.h>
#include <nget.h>
#include <stdio.h>
#include <stdlib.h>
    
#define T SOTentFunction

void SOTentFunction_EvalMth(T *f, double *p, double *fp);
void SOTentFunction_GradMth(T *f, double *p, double *dp);
void SOTentFunction_HessMth(T *f, double *p, double *ddp);
void SOTentFunction_WriteMth(T *f, FILE *wr);
SOTentFunction *SOTentFunction_CopyMth(T *f);

SOTentFunction_Methods *SOTentFunction_Methods_New(void);
SOTentFunction_Data *SOTentFunction_Data_New(void);
SOTentFunction *SOTentFunction_New(void);

SOTentFunction *SOTentFunction_FullNew(void);
  /* Allocates a new {SOTentFunction} object {f}
    (and its data record), with the proper {type}
    fields and methods record.*/

SOTentFunction *SOTentFunction_Cast(OBJ *f)
  { SOFunction *ff = (SOFunction *)f;
    if ((f != NULL) && isprefix(SOTentFunction_TypeId, ff->type))
      { return (SOTentFunction *)f; }
    else
      { return NULL; }
  }

#define SOTentFunction_FileFormat "2003-04-01"

/* OVERRIDES FOR PARENT CLASS METHODS */
  
void SOTentFunction_EvalMth(T *f, double *p, double *fp) 
  { interval_t box[f->d->fn.pDim];
    double tpoint[f->d->fn.pDim];
    dg_dim_t d = (dg_dim_t)f->d->fn.pDim;
    int i;
  
    affirm(isprefix(SOTentFunction_TypeId, f->type), "type mismatch");

    dg_cell_box_root_relative(d, f->d->index, box);

    for(i = 0; i < d; i++) 
      tpoint[i] = (double)((p[i] - LO(box[i])) / (HI(box[i]) - LO(box[i])));

    *fp = SOTent_eval(d, f->d->rank, tpoint);
  }
  /* Returns the SOTentFunction value relative to the (root cell) 
     considering its coordinates [0..1].*/

void SOTentFunction_GradMth(T *f, double *p, double *dp) 
  { affirm(isprefix(SOTentFunction_TypeId, f->type), "type mismatch");
    //dg_tent_grad(t, p, fp); // problema com coordenadas relativas!
    affirm(FALSE, "Grad not implemented");
  }
  
void SOTentFunction_HessMth(T *f, double *p, double *ddp) 
  { affirm(FALSE, "Hess not implemented"); }
  
T *SOTentFunction_CopyMth(T *f)
  { affirm(isprefix(SOTentFunction_TypeId, f->type), "type/method bug");
    SOTentFunction *fnew = SOTentFunction_FullNew();
    fnew->d->fn.pDim = f->d->fn.pDim;
    fnew->d->fn.fDim = f->d->fn.fDim;
    affirm(fnew->d->fn.fDim == 1, "tent functions should be scalar");
    fnew->d->index = f->d->index;
    fnew->d->rank = f->d->rank;
    return fnew;
  }

/* CLASS-SPECIFIC METHODS */
  
void SOTentFunction_WriteMth(T *f, FILE *wr)
  { affirm(isprefix(SOTentFunction_TypeId, f->type), "type mismatch");
    
    filefmt_write_header(wr, "SOTentFunction", SOTentFunction_FileFormat);

    fprintf(wr, "index = %llu", f->d->index); fputc('\n', wr);
    fprintf(wr, "rank = %d", f->d->rank); fputc('\n', wr);

    filefmt_write_footer(wr, "SOTentFunction");
  }
    
/* OTHER PROCS */
  
SOTentFunction_Methods *SOTentFunction_Methods_New(void)
  { void *v = malloc(sizeof(SOTentFunction_Methods));
    return (SOTentFunction_Methods *)notnull(v, "no mem for SOTentFunction_Methods");
  }

SOTentFunction_Data *SOTentFunction_Data_New(void)
  { void *v = malloc(sizeof(SOTentFunction_Data));
    return (SOTentFunction_Data *)notnull(v, "no mem for SOTentFunction_Data");
  }

SOTentFunction *SOTentFunction_New(void)
  { void *v = malloc(sizeof(SOTentFunction));
    return (SOTentFunction *)notnull(v, "no mem for SOTentFunction");
  }

static SOTentFunction_Methods *TentMths = NULL;

SOTentFunction *SOTentFunction_FullNew(void)
  { SOTentFunction *f = SOTentFunction_New();
    f->type = SOTentFunction_TypeId;
    f->d = SOTentFunction_Data_New();
    if (TentMths == NULL)
      { TentMths = SOTentFunction_Methods_New();
        TentMths->fn.eval = (EvalMth *)&SOTentFunction_EvalMth;
        TentMths->fn.grad = (GradMth *)&SOTentFunction_GradMth;
        TentMths->fn.hess = (HessMth *)&SOTentFunction_HessMth;
	TentMths->fn.write = (WriteMth *)&SOFunction_WriteMth; /* NOT SOTent! */
        TentMths->fn.copy = (CopyMth *)&SOTentFunction_CopyMth;
        /* Class-specific methods */
        TentMths->write = (WriteMth *)&SOTentFunction_WriteMth;
      }
    f->m = TentMths;
    return f;
  }

SOTentFunction *SOTentFunction_Read(FILE *rd, dg_dim_t dDim)
  { 
    SOTentFunction *f = SOTentFunction_FullNew();
    
    filefmt_read_header(rd, "SOTentFunction", SOTentFunction_FileFormat);
    
    f->d->fn.pDim = dDim;
    f->d->fn.fDim = 1;
  
    f->d->index = nget_int(rd, "index"); fget_eol(rd);
    f->d->rank = nget_int(rd, "rank"); fget_eol(rd);

    filefmt_read_footer(rd, "SOTentFunction");
    return f;
  }
  
SOTentFunction *SOTentFunction_FromHighBrick(dg_dim_t pDim, dg_cell_index_t k)
  { 
    SOTentFunction *f = SOTentFunction_FullNew();
    f->d->fn.pDim = pDim;
    f->d->fn.fDim = 1;
    f->d->index = k;
    f->d->rank = dg_cell_rank(k);
    return f;
  }
  
double SOTentFunction_DotBoth(SOTentFunction *ftf, SOTentFunction *gtf)
  { SOTent_Pair tp;   
    dg_dim_t d;             /* Dimension of domain. */
    double sum = 0;       /* Accumulator for dot product. */
    double corr = 0;      /* Low-order bits of {*sum}. */

    tp.cx[0] = ftf->d->index; 
    tp.cx[1] = gtf->d->index;  

    tp.r[0] = ftf->d->rank;
    tp.r[1] = gtf->d->rank;

    affirm(ftf->d->fn.pDim == gtf->d->fn.pDim, "Tent Functions domain dimension must be equal");

    d = ftf->d->fn.pDim;

    SOTent_tt_dot(&tp, d, &sum, &corr);

    return(sum);
  }

   
void SOTentFunction_DotSingle(SOTentFunction *tf, SOFunction *f, FuncMap *FMap, double *sum)
  { dg_dim_t pd, fd;                /* Dimension of domain and general function range */
    double corr[f->d->fDim];      /* Low-order bits of {*sum}. */
    interval_t box[f->d->pDim];  /* Coordinates of {tf's support} rel. root-cell. */
    double mink[f->d->pDim], maxk[f->d->pDim], diff[f->d->pDim];
    int i, k, max_cells = (1 << f->d->pDim);

    pd = f->d->pDim;
    if(FMap->map == NULL) { fd = f->d->fDim; } else { fd = FMap->vDim; }
    for(i = 0; i < fd; i++) { sum[i] = 0.0; corr[i] = 0.0; }

    affirm(pd == tf->d->fn.pDim, "domain dims mismatch");
   
    auto void integrand(double *var, double *fvar);
    /* Integrand of dot product: computes {h(x) = f(x)*t(x)},
       where {t} is the tent function, and {x} is a point within
       the root cell (defined as [0..1] along each coordinate). */

    void integrand(double *x, double *fx)
    { double tx;
      int j;
    
      /* Evaluate the function {f} at {var}: */
      if(FMap->map == NULL)
        { f->m->eval((SOFunction *)f, x, fx); } 
      else
        { double val[fd = f->d->fDim];
          f->m->eval((SOFunction *)f, x, val);
          FMap->map(val, x, fx);  
        }

      /* Evaluate the tent function {tf} at {var}: */
      tf->m->fn.eval((SOTentFunction *)tf, x, &tx);

      /* Multiply fvar by the tent value: */
      for(j = 0; j < fd; j++) { fx[j] *= tx; }
    }

    /* Find limits for the high cell */
    dg_cell_box_root_relative(pd, tf->d->index, box); 
    
    /* Get sizes of the high cell: */
    for(i = 0; i < pd; i++) { diff[i] = HI(box[i]) - LO(box[i]); }
  
    for(k = 0; k < max_cells; k++)
      {
        for(i = 0; i < pd; i++)
          {
            double shift = (double)(((k >> i)&1) - 1)*diff[i];
            mink[i] = LO(box[i]) + shift;    
            maxk[i] = HI(box[i]) + shift;
          }
        SOIntegral_Gauss_Limits(integrand, pd, fd, sum, corr, mink, maxk);
      }
    
  }

double SOTentFunction_GradDotBoth(SOTentFunction *ftf, SOTentFunction *gtf)
  { SOTent_Pair tp;   
    dg_dim_t d;         /* Dimension of domain. */
    double sum = 0;       /* Accumulator for dot product. */
    double corr = 0;      /* Low-order bits of {*sum}. */

    tp.cx[0] = ftf->d->index; 
    tp.cx[1] = gtf->d->index;  

    tp.r[0] = ftf->d->rank;
    tp.r[1] = gtf->d->rank;

    affirm(ftf->d->fn.pDim == gtf->d->fn.pDim, "domain dims mismatch");

    d = ftf->d->fn.pDim;

    SOTent_tt_grad_dot(&tp, d, &sum, &corr);

    return(sum);
  }
