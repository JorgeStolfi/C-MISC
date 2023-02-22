/* See SOShiftedFunction.h */
/* Last edited on 2023-02-12 07:50:55 by stolfi */

#include <SOShiftedFunction.h>
#include <SOFunction.h>

#include <dg_grid.h>

#include <filefmt.h>
#include <fget.h>
#include <nget.h>
#include <stdio.h>
#include <stdlib.h>
    
#define T SOShiftedFunction

void SOShiftedFunction_EvalMth(T *f, double *p, double *fp);
void SOShiftedFunction_GradMth(T *f, double *p, double *dp);
void SOShiftedFunction_HessMth(T *f, double *p, double *ddp);
void SOShiftedFunction_WriteMth(T *f, FILE *wr);
SOShiftedFunction *SOShiftedFunction_CopyMth(T *f);

SOShiftedFunction_Methods *SOShiftedFunction_Methods_New(void);
SOShiftedFunction_Data *SOShiftedFunction_Data_New(void);
SOShiftedFunction *SOShiftedFunction_New(void);

SOShiftedFunction *SOShiftedFunction_FullNew(void);
  /* Allocates a new {SOShiftedFunction} object {f}
    (and its data record), with the proper {type}
    fields and methods record.*/

SOShiftedFunction *SOShiftedFunction_Cast(OBJ *f)
  { SOFunction *ff = (SOFunction *)f;
    if ((f != NULL) && isprefix(SOShiftedFunction_TypeId, ff->type))
      { return (SOShiftedFunction *)f; }
    else
      { return NULL; }
  }

#define SOShiftedFunction_FileFormat "2004-06-19"

/* OVERRIDES FOR PARENT CLASS METHODS */
  
void SOShiftedFunction_EvalMth(T *f, double *p, double *fp) 
  { dg_dim_t pDim = (dg_dim_t)f->d->fn.pDim;
    double q[pDim];
    interval_t *B = &(f->d->B[0]);
    SOFunction *g = f->d->g;
  
    affirm(isprefix(SOShiftedFunction_TypeId, f->type), "type mismatch");
    box_point_unmap(pDim, p, B, q);
    g->m->eval(g, q, fp);
  }

void SOShiftedFunction_GradMth(T *f, double *p, double *dp) 
  { dg_dim_t pDim = (dg_dim_t)f->d->fn.pDim;
    double q[pDim];
    interval_t *B = &(f->d->B[0]);
    SOFunction *g = f->d->g;
    int i;
  
    affirm(isprefix(SOShiftedFunction_TypeId, f->type), "type mismatch");
    box_point_unmap(pDim, p, B, q);
    g->m->grad(g, q, dp);
    /* Adjust gradient for domain compression: */
    for (i = 0; i < pDim; i++) 
      { dp[i] /= (HI(B[i]) - LO(B[i])); }
  }
  
void SOShiftedFunction_HessMth(T *f, double *p, double *ddp) 
  { affirm(FALSE, "Hess not implemented"); }
  
T *SOShiftedFunction_CopyMth(T *f)
  { affirm(isprefix(SOShiftedFunction_TypeId, f->type), "type/method bug");
    SOShiftedFunction *fnew = SOShiftedFunction_FullNew();
    fnew->d->fn.pDim = f->d->fn.pDim;
    fnew->d->fn.fDim = f->d->fn.fDim;
    fnew->d->g = f->d->g;
    int i;
    for (i = 0; i < f->d->fn.pDim; i++) { fnew->d->B[i] = f->d->B[i]; }
    return fnew;
  }

/* CLASS-SPECIFIC METHODS */
  
void SOShiftedFunction_WriteMth(T *f, FILE *wr)
  { affirm(isprefix(SOShiftedFunction_TypeId, f->type), "type mismatch");
    interval_t *B = &(f->d->B[0]);
    SOFunction *g = f->d->g;
    
    filefmt_write_header(wr, "SOShiftedFunction", SOShiftedFunction_FileFormat);

    g->m->write(g, wr);
    
    int i;
    for (i = 0; i < f->d->fn.pDim; i++) 
      { fprintf(wr, "%d %22.15e %22.15e", i, LO(B[i]), HI(B[i])); 
        fputc('\n', wr);
      }

    filefmt_write_footer(wr, "SOShiftedFunction");
  }
    
/* OTHER PROCS */
  
SOShiftedFunction_Methods *SOShiftedFunction_Methods_New(void)
  { void *v = malloc(sizeof(SOShiftedFunction_Methods));
    return (SOShiftedFunction_Methods *)notnull(v, "no mem for SOShiftedFunction_Methods");
  }

SOShiftedFunction_Data *SOShiftedFunction_Data_New(void)
  { void *v = malloc(sizeof(SOShiftedFunction_Data));
    return (SOShiftedFunction_Data *)notnull(v, "no mem for SOShiftedFunction_Data");
  }

SOShiftedFunction *SOShiftedFunction_New(void)
  { void *v = malloc(sizeof(SOShiftedFunction));
    return (SOShiftedFunction *)notnull(v, "no mem for SOShiftedFunction");
  }

static SOShiftedFunction_Methods *ShiftedMths = NULL;

SOShiftedFunction *SOShiftedFunction_FullNew(void)
  { SOShiftedFunction *f = SOShiftedFunction_New();
    f->type = SOShiftedFunction_TypeId;
    f->d = SOShiftedFunction_Data_New();
    if (ShiftedMths == NULL)
      { ShiftedMths = SOShiftedFunction_Methods_New();
        ShiftedMths->fn.eval = (EvalMth *)&SOShiftedFunction_EvalMth;
        ShiftedMths->fn.grad = (GradMth *)&SOShiftedFunction_GradMth;
        ShiftedMths->fn.hess = (HessMth *)&SOShiftedFunction_HessMth;
	ShiftedMths->fn.write = (WriteMth *)&SOFunction_WriteMth; /* NOT SOShifted! */
        ShiftedMths->fn.copy = (CopyMth *)&SOShiftedFunction_CopyMth;
        /* Class-specific methods */
        ShiftedMths->write = (WriteMth *)&SOShiftedFunction_WriteMth;
      }
    f->m = ShiftedMths;
    return f;
  }

SOShiftedFunction *SOShiftedFunction_Read(FILE *rd)
  { SOShiftedFunction *f = SOShiftedFunction_FullNew();
    interval_t B[MAX_PDIM];
    SOFunction *g;
    
    filefmt_read_header(rd, "SOShiftedFunction", SOShiftedFunction_FileFormat);
    
    g = SOFunction_Read(rd);
    
    int i;
    for (i = 0; i < MAX_PDIM; i++) 
      { if (i < f->d->fn.pDim)
          { int ir = fget_int32(rd); 
            affirm(i == ir, "bad seq in file");
            LO(B[i]) = fget_double(rd);
            HI(B[i]) = fget_double(rd);
          }
        else
          { LO(B[i]) = 0.0; HI(B[i]) = 1.0; }
        fget_eol(rd);
      }

    filefmt_read_footer(rd, "SOShiftedFunction");
    return SOShiftedFunction_Make(g, B);
  }
  
SOShiftedFunction *SOShiftedFunction_Make(SOFunction *g, interval_t B[])
  { SOShiftedFunction *f = SOShiftedFunction_FullNew();
    f->d->g = g;
    
    int i;
    for (i = 0; i < MAX_PDIM; i++) 
      { if (i < g->d->pDim)
          { f->d->B[i] = B[i]; }
        else
          { LO(f->d->B[i]) = 0.0; HI(f->d->B[i]) = 1.0; }
      }

    return f;
  }
