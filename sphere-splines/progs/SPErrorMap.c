/* See SPErrorMap.h */
/* Last edited on 2005-10-27 08:40:23 by stolfi */

#include <SPErrorMap.h>
#include <SPFunction.h>
#include <SPBasic.h>
#include <fget.h>
#include <nget.h>
#include <filefmt.h>
#include <affirm.h>
#include <nat.h>
#include <stdio.h>
#include <stdlib.h>

#define T SPErrorMap

double SPErrorMap_M_Eval(T *f, R3Point *p);
R3Gradient SPErrorMap_M_Grad(T *f, R3Point *p);
R3Hessian SPErrorMap_M_Hess(T *f, R3Point *p);
void SPErrorMap_M_Write(T *f, FILE *wr);
void SPErrorMap_M_Add(T *f, double a, T *h);
void SPErrorMap_M_Scale(T *f, double a);
void SPErrorMap_M_Maple(T *f, FILE *wr);
SPErrorMap *SPErrorMap_M_Copy(T *f);
void SPErrorMap_M_Free(T *f);

SPErrorMap_Methods *SPErrorMap_Methods_New(void);
SPErrorMap_Data *SPErrorMap_Data_New(void);
SPErrorMap *SPErrorMap_New(void);

SPErrorMap *SPErrorMap_FullNew(void);
  /* Allocates a new {SPErrorMap} object {f} (and its data
    record), with the proper {type} fields and methods record. The
    arrays {f->d->F} and {f->d->G} are left undefined. */

SPErrorMap *SPErrorMap_Cast(OBJ *f)
  { SPFunction *ff = (SPFunction *)f;
    if ((f != NULL) && isprefix(SPErrorMap_TypeId, ff->type))
      { return (SPErrorMap *)f; }
    else
      { return NULL; }
  }

#define SPErrorMap_FileFormat "2002-11-18"

/* OVERRIDES FOR PARENT CLASS METHODS */
  
double SPErrorMap_M_Eval(T *f, R3Point *p) 
  { double sum = 0.0, corr = 0.0;
    Basis *F = &(f->d->F);
    Basis *G = &(f->d->G);
    int i;
    for (i = 0; i < F->ne;  i++)
      { SPFunction *Fi = F->e[i];
        SPFunction *Gi = G->e[i];
        double ei = Fi->m->eval(Fi, p) - Gi->m->eval(Gi, p);
        double term = ei*ei;
        /* Kahan's summation formula: */
        double tcorr = term - corr;
        double newSum = sum + tcorr;
        corr = (newSum - sum) - tcorr;
        sum = newSum;
      }
    return sum;
  }

R3Gradient SPErrorMap_M_Grad(T *f, R3Point *p) 
  { r3_t sum = (r3_t){{0.0, 0.0, 0.0}};
    r3_t corr = (r3_t){{0.0, 0.0, 0.0}};
    Basis *F = &(f->d->F);
    Basis *G = &(f->d->G);
    int i;
    for (i = 0; i < F->ne;  i++)
      { SPFunction *Fi = F->e[i];
        SPFunction *Gi = G->e[i];
        double ei = Fi->m->eval(Fi, p) - Gi->m->eval(Gi, p);
        r3_t DFi = Fi->m->grad(Fi, p);
        r3_t DGi = Gi->m->grad(Gi, p);
        r3_t Dei;
        r3_sub(&DFi, &DGi, &Dei);
        int r;
        for (r = 0; r < 3; r++)
          { double term = 2.0*ei*Dei.c[r];
            /* Kahan's summation formula: */
            double tcorr = term - corr.c[r];
            double newSum = sum.c[r] + tcorr;
            corr.c[r] = (newSum - sum.c[r]) - tcorr;
            sum.c[r] = newSum;
          }
      }
    return sum;
  }

R3Hessian SPErrorMap_M_Hess(T *f, R3Point *p) 
  { r6_t sum = (r6_t){{0.0, 0.0, 0.0, 0.0, 0.0, 0.0}};
    r6_t corr = (r6_t){{0.0, 0.0, 0.0, 0.0, 0.0, 0.0}};
    Basis *F = &(f->d->F);
    Basis *G = &(f->d->G);
    int i;
    for (i = 0; i < F->ne;  i++)
      { SPFunction *Fi = F->e[i];
        SPFunction *Gi = G->e[i];
        double ei = Fi->m->eval(Fi, p) - Gi->m->eval(Gi, p);
        r3_t DFi = Fi->m->grad(Fi, p);
        r3_t DGi = Gi->m->grad(Gi, p);
        r3_t Dei;
        r3_sub(&DFi, &DGi, &Dei);
        r6_t DDFi = Fi->m->hess(Fi, p);
        r6_t DDGi = Gi->m->hess(Gi, p);
        r6_t DDei;
        r6_sub(&DDFi, &DDGi, &DDei);
        int t = 0;
        int r, s;
        for (s = 0; s < 3; s++)
          for (r = 0; r <= r; r++)
            { double term = 2.0*(Dei.c[s]*Dei.c[r] + ei*DDei.c[t]);
              /* Kahan's summation formula: */
              double tcorr = term - corr.c[t];
              double newSum = sum.c[t] + tcorr;
              corr.c[t] = (newSum - sum.c[t]) - tcorr;
              sum.c[t] = newSum;
              
              t++;
            }
      }
    return sum;
  }

void SPErrorMap_M_Maple(T *f, FILE *wr)
  { affirm(FALSE , "maple method not implemented yet"); }

void SPErrorMap_M_Add(T *f, double a, T *h)
  { affirm(FALSE , "add method is undefined"); }

void SPErrorMap_M_Scale(T *f, double a)
  { Basis *G = &(f->d->G);
    Basis *F = &(f->d->F);
    int i;
    for (i = 0; i < F->ne; i++)
      { SPFunction *Fi = F->e[i];
        SPFunction *Gi = G->e[i];
        Fi->m->scale(Fi, a);
        Gi->m->scale(Gi, a); 
      }
  }

T *SPErrorMap_M_Copy(T *f)
  { affirm(isprefix(SPErrorMap_TypeId, f->type), "type/method bug");
    { SPErrorMap *g = SPErrorMap_FullNew();
      g->d->F = SPFunction_CopyBasis(f->d->F);
      g->d->G = SPFunction_CopyBasis(f->d->G);
      return g;
    }
  }

void SPErrorMap_M_Free(T *f)
  { affirm(isprefix(SPErrorMap_TypeId, f->type), "type/method bug");
    { Basis F = f->d->F;
      Basis G = f->d->G;
      int i;
      affirm(F.ne == G.ne, "size mismatch");
      for (i = 0; i < F.ne; i++) 
        { SPFunction *f = F.e[i], *g = G.e[i];
          if (f != NULL) { f->m->free(f); }
          if (g != NULL) { g->m->free(g); }
        }
      free(f->d->F.e);
      free(f->d->G.e);
      free(f->d);
      free(f);
    }
  }

/* CLASS-SPECIFIC METHODS */
  
void SPErrorMap_M_Write(T *f, FILE *wr)
  { filefmt_write_header(wr, "SPErrorMap", SPErrorMap_FileFormat);
    fprintf(wr, "GaugeBasis:\n");
    SPFunction_WriteBasis(wr, f->d->F);
    fprintf(wr, "\n");
    fprintf(wr, "GaugeApproxBasis:\n");
    SPFunction_WriteBasis(wr, f->d->G);
    fprintf(wr, "\n");
    filefmt_write_footer(wr, "SPErrorMap");
    fflush(wr);
  }
    
/* OTHER PROCS */
  
SPErrorMap_Methods *SPErrorMap_Methods_New(void)
  { void *v = malloc(sizeof(SPErrorMap_Methods));
    return (SPErrorMap_Methods *)notnull(v, "no mem for SPErrorMap_Methods");
  }

SPErrorMap_Data *SPErrorMap_Data_New(void)
  { void *v = malloc(sizeof(SPErrorMap_Data));
    return (SPErrorMap_Data *)notnull(v, "no mem for SPErrorMap_Data");
  }

SPErrorMap *SPErrorMap_New(void)
  { void *v = malloc(sizeof(SPErrorMap));
    return (SPErrorMap *)notnull(v, "no mem for SPErrorMap");
  }

static SPErrorMap_Methods *ErrorMapMths = NULL;

SPErrorMap *SPErrorMap_FullNew(void)
  { SPErrorMap *f = SPErrorMap_New();
    f->type = SPErrorMap_TypeId;
    f->d = SPErrorMap_Data_New();
    if (ErrorMapMths == NULL)
      { ErrorMapMths = SPErrorMap_Methods_New();
        ErrorMapMths->fn.eval = (SPFunction_EvalMth *)&SPErrorMap_M_Eval;
        ErrorMapMths->fn.grad = (SPFunction_GradMth *)&SPErrorMap_M_Grad;
        ErrorMapMths->fn.hess = (SPFunction_HessMth *)&SPErrorMap_M_Hess;
        ErrorMapMths->fn.maple = (SPFunction_MapleMth *)&SPErrorMap_M_Maple;
        /* Note: the {fn.write} method is inherited from {SPFunction}! */
        ErrorMapMths->fn.write = (SPFunction_WriteMth *)&SPFunction_M_Write;
        ErrorMapMths->fn.scale = (SPFunction_ScaleMth *)&SPErrorMap_M_Scale;
        ErrorMapMths->fn.add = (SPFunction_AddMth *)&SPErrorMap_M_Add;
        ErrorMapMths->fn.copy = (SPFunction_CopyMth *)&SPErrorMap_M_Copy;
        ErrorMapMths->fn.free = (SPFunction_FreeMth *)&SPErrorMap_M_Free;
        /* Class-specific methods */
        ErrorMapMths->write = (SPFunction_WriteMth *)&SPErrorMap_M_Write;
      }
    f->m = ErrorMapMths;
    return f;
  }

SPErrorMap *SPErrorMap_FromBases(Basis F, Basis G)
  { SPErrorMap *f = SPErrorMap_FullNew();
    affirm(F.ne == G.ne, "basis size mismatch");
    f->d->F = F;
    f->d->G = G;
    return f;
  }
    
SPErrorMap *SPErrorMap_Read(FILE *rd)
  { SPErrorMap *f = SPErrorMap_FullNew();
    filefmt_read_header(rd, "SPErrorMap", SPErrorMap_FileFormat);
    fget_skip_formatting_chars(rd);
    
    fget_match(rd, "GaugeBasis"); fget_skip_spaces(rd);
    fget_match(rd, ":"); fget_eol(rd);
    f->d->F = SPFunction_ReadBasis(rd);
    fget_skip_formatting_chars(rd);
    
    fget_match(rd, "GaugeApproxBasis"); fget_skip_spaces(rd);
    fget_match(rd, ":"); fget_eol(rd);
    f->d->G = SPFunction_ReadBasis(rd);
    fget_skip_formatting_chars(rd);
    
    filefmt_read_footer(rd, "SPErrorMap");
    return f;
  }
