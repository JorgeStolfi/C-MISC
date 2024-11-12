/* See SPHBezFunction.h */
/* Last edited on 2023-02-12 07:50:29 by stolfi */

#include <SPHBezFunction.h>
#include <SPFunction.h>
#include <SPDeCasteljau.h>
#include <SPQuad.h>

#include <filefmt.h>
#include <fget.h>
#include <nget.h>
#include <r3.h>
#include <affirm.h>
#include <nat.h>

#include <stdio.h>
#include <stdlib.h>
    
#define T SPHBezFunction

double SPHBezFunction_M_Eval(T *f, R3Point *p);
R3Gradient SPHBezFunction_M_Grad(T *f, R3Point *p);
R3Hessian SPHBezFunction_M_Hess(T *f, R3Point *p);
void SPHBezFunction_M_Write(T *f, FILE *wr);
void SPHBezFunction_M_Add(T *f, double a, T *h);
void SPHBezFunction_M_Scale(T *f, double a);
void SPHBezFunction_M_Maple(T *f, FILE *wr);
SPHBezFunction *SPHBezFunction_M_Copy(T *f);
void SPHBezFunction_M_Free(T *f);

SPHBezFunction_Methods *SPHBezFunction_Methods_New(void);
SPHBezFunction_Data *SPHBezFunction_Data_New(void);
SPHBezFunction *SPHBezFunction_New(void);

SPHBezFunction *SPHBezFunction_FullNew(void);
  /* Allocates a new {SPHBezFunction} object {f}
    (and its data record), with the proper {type}
    fields and methods record. The array {f->d->c}
    is left undefined. */

SPHBezFunction *SPHBezFunction_Cast(OBJ *f)
  { SPFunction *ff = (SPFunction *)f;
    if ((f != NULL) && isprefix(SPHBezFunction_TypeId, ff->type))
      { return (SPHBezFunction *)f; }
    else
      { return NULL; }
  }

#define SPHBezFunction_FileFormat "2002-11-18"

/* OVERRIDES FOR PARENT CLASS METHODS */
  
double SPHBezFunction_M_Eval(T *f, R3Point *p) 
  { BezCoeff_vec_t coeff = f->d->c;
    double b;
    int deg = f->d->deg;
    switch (deg)
      { case 0: SPDeCasteljau_Eval0(coeff, p, &b); break;
        case 1: SPDeCasteljau_Eval1(coeff, p, &b); break;
        case 2: SPDeCasteljau_Eval2(coeff, p, &b); break;
        case 3: SPDeCasteljau_Eval3(coeff, p, &b); break;
        case 4: SPDeCasteljau_Eval4(coeff, p, &b); break;
        case 5: SPDeCasteljau_Eval5(coeff, p, &b); break;
        case 6: SPDeCasteljau_Eval6(coeff, p, &b); break;
        case 7: SPDeCasteljau_Eval7(coeff, p, &b); break;
        default: 
          affirm(deg >=0 , "invalid degree");
          SPDeCasteljau_EvalGen(coeff, deg, p, &b);
      }
    return b;
  }

R3Gradient SPHBezFunction_M_Grad(T *f, R3Point *p) 
  { BezCoeff_vec_t coeff = f->d->c;
    double b;
    r3_t db;
    int deg = f->d->deg;
    switch (deg)
      { case 0: SPDeCasteljau_Grad0(coeff, p, &b, &db); break;
        case 1: SPDeCasteljau_Grad1(coeff, p, &b, &db); break;
        case 2: SPDeCasteljau_Grad2(coeff, p, &b, &db); break;
        case 3: SPDeCasteljau_Grad3(coeff, p, &b, &db); break;
        case 4: SPDeCasteljau_Grad4(coeff, p, &b, &db); break;
        case 5: SPDeCasteljau_Grad5(coeff, p, &b, &db); break;
        case 6: SPDeCasteljau_Grad6(coeff, p, &b, &db); break;
        case 7: SPDeCasteljau_Grad7(coeff, p, &b, &db); break;
        default:
          affirm(deg >=0 , "invalid degree");
          SPDeCasteljau_GradGen(coeff, deg, p, &b, &db);
      }
    return db;
  }
  
R3Hessian SPHBezFunction_M_Hess(T *f, R3Point *p) 
  { BezCoeff_vec_t coeff = f->d->c;
    double b;
    r3_t db;
    r6_t ddb;
    int deg = f->d->deg;
    affirm(deg >=0 , "invalid degree");
    SPDeCasteljau_HessGen(coeff, deg, p, &b, &db, &ddb);
    return ddb;
  }

void SPHBezFunction_M_Maple(T *f, FILE *wr)
  { BezCoeff_vec_t coeff = f->d->c;
    int deg = f->d->deg;
    int nc; int j, k;
    nc = 0;
    for (j = 0; j <= deg; j++)
      { for (k = 0; k <= j; k++)
          { if (nc > 0) { fprintf(wr, " + "); } 
            fprintf(wr, "%.16g", coeff.e[nc]);
            fprintf(wr, "*B(%d,%d,%d)", deg-j, j-k, k);
            nc++;
          }
        fputc('\n', wr);
      }
    fputc('\n', wr);
    affirm (nc == coeff.ne, "inconsistent coeff count");
    fflush(wr);
  }
  
void SPHBezFunction_M_Add(T *f, double a, T *h)
  { affirm(isprefix(SPHBezFunction_TypeId, f->type), "type/method bug");
    affirm(isprefix(SPHBezFunction_TypeId, h->type), "wrong operand type");
    { BezCoeff_vec_t fc = f->d->c;
      BezCoeff_vec_t hc = h->d->c;
      int i;
      for (i = 0; i < fc.ne; i++) { fc.e[i] = fc.e[i] + a*hc.e[i]; }
    }
  }
  
void SPHBezFunction_M_Scale(T *f, double a)
  { affirm(isprefix(SPHBezFunction_TypeId, f->type), "type/method bug");
    { BezCoeff_vec_t fc = f->d->c;
      int i;
      for (i = 0; i < fc.ne; i++) { fc.e[i] *= a; }
    }
  }
  
T *SPHBezFunction_M_Copy(T *f)
  { affirm(isprefix(SPHBezFunction_TypeId, f->type), "type/method bug");
    { int i;
      BezCoeff_vec_t fc = f->d->c;
      SPHBezFunction *g = SPHBezFunction_FullNew();
      BezCoeff_vec_t gc = BezCoeff_vec_new(fc.ne);
      g->d->deg = f->d->deg;
      for (i = 0; i < fc.ne; i++) { gc.e[i] = fc.e[i]; }
      g->d->c = gc;
      return g;
    }
  }

void SPHBezFunction_M_Free(T *f)
  { affirm(isprefix(SPHBezFunction_TypeId, f->type), "type/method bug");
    free(f->d->c.e);
    free(f->d);
    free(f);
  }

/* CLASS-SPECIFIC METHODS */
  
void SPHBezFunction_M_Write(T *f, FILE *wr)
  { BezCoeff_vec_t coeff = f->d->c;
    int deg = f->d->deg;
    int nc; int j, k;

    filefmt_write_header(wr, "SPHBezFunction", SPHBezFunction_FileFormat);
    fprintf(wr, "degree = %d", deg); fputc('\n', wr);
    nc = 0;
    for (j = 0; j <= deg; j++)
      { for (k = 0; k <= j; k++)
          { if (k > 0) { fputc(' ', wr); } 
            fprintf(wr, "%.16g", coeff.e[nc]);
            nc++;
          }
        fputc('\n', wr);
      }
    affirm (nc == coeff.ne, "inconsistent coeff count");
    filefmt_write_footer(wr, "SPHBezFunction");
    fflush(wr);
  }
    
/* OTHER PROCS */
  
SPHBezFunction_Methods *SPHBezFunction_Methods_New(void)
  { void *v = malloc(sizeof(SPHBezFunction_Methods));
    return (SPHBezFunction_Methods *)notnull(v, "no mem for SPHBezFunction_Methods");
  }

SPHBezFunction_Data *SPHBezFunction_Data_New(void)
  { void *v = malloc(sizeof(SPHBezFunction_Data));
    return (SPHBezFunction_Data *)notnull(v, "no mem for SPHBezFunction_Data");
  }

SPHBezFunction *SPHBezFunction_New(void)
  { void *v = malloc(sizeof(SPHBezFunction));
    return (SPHBezFunction *)notnull(v, "no mem for SPHBezFunction");
  }

static SPHBezFunction_Methods *HBezMths = NULL;

SPHBezFunction *SPHBezFunction_FullNew(void)
  { SPHBezFunction *f = SPHBezFunction_New();
    f->type = SPHBezFunction_TypeId;
    f->d = SPHBezFunction_Data_New();
    if (HBezMths == NULL)
      { HBezMths = SPHBezFunction_Methods_New();
        HBezMths->fn.eval = (SPFunction_EvalMth *)&SPHBezFunction_M_Eval;
        HBezMths->fn.grad = (SPFunction_GradMth *)&SPHBezFunction_M_Grad;
        HBezMths->fn.hess = (SPFunction_HessMth *)&SPHBezFunction_M_Hess;
        HBezMths->fn.maple = (SPFunction_MapleMth *)&SPHBezFunction_M_Maple;
        /* Note: the {fn.write} method is inherited from {SPFunction}! */
        HBezMths->fn.write = (SPFunction_WriteMth *)&SPFunction_M_Write;
        HBezMths->fn.scale = (SPFunction_ScaleMth *)&SPHBezFunction_M_Scale;
        HBezMths->fn.add = (SPFunction_AddMth *)&SPHBezFunction_M_Add;
        HBezMths->fn.copy = (SPFunction_CopyMth *)&SPHBezFunction_M_Copy;
        HBezMths->fn.free = (SPFunction_FreeMth *)&SPHBezFunction_M_Free;
        /* Class-specific methods */
        HBezMths->write = (SPFunction_WriteMth *)&SPHBezFunction_M_Write;
      }
    f->m = HBezMths;
    return f;
  }

SPHBezFunction *SPHBezFunction_FromCoeffs(int deg, BezCoeff_vec_t coeff)
  { SPHBezFunction *f = SPHBezFunction_FullNew();
    affirm (coeff.ne == SPHBezFunction_NumCoeffs(deg), "inconsistent coeff count");
    f->d->deg = deg;
    f->d->c = coeff;
    return f;
  }
    
SPHBezFunction *SPHBezFunction_Read(FILE *rd)
  { int deg;
    SPHBezFunction *f = SPHBezFunction_FullNew();
    filefmt_read_header(rd, "SPHBezFunction", SPHBezFunction_FileFormat);
    deg = nget_int32(rd, "degree"); fget_eol(rd);
    affirm(deg >= 0, "bad degree");
    f->d->deg = deg;
    /* Read coefficients: */
    { int NC = SPHBezFunction_NumCoeffs(deg);
      BezCoeff_vec_t coeff = BezCoeff_vec_new(NC);
      int j, k; int nc = 0;
      for (j = 0; j <= deg; j++)
        { for (k = 0; k <= j; k++)
            { fget_skip_formatting_chars(rd);
              coeff.e[nc] = fget_double(rd);
              nc++;
            }
          fget_eol(rd);
        }
      f->d->c = coeff;
      affirm (nc == coeff.ne, "inconsistent coeff count");
    }
    filefmt_read_footer(rd, "SPHBezFunction");
    return f;
  }

nat_t SPHBezFunction_NumCoeffs(int deg) 
  { return (deg+1)*(deg+2)/2; }
