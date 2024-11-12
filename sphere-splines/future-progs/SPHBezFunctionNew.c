/* See SPHBezFunctionNew.h */
/* Last edited on 2023-02-12 07:56:02 by stolfi */

#include <SPHBezFunctionNew.h>
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

double SPHBezFunction_EvalMth(T *f, R3Point *p);
R3Gradient SPHBezFunction_GradMth(T *f, R3Point *p);
R3Hessian SPHBezFunction_HessMth(T *f, R3Point *p);
void SPHBezFunction_WriteMth(T *f, FILE *wr);
void SPHBezFunction_AddMth(T *f, double a, T *h);
void SPHBezFunction_ScaleMth(T *f, double a);
void SPHBezFunction_MapleMth(T *f, FILE *wr);
SPHBezFunction *SPHBezFunction_CopyMth(T *f);

SPHBezFunction_Methods *SPHBezFunction_Methods_New(void);
SPHBezFunction_Data *SPHBezFunction_Data_New(void);
SPHBezFunction *SPHBezFunction_New(void);

SPHBezFunction *SPHBezFunction_FullNew(void);
  /* Allocates a new {SPHBezFunction} object {f} (and its data
    record), with the proper {type} fields and methods record. The
    array {f->d->c} is left undefined. */

SPHBezFunction *SPHBezFunction_Cast(OBJ *f)
  { SPFunction *ff = (SPFunction *)f;
    if ((f != NULL) && isprefix(SPHBezFunction_TypeId, ff->type))
      { return (SPHBezFunction *)f; }
    else
      { return NULL; }
  }

#define SPHBezFunction_FileFormat "2002-11-18"

/* OVERRIDES FOR PARENT CLASS METHODS */
  
double SPHBezFunction_EvalMth(T *f, R3Point *p) 
  { BezCoeff_vec_t coeff = f->d->c;
    int *fex = f->d->fex;
    int tfex = fex[0]+fex[1]+fex[2];
    double h = (tfex == 0 ? 1.0 : SPHBezFunction_EvalMonomial(fex, p)), b;
    int degb = f->d->deg - tfex;
    switch (degb)
      { case 0: b = SPDeCasteljau_Eval0(coeff, p); break;
        case 1: b = SPDeCasteljau_Eval1(coeff, p); break;
        case 2: b = SPDeCasteljau_Eval2(coeff, p); break;
        case 3: b = SPDeCasteljau_Eval3(coeff, p); break;
        case 4: b = SPDeCasteljau_Eval4(coeff, p); break;
        case 5: b = SPDeCasteljau_Eval5(coeff, p); break;
        case 6: b = SPDeCasteljau_Eval6(coeff, p); break;
        case 7: b = SPDeCasteljau_Eval7(coeff, p); break;
        default: 
          affirm(degb >= 0, "invalid degree");
          b = SPDeCasteljau_EvalGen(coeff, degb, p);
      }
    return h*b;
  }

R3Gradient SPHBezFunction_GradMth(T *f, R3Point *p) 
  { BezCoeff_vec_t coeff = f->d->c;
    int fex[3] = f->d->fex, tfex = fex[0]+fex[1]+fex[2];
    int degb = f->d->deg - tfex
    if (tfex == 0)
      { switch (degb)
          { case 0: return SPDeCasteljau_EvalGrad0(coeff, p); break;
            case 1: return SPDeCasteljau_EvalGrad1(coeff, p); break;
            case 2: return SPDeCasteljau_EvalGrad2(coeff, p); break;
            case 3: return SPDeCasteljau_EvalGrad3(coeff, p); break;
            case 4: return SPDeCasteljau_EvalGrad4(coeff, p); break;
            case 5: return SPDeCasteljau_EvalGrad5(coeff, p); break;
            case 6: return SPDeCasteljau_EvalGrad6(coeff, p); break;
            case 7: return SPDeCasteljau_EvalGrad7(coeff, p); break;
            default: 
              affirm(degb >= 0, "invalid degree");
              return SPDeCasteljau_EvalGradGen(coeff, degb, p);
          }
      }
    else
      { double h = SPHBezFunction_EvalMonomial(fex, p), b;
        R3Gradient dh = SPHBezFunction_GradMonomial(fex, p), db, df;
        switch (degb)
          { case 0: 
              b = SPDeCasteljau_Eval0(coeff, p);
              db = SPDeCasteljau_EvalGrad0(coeff, p); 
              break;
            case 1: 
              b = SPDeCasteljau_Eval1(coeff, p);
              db = SPDeCasteljau_EvalGrad1(coeff, p); 
              break;
            case 2: 
              b = SPDeCasteljau_Eval2(coeff, p);
              db = SPDeCasteljau_EvalGrad2(coeff, p); 
              break;
            case 3: 
              b = SPDeCasteljau_Eval3(coeff, p);
              db = SPDeCasteljau_EvalGrad3(coeff, p); 
              break;
            case 4: 
              b = SPDeCasteljau_Eval4(coeff, p);
              db = SPDeCasteljau_EvalGrad4(coeff, p); 
              break;
            case 5: 
              b = SPDeCasteljau_Eval5(coeff, p);
              db = SPDeCasteljau_EvalGrad5(coeff, p); 
              break;
            case 6: 
              b = SPDeCasteljau_Eval6(coeff, p);
              db = SPDeCasteljau_EvalGrad6(coeff, p); 
              break;
            case 7: 
              b = SPDeCasteljau_Eval7(coeff, p);
              db = SPDeCasteljau_EvalGrad7(coeff, p); 
              break;
            default: 
              affirm(degb >= 0, "invalid degree"); 
              b = SPDeCasteljau_EvalGen(coeff, degb, p);
              db = SPDeCasteljau_EvalGradGen(coeff, degb, p); 
          }
        r3_mix(h, &db, b, &dh, &df);
        return df;
      }
  }
  
R3Hessian SPHBezFunction_HessMth(T *f, R3Point *p) 
  { BezCoeff_vec_t coeff = f->d->c;
    int fex[3] = f->d->fex, tfex = fex[0]+fex[1]+fex[2];
    int degb = f->d->deg - tfex
    if (tfex == 0)
      { switch (degb)
          { case 0: return SPDeCasteljau_EvalHess0(coeff, p); break;
            case 1: return SPDeCasteljau_EvalHess1(coeff, p); break;
            case 2: return SPDeCasteljau_EvalHess2(coeff, p); break;
            case 3: return SPDeCasteljau_EvalHess3(coeff, p); break;
            case 4: return SPDeCasteljau_EvalHess4(coeff, p); break;
            case 5: return SPDeCasteljau_EvalHess5(coeff, p); break;
            case 6: return SPDeCasteljau_EvalHess6(coeff, p); break;
            case 7: return SPDeCasteljau_EvalHess7(coeff, p); break;
            default: 
              affirm(degb >= 0, "invalid degree");
              return SPDeCasteljau_EvalHessGen(coeff, degb, p);
          }
      }
    else
      { double h = SPHBezFunction_EvalMonomial(fex, p), b;
        R3Gradient dh = SPHBezFunction_GradMonomial(fex, p), db, df;
        R3Hessian ddh = SPHBezFunction_HessMonomial(fex, p), ddb, ddf;
        switch (degb)
          { case 0: 
              b = SPDeCasteljau_Eval0(coeff, p);
              db = SPDeCasteljau_EvalGrad0(coeff, p); 
              break;
            case 1: 
              b = SPDeCasteljau_Eval1(coeff, p);
              db = SPDeCasteljau_EvalGrad1(coeff, p); 
              break;
            case 2: 
              b = SPDeCasteljau_Eval2(coeff, p);
              db = SPDeCasteljau_EvalGrad2(coeff, p); 
              break;
            case 3: 
              b = SPDeCasteljau_Eval3(coeff, p);
              db = SPDeCasteljau_EvalGrad3(coeff, p); 
              break;
            case 4: 
              b = SPDeCasteljau_Eval4(coeff, p);
              db = SPDeCasteljau_EvalGrad4(coeff, p); 
              break;
            case 5: 
              b = SPDeCasteljau_Eval5(coeff, p);
              db = SPDeCasteljau_EvalGrad5(coeff, p); 
              break;
            case 6: 
              b = SPDeCasteljau_Eval6(coeff, p);
              db = SPDeCasteljau_EvalGrad6(coeff, p); 
              break;
            case 7: 
              b = SPDeCasteljau_Eval7(coeff, p);
              db = SPDeCasteljau_EvalGrad7(coeff, p); 
              break;
            default: 
              affirm(deg >= 0, "invalid degree"); 
              b = SPDeCasteljau_EvalGen(coeff, deg, p);
              db = SPDeCasteljau_EvalGradGen(coeff, deg, p); 
          }
        ddb = SPDeCasteljau_EvalHessGen(coeff, deg, p); 
        r3_mix(h, &db, b, &dh, &df);
        return df;
      }
  
  
    affirm(FALSE , "hess method not implemented yet");
    return (R3Hessian){{0.0, 0.0, 0.0, 0.0, 0.0, 0.0}};
  }

void SPHBezFunction_MapleMth(T *f, FILE *wr)
  { BezCoeff_vec_t coeff = f->d->c;
    int deg = f->d->deg;
    nat nc; int j, k;
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
  
void SPHBezFunction_AddMth(T *f, double a, T *h)
  { affirm(isprefix(SPHBezFunction_TypeId, f->type), "type/method bug");
    affirm(isprefix(SPHBezFunction_TypeId, h->type), "wrong operand type");
    { BezCoeff_vec_t fc = f->d->c;
      BezCoeff_vec_t hc = h->d->c;
      int i;
      for (i = 0; i < fc.ne; i++) { fc.e[i] = fc.e[i] + a*hc.e[i]; }
    }
  }
  
void SPHBezFunction_ScaleMth(T *f, double a)
  { affirm(isprefix(SPHBezFunction_TypeId, f->type), "type/method bug");
    { BezCoeff_vec_t fc = f->d->c;
      int i;
      for (i = 0; i < fc.ne; i++) { fc.e[i] *= a; }
    }
  }
  
T *SPHBezFunction_CopyMth(T *f)
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

/* CLASS-SPECIFIC METHODS */
  
void SPHBezFunction_WriteMth(T *f, FILE *wr)
  { BezCoeff_vec_t coeff = f->d->c;
    int deg = f->d->deg;
    nat nc; int j, k;

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

static SPHBezFunction_Methods *HBezMths;

SPHBezFunction *SPHBezFunction_FullNew(void)
  { SPHBezFunction *f = SPHBezFunction_New();
    f->type = SPHBezFunction_TypeId;
    f->d = SPHBezFunction_Data_New();
    if (HBezMths == NULL)
      { HBezMths = SPHBezFunction_Methods_New();
        HBezMths->fn.eval = (EvalMth *)&SPHBezFunction_EvalMth;
        HBezMths->fn.grad = (GradMth *)&SPHBezFunction_GradMth;
        HBezMths->fn.hess = (HessMth *)&SPHBezFunction_HessMth;
        HBezMths->fn.maple = (MapleMth *)&SPHBezFunction_MapleMth;
        /* Note: the {fn.write} method is inherited from {SPFunction}! */
        HBezMths->fn.write = (WriteMth *)&SPFunction_WriteMth;
        HBezMths->fn.scale = (ScaleMth *)&SPHBezFunction_ScaleMth;
        HBezMths->fn.add = (AddMth *)&SPHBezFunction_AddMth;
        HBezMths->fn.copy = (CopyMth *)&SPHBezFunction_CopyMth;
        /* Class-specific methods */
        HBezMths->write = (WriteMth *)&SPHBezFunction_WriteMth;
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
    { nat NC = SPHBezFunction_NumCoeffs(deg);
      BezCoeff_vec_t coeff = BezCoeff_vec_new(NC);
      int j, k; nat nc = 0;
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
