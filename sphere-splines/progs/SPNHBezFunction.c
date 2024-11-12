/* see SPNHBezFunction.h */
/* Last edited on 2023-02-12 07:55:34 by stolfi */

#include <SPNHBezFunction.h>
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
    
#define T SPNHBezFunction

double SPNHBezFunction_M_Eval(T *f, R3Point *p);
R3Gradient SPNHBezFunction_M_Grad(T *f, R3Point *p);
R3Hessian SPNHBezFunction_M_Hess(T *f, R3Point *p);
void SPNHBezFunction_M_Write(T *f, FILE *wr);
void SPNHBezFunction_M_Add(T *f, double a, T *h);
void SPNHBezFunction_M_Scale(T *f, double a);
void SPNHBezFunction_M_Maple(T *f, FILE *wr);
SPNHBezFunction *SPNHBezFunction_M_Copy(T *f);
void SPNHBezFunction_M_Free(T *f);

SPNHBezFunction_Methods *SPNHBezFunction_Methods_New(void);
SPNHBezFunction_Data *SPNHBezFunction_Data_New(void);
SPNHBezFunction *SPNHBezFunction_New(void);

SPNHBezFunction *SPNHBezFunction_FullNew(void);
  /* Allocates a new {SPNHBezFunction} object {f}
    (and its data record), with the proper {type}
    fields and methods record. The array {f->d->c}
    is left undefined. */

SPNHBezFunction *SPNHBezFunction_Cast(OBJ *f)
  { SPFunction *ff = (SPFunction *)f;
    if ((f != NULL) && isprefix(SPNHBezFunction_TypeId, ff->type))
      { return (SPNHBezFunction *)f; }
    else
      { return NULL; }
  }

#define SPNHBezFunction_FileFormat "2002-11-18"

/* OVERRIDES FOR PARENT CLASS METHODS */
  
double SPNHBezFunction_M_Eval(T *f, R3Point *p) 
  { BezCoeff_vec_t c0 = f->d->c0;
    BezCoeff_vec_t c1 = f->d->c1;
    double b0, b1;
    int deg = f->d->deg;
    int deg0 = f->d->deg;
    int deg1 = deg0 - 1;
    switch (deg)
      { case 0: 
          SPDeCasteljau_Eval0(c0, p, &b0); 
          b1 = 0.0; break;
        case 1: 
          SPDeCasteljau_Eval1(c0, p, &b0); 
          SPDeCasteljau_Eval0(c1, p, &b1); 
          break;
        case 2: 
          SPDeCasteljau_Eval2(c0, p, &b0); 
          SPDeCasteljau_Eval1(c1, p, &b1); 
          break;
        case 3: 
          SPDeCasteljau_Eval3(c0, p, &b0); 
          SPDeCasteljau_Eval2(c1, p, &b1); 
          break;
        case 4: 
          SPDeCasteljau_Eval4(c0, p, &b0); 
          SPDeCasteljau_Eval3(c1, p, &b1); 
          break;
        case 5: 
          SPDeCasteljau_Eval5(c0, p, &b0); 
          SPDeCasteljau_Eval4(c1, p, &b1); 
          break;
        case 6: 
          SPDeCasteljau_Eval6(c0, p, &b0); 
          SPDeCasteljau_Eval5(c1, p, &b1); 
          break;
        case 7: 
          SPDeCasteljau_Eval7(c0, p, &b0); 
          SPDeCasteljau_Eval6(c1, p, &b1); 
          break;
        default: 
          affirm(deg0 >=0, "invalid degree");
          SPDeCasteljau_EvalGen(c0, deg0, p, &b0);
          SPDeCasteljau_EvalGen(c1, deg1, p, &b1);
      }
    return b0 + b1;
  }

R3Gradient SPNHBezFunction_M_Grad(T *f, R3Point *p) 
  { BezCoeff_vec_t c0 = f->d->c0;
    BezCoeff_vec_t c1 = f->d->c1;
    double b0, b1;
    r3_t db0, db1, db;
    int deg0 = f->d->deg;
    int deg1 = deg0 - 1;
    switch (deg0)
      { case 0:  
          SPDeCasteljau_Grad0(c0, p, &b0, &db0); 
          b1 = 0.0; db1 = (r3_t){{0.0, 0.0, 0.0}};  
          break;
        case 1:  
          SPDeCasteljau_Grad1(c0, p, &b0, &db0); 
          SPDeCasteljau_Grad0(c1, p, &b1, &db1);  
          break;
        case 2:  
          SPDeCasteljau_Grad2(c0, p, &b0, &db0); 
          SPDeCasteljau_Grad1(c1, p, &b1, &db1);  
          break;
        case 3:  
          SPDeCasteljau_Grad3(c0, p, &b0, &db0); 
          SPDeCasteljau_Grad2(c1, p, &b1, &db1);  
          break;
        case 4:  
          SPDeCasteljau_Grad4(c0, p, &b0, &db0); 
          SPDeCasteljau_Grad3(c1, p, &b1, &db1);  
          break;
        case 5:  
          SPDeCasteljau_Grad5(c0, p, &b0, &db0);  
          SPDeCasteljau_Grad4(c1, p, &b1, &db1);  
          break;
        case 6:  
          SPDeCasteljau_Grad6(c0, p, &b0, &db0);  
          SPDeCasteljau_Grad5(c1, p, &b1, &db1);  
          break;
        case 7:  
          SPDeCasteljau_Grad7(c0, p, &b0, &db0);  
          SPDeCasteljau_Grad6(c1, p, &b1, &db1);  
          break;
        default: 
          affirm(deg0 >= 0, "invalid degree"); 
          SPDeCasteljau_GradGen(c0, deg0, p, &b0, &db0); 
          SPDeCasteljau_GradGen(c1, deg1, p, &b1, &db1);
      }
    r3_add(&db0, &db1, &db);
    return db;
  }
  
R3Hessian SPNHBezFunction_M_Hess(T *f, R3Point *p) 
  { BezCoeff_vec_t c0 = f->d->c0;
    BezCoeff_vec_t c1 = f->d->c1;
    double b0, b1;
    r3_t db0, db1;
    r6_t ddb0, ddb1, ddb;
    int deg0 = f->d->deg;
    int deg1 = deg0 - 1;
    affirm(deg0 >= 0, "invalid degree"); 
    SPDeCasteljau_HessGen(c0, deg0, p, &b0, &db0, &ddb0); 
    SPDeCasteljau_HessGen(c1, deg1, p, &b1, &db1, &ddb1);
    r6_add(&ddb0, &ddb1, &ddb);
    return ddb;
  }

void SPNHBezFunction_M_Maple(T *f, FILE *wr)
  {
    BezCoeff_vec_t c0 = f->d->c0;
    BezCoeff_vec_t c1 = f->d->c1;
    int deg = f->d->deg;
    int nc;
    int j, k;
    
    /* Coefficients of {deg} component: */
    nc = 0;
    for (j = 0; j <= deg; j++)
      { for (k = 0; k <= j; k++)
          { if (nc > 0) { fprintf(wr, " + "); } 
            fprintf(wr, "%.16g", c0.e[nc]);
            fprintf(wr, "*B(%d,%d,%d)", deg-j, j-k, k);
            nc++;
          }
        fputc('\n', wr);
      }
    fputc('\n', wr);
    affirm (nc == c0.ne, "inconsistent coeff count");
    
    /* Coefficients of {deg-1} component: */
    if (deg > 0)
      { nc = 0;
        for (j = 0; j <= deg-1; j++)
          { for (k = 0; k <= j; k++)
              { fprintf(wr, " + ");
                fprintf(wr, "%.16g", c1.e[nc]);
                fprintf(wr, "*B(%d,%d,%d)", deg-1-j, j-k, k);
                nc++;
              }
            fputc('\n', wr);
          }
        fputc('\n', wr);
        affirm (nc == c1.ne, "inconsistent coeff count");
      }
    fflush(wr);
  }
  
void SPNHBezFunction_M_Add(T *f, double a, T *h)
  { affirm(isprefix(SPNHBezFunction_TypeId, f->type), "type/method bug");
    affirm(isprefix(SPNHBezFunction_TypeId, h->type), "wrong operand type");
    { BezCoeff_vec_t fc = f->d->c0;
      BezCoeff_vec_t hc = h->d->c0;
      int i;
      for (i = 0; i < fc.ne; i++) { fc.e[i] = fc.e[i] + a*hc.e[i]; }
    }
    { BezCoeff_vec_t fc = f->d->c1;
      BezCoeff_vec_t hc = h->d->c1;
      int i;
      for (i = 0; i < fc.ne; i++) { fc.e[i] = fc.e[i] + a*hc.e[i]; }
    }
  }
  
void SPNHBezFunction_M_Scale(T *f, double a)
  { affirm(isprefix(SPNHBezFunction_TypeId, f->type), "type/method bug");
    { BezCoeff_vec_t fc = f->d->c0;
      int i;
      for (i = 0; i < fc.ne; i++) { fc.e[i] *= a; }
    }
    { BezCoeff_vec_t fc = f->d->c1;
      int i;
      for (i = 0; i < fc.ne; i++) { fc.e[i] *= a; }
    }
  }
  
T *SPNHBezFunction_M_Copy(T *f)
  { affirm(isprefix(SPNHBezFunction_TypeId, f->type), "type/method bug");
    { int i;
      SPNHBezFunction *g = SPNHBezFunction_FullNew();
      g->d->deg = f->d->deg;
      { BezCoeff_vec_t fc = f->d->c0;
        BezCoeff_vec_t gc = double_vec_new(fc.ne);
        for (i = 0; i < fc.ne; i++) { gc.e[i] = fc.e[i]; }
        g->d->c0 = gc;
      }
      { BezCoeff_vec_t fc = f->d->c1;
        BezCoeff_vec_t gc = double_vec_new(fc.ne);
        for (i = 0; i < fc.ne; i++) { gc.e[i] = fc.e[i]; }
        g->d->c1 = gc;
      }
      return g;
    }
  }
 
void SPNHBezFunction_M_Free(T *f)
  { affirm(isprefix(SPNHBezFunction_TypeId, f->type), "type/method bug");
    free(f->d->c0.e);
    free(f->d->c1.e);
    free(f->d);
    free(f);
  }

/* CLASS-SPECIFIC METHODS */
  
void SPNHBezFunction_M_Write(T *f, FILE *wr)
  { BezCoeff_vec_t c0 = f->d->c0;
    BezCoeff_vec_t c1 = f->d->c1;
    int deg = f->d->deg;
    int nc;
    int j, k;

    filefmt_write_header(wr, "SPNHBezFunction", SPNHBezFunction_FileFormat);
    fprintf(wr, "degree = %d", deg); fputc('\n', wr);

    /* Coefficients of {deg} component: */
    nc = 0;
    for (j = 0; j <= deg; j++)
      { for (k = 0; k <= j; k++)
          { if (k > 0) { fputc(' ', wr); } 
            fprintf(wr, "%.16g", c0.e[nc]);
            nc++;
          }
        fputc('\n', wr);
      }
    affirm (nc == c0.ne, "inconsistent coeff count");

    /* Coefficients of {deg-1} component: */
    if (deg > 0)
      { nc = 0;
        for (j = 0; j <= deg-1; j++)
          { for (k = 0; k <= j; k++)
              { if (k > 0) { fputc(' ', wr); } 
                fprintf(wr, "%.16g", c1.e[nc]);
                nc++;
              }
            fputc('\n', wr);
          }
        fputc('\n', wr);
        affirm (nc == c1.ne, "inconsistent coeff count");
      }

    filefmt_write_footer(wr, "SPNHBezFunction");
    fflush(wr);
  }
    
/* OTHER PROCS */
  
SPNHBezFunction *SPNHBezFunction_FromHomo(SPHBezFunction *f, int deg)
  { int nc0 = SPHBezFunction_NumCoeffs(deg);
    int nc1 = SPHBezFunction_NumCoeffs(deg-1);
    BezCoeff_vec_t fc = f->d->c;
    BezCoeff_vec_t gc0 = double_vec_new(nc0);
    BezCoeff_vec_t gc1 = double_vec_new(nc1);
    SPNHBezFunction *g = SPNHBezFunction_FullNew();
    int i;
    
    if (f->d->deg == deg)
      { for (i = 0; i < gc0.ne; i++) { gc0.e[i] = fc.e[i]; }
        for (i = 0; i < gc1.ne; i++) { gc1.e[i] = 0.0; }
      }
    else if (f->d->deg == deg-1)
      { for (i = 0; i < gc0.ne; i++) { gc0.e[i] = 0.0; } 
        for (i = 0; i < gc1.ne; i++) { gc1.e[i] = fc.e[i]; } 
      }
    else
      { affirm(FALSE , "invalid homo degree"); }
    
    g->d->deg = deg;
    g->d->c0 = gc0;
    g->d->c1 = gc1;
    
    return g;
  }

SPNHBezFunction_Methods *SPNHBezFunction_Methods_New(void)
  { void *v = malloc(sizeof(SPNHBezFunction_Methods));
    return (SPNHBezFunction_Methods *)notnull(v, "no mem for SPNHBezFunction_Methods");
  }

SPNHBezFunction_Data *SPNHBezFunction_Data_New(void)
  { void *v = malloc(sizeof(SPNHBezFunction_Data));
    return (SPNHBezFunction_Data *)notnull(v, "no mem for SPNHBezFunction_Data");
  }

SPNHBezFunction *SPNHBezFunction_New(void)
  { void *v = malloc(sizeof(SPNHBezFunction));
    return (SPNHBezFunction *)notnull(v, "no mem for SPNHBezFunction");
  }

static SPNHBezFunction_Methods *NHBezMths = NULL;

SPNHBezFunction *SPNHBezFunction_FullNew(void)
  { SPNHBezFunction *f = SPNHBezFunction_New();
    f->type = SPNHBezFunction_TypeId;
    f->d = SPNHBezFunction_Data_New();
    if (NHBezMths == NULL)
      { NHBezMths = SPNHBezFunction_Methods_New();
        NHBezMths->fn.eval = (SPFunction_EvalMth *)&SPNHBezFunction_M_Eval;
        NHBezMths->fn.grad = (SPFunction_GradMth *)&SPNHBezFunction_M_Grad;
        NHBezMths->fn.hess = (SPFunction_HessMth *)&SPNHBezFunction_M_Hess;
        NHBezMths->fn.maple = (SPFunction_MapleMth *)&SPNHBezFunction_M_Maple;
        /* Note: the {fn.write} method is inherited from {SPFunction}! */
        NHBezMths->fn.write = (SPFunction_WriteMth *)&SPFunction_M_Write;
        NHBezMths->fn.scale = (SPFunction_ScaleMth *)&SPNHBezFunction_M_Scale;
        NHBezMths->fn.add = (SPFunction_AddMth *)&SPNHBezFunction_M_Add;
        NHBezMths->fn.copy = (SPFunction_CopyMth *)&SPNHBezFunction_M_Copy;
        NHBezMths->fn.free = (SPFunction_FreeMth *)&SPNHBezFunction_M_Free;
        /* Class-specific methods */
        NHBezMths->write = (SPFunction_WriteMth *)&SPNHBezFunction_M_Write;
      }
    f->m = NHBezMths;
    return f;
  }

SPNHBezFunction *SPNHBezFunction_Read(FILE *rd)
  { int deg;
    SPNHBezFunction *f = SPNHBezFunction_FullNew();
    filefmt_read_header(rd, "SPNHBezFunction", SPNHBezFunction_FileFormat);
    deg = nget_int32(rd, "degree"); fget_eol(rd);
    f->d->deg = deg;
    
    /* Read coefficients of {deg} component: */
    { int nc = SPHBezFunction_NumCoeffs(deg);
      BezCoeff_vec_t c = double_vec_new(nc);
      int j, k;
      nc = 0;
      for (j = 0; j <= deg; j++)
        { for (k = 0; k <= j; k++)
            { fget_skip_formatting_chars(rd);
              c.e[nc] = fget_double(rd);
              nc++;
            }
          fget_eol(rd);
        }
      f->d->c0 = c;
    }
    
    /* Read coefficients of {deg-1} component: */
    { int nc = SPHBezFunction_NumCoeffs(deg-1);
      BezCoeff_vec_t c = double_vec_new(nc);
      if (deg > 0)
        { int j, k;
          nc = 0;
          for (j = 0; j <= deg-1; j++)
            { for (k = 0; k <= j; k++)
                { fget_skip_formatting_chars(rd);
                  c.e[nc] = fget_double(rd);
                  nc++;
                }
              fget_eol(rd);
            }
        }
      f->d->c1 = c;
    }
    
    filefmt_read_footer(rd, "SPNHBezFunction");
    return f;
  }

SPNHBezFunction *SPNHBezFunction_FromCoeffs
  ( int deg, 
    BezCoeff_vec_t c0,
    BezCoeff_vec_t c1
  )
  { SPNHBezFunction *f = SPNHBezFunction_FullNew();
    affirm (c0.ne == SPHBezFunction_NumCoeffs(deg), "inconsistent coeff count");
    affirm (c1.ne == SPHBezFunction_NumCoeffs(deg-1), "inconsistent coeff count");
    f->d->deg = deg;
    f->d->c0 = c0;
    f->d->c1 = c1;
    return f;
  }
    
