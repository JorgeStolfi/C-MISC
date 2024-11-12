/* See SPTimeOnlyProcFunction.h */
/* Last edited on 2005-08-21 16:38:06 by stolfi */

#include <SPTimeOnlyProcFunction.h>
#include <filefmt.h>
#include <nget.h>
#include <fget.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <affirm.h>
#include <nat.h>

#define T SPTimeOnlyProcFunction 
    
void SPTimeOnlyProcFunction_M_Write(T *f, FILE *wr);
SPTimeOnlyProcFunction *SPTimeOnlyProcFunction_M_Copy(T *f);
void SPTimeOnlyProcFunction_M_Free(T *f);

SPTimeOnlyProcFunction_Methods *SPTimeOnlyProcFunction_Methods_New(void);
SPTimeOnlyProcFunction_Data *SPTimeOnlyProcFunction_Data_New(void);
SPTimeOnlyProcFunction *SPTimeOnlyProcFunction_New(void);

SPTimeOnlyProcFunction *SPTimeOnlyProcFunction_Cast(OBJ *f)
  { SPTimeOnlyFunction *ff = (SPTimeOnlyFunction *)f;
    if ((f != NULL) && isprefix(SPTimeOnlyProcFunction_TypeId, ff->type))
      { return (SPTimeOnlyProcFunction *)f; }
    else
      { return NULL; }
  }

SPTimeOnlyProcFunction *SPTimeOnlyProcFunction_Make
  ( char *type,
    char *descr,
    SPTimeOnlyProcFunction_Methods **mp,
    double (*eval)(T *f, double t)
  );
  /* Creates a new {SPTimeOnlyProcFunction} object (and its data record),
    with given {type} and {descr} fields. 

    If {*mp != NULL}, uses {*mp} as the methods record, ignoring the
    other arguments. If {*mp == NULL}, allocates a new
    {SPTimeOnlyProcFunction_Methods} record, with the given {eval},
    {grad}, and {hess} methods, stores its address in {*mp}, and uses
    that for the new object. */

/*** {unit(t)} **********************************/

static char *SPTOPF_UnitDesc = "1";
static char *SPTOPF_UnitType = (SPTimeOnlyProcFunction_TypeId "unit.");
static SPTimeOnlyProcFunction_Methods *SPTOPF_UnitMths = NULL;

double SPTOPF_UnitEval(T *f, double t);

double SPTOPF_UnitEval(T *f, double t)
  { double s = f->d->scale;
    return s;
  }

/*** {lint(t)} **********************************/

static char *SPTOPF_LintDesc = "t";
static char *SPTOPF_LintType = SPTimeOnlyProcFunction_TypeId "lint.";
static SPTimeOnlyProcFunction_Methods *SPTOPF_LintMths = NULL;

double SPTOPF_LintEval(T *f, double t);

double SPTOPF_LintEval(T *f, double t)
  { double s = f->d->scale;
    return s*t;
  }
    
/*** {htsq(t)} **********************************/

/* Solution of the oscillator equation with 
   {M=1}, {R=0}, {K=0}, {TMap(f,t)=1}, {a0=0}, {a1=0}. */

static char *SPTOPF_HtsqDesc = "t^2/2";
static char *SPTOPF_HtsqType = SPTimeOnlyProcFunction_TypeId "htsq.";
static SPTimeOnlyProcFunction_Methods *SPTOPF_HtsqMths = NULL;

double SPTOPF_HtsqEval(T *f, double t);

double SPTOPF_HtsqEval(T *f, double t)
  { double s = f->d->scale;
    return s*t*t/2;
  }
    
/*** {cost(t)} **********************************/

/* Solution of the oscillator equation with 
   {M=1}, {R=0}, {K=1}, {TMap(f,t)=0}, {a0=1}, {a1=0}. */

static char *SPTOPF_CostDesc = "cos(t)";
static char *SPTOPF_CostType = SPTimeOnlyProcFunction_TypeId "cost.";
static SPTimeOnlyProcFunction_Methods *SPTOPF_CostMths = NULL;

double SPTOPF_CostEval(T *f, double t);

double SPTOPF_CostEval(T *f, double t)
  { double s = f->d->scale;
    return s*cos(t);
  }

/*** {exnt(t)} **********************************/

/* Solution of the oscillator equation with 
   {M=0}, {R=1}, {K=1}, {TMap(f,t)=0}, {a0=1}, {a1=-1};
   {M=0}, {R=1}, {K=0}, {TMap(f,t)=-f}, {a0=1}, {a1=-1}; 
   {M=1}, {R=0}, {K=0}, {TMap(f,t)=f}, {a0=1}, {a1=-1}; 
   {M=1}, {R=0}, {K=-1}, {TMap(f,t)=0}, {a0=1}, {a1=-1}; 
   */

static char *SPTOPF_ExntDesc = "exp(-t)";
static char *SPTOPF_ExntType = SPTimeOnlyProcFunction_TypeId "exnt.";
static SPTimeOnlyProcFunction_Methods *SPTOPF_ExntMths = NULL;

double SPTOPF_ExntEval(T *f, double t);

double SPTOPF_ExntEval(T *f, double t)
  { double s = f->d->scale;
    return s*exp(-t);
  }

/*** {amor(t)} **********************************/

/* Solution of the oscillator equation with 
   {M=1}, {R=2}, {K=2}, {TMap(f,t)=0}, {a0=1}, {a1=-1};
   */

static char *SPTOPF_AmorDesc = "cos(t)*exp(-t)";
static char *SPTOPF_AmorType = SPTimeOnlyProcFunction_TypeId "amor.";
static SPTimeOnlyProcFunction_Methods *SPTOPF_AmorMths = NULL;

double SPTOPF_AmorEval(T *f, double t);

double SPTOPF_AmorEval(T *f, double t)
  { double s = f->d->scale;
    return s*cos(t)*exp(-t);
  }
  
/*************************************************/ 

#define SPTimeOnlyProcFunction_FileFormat "2005-08-17"

/* OVERRIDES FOR PARENT CLASS METHODS */
  
SPTimeOnlyProcFunction *SPTimeOnlyProcFunction_M_Copy(T *f)
  { SPTimeOnlyProcFunction *fnew = SPTimeOnlyProcFunction_New();
    fnew->d = SPTimeOnlyProcFunction_Data_New();
    *(fnew->d) = *(f->d);
    fnew->type = f->type;
    fnew->m = f->m;
    return fnew;
  }
 
void SPTimeOnlyProcFunction_M_Free(T *f)
  { affirm(isprefix(SPTimeOnlyProcFunction_TypeId, f->type), "type/method bug");
    free(f->d);
    free(f);
  }

/* CLASS-SPECIFIC METHODS */
  
void SPTimeOnlyProcFunction_M_Write(T *f, FILE *wr)
  { filefmt_write_header(wr, "SPTimeOnlyProcFunction", SPTimeOnlyProcFunction_FileFormat);
    fprintf(wr, "type = %s\n",  f->type);
    fprintf(wr, "descr = %s\n",  f->d->descr);
    fprintf(wr, "scale = %.16g\n",  f->d->scale);
    filefmt_write_footer(wr, "SPTimeOnlyProcFunction");
    fflush(wr);
  }

/* OTHER PROCS */
  
SPTimeOnlyProcFunction_Methods *SPTimeOnlyProcFunction_Methods_New(void)
  { void *v = malloc(sizeof(SPTimeOnlyProcFunction_Methods));
    return (SPTimeOnlyProcFunction_Methods *)notnull(v, "out of mem for SPTimeOnlyProcFunction_Methods");
  }

SPTimeOnlyProcFunction_Data *SPTimeOnlyProcFunction_Data_New(void)
  { void *v = malloc(sizeof(SPTimeOnlyProcFunction_Data));
    return (SPTimeOnlyProcFunction_Data *)notnull(v, "out of mem for SPTimeOnlyProcFunction_Data");
  }

SPTimeOnlyProcFunction *SPTimeOnlyProcFunction_New(void)
  { void *v = malloc(sizeof(SPTimeOnlyProcFunction));
    return (SPTimeOnlyProcFunction *)notnull(v, "no mem for SPTimeOnlyProcFunction");
  }

SPTimeOnlyProcFunction *SPTimeOnlyProcFunction_Make
  ( char *type,
    char *descr,
    SPTimeOnlyProcFunction_Methods **mp,
    double (*eval)(T *f, double t)
  )
  {
    SPTimeOnlyProcFunction *f = SPTimeOnlyProcFunction_New();
    f->type = type;
    if ((*mp) == NULL)
      { SPTimeOnlyProcFunction_Methods *m = SPTimeOnlyProcFunction_Methods_New();
        /* Superclass methods: */
        m->fn.eval = (SPTimeOnlyFunction_EvalMth *)eval;
        /* Note: the {fn.write} method is inherited from {SPTimeOnlyFunction}! */
        m->fn.write = (SPTimeOnlyFunction_WriteMth *)&SPTimeOnlyFunction_M_Write;
        m->fn.free = (SPTimeOnlyFunction_FreeMth *)&SPTimeOnlyProcFunction_M_Free;
        /* Class-specific methods */
        m->write = (SPTimeOnlyFunction_WriteMth *)&SPTimeOnlyProcFunction_M_Write;
        (*mp) = m;
      }
    f->m = (*mp);
    f->d = SPTimeOnlyProcFunction_Data_New();
    f->d->scale = 1.0;
    f->d->descr = descr;
    f->type = type;
    return f;
  }

#define MkFn SPTimeOnlyProcFunction_Make

SPTimeOnlyProcFunction *SPTimeOnlyProcFunction_FromName(char *name)
  { SPTimeOnlyProcFunction *f = SPTimeOnlyProcFunction_New();
    if (strcmp(name, "unit") == 0)
      { f = MkFn(SPTOPF_UnitType, SPTOPF_UnitDesc, &SPTOPF_UnitMths, SPTOPF_UnitEval); }
    else if (strcmp(name, "lint") == 0)
      { f = MkFn(SPTOPF_LintType, SPTOPF_LintDesc, &SPTOPF_LintMths, SPTOPF_LintEval); }
    else if (strcmp(name, "htsq") == 0)
      { f = MkFn(SPTOPF_HtsqType, SPTOPF_HtsqDesc, &SPTOPF_HtsqMths, SPTOPF_HtsqEval); }
    else if (strcmp(name, "cost") == 0)
      { f = MkFn(SPTOPF_CostType, SPTOPF_CostDesc, &SPTOPF_CostMths, SPTOPF_CostEval); }
    else if (strcmp(name, "exnt") == 0)
      { f = MkFn(SPTOPF_ExntType, SPTOPF_ExntDesc, &SPTOPF_ExntMths, SPTOPF_ExntEval); }
    else if (strcmp(name, "amor") == 0)
      { f = MkFn(SPTOPF_AmorType, SPTOPF_AmorDesc, &SPTOPF_AmorMths, SPTOPF_AmorEval); }
    else 
      { fprintf (stderr, "bad SPTimeOnlyProcFunction function name = %s\n", name);
        affirm(FALSE, "aborted");
      };
    return f;
  }

SPTimeOnlyProcFunction *SPTimeOnlyProcFunction_Read(FILE *rd)
  { char *t, *descr;
    double scale;
    SPTimeOnlyProcFunction *f;
    filefmt_read_header(rd, "SPTimeOnlyProcFunction", SPTimeOnlyProcFunction_FileFormat);
    t = nget_string(rd, "type"); fget_eol(rd);
    if (! isprefix(SPTimeOnlyProcFunction_TypeId, t))
      { fprintf (stderr, "bad SPTimeOnlyProcFunction type = %7s\n", t);
        affirm(FALSE, "expected \"" SPTimeOnlyProcFunction_TypeId "\"");
      }
    /* Remove final "." from subtype, and look it up: */
    { char *funcName = txtcat(t + 8, "#");
      affirm(strlen(funcName) > 0, "missing function name");
      affirm(funcName[strlen(funcName)-2] == '.', "missing \".\" in func type");
      funcName[strlen(funcName)-2] = '\000';
      f = SPTimeOnlyProcFunction_FromName(funcName);
      free(funcName);
    }
    descr = nget_string(rd, "descr"); fget_eol(rd);
    if (strcmp(f->d->descr, descr) != 0)
      { fprintf (stderr, "description mismatch: \"%s\", \"%s\"\n", 
          f->d->descr, descr);
        affirm(FALSE, "aborted");
      }
    scale = nget_double(rd, "scale"); fget_eol(rd);
    f->d->scale = scale;
    filefmt_read_footer(rd, "SPTimeOnlyProcFunction");
    free(t);
    free(descr);
    return f;
  }

