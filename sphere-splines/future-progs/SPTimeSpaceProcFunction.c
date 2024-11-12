/* See SPTimeSpaceProcFunction.h */
/* Last edited on 2005-10-01 17:06:22 by stolfi */

#include <SPTimeSpaceProcFunction.h>
#include <filefmt.h>
#include <nget.h>
#include <fget.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <r3.h>
#include <r6.h>
#include <affirm.h>
#include <nat.h>

#define T SPTimeSpaceProcFunction 
    
void SPTimeSpaceProcFunction_M_Write(T *f, FILE *wr);
SPTimeSpaceProcFunction *SPTimeSpaceProcFunction_M_Copy(T *f);
void SPTimeSpaceProcFunction_M_Free(T *f);

SPTimeSpaceProcFunction_Methods *SPTimeSpaceProcFunction_Methods_New(void);
SPTimeSpaceProcFunction_Data *SPTimeSpaceProcFunction_Data_New(void);
SPTimeSpaceProcFunction *SPTimeSpaceProcFunction_New(void);

SPTimeSpaceProcFunction *SPTimeSpaceProcFunction_Cast(OBJ *f)
  { SPTimeSpaceFunction *ff = (SPTimeSpaceFunction *)f;
    if ((f != NULL) && isprefix(SPTimeSpaceProcFunction_TypeId, ff->type))
      { return (SPTimeSpaceProcFunction *)f; }
    else
      { return NULL; }
  }

SPTimeSpaceProcFunction *SPTimeSpaceProcFunction_Make
  ( char *type,
    char *descr,
    SPTimeSpaceProcFunction_Methods **mp,
    double (*eval)(T *f, R3Point *p, double t)
  );
  /* Creates a new {SPTimeSpaceProcFunction} object (and its data record),
    with given {type} and {descr} fields. 

    If {*mp != NULL}, uses {*mp} as the methods record, ignoring the
    other arguments. If {*mp == NULL}, allocates a new
    {SPTimeSpaceProcFunction_Methods} record, with the given {eval}, {grad}, and
    {hess} methods, stores its address in {*mp}, and uses that for the
    new object. */

/*** {unit(p,t)} **********************************/

static char *SPSTPF_UnitDesc = "1";
static char *SPSTPF_UnitType = "STF.Proc.unit.";
static SPTimeSpaceProcFunction_Methods *SPSTPF_UnitMths = NULL;

double SPSTPF_UnitEval(T *f, R3Point *p, double t);

double SPSTPF_UnitEval(T *f, R3Point *p, double t)
  { double s = f->d->scale;
    return s;
  }

/*** {linx(p,t)} **********************************/

static char *SPSTPF_LinxDesc = "x";
static char *SPSTPF_LinxType = "STF.Proc.linx.";
static SPTimeSpaceProcFunction_Methods *SPSTPF_LinxMths = NULL;

double SPSTPF_LinxEval(T *f, R3Point *p, double t);

double SPSTPF_LinxEval(T *f, R3Point *p, double t)
  { double s = f->d->scale, x = p->c[0];
    return s*x;
  }
    
/*** {liny(p,t)} **********************************/

static char *SPSTPF_LinyDesc = "y";
static char *SPSTPF_LinyType = "STF.Proc.liny.";
static SPTimeSpaceProcFunction_Methods *SPSTPF_LinyMths = NULL;

double SPSTPF_LinyEval(T *f, R3Point *p, double t);

double SPSTPF_LinyEval(T *f, R3Point *p, double t)
  { double s = f->d->scale, y = p->c[1];
    return s*y;
  }
    
/*** {linz(p,t)} **********************************/

static char *SPSTPF_LinzDesc = "z";
static char *SPSTPF_LinzType = "STF.Proc.linz.";
static SPTimeSpaceProcFunction_Methods *SPSTPF_LinzMths = NULL;

double SPSTPF_LinzEval(T *f, R3Point *p, double t);

double SPSTPF_LinzEval(T *f, R3Point *p, double t)
  { double s = f->d->scale, z = p->c[2];
    return s*z;
  }
    
/*** {lint(p,t)} **********************************/

static char *SPSTPF_LintDesc = "t";
static char *SPSTPF_LintType = "STF.Proc.lint.";
static SPTimeSpaceProcFunction_Methods *SPSTPF_LintMths = NULL;

double SPSTPF_LintEval(T *f, R3Point *p, double t);

double SPSTPF_LintEval(T *f, R3Point *p, double t)
  { double s = f->d->scale;
    return s*t;
  }
    
/*** {sqrx(p,t)} **********************************/

static char *SPSTPF_SqrxDesc = "x^2";
static char *SPSTPF_SqrxType = "STF.Proc.sqrx.";
static SPTimeSpaceProcFunction_Methods *SPSTPF_SqrxMths = NULL;

double SPSTPF_SqrxEval(T *f, R3Point *p, double t);

double SPSTPF_SqrxEval(T *f, R3Point *p, double t)
  { double s = f->d->scale, x = p->c[0];
    return s * x*x;
  }

/*** {sqry(p,t)} **********************************/

static char *SPSTPF_SqryDesc = "y^2";
static char *SPSTPF_SqryType = "STF.Proc.sqry.";
static SPTimeSpaceProcFunction_Methods *SPSTPF_SqryMths = NULL;

double SPSTPF_SqryEval(T *f, R3Point *p, double t);

double SPSTPF_SqryEval(T *f, R3Point *p, double t)
  { double s = f->d->scale, y = p->c[1];
    return s*y*y;
  }

/*** {sqtm(p,t)} **********************************/

static char *SPSTPF_SqtmDesc = "t^2";
static char *SPSTPF_SqtmType = "STF.Proc.sqtm.";
static SPTimeSpaceProcFunction_Methods *SPSTPF_SqtmMths = NULL;

double SPSTPF_SqtmEval(T *f, R3Point *p, double t);

double SPSTPF_SqtmEval(T *f, R3Point *p, double t)
  { double s = f->d->scale;
    return s*t*t;
  }

/*** {cost(p,t)} **********************************/

static char *SPSTPF_CostDesc = "cos(t)";
static char *SPSTPF_CostType = "STF.Proc.cost.";
static SPTimeSpaceProcFunction_Methods *SPSTPF_CostMths = NULL;

double SPSTPF_CostEval(T *f, R3Point *p, double t);

double SPSTPF_CostEval(T *f, R3Point *p, double t)
  { double s = f->d->scale;
    return s*cos(t);
  }

/*** {hent(p,t)} **********************************/

/* Solution of the heat diffusion equation on the 
  sphere with {K=1/6}, {R=0}, {L=0}, {f(p,0)=hrm2(p)}. */

static char *SPSTPF_HentDesc = "(x^2-1/3)*exp(-t)";
static char *SPSTPF_HentType = "STF.Proc.hent.";
static SPTimeSpaceProcFunction_Methods *SPSTPF_HentMths = NULL;

double SPSTPF_HentEval(T *f, R3Point *p, double t);

double SPSTPF_HentEval(T *f, R3Point *p, double t)
  { double s = f->d->scale, x = p->c[0];
    return s*(x*x - 1.0/3.0)*exp(-t);
  }
  
/*** {hspn(p,t)} **********************************/

/* Solution of the heat diffusion equation on the 
  rotating sphere with {K=1/6}, {R=1/32}, {L=0}, {f(p,0)=hrm2(p)}. */

static char *SPSTPF_HspnDesc = "m=x*cos(t/32)+y*sin(t/32);(m^2-1/3)*exp(-t)";
static char *SPSTPF_HspnType = "STF.Proc.hspn.";
static SPTimeSpaceProcFunction_Methods *SPSTPF_HspnMths = NULL;

double SPSTPF_HspnEval(T *f, R3Point *p, double t);

double SPSTPF_HspnEval(T *f, R3Point *p, double t)
  { double s = f->d->scale, x = p->c[0], y = p->c[1];
    double w = t/32;
    double m = x*cos(w) + y*sin(w);
    return s*(m*m - 1.0/3.0)*exp(-t);
  }
  
/*** {xent(p,t)} **********************************/

static char *SPSTPF_XentDesc = "x*exp(-t)";
static char *SPSTPF_XentType = "STF.Proc.xent.";
static SPTimeSpaceProcFunction_Methods *SPSTPF_XentMths = NULL;

double SPSTPF_XentEval(T *f, R3Point *p, double t);

double SPSTPF_XentEval(T *f, R3Point *p, double t)
  { double s = f->d->scale, x = p->c[0];
    return s*x*exp(-t);
  }
  
/*** {xsnt(p,t)} **********************************/

static char *SPSTPF_XsntDesc = "x*sin(t)";
static char *SPSTPF_XsntType = "STF.Proc.xsnt.";
static SPTimeSpaceProcFunction_Methods *SPSTPF_XsntMths = NULL;

double SPSTPF_XsntEval(T *f, R3Point *p, double t);

double SPSTPF_XsntEval(T *f, R3Point *p, double t)
  { double s = f->d->scale, x = p->c[0];
    return s*x*sin(t);
  }

/*** {turn(p,t)} **********************************/

static char *SPSTPF_TurnDesc = "x*cos(t)+y*sin(t)";
static char *SPSTPF_TurnType = "STF.Proc.turn.";
static SPTimeSpaceProcFunction_Methods *SPSTPF_TurnMths = NULL;

double SPSTPF_TurnEval(T *f, R3Point *p, double t);

double SPSTPF_TurnEval(T *f, R3Point *p, double t)
  { double s = f->d->scale, x = p->c[0], y = p->c[1];
    return s*(x*cos(t) + y*sin(t));
  }

/*** {ripl(p,t)} **********************************/

static char *SPSTPF_RiplDesc = "sin(6*PI*(x-t))";
static char *SPSTPF_RiplType = "STF.Proc.ripl.";
static SPTimeSpaceProcFunction_Methods *SPSTPF_RiplMths = NULL;

#define RiplBands 8

double SPSTPF_RiplEval(T *f, R3Point *p, double t);

double SPSTPF_RiplEval(T *f, R3Point *p, double t)
  { double s = f->d->scale, x = p->c[0];
    return s*sin(RiplBands*PI*(x-t));
  }

/*** {spir(p,t)} **********************************/

static char *SPSTPF_SpirDesc = "SpiralFn(p,t)";
static char *SPSTPF_SpirType = "STF.Proc.spir.";
static SPTimeSpaceProcFunction_Methods *SPSTPF_SpirMths = NULL;

double SPSTPF_SpirEval(T *f, R3Point *p, double t);

#define SpirC (2*PI/16)
#define SpirU 1.0625

double SPSTPF_SpirEval(T *f, R3Point *p, double t)
  { double s = f->d->scale;
    double x = p->c[0], y = p->c[1], z = p->c[2];
    double ct = cos(t);
    double st = sin(t);
    double yy = ct*y - st*z;
    double zz = st*y + ct*z;
    double h = x/SpirU, h2 = h*h, h4 = h2*h2, h8 = h4*h4, h10=h8*h2;
    double m = 1.0 - h2;
    double q = 1.0 - h10;
    double w = 2*h/m;
    double cw = (fabs(w) > 1.0e6 ? 0.0 : cos(SpirC*w));
    double sw = (fabs(w) > 1.0e6 ? 0.0 : sin(SpirC*w));
    double fv = q*(cw*yy - sw*zz);
    return s*fv;
  }
  
/*** {rain(p,t)} **********************************/

static char *SPSTPF_RainDesc = "Raindrop(p,t)";
static char *SPSTPF_RainType = "STF.Proc.rain.";
static SPTimeSpaceProcFunction_Methods *SPSTPF_RainMths = NULL;

double SPSTPF_RainEval(T *f, R3Point *p, double t);

#define RainShort 0.95
#define RainBands 8

double SPSTPF_RainEval(T *f, R3Point *p, double t)
  { double s = f->d->scale, x = p->c[0], y = p->c[1], z = p->c[2];
    double ct = RainShort*(x + y + z)/SQRT3;
    if (ct > +1.0) { ct = +1.0; }
    if (ct < -1.0) { ct = -1.0; }
    double st = sqrt(1.0 - ct*ct);
    double g = acos(ct);
    double m = 1/st;
    double w = RainBands*PI*(g-t);
    double fv = m*cos(w);
    return s*fv;
  }
  
/*************************************************/ 

#define SPTimeSpaceProcFunction_FileFormat "2002-11-12"

/* OVERRIDES FOR PARENT CLASS METHODS */
  
SPTimeSpaceProcFunction *SPTimeSpaceProcFunction_M_Copy(T *f)
  { SPTimeSpaceProcFunction *fnew = SPTimeSpaceProcFunction_New();
    fnew->d = SPTimeSpaceProcFunction_Data_New();
    *(fnew->d) = *(f->d);
    fnew->type = f->type;
    fnew->m = f->m;
    return fnew;
  }
 
void SPTimeSpaceProcFunction_M_Free(T *f)
  { affirm(isprefix(SPTimeSpaceProcFunction_TypeId, f->type), "type/method bug");
    free(f->d);
    free(f);
  }

/* CLASS-SPECIFIC METHODS */
  
void SPTimeSpaceProcFunction_M_Write(T *f, FILE *wr)
  { filefmt_write_header(wr, "SPTimeSpaceProcFunction", SPTimeSpaceProcFunction_FileFormat);
    fprintf(wr, "type = %s\n",  f->type);
    fprintf(wr, "descr = %s\n",  f->d->descr);
    fprintf(wr, "scale = %.16g\n",  f->d->scale);
    filefmt_write_footer(wr, "SPTimeSpaceProcFunction");
    fflush(wr);
  }

/* OTHER PROCS */
  
SPTimeSpaceProcFunction_Methods *SPTimeSpaceProcFunction_Methods_New(void)
  { void *v = malloc(sizeof(SPTimeSpaceProcFunction_Methods));
    return (SPTimeSpaceProcFunction_Methods *)notnull(v, "out of mem for SPTimeSpaceProcFunction_Methods");
  }

SPTimeSpaceProcFunction_Data *SPTimeSpaceProcFunction_Data_New(void)
  { void *v = malloc(sizeof(SPTimeSpaceProcFunction_Data));
    return (SPTimeSpaceProcFunction_Data *)notnull(v, "out of mem for SPTimeSpaceProcFunction_Data");
  }

SPTimeSpaceProcFunction *SPTimeSpaceProcFunction_New(void)
  { void *v = malloc(sizeof(SPTimeSpaceProcFunction));
    return (SPTimeSpaceProcFunction *)notnull(v, "no mem for SPTimeSpaceProcFunction");
  }

SPTimeSpaceProcFunction *SPTimeSpaceProcFunction_Make
  ( char *type,
    char *descr,
    SPTimeSpaceProcFunction_Methods **mp,
    double (*eval)(T *f, R3Point *p, double t)
  )
  {
    SPTimeSpaceProcFunction *f = SPTimeSpaceProcFunction_New();
    f->type = type;
    if ((*mp) == NULL)
      { SPTimeSpaceProcFunction_Methods *m = SPTimeSpaceProcFunction_Methods_New();
        /* Superclass methods: */
        m->fn.eval = (SPTimeSpaceFunction_EvalMth *)eval;
        /* Note: the {fn.write} method is inherited from {SPTimeSpaceFunction}! */
        m->fn.write = (SPTimeSpaceFunction_WriteMth *)&SPTimeSpaceFunction_M_Write;
        m->fn.free = (SPTimeSpaceFunction_FreeMth *)&SPTimeSpaceProcFunction_M_Free;
        /* Class-specific methods */
        m->write = (SPTimeSpaceFunction_WriteMth *)&SPTimeSpaceProcFunction_M_Write;
        (*mp) = m;
      }
    f->m = (*mp);
    f->d = SPTimeSpaceProcFunction_Data_New();
    f->d->scale = 1.0;
    f->d->descr = descr;
    f->type = type;
    return f;
  }

#define MkFn SPTimeSpaceProcFunction_Make

SPTimeSpaceProcFunction *SPTimeSpaceProcFunction_FromName(char *name)
  { SPTimeSpaceProcFunction *f = SPTimeSpaceProcFunction_New();
    if (strcmp(name, "unit") == 0)
      { f = MkFn(SPSTPF_UnitType, SPSTPF_UnitDesc, &SPSTPF_UnitMths, SPSTPF_UnitEval); }
    else if (strcmp(name, "linx") == 0)
      { f = MkFn(SPSTPF_LinxType, SPSTPF_LinxDesc, &SPSTPF_LinxMths, SPSTPF_LinxEval); }
    else if (strcmp(name, "liny") == 0)
      { f = MkFn(SPSTPF_LinyType, SPSTPF_LinyDesc, &SPSTPF_LinyMths, SPSTPF_LinyEval); }
    else if (strcmp(name, "linz") == 0)
      { f = MkFn(SPSTPF_LinzType, SPSTPF_LinzDesc, &SPSTPF_LinzMths, SPSTPF_LinzEval); }
    else if (strcmp(name, "lint") == 0)
      { f = MkFn(SPSTPF_LintType, SPSTPF_LintDesc, &SPSTPF_LintMths, SPSTPF_LintEval); }
    else if (strcmp(name, "sqrx") == 0)
      { f = MkFn(SPSTPF_SqrxType, SPSTPF_SqrxDesc, &SPSTPF_SqrxMths, SPSTPF_SqrxEval); }
    else if (strcmp(name, "sqry") == 0)
      { f = MkFn(SPSTPF_SqryType, SPSTPF_SqryDesc, &SPSTPF_SqryMths, SPSTPF_SqryEval); }
    else if (strcmp(name, "sqtm") == 0)
      { f = MkFn(SPSTPF_SqtmType, SPSTPF_SqtmDesc, &SPSTPF_SqtmMths, SPSTPF_SqtmEval); }
    else if (strcmp(name, "cost") == 0)
      { f = MkFn(SPSTPF_CostType, SPSTPF_CostDesc, &SPSTPF_CostMths, SPSTPF_CostEval); }
    else if (strcmp(name, "hent") == 0)
      { f = MkFn(SPSTPF_HentType, SPSTPF_HentDesc, &SPSTPF_HentMths, SPSTPF_HentEval); }
    else if (strcmp(name, "hspn") == 0)
      { f = MkFn(SPSTPF_HspnType, SPSTPF_HspnDesc, &SPSTPF_HspnMths, SPSTPF_HspnEval); }
    else if (strcmp(name, "xent") == 0)
      { f = MkFn(SPSTPF_XentType, SPSTPF_XentDesc, &SPSTPF_XentMths, SPSTPF_XentEval); }
    else if (strcmp(name, "turn") == 0)
      { f = MkFn(SPSTPF_TurnType, SPSTPF_TurnDesc, &SPSTPF_TurnMths, SPSTPF_TurnEval); }
    else if (strcmp(name, "xsnt") == 0)
      { f = MkFn(SPSTPF_XsntType, SPSTPF_XsntDesc, &SPSTPF_XsntMths, SPSTPF_XsntEval); }
    else if (strcmp(name, "ripl") == 0)
      { f = MkFn(SPSTPF_RiplType, SPSTPF_RiplDesc, &SPSTPF_RiplMths, SPSTPF_RiplEval); }
    else if (strcmp(name, "spir") == 0)
      { f = MkFn(SPSTPF_SpirType, SPSTPF_SpirDesc, &SPSTPF_SpirMths, SPSTPF_SpirEval); }
    else if (strcmp(name, "rain") == 0)
      { f = MkFn(SPSTPF_RainType, SPSTPF_RainDesc, &SPSTPF_RainMths, SPSTPF_RainEval); }
    else 
      { fprintf (stderr, "bad SPTimeSpaceProcFunction function name = %s\n", name);
        affirm(FALSE, "aborted");
      };
    return f;
  }

SPTimeSpaceProcFunction *SPTimeSpaceProcFunction_Read(FILE *rd)
  { char *t, *descr;
    double scale;
    SPTimeSpaceProcFunction *f;
    filefmt_read_header(rd, "SPTimeSpaceProcFunction", SPTimeSpaceProcFunction_FileFormat);
    t = nget_string(rd, "type"); fget_eol(rd);
    if (! isprefix("STF.Proc.", t))
      { fprintf (stderr, "bad SPTimeSpaceProcFunction type = %7s\n", t);
        affirm(FALSE, "expected \"STF.Proc.\"");
      }
    /* Remove final "." from subtype, and look it up: */
    { char *funcName = txtcat(t + 9, "#");
      affirm(strlen(funcName) > 0, "missing function name");
      affirm(funcName[strlen(funcName)-2] == '.', "missing \".\" in func type");
      funcName[strlen(funcName)-2] = '\000';
      f = SPTimeSpaceProcFunction_FromName(funcName);
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
    filefmt_read_footer(rd, "SPTimeSpaceProcFunction");
    free(t);
    free(descr);
    return f;
  }

