/* See SPProcFunction.h */
/* Last edited on 2005-10-27 15:55:32 by stolfi */

#include <SPProcFunction.h>
#include <SPHarmonic.h>

#include <filefmt.h>
#include <nget.h>
#include <fget.h>
#include <r3.h>
#include <r6.h>
#include <affirm.h>
#include <nat.h>

#include <stdlib.h>
#include <string.h>
#include <math.h>

#define T SPProcFunction 
    
void SPProcFunction_M_Write(T *f, FILE *wr);
void SPProcFunction_M_Add(T *f, double a, T *h);
void SPProcFunction_M_Scale(T *f, double a);
void SPProcFunction_M_Maple(T *f, FILE *wr);
SPProcFunction *SPProcFunction_M_Copy(T *f);
void SPProcFunction_M_Free(T *f);

SPProcFunction_Methods *SPProcFunction_Methods_New(void);
SPProcFunction_Data *SPProcFunction_Data_New(void);
SPProcFunction *SPProcFunction_New(void);

SPProcFunction *SPProcFunction_Cast(OBJ *f)
  { SPFunction *ff = (SPFunction *)f;
    if ((f != NULL) && isprefix(SPProcFunction_TypeId, ff->type))
      { return (SPProcFunction *)f; }
    else
      { return NULL; }
  }

SPProcFunction *SPProcFunction_Make
  ( char *type,
    char *descr,
    SPProcFunction_Methods **mp,
    double (*eval)(T *f, R3Point *p),
    R3Gradient (*grad)(T *f, R3Point *p),
    R3Hessian (*hess)(T *f, R3Point *p)
  );
  /* Creates a new {SPProcFunction} object (and its data record),
    with given {type} and {descr} fields. 

    If {*mp != NULL}, uses {*mp} as the methods record, ignoring the
    other arguments. If {*mp == NULL}, allocates a new
    {SPProcFunction_Methods} record, with the given {eval}, {grad}, and
    {hess} methods, stores its address in {*mp}, and uses that for the
    new object. */

/*** {unit(p)} **********************************/

static char *SPPF_UnitDesc = "1";
static char *SPPF_UnitType = "SF.Proc.unit.";
static SPProcFunction_Methods *SPPF_UnitMths = NULL;

double SPPF_UnitEval(T *f, R3Point *p);
R3Gradient UnitGrad(T *f, R3Point *p);
R3Hessian SPPF_UnitHess(T *f, R3Point *p);

double SPPF_UnitEval(T *f, R3Point *p)
  { double s = f->d->scale;
    return s;
  }
  
R3Gradient UnitGrad(T *f, R3Point *p)
  { return (r3_t){{0.0, 0.0, 0.0}};
  }

R3Hessian SPPF_UnitHess(T *f, R3Point *p)
  { return (R3Hessian){{0.0, 0.0, 0.0, 0.0, 0.0, 0.0}};
  }

/*** {linx(p)} **********************************/

static char *SPPF_LinxDesc = "x";
static char *SPPF_LinxType = "SF.Proc.linx.";
static SPProcFunction_Methods *SPPF_LinxMths = NULL;

double SPPF_LinxEval(T *f, R3Point *p);
R3Gradient LinxGrad(T *f, R3Point *p);
R3Hessian SPPF_LinxHess(T *f, R3Point *p);

double SPPF_LinxEval(T *f, R3Point *p)
  { double s = f->d->scale, x = p->c[0];
    return s*x;
  }
  
R3Gradient LinxGrad(T *f, R3Point *p)
  { double s = f->d->scale;
    return (r3_t){{s, 0.0, 0.0}};
  }

R3Hessian SPPF_LinxHess(T *f, R3Point *p)
  { return (R3Hessian){{0.0, 0.0, 0.0, 0.0, 0.0, 0.0}};
  }
    
/*** {liny(p)} **********************************/

static char *SPPF_LinyDesc = "y";
static char *SPPF_LinyType = "SF.Proc.liny.";
static SPProcFunction_Methods *SPPF_LinyMths = NULL;

double SPPF_LinyEval(T *f, R3Point *p);
R3Gradient LinyGrad(T *f, R3Point *p);
R3Hessian SPPF_LinyHess(T *f, R3Point *p);

double SPPF_LinyEval(T *f, R3Point *p)
  { double s = f->d->scale, y = p->c[1];
    return s*y;
  }
  
R3Gradient LinyGrad(T *f, R3Point *p)
  { double s = f->d->scale;
    return (r3_t){{0.0, s, 0.0}};
  }

R3Hessian SPPF_LinyHess(T *f, R3Point *p)
  { return (R3Hessian){{0.0, 0.0, 0.0, 0.0, 0.0, 0.0}};
  }
    
/*** {linz(p)} **********************************/

static char *SPPF_LinzDesc = "z";
static char *SPPF_LinzType = "SF.Proc.linz.";
static SPProcFunction_Methods *SPPF_LinzMths = NULL;

double SPPF_LinzEval(T *f, R3Point *p);
R3Gradient LinzGrad(T *f, R3Point *p);
R3Hessian SPPF_LinzHess(T *f, R3Point *p);

double SPPF_LinzEval(T *f, R3Point *p)
  { double s = f->d->scale, z = p->c[2];
    return s*z;
  }
  
R3Gradient LinzGrad(T *f, R3Point *p)
  { double s = f->d->scale;
    return (r3_t){{0.0, 0.0, s}};
  }

R3Hessian SPPF_LinzHess(T *f, R3Point *p)
  { return (R3Hessian){{0.0, 0.0, 0.0, 0.0, 0.0, 0.0}};
  }
    
/*** {sqrx(p)} **********************************/

static char *SPPF_SqrxDesc = "x^2";
static char *SPPF_SqrxType = "SF.Proc.sqrx.";
static SPProcFunction_Methods *SPPF_SqrxMths = NULL;

double SPPF_SqrxEval(T *f, R3Point *p);
R3Gradient SqrxGrad(T *f, R3Point *p);
R3Hessian SPPF_SqrxHess(T *f, R3Point *p);

double SPPF_SqrxEval(T *f, R3Point *p)
  { double s = f->d->scale, x = p->c[0];
    return s * x*x;
  }
  
R3Gradient SqrxGrad(T *f, R3Point *p)
  { double s = f->d->scale, x = p->c[0];
    return (r3_t){{2.0*s*x, 0.0, 0.0}};
  } 

R3Hessian SPPF_SqrxHess(T *f, R3Point *p)
  { double s = f->d->scale;
    return (R3Hessian){{s*2.0, 0.0, 0.0, 0.0, 0.0, 0.0}};
  }

/*** {sqry(p)} **********************************/

static char *SPPF_SqryDesc = "y^2";
static char *SPPF_SqryType = "SF.Proc.sqry.";
static SPProcFunction_Methods *SPPF_SqryMths = NULL;

double SPPF_SqryEval(T *f, R3Point *p);
R3Gradient SqryGrad(T *f, R3Point *p);
R3Hessian SPPF_SqryHess(T *f, R3Point *p);

double SPPF_SqryEval(T *f, R3Point *p)
  { double s = f->d->scale, y = p->c[1];
    return s*y*y;
  }
  
R3Gradient SqryGrad(T *f, R3Point *p)
  { double s = f->d->scale, y = p->c[1];
    return (r3_t){{ 0.0, 2.0*s*y, 0.0}};
  }

R3Hessian SPPF_SqryHess(T *f, R3Point *p)
  { double s = f->d->scale;
    return (R3Hessian){{0.0, 0.0, s*2.0, 0.0, 0.0, 0.0}};
  }  

/*** {cubx(p)} **********************************/

static char *SPPF_CubxDesc = "x^3";
static char *SPPF_CubxType = "SF.Proc.cubx.";
static SPProcFunction_Methods *SPPF_CubxMths = NULL;

double SPPF_CubxEval(T *f, R3Point *p);
R3Gradient CubxGrad(T *f, R3Point *p);
R3Hessian SPPF_CubxHess(T *f, R3Point *p);

double SPPF_CubxEval(T *f, R3Point *p)
  { double s = f->d->scale, x = p->c[0];
    return s*x*x*x;
  }
  
R3Gradient CubxGrad(T *f, R3Point *p)
  { double s = f->d->scale, x = p->c[0];
    return (r3_t){{3.0*s*x*x, 0.0, 0.0}};
  }

R3Hessian SPPF_CubxHess(T *f, R3Point *p)
  { double s = f->d->scale, x = p->c[0];
    return (R3Hessian){{s*6.0*x, 0.0, 0.0, 0.0, 0.0, 0.0}};
  }

/*** {qrtx(p)} **********************************/

static char *SPPF_QrtxDesc = "x^4";
static char *SPPF_QrtxType = "SF.Proc.qrtx.";
static SPProcFunction_Methods *SPPF_QrtxMths = NULL;

double SPPF_QrtxEval(T *f, R3Point *p);
R3Gradient QrtxGrad(T *f, R3Point *p);
R3Hessian SPPF_QrtxHess(T *f, R3Point *p);

double SPPF_QrtxEval(T *f, R3Point *p)
  { double s = f->d->scale, x = p->c[0];
    double x2 = x*x;
    return s*x2*x2;
  }
  
R3Gradient QrtxGrad(T *f, R3Point *p)
  { double s = f->d->scale, x = p->c[0];
    return (r3_t){{4.0*s*x*x*x, 0.0, 0.0}};
  }

R3Hessian SPPF_QrtxHess(T *f, R3Point *p)
  { double s = f->d->scale, x = p->c[0];
    return (R3Hessian){{s*12.0*x*x, 0.0, 0.0, 0.0, 0.0, 0.0}};
  }
  
/*** {qrtz(p)} **********************************/

static char *SPPF_QrtzDesc = "z^4";
static char *SPPF_QrtzType = "SF.Proc.qrtz.";
static SPProcFunction_Methods *SPPF_QrtzMths = NULL;

double SPPF_QrtzEval(T *f, R3Point *p);
R3Gradient QrtzGrad(T *f, R3Point *p);
R3Hessian SPPF_QrtzHess(T *f, R3Point *p);

double SPPF_QrtzEval(T *f, R3Point *p)
  { double s = f->d->scale, z = p->c[2];
    double z2 = z*z;
    return s*z2*z2;
  }
  
R3Gradient QrtzGrad(T *f, R3Point *p)
  { double s = f->d->scale, z = p->c[2];
    return (r3_t){{0.0, 0.0,4.0*s*z*z*z}};
  }

R3Hessian SPPF_QrtzHess(T *f, R3Point *p)
  { double s = f->d->scale, z = p->c[2];
    return (R3Hessian){{0.0, 0.0, 0.0, 0.0, 0.0, s*12.0*z*z}};
  }

/*** {quix(p)} **********************************/

static char *SPPF_QuixDesc = "x^5";
static char *SPPF_QuixType = "SF.Proc.quix.";
static SPProcFunction_Methods *SPPF_QuixMths = NULL;

double SPPF_QuixEval(T *f, R3Point *p);
R3Gradient QuixGrad(T *f, R3Point *p);
R3Hessian SPPF_QuixHess(T *f, R3Point *p);

double SPPF_QuixEval(T *f, R3Point *p)
  { double s = f->d->scale, x = p->c[0];
    double x2 = x*x, x3 = x*x2;
    return s*x2*x3;
  }
  
R3Gradient QuixGrad(T *f, R3Point *p)
  { double s = f->d->scale, x = p->c[0];
    double x2 = x*x;
    return (r3_t){{5.0*s*x2*x2, 0.0, 0.0}};
  }

R3Hessian SPPF_QuixHess(T *f, R3Point *p)
  { double s = f->d->scale, x = p->c[0];
    double x2 = x*x, x3 = x*x2;
    return (R3Hessian){{s*20.0*x3, 0.0, 0.0, 0.0, 0.0, 0.0}};
  }

/*** {sexx(p)} **********************************/

static char *SPPF_SexxDesc = "x^6";
static char *SPPF_SexxType = "SF.Proc.sexx.";
static SPProcFunction_Methods *SPPF_SexxMths = NULL;

double SPPF_SexxEval(T *f, R3Point *p);
R3Gradient SexxGrad(T *f, R3Point *p);
R3Hessian SPPF_SexxHess(T *f, R3Point *p);

double SPPF_SexxEval(T *f, R3Point *p)
  { double s = f->d->scale, x = p->c[0];
    double x2 = x*x, x3 = x*x2;
    return s*x3*x3;
  }
  
R3Gradient SexxGrad(T *f, R3Point *p)
  { double s = f->d->scale, x = p->c[0];
    double x2 = x*x, x3 = x*x2;
    return (r3_t){{6.0*s*x2*x3, 0.0, 0.0}};
  }

R3Hessian SPPF_SexxHess(T *f, R3Point *p)
  { double s = f->d->scale, x = p->c[0];
    double x2 = x*x, x4 = x2*x2;
    return (R3Hessian){{s*30.0*x4, 0.0, 0.0, 0.0, 0.0, 0.0}};
  }

/*** {sepx(p)} **********************************/

static char *SPPF_SepxDesc = "x^7";
static char *SPPF_SepxType = "SF.Proc.sepx.";
static SPProcFunction_Methods *SPPF_SepxMths = NULL;

double SPPF_SepxEval(T *f, R3Point *p);
R3Gradient SepxGrad(T *f, R3Point *p);
R3Hessian SPPF_SepxHess(T *f, R3Point *p);

double SPPF_SepxEval(T *f, R3Point *p)
  { double s = f->d->scale, x = p->c[0];
    double x2 = x*x, x3 = x*x2, x4 = x2*x2;
    return s*x3*x4;
  }
  
R3Gradient SepxGrad(T *f, R3Point *p)
  { double s = f->d->scale, x = p->c[0];
    double  x2 = x*x, x3 = x*x2;
    return (r3_t){{7.0*s*x3*x3, 0.0, 0.0}};
  }

R3Hessian SPPF_SepxHess(T *f, R3Point *p)
  { double s = f->d->scale, x = p->c[0];
    double x2 = x*x, x3 = x*x2, x5 = x3*x2;
    return (R3Hessian){{s*42.0*x5, 0.0, 0.0, 0.0, 0.0, 0.0}};
  }

/*** {octx(p)} **********************************/

static char *SPPF_OctxDesc = "x^8";
static char *SPPF_OctxType = "SF.Proc.octx.";
static SPProcFunction_Methods *SPPF_OctxMths = NULL;

double SPPF_OctxEval(T *f, R3Point *p);
R3Gradient OctxGrad(T *f, R3Point *p);
R3Hessian SPPF_OctxHess(T *f, R3Point *p);

double SPPF_OctxEval(T *f, R3Point *p)
  { double s = f->d->scale, x = p->c[0];
    double x2 = x*x, x4 = x2*x2;
    return s*x4*x4;
  }
  
R3Gradient OctxGrad(T *f, R3Point *p)
  { double s = f->d->scale, x = p->c[0];
    double x2 = x*x, x3 = x*x2, x4 = x2*x2;
    return (r3_t){{8.0*s*x3*x4, 0.0, 0.0}};
  }

R3Hessian SPPF_OctxHess(T *f, R3Point *p)
  { double s = f->d->scale, x = p->c[0];
    double x2 = x*x, x3 = x*x2, x6 = x3*x3;
    return (R3Hessian){{s*56.0*x6, 0.0, 0.0, 0.0, 0.0, 0.0}};
  }
  
/*** {band(p)} **********************************/

static char *SPPF_BandDesc = "1-x^2";
static char *SPPF_BandType = "SF.Proc.band.";
static SPProcFunction_Methods *SPPF_BandMths = NULL;

double SPPF_BandEval(T *f, R3Point *p);
R3Gradient BandGrad(T *f, R3Point *p);
R3Hessian SPPF_BandHess(T *f, R3Point *p);

double SPPF_BandEval(T *f, R3Point *p)
  { double s = f->d->scale, x = p->c[0];
    return s*(1.0 - x*x);
  }
  
R3Gradient BandGrad(T *f, R3Point *p)
  { double s = f->d->scale, x = p->c[0];
    return (r3_t){{-2.0*s*x, 0.0, 0.0}};
  }

R3Hessian SPPF_BandHess(T *f, R3Point *p)
  { double s = f->d->scale;
    return (R3Hessian){{-s*2.0, 0.0, 0.0, 0.0, 0.0, 0.0}};
  }

/*** {hrm2(p)} **********************************/

/* This is an eigenfunction of the spherical laplacian
 operator, with eigenvalue -6. */

static char *SPPF_Hrm2Desc = "x^2-1/3";
static char *SPPF_Hrm2Type = "SF.Proc.hrm2.";
static SPProcFunction_Methods *SPPF_Hrm2Mths = NULL;

double SPPF_Hrm2Eval(T *f, R3Point *p);
R3Gradient Hrm2Grad(T *f, R3Point *p);
R3Hessian SPPF_Hrm2Hess(T *f, R3Point *p);

double SPPF_Hrm2Eval(T *f, R3Point *p)
  { double s = f->d->scale, x = p->c[0];
    return s*(x*x - 1.0/3.0);
  }
  
R3Gradient Hrm2Grad(T *f, R3Point *p)
  { double s = f->d->scale, x = p->c[0];
    return (r3_t){{2.0*s*x, 0.0, 0.0}};
  }

R3Hessian SPPF_Hrm2Hess(T *f, R3Point *p)
  { double s = f->d->scale;
    return (R3Hessian){{s*2.0, 0.0, 0.0, 0.0, 0.0, 0.0}};
  }

/*** {hr75(p)} **********************************/

/* This is an harmonic term with {d = 7} and {m = 5}.
  It is an eigenfunction of the spherical laplacian
  operator, with eigenvalue ??.*/

static char *SPPF_Hr75Desc = "x^2-1/3";
static char *SPPF_Hr75Type = "SF.Proc.hr75.";
static SPProcFunction_Methods *SPPF_Hr75Mths = NULL;

#define SPPF_Hr75_d 7
#define SPPF_Hr75_m 5

static 

double SPPF_Hr75Eval(T *f, R3Point *p);
R3Gradient Hr75Grad(T *f, R3Point *p);
R3Hessian SPPF_Hr75Hess(T *f, R3Point *p);

double SPPF_Hr75Eval(T *f, R3Point *p)
  { double s = f->d->scale, x = p->c[0], y = p->c[1], z = p->c[2];
    HarmonicTerm t = (HarmonicTerm){ SPPF_Hr75_d, SPPF_Hr75_m, s };
    double clon = x, slon = y;
    double r = hypot(clon, slon);
    if (r > 0) { clon /= r; slon /= r; }

    double slat = z;
    double R = hypot(r, slat);
    if (R > 0) { slat /= R; }
    double v;
    SPHarmonic_EvalTerm(&t, clon, slon, slat, R, 0, &v, NULL, NULL);
    return v;
  }
  
R3Gradient Hr75Grad(T *f, R3Point *p)
  { double s = f->d->scale, x = p->c[0], y = p->c[1], z = p->c[2];
    HarmonicTerm t = (HarmonicTerm){ SPPF_Hr75_d, SPPF_Hr75_m, s };
    double clon = x, slon = y;
    double r = hypot(clon, slon);
    if (r > 0) { clon /= r; slon /= r; }

    double slat = z;
    double R = hypot(r, slat);
    if (R > 0) { slat /= R; }
    double v;
    r3_t dv;
    SPHarmonic_EvalTerm(&t, clon, slon, slat, R, 1, &v, &dv, NULL);
    return dv;
  }

R3Hessian SPPF_Hr75Hess(T *f, R3Point *p)
  { double s = f->d->scale, x = p->c[0], y = p->c[1], z = p->c[2];
    HarmonicTerm t = (HarmonicTerm){ SPPF_Hr75_d, SPPF_Hr75_m, s };
    double clon = x, slon = y;
    double r = hypot(clon, slon);
    if (r > 0) { clon /= r; slon /= r; }

    double slat = z;
    double R = hypot(r, slat);
    if (R > 0) { slat /= R; }
    
    double v;
    r3_t dv;
    r6_t ddv;
    SPHarmonic_EvalTerm(&t, clon, slon, slat, R, 1, &v, &dv, &ddv);
    return ddv;
  }

/*** {loca(p)} **********************************/

static char *SPPF_LocaDesc = "x*(2x^2-1)*(3x^2-2)";
static char *SPPF_LocaType = "SF.Proc.loca.";
static SPProcFunction_Methods *SPPF_LocaMths = NULL;

double SPPF_LocaEval(T *f, R3Point *p);
R3Gradient LocaGrad(T *f, R3Point *p);
R3Hessian SPPF_LocaHess(T *f, R3Point *p);

double SPPF_LocaEval(T *f, R3Point *p)
  { double s = f->d->scale, x = p->c[0];
    double x2 = x*x;
    return s*x*(2.0*x2 - 1.0)*(3.0*x2 - 2.0);
  }
  
R3Gradient LocaGrad(T *f, R3Point *p)
  { double s = f->d->scale, x = p->c[0];
    double x2 = x*x, x4 = x2*x2;
    return (r3_t){{s*(30.0*x4 - 21.0*x2 + 2.0), 0.0, 0.0}};
  }

R3Hessian SPPF_LocaHess(T *f, R3Point *p)
  { double s = f->d->scale, x = p->c[0];
    double x3 = x*x*x;
    return (R3Hessian){{s*(120.0*x3 - 42.0*x), 0.0, 0.0, 0.0, 0.0, 0.0}};
  }

/*** {fxyz(p)} **********************************/

static char *SPPF_FxyzDesc = "x*y*z";
static char *SPPF_FxyzType = "SF.Proc.fxyz.";
static SPProcFunction_Methods *SPPF_FxyzMths = NULL;

double SPPF_FxyzEval(T *f, R3Point *p);
R3Gradient FxyzGrad(T *f, R3Point *p);
R3Hessian SPPF_FxyzHess(T *f, R3Point *p);

double SPPF_FxyzEval(T *f, R3Point *p)
  { double s = f->d->scale, x = p->c[0], y = p->c[1], z = p->c[2];
    return s*x*y*z;
  }
  
R3Gradient FxyzGrad(T *f, R3Point *p)
  { double s = f->d->scale, x = p->c[0], y = p->c[1], z = p->c[2];
    return (r3_t){{s*y*z, s*x*z, s*x*y}};
  }

R3Hessian SPPF_FxyzHess(T *f, R3Point *p)
  { double s = f->d->scale, x = p->c[0], y = p->c[1], z = p->c[2];
    return (R3Hessian){{0.0, s*z, 0.0, s*y, s*x, 0.0}};
  }

/*** {expx(p)} **********************************/

static char *SPPF_ExpxDesc = "exp(x)";
static char *SPPF_ExpxType = "SF.Proc.expx.";
static SPProcFunction_Methods *SPPF_ExpxMths = NULL;

double SPPF_ExpxEval(T *f, R3Point *p);
R3Gradient ExpxGrad(T *f, R3Point *p);
R3Hessian SPPF_ExpxHess(T *f, R3Point *p);

double SPPF_ExpxEval(T *f, R3Point *p)
  { double s = f->d->scale, x = p->c[0];
    return s*exp(x);
  }
  
R3Gradient ExpxGrad(T *f, R3Point *p)
  { double s = f->d->scale, x = p->c[0];
    return (r3_t){{s*exp(x), 0.0, 0.0}};
  }

R3Hessian SPPF_ExpxHess(T *f, R3Point *p)
  { double s = f->d->scale, x = p->c[0];
    return (R3Hessian){{s*exp(x), 0.0, 0.0, 0.0, 0.0, 0.0}};
  }
  
/*** {expz(p)} **********************************/

static char *SPPF_ExpzDesc = "exp(z)";
static char *SPPF_ExpzType = "SF.Proc.expz.";
static SPProcFunction_Methods *SPPF_ExpzMths = NULL;

double SPPF_ExpzEval(T *f, R3Point *p);
R3Gradient ExpzGrad(T *f, R3Point *p);
R3Hessian SPPF_ExpzHess(T *f, R3Point *p);

double SPPF_ExpzEval(T *f, R3Point *p)
  { double s = f->d->scale,z = p->c[2];
    return s*exp(z);
  }
  
R3Gradient ExpzGrad(T *f, R3Point *p)
  { double s = f->d->scale, z = p->c[2];
    return (r3_t){{0.0, 0.0, s*exp(z)}};
  }

R3Hessian SPPF_ExpzHess(T *f, R3Point *p)
  { double s = f->d->scale, z = p->c[2];
    return (R3Hessian){{0.0, 0.0, 0.0,  0.0, 0.0, s*exp(z)}};
  }  

/*** {sinx(p)} **********************************/

static char *SPPF_SinxDesc = "sin(x)";
static char *SPPF_SinxType = "SF.Proc.sinx.";
static SPProcFunction_Methods *SPPF_SinxMths = NULL;

double SPPF_SinxEval(T *f, R3Point *p);
R3Gradient SinxGrad(T *f, R3Point *p);
R3Hessian SPPF_SinxHess(T *f, R3Point *p);

double SPPF_SinxEval(T *f, R3Point *p)
  { double s = f->d->scale, x = p->c[0];
    return s*sin(x);
  }
  
R3Gradient SinxGrad(T *f, R3Point *p)
  { double s = f->d->scale, x = p->c[0];
    return (r3_t){{s*cos(x), 0.0, 0.0}};
  }

R3Hessian SPPF_SinxHess(T *f, R3Point *p)
  { double s = f->d->scale, x = p->c[0];
    return (R3Hessian){{-s*sin(x), 0.0, 0.0, 0.0, 0.0, 0.0}};
  }

/*** {cosx(p)} **********************************/

static char *SPPF_CosxDesc = "cos(x)";
static char *SPPF_CosxType = "SF.Proc.cosx.";
static SPProcFunction_Methods *SPPF_CosxMths = NULL;

double SPPF_CosxEval(T *f, R3Point *p);
R3Gradient CosxGrad(T *f, R3Point *p);
R3Hessian SPPF_CosxHess(T *f, R3Point *p);

double SPPF_CosxEval(T *f, R3Point *p)
  { double s = f->d->scale, x = p->c[0];
    return s*cos(x);
  }
  
R3Gradient CosxGrad(T *f, R3Point *p)
  { double s = f->d->scale, x = p->c[0];
    return (r3_t){{-s*sin(x), 0.0, 0.0}};
  }

R3Hessian SPPF_CosxHess(T *f, R3Point *p)
  { double s = f->d->scale, x = p->c[0];
    return (R3Hessian){{-s*cos(x), 0.0, 0.0, 0.0, 0.0, 0.0}};
  }  

/*** {lrry(p)} **********************************/

static char *SPPF_LrryDesc = "LarrysFn(p)";
static char *SPPF_LrryType = "SF.Proc.lrry.";
static SPProcFunction_Methods *SPPF_LrryMths = NULL;

double SPPF_LrryEval(T *f, R3Point *p);
R3Gradient LrryGrad(T *f, R3Point *p);
R3Hessian SPPF_LrryHess(T *f, R3Point *p);

/* Scale Larry's example func so that its range is inside [0 _ 1]: */
/* #define LrryMax 0.389056 */
#define LrrySup (0.4)

double SPPF_LrryEval(T *f, R3Point *p)
  { double s = f->d->scale/LrrySup, x = p->c[0], y = p->c[1], z = p->c[2];
    double x2 = x*x, x4 = x2*x2, x8 = x4*x4;
    double y3 = y*y*y;
    double z2 = z*z;
    double fv = 1.0 + x8 + exp(2.0*y3) + exp(2.0*z2) + 10.0*x*y*z;
    return s*fv;
  }
  
R3Gradient LrryGrad(T *f, R3Point *p)
  { double s = f->d->scale, x = p->c[0], y = p->c[1], z = p->c[2];
    double x3 = x*x*x, x7 = x3*x3*x;
    double y2 = y*y, y3 = y2*y;
    double z2 = z*z;
    s /= LrrySup;
    return (r3_t)
      { { s*(8.0*x7 + 10.0*y*z),
          s*(6.0*y2*exp(2.0*y3) + 10.0*x*z),
          s*(4.0*z*exp(2.0*z2) + 10.0*x*y)
        }
      };
  }

R3Hessian SPPF_LrryHess(T *f, R3Point *p)
  { double s = f->d->scale/LrrySup, x = p->c[0], y = p->c[1], z = p->c[2];
    double x3 = x*x*x, x6 = x3*x3;
    double y3 = y*y*y, y4 = y3*y;
    double z2 = z*z;
    return (R3Hessian)
      { { s*56.0*x6, 
          s*10.0*z, s*(12.0*y + 36.0*y4)*exp(2.0*y3),
          s*10.0*y, s*10.0*x, s*(4.0 + 4.0*z2)*exp(2.0*z2)
        }
      };
  }  
  
/*** {spir(p)} **********************************/

static char *SPPF_SpirDesc = "SpiralFn(p)";
static char *SPPF_SpirType = "SF.Proc.spir.";
static SPProcFunction_Methods *SPPF_SpirMths = NULL;

double SPPF_SpirEval(T *f, R3Point *p);
R3Gradient SpirGrad(T *f, R3Point *p);
R3Hessian SPPF_SpirHess(T *f, R3Point *p);

#define SpirC (2*PI/16)
#define SpirU 1.0625

double SPPF_SpirEval(T *f, R3Point *p)
  { double s = f->d->scale, x = p->c[0], y = p->c[1], z = p->c[2];
    double h = x/SpirU, h2 = h*h, h4 = h2*h2, h8 = h4*h4, h10=h8*h2;
    double m = 1.0 - h2;
    double q = 1.0 - h10;
    double t = 2*h/m;
    double ct = (fabs(t) > 1.0e6 ? 0.0 : cos(SpirC*t));
    double st = (fabs(t) > 1.0e6 ? 0.0 : sin(SpirC*t));
    double fv = q*(ct*y - st*z);
    return s*fv;
  }
  
R3Gradient SpirGrad(T *f, R3Point *p)
  { double s = f->d->scale, x = p->c[0], y = p->c[1], z = p->c[2];
    double h = x/SpirU, h2 = h*h, h4 = h2*h2, h8 = h4*h4, h10=h8*h2;
    double m = 1.0 - h2;
    double q = 1.0 - h10;
    double t = 2*h/m;
    double ct = (fabs(t) > 1.0e6 ? 0.0 : cos(SpirC*t));
    double st = (fabs(t) > 1.0e6 ? 0.0 : sin(SpirC*t));
    /* double fv = q*(ct*y - st*z); */
    double DfvDst = -q*z;
    double DfvDct = q*y;
    double DfvDt = (fabs(t) > 1.0e6 ? 0.0 : DfvDst*ct - DfvDct*st)*SpirC;
    double DfvDq = ct*y - st*z;
    double DfvDm = DfvDt*(-2*h/(m*m));
    double DfvDh = -h*(DfvDm*2.0 + DfvDq*10.0*h8) + DfvDt*2.0/m;
    double DfvDx = DfvDh/SpirU;
    double DfvDy = q*ct;
    double DfvDz = -q*st;
    return (r3_t){{ s*DfvDx, s*DfvDy, s*DfvDz }};
  }

R3Hessian SPPF_SpirHess(T *f, R3Point *p)
  { affirm(FALSE, "hess method not implemented yet"); 
    return (R3Hessian){{0.0, 0.0, 0.0, 0.0, 0.0, 0.0}};
  }  
  
/*** {rain(p)} **********************************/

static char *SPPF_RainDesc = "Raindrop(p)";
static char *SPPF_RainType = "SF.Proc.rain.";
static SPProcFunction_Methods *SPPF_RainMths = NULL;

double SPPF_RainEval(T *f, R3Point *p);
R3Gradient RainGrad(T *f, R3Point *p);
R3Hessian SPPF_RainHess(T *f, R3Point *p);

#define RainFmax 16.0
#define RainFmin 1.0

double SPPF_RainEval(T *f, R3Point *p)
  { double s = f->d->scale, x = p->c[0], y = p->c[1], z = p->c[2];
    double ct = (x + y + z)/SQRT3;
    double g = fabs(1.0 - ct)/2.0;
    double h = sqrt(g);
    double u = -12*h*h;
    double m = exp(u);
    double w = 8*PI*h;
    double fv = m*cos(w);
    return s*fv;
  }
  
R3Gradient RainGrad(T *f, R3Point *p)
  { double s = f->d->scale, x = p->c[0], y = p->c[1], z = p->c[2];
    double ct = (x + y + z)/SQRT3;
    double g = fabs(1.0 - ct)/2.0;
    double h = sqrt(g);
    double u = -12*h*h;
    double m = exp(u);
    double w = 8*PI*h;
    /* double fv = m*cos(w); */
    double DfvDm = cos(w);
    double DfvDw = -m*sin(w);
    double DfvDu = DfvDm*m;
    double DfvDh = 8*PI*DfvDw -24*h*DfvDu;
    double DfvDg = DfvDh*0.5/h;
    double DfvDct = -0.5*DfvDg;
    double DfvDx = DfvDct/SQRT3;
    double DfvDy = DfvDct/SQRT3;
    double DfvDz = DfvDct/SQRT3;
    return (r3_t){{ s*DfvDx, s*DfvDy, s*DfvDz }};
  }

R3Hessian SPPF_RainHess(T *f, R3Point *p)
  { affirm(FALSE, "hess method not implemented yet"); 
    return (R3Hessian){{0.0, 0.0, 0.0, 0.0, 0.0, 0.0}};
  }  
  
/*************************************************/ 

#define SPProcFunction_FileFormat "2002-11-12"

/* OVERRIDES FOR PARENT CLASS METHODS */
  
void SPProcFunction_M_Add(T *f, double a, T *h)
  { if (strcmp(f->type, h->type) != 0)
      { fprintf (stderr, "copy: type mismatch: \"%s\", \"%s\"\n", 
          f->type, h->type);
        affirm(FALSE, "aborted");
      }
    f->d->scale += h->d->scale;
  }

void SPProcFunction_M_Scale(T *f, double a)
  { f->d->scale *= a;
  }

void SPProcFunction_M_Maple(T *f, FILE *wr)
  { affirm(FALSE , "maple method not implemented yet");
  }

SPProcFunction *SPProcFunction_M_Copy(T *f)
  { SPProcFunction *fnew = SPProcFunction_New();
    fnew->d = SPProcFunction_Data_New();
    *(fnew->d) = *(f->d);
    fnew->type = f->type;
    fnew->m = f->m;
    return fnew;
  }
 
void SPProcFunction_M_Free(T *f)
  { affirm(isprefix(SPProcFunction_TypeId, f->type), "type/method bug");
    free(f->d);
    free(f);
  }

/* CLASS-SPECIFIC METHODS */
  
void SPProcFunction_M_Write(T *f, FILE *wr)
  { filefmt_write_header(wr, "SPProcFunction", SPProcFunction_FileFormat);
    fprintf(wr, "type = %s\n",  f->type);
    fprintf(wr, "descr = %s\n",  f->d->descr);
    fprintf(wr, "scale = %.16g\n",  f->d->scale);
    filefmt_write_footer(wr, "SPProcFunction");
    fflush(wr);
  }

/* OTHER PROCS */
  
SPProcFunction_Methods *SPProcFunction_Methods_New(void)
  { void *v = malloc(sizeof(SPProcFunction_Methods));
    return (SPProcFunction_Methods *)notnull(v, "out of mem for SPProcFunction_Methods");
  }

SPProcFunction_Data *SPProcFunction_Data_New(void)
  { void *v = malloc(sizeof(SPProcFunction_Data));
    return (SPProcFunction_Data *)notnull(v, "out of mem for SPProcFunction_Data");
  }

SPProcFunction *SPProcFunction_New(void)
  { void *v = malloc(sizeof(SPProcFunction));
    return (SPProcFunction *)notnull(v, "no mem for SPProcFunction");
  }

SPProcFunction *SPProcFunction_Make(
    char *type,
    char *descr,
    SPProcFunction_Methods **mp,
    double (*eval)(T *f, R3Point *p),
    R3Gradient (*grad)(T *f, R3Point *p),
    R3Hessian (*hess)(T *f, R3Point *p)
)
{
  SPProcFunction *f = SPProcFunction_New();
  f->type = type;
  if ((*mp) == NULL)
    { SPProcFunction_Methods *m = SPProcFunction_Methods_New();
      /* Superclass methods: */
      m->fn.eval = (SPFunction_EvalMth *)eval;
      m->fn.grad = (SPFunction_GradMth *)grad;
      m->fn.hess = (SPFunction_HessMth *)hess;
      m->fn.maple = (SPFunction_MapleMth *)&SPProcFunction_M_Maple;
      /* Note: the {fn.write} method is inherited from {SPFunction}! */
      m->fn.write = (SPFunction_WriteMth *)&SPFunction_M_Write;
      m->fn.scale = (SPFunction_ScaleMth *)&SPProcFunction_M_Scale;
      m->fn.add = (SPFunction_AddMth *)&SPProcFunction_M_Add;
      m->fn.copy = (SPFunction_CopyMth *)&SPProcFunction_M_Copy;
      m->fn.free = (SPFunction_FreeMth *)&SPProcFunction_M_Free;
      /* Class-specific methods */
      m->write = (SPFunction_WriteMth *)&SPProcFunction_M_Write;
      (*mp) = m;
    }
  f->m = (*mp);
  f->d = SPProcFunction_Data_New();
  f->d->scale = 1.0;
  f->d->descr = descr;
  f->type = type;
  return f;
}

#define MkFn SPProcFunction_Make

SPProcFunction *SPProcFunction_FromName(char *name)
  { SPProcFunction *f = SPProcFunction_New();
    if (strcmp(name, "unit") == 0)
      { f = MkFn(SPPF_UnitType, SPPF_UnitDesc, &SPPF_UnitMths, SPPF_UnitEval, UnitGrad, SPPF_UnitHess); }
    else if (strcmp(name, "linx") == 0)
      { f = MkFn(SPPF_LinxType, SPPF_LinxDesc, &SPPF_LinxMths, SPPF_LinxEval, LinxGrad, SPPF_LinxHess); }
    else if (strcmp(name, "liny") == 0)
      { f = MkFn(SPPF_LinyType, SPPF_LinyDesc, &SPPF_LinyMths, SPPF_LinyEval, LinyGrad, SPPF_LinyHess); }
    else if (strcmp(name, "linz") == 0)
      { f = MkFn(SPPF_LinzType, SPPF_LinzDesc, &SPPF_LinzMths, SPPF_LinzEval, LinzGrad, SPPF_LinzHess); }
    else if (strcmp(name, "sqrx") == 0)
      { f = MkFn(SPPF_SqrxType, SPPF_SqrxDesc, &SPPF_SqrxMths, SPPF_SqrxEval, SqrxGrad, SPPF_SqrxHess); }
    else if (strcmp(name, "sqry") == 0)
      { f = MkFn(SPPF_SqryType, SPPF_SqryDesc, &SPPF_SqryMths, SPPF_SqryEval, SqryGrad, SPPF_SqryHess); }
    else if (strcmp(name, "cubx") == 0)
      { f = MkFn(SPPF_CubxType, SPPF_CubxDesc, &SPPF_CubxMths, SPPF_CubxEval, CubxGrad, SPPF_CubxHess); }
    else if (strcmp(name, "qrtx") == 0)
      { f = MkFn(SPPF_QrtxType, SPPF_QrtxDesc, &SPPF_QrtxMths, SPPF_QrtxEval, QrtxGrad, SPPF_QrtxHess); }
    else if (strcmp(name, "qrtz") == 0)
      { f = MkFn(SPPF_QrtzType, SPPF_QrtzDesc, &SPPF_QrtzMths, SPPF_QrtzEval, QrtzGrad, SPPF_QrtzHess); }
    else if (strcmp(name, "quix") == 0)
      { f = MkFn(SPPF_QuixType, SPPF_QuixDesc, &SPPF_QuixMths, SPPF_QuixEval, QuixGrad, SPPF_QuixHess); }
    else if (strcmp(name, "sexx") == 0)
      { f = MkFn(SPPF_SexxType, SPPF_SexxDesc, &SPPF_SexxMths, SPPF_SexxEval, SexxGrad, SPPF_SexxHess); }
    else if (strcmp(name, "sepx") == 0)
      { f = MkFn(SPPF_SepxType, SPPF_SepxDesc, &SPPF_SepxMths, SPPF_SepxEval, SepxGrad, SPPF_SepxHess); }
    else if (strcmp(name, "octx") == 0)
      { f = MkFn(SPPF_OctxType, SPPF_OctxDesc, &SPPF_OctxMths, SPPF_OctxEval, OctxGrad, SPPF_OctxHess); }
    else if (strcmp(name, "expx") == 0)
      { f = MkFn(SPPF_ExpxType, SPPF_ExpxDesc, &SPPF_ExpxMths, SPPF_ExpxEval, ExpxGrad, SPPF_ExpxHess); }
    else if (strcmp(name, "expz") == 0)
      { f = MkFn(SPPF_ExpzType, SPPF_ExpzDesc, &SPPF_ExpzMths, SPPF_ExpzEval, ExpzGrad, SPPF_ExpzHess); }
    else if (strcmp(name, "band") == 0)
      { f = MkFn(SPPF_BandType, SPPF_BandDesc, &SPPF_BandMths, SPPF_BandEval, BandGrad, SPPF_BandHess); }
    else if (strcmp(name, "hrm2") == 0)
      { f = MkFn(SPPF_Hrm2Type, SPPF_Hrm2Desc, &SPPF_Hrm2Mths, SPPF_Hrm2Eval, Hrm2Grad, SPPF_Hrm2Hess); }
    else if (strcmp(name, "hr75") == 0)
      { f = MkFn(SPPF_Hr75Type, SPPF_Hr75Desc, &SPPF_Hr75Mths, SPPF_Hr75Eval, Hr75Grad, SPPF_Hr75Hess); }
    else if (strcmp(name, "loca") == 0)
      { f = MkFn(SPPF_LocaType, SPPF_LocaDesc, &SPPF_LocaMths, SPPF_LocaEval, LocaGrad, SPPF_LocaHess); }
    else if (strcmp(name, "fxyz") == 0)
      { f = MkFn(SPPF_FxyzType, SPPF_FxyzDesc, &SPPF_FxyzMths, SPPF_FxyzEval, FxyzGrad, SPPF_FxyzHess); }
    else if (strcmp(name, "sinx") == 0)
      { f = MkFn(SPPF_SinxType, SPPF_SinxDesc, &SPPF_SinxMths, SPPF_SinxEval, SinxGrad, SPPF_SinxHess); }
    else if (strcmp(name, "cosx") == 0)
      { f = MkFn(SPPF_CosxType, SPPF_CosxDesc, &SPPF_CosxMths, SPPF_CosxEval, CosxGrad, SPPF_CosxHess); }
    else if (strcmp(name, "lrry") == 0)
      { f = MkFn(SPPF_LrryType, SPPF_LrryDesc, &SPPF_LrryMths, SPPF_LrryEval, LrryGrad, SPPF_LrryHess); }
    else if (strcmp(name, "spir") == 0)
      { f = MkFn(SPPF_SpirType, SPPF_SpirDesc, &SPPF_SpirMths, SPPF_SpirEval, SpirGrad, SPPF_SpirHess); }
    else if (strcmp(name, "rain") == 0)
      { f = MkFn(SPPF_RainType, SPPF_RainDesc, &SPPF_RainMths, SPPF_RainEval, RainGrad, SPPF_RainHess); }
    else 
      { fprintf (stderr, "bad SPProcFunction function name = %s\n", name);
        affirm(FALSE, "aborted");
      };
    return f;
  }

SPProcFunction *SPProcFunction_Read(FILE *rd)
  { char *t, *descr;
    double scale;
    SPProcFunction *f;
    filefmt_read_header(rd, "SPProcFunction", SPProcFunction_FileFormat);
    t = nget_string(rd, "type"); fget_eol(rd);
    if (! isprefix("SF.Proc.", t))
      { fprintf (stderr, "bad SPProcFunction type = %7s\n", t);
        affirm(FALSE, "expected \"SF.Proc.\"");
      }
    /* Remove final "." from subtype, and look it up: */
    { char *funcName = txtcat(t + 8, "#");
      affirm(strlen(funcName) > 0, "missing function name");
      affirm(funcName[strlen(funcName)-2] == '.', "missing \".\" in func type");
      funcName[strlen(funcName)-2] = '\000';
      f = SPProcFunction_FromName(funcName);
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
    filefmt_read_footer(rd, "SPProcFunction");
    free(t);
    free(descr);
    return f;
  }

