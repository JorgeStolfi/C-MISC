/* See SOProcFunction.h */
/* Last edited on 2005-06-05 21:05:05 by stolfi */

#include <SOProcFunction.h>
#include <SOFunction.h>

#include <filefmt.h>
#include <fget.h>
#include <nget.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <affirm.h>
#include <nat.h>

#define T SOProcFunction 
    
void SOProcFunction_WriteMth(T *f, FILE *wr);
SOProcFunction *SOProcFunction_CopyMth(T *f);

SOProcFunction_Methods *SOProcFunction_Methods_New(void);
SOProcFunction_Data *SOProcFunction_Data_New(void);
SOProcFunction *SOProcFunction_New(void);

SOProcFunction *SOProcFunction_Cast(OBJ *f)
  { SOFunction *ff = (SOFunction *)f;
    if ((f != NULL) && isprefix(SOProcFunction_TypeId, ff->type))
      { return (SOProcFunction *)f; }
    else
      { return NULL; }
  }

SOProcFunction *SOProcFunction_Make
  ( char *type,
    char *descr,
    SOProcFunction_Methods **mp,
    void (*eval)(T *f, double *p, double *fp),
    void (*grad)(T *f, double *p, double *dfp),
    void (*hess)(T *f, double *p, double *ddfp)
  );
  /* Creates a new {SOProcFunction} object (and its data record),
    with given {type} and {descr} fields. 

    If {*mp != NULL}, uses {*mp} as the methods record, ignoring the
    other arguments. If {*mp == NULL}, allocates a new
    {SOProcFunction_Methods} record, with the given {eval}, {grad}, and
    {hess} methods, stores its address in {*mp}, and uses that for the
    new object. */

/* SPECIFIC FUNCTIONS */

/* Unless said otherwise, these functions set all components of the
  result vector to the same value.
  
  Unless said otherwise, the descriptions assume root-relative
  coordinates (where the root cell is the square {[0 _ 1] × [0 _ 1]}.
  A few functions use world-relative coordinates, where the root cell
  has length {(1/2)^{i/d}} along axis {i}. */

/*** unit(p) **********************************/

/* The constant "1" as a function. */

static char *UnitDesc = "1";
static char *UnitType = "SOF.Proc.unit.";
static SOProcFunction_Methods *UnitMths = NULL;

void UnitEval(T *f, double *p, double *fp);
void UnitGrad(T *f, double *p, double *dfp);
void UnitHess(T *f, double *p, double *ddfp);

void UnitEval(T *f, double *p, double *fp)
  { nat_t fDim = f->d->fn.fDim;
    int i;
    for (i = 0; i < fDim; i++) { fp[i] = 1.0; }
  }
  
void UnitGrad(T *f, double *p, double *dfp)
  { affirm(FALSE, "not implemented yet"); }

void UnitHess(T *f, double *p, double *ddfp)
  { affirm(FALSE, "not implemented yet"); }

/*** linx(p) **********************************/

/* A linear ramp in the {X} direction, with unit slope. */

static char *LinxDesc = "x";
static char *LinxType = "SOF.Proc.linx.";
static SOProcFunction_Methods *LinxMths = NULL;

void LinxEval(T *f, double *p, double *fp);
void LinxGrad(T *f, double *p, double *dfp);
void LinxHess(T *f, double *p, double *ddfp);

void LinxEval(T *f, double *p, double *fp)
  { nat_t fDim = f->d->fn.fDim;
    int i;
    for (i = 0; i < fDim; i++) { fp[i] = p[0]; }
  }
  
void LinxGrad(T *f, double *p, double *dfp)
  { affirm(FALSE, "not implemented yet"); }

void LinxHess(T *f, double *p, double *ddfp)
  { affirm(FALSE, "not implemented yet"); }
    
/*** liny(p) **********************************/

/* A linear ramp in the {Y} direction, with unit slope. */

static char *LinyDesc = "y";
static char *LinyType = "SOF.Proc.liny.";
static SOProcFunction_Methods *LinyMths = NULL;

void LinyEval(T *f, double *p, double *fp);
void LinyGrad(T *f, double *p, double *dfp);
void LinyHess(T *f, double *p, double *ddfp);

void LinyEval(T *f, double *p, double *fp)
  { nat_t pDim = f->d->fn.pDim;
    nat_t fDim = f->d->fn.fDim;
    int i;
    affirm(pDim >= 2, "domain dimension too small");
    for (i = 0; i < fDim; i++) { fp[i] = p[1]; }
  }
  
void LinyGrad(T *f, double *p, double *dfp)
  { affirm(FALSE, "not implemented yet"); }

void LinyHess(T *f, double *p, double *ddfp)
  { affirm(FALSE, "not implemented yet"); }

/*** linz(p) **********************************/

/* A linear ramp in the {Z} direction, with unit slope. */

static char *LinzDesc = "z";
static char *LinzType = "SOF.Proc.linz.";
static SOProcFunction_Methods *LinzMths = NULL;

void LinzEval(T *f, double *p, double *fp);
void LinzGrad(T *f, double *p, double *dfp);
void LinzHess(T *f, double *p, double *ddfp);

void LinzEval(T *f, double *p, double *fp)
  { nat_t pDim = f->d->fn.pDim;
    nat_t fDim = f->d->fn.fDim;
    int i;
    affirm(pDim >= 3, "domain dimension too small");
    for (i = 0; i < fDim; i++) { fp[i] = p[2]; }
  }
  
void LinzGrad(T *f, double *p, double *dfp)
  { affirm(FALSE, "not implemented yet"); }

void LinzHess(T *f, double *p, double *ddfp)
  { affirm(FALSE, "not implemented yet"); }
    
/*** lind(p) **********************************/

/* A linear ramp in the (1,1) direction (world), with unit slope. */

static char *LindDesc = "(x/sqrt(2)+y/2)";
static char *LindType = "SOF.Proc.lind.";
static SOProcFunction_Methods *LindMths = NULL;

void LindEval(T *f, double *p, double *fp);
void LindGrad(T *f, double *p, double *dfp);
void LindHess(T *f, double *p, double *ddfp);

void LindEval(T *f, double *p, double *fp)
  { nat_t pDim = f->d->fn.pDim;
    nat_t fDim = f->d->fn.fDim;
    affirm(pDim >= 2, "domain dimension too small");
    double x = p[0]; 
    double y = p[1];
    double d = (x*SQRTHALF + y*0.5);
    int i;
    for (i = 0; i < fDim; i++) { fp[i] = d; }
  }
  
void LindGrad(T *f, double *p, double *dfp)
  { affirm(FALSE, "not implemented yet"); }

void LindHess(T *f, double *p, double *ddfp)
  { affirm(FALSE, "not implemented yet"); }
    
/*** cbsm(p) **********************************/

/* A cubic function of {x,y} for ODE testing */

static char *CbsmDesc = "x^3+y^3+6*x+6*y";
static char *CbsmType = "SOF.Proc.cbsm.";
static SOProcFunction_Methods *CbsmMths = NULL;

void CbsmEval(T *f, double *p, double *fp);
void CbsmGrad(T *f, double *p, double *dfp);
void CbsmHess(T *f, double *p, double *ddfp);

void CbsmEval(T *f, double *p, double *fp)
  { nat_t pDim = f->d->fn.pDim;
    nat_t fDim = f->d->fn.fDim;
    int i;
    affirm(pDim >= 2, "domain dimension too small");
    for (i = 0; i < fDim; i++) 
      { fp[i] = (double) p[0]*p[0]*p[0] + (double) p[1]*p[1]*p[1]; }
  }
  
void CbsmGrad(T *f, double *p, double *dfp)
  { affirm(FALSE, "not implemented yet"); }

void CbsmHess(T *f, double *p, double *ddfp)
  { affirm(FALSE, "not implemented yet"); }
    
/*** sqrx(p) **********************************/

/* A parabolic ramp in the {X} direction. */

static char *SqrxDesc = "x^2";
static char *SqrxType = "SOF.Proc.sqrx.";
static SOProcFunction_Methods *SqrxMths = NULL;

void SqrxEval(T *f, double *p, double *fp);
void SqrxGrad(T *f, double *p, double *dfp);
void SqrxHess(T *f, double *p, double *ddfp);

void SqrxEval(T *f, double *p, double *fp)
  { nat_t fDim = f->d->fn.fDim;
    int i;
    for (i = 0; i < fDim; i++) { fp[i] = p[0]*p[0]; }
  }
  
void SqrxGrad(T *f, double *p, double *dfp)
  { affirm(FALSE, "not implemented yet"); } 

void SqrxHess(T *f, double *p, double *ddfp)
  { affirm(FALSE, "not implemented yet"); }

/*** para(p) **********************************/

/* An isotropic paraboloid centered at {(1/2,1/2,...)}: */

static char *ParaDesc = "|p-CTR|^2";
static char *ParaType = "SOF.Proc.para.";
static SOProcFunction_Methods *ParaMths = NULL;

void ParaEval(T *f, double *p, double *fp);
void ParaGrad(T *f, double *p, double *dfp);
void ParaHess(T *f, double *p, double *ddfp);

void ParaEval(T *f, double *p, double *fp)
  { nat_t pDim = f->d->fn.pDim;
    nat_t fDim = f->d->fn.fDim;
    int i;
    double r2 = 0.0;
    for (i = 0; i < pDim; i++) { double dxi = p[i] - 0.5; r2 += dxi*dxi; }
    for (i = 0; i < fDim; i++) { fp[i] = r2; }
  }
  
void ParaGrad(T *f, double *p, double *dfp)
  { affirm(FALSE, "not implemented yet"); } 

void ParaHess(T *f, double *p, double *ddfp)
  { affirm(FALSE, "not implemented yet"); }

/*** gaus(p) **********************************/

/* A  Gaussian hump with unit HEIGHT, centered at {(1/2,1/2,...)}: */

/* Standard deviation of Gaussian, projected on each axis: */
#define GAUSS_SIGMA 0.20

static char *GausDesc = "exp(-|p-CTR|^2/SGM^2)";
static char *GausType = "SOF.Proc.gaus.";
static SOProcFunction_Methods *GausMths = NULL;

void GausEval(T *f, double *p, double *fp);
void GausGrad(T *f, double *p, double *dfp);
void GausHess(T *f, double *p, double *ddfp);

void GausEval(T *f, double *p, double *fp)
  { nat_t pDim = f->d->fn.pDim;
    nat_t fDim = f->d->fn.fDim;
    int i;
    double s2 = GAUSS_SIGMA*GAUSS_SIGMA;
    double r2 = 0.0;
    for (i = 0; i < pDim; i++) 
      { double dxi = p[i] - 0.5; r2 += dxi*dxi; }
    double e2 = r2/s2;
    double ga = exp(-e2/2);
    for (i = 0; i < fDim; i++) { fp[i] = ga; }
  }
  
void GausGrad(T *f, double *p, double *dfp)
  { affirm(FALSE, "not implemented yet"); } 

void GausHess(T *f, double *p, double *ddfp)
  { affirm(FALSE, "not implemented yet"); }

/*** dgdx(p) **********************************/

/* The {x}-derivative of {gaus}: */

static char *DgDxDesc = "-2*(x-0.5)/SGM^2*exp(-|p-CTR|^2/SGM^2)";
static char *DgDxType = "SOF.Proc.dgdx.";
static SOProcFunction_Methods *DgDxMths = NULL;

void DgDxEval(T *f, double *p, double *fp);
void DgDxGrad(T *f, double *p, double *dfp);
void DgDxHess(T *f, double *p, double *ddfp);

void DgDxEval(T *f, double *p, double *fp)
  { nat_t pDim = f->d->fn.pDim;
    nat_t fDim = f->d->fn.fDim;
    int i;
    double s2 = GAUSS_SIGMA*GAUSS_SIGMA;
    double r2 = 0.0;
    for (i = 0; i < pDim; i++) 
      { double dxi = p[i] - 0.5; r2 += dxi*dxi; }
    double e2 = r2/s2;
    double dg = -(p[0]-0.5)*exp(-e2/2)/s2;
    for (i = 0; i < fDim; i++) { fp[i] = dg; }
  }
  
void DgDxGrad(T *f, double *p, double *dfp)
  { affirm(FALSE, "not implemented yet"); } 

void DgDxHess(T *f, double *p, double *ddfp)
  { affirm(FALSE, "not implemented yet"); }

/*** mhat(p) **********************************/

/* The Laplacian of {gaus}, scaled to curvature -3 at center. */

static char *MHatDesc = "(Lapl(gaus))(p)";
static char *MHatType = "SOF.Proc.mhat.";
static SOProcFunction_Methods *MHatMths = NULL;

void MHatEval(T *f, double *p, double *fp);
void MHatGrad(T *f, double *p, double *dfp);
void MHatHess(T *f, double *p, double *ddfp);

void MHatEval(T *f, double *p, double *fp)
  { nat_t pDim = f->d->fn.pDim;
    nat_t fDim = f->d->fn.fDim;
    int i;
    double s2 = GAUSS_SIGMA*GAUSS_SIGMA;
    double r2 = 0.0;
    for (i = 0; i < pDim; i++)
      { double dxi = p[i] - 0.5; r2 += dxi*dxi; }
    double e2 = r2/s2;
    double mh = exp(-e2)*(r2 - s2*((double)pDim));
    for (i = 0; i < fDim; i++) { fp[i] = mh; }
  }
  
void MHatGrad(T *f, double *p, double *dfp)
  { affirm(FALSE, "not implemented yet"); } 

void MHatHess(T *f, double *p, double *ddfp)
  { affirm(FALSE, "not implemented yet"); }

/*** hump ********************************/

/* A multi-quadratic hump centered at {(1/2,1/2,...)} */

static char *HumpDesc = "2*x*(1-x)*y*(1-y)*···";
static char *HumpType = "SOF.Proc.hump.";
static SOProcFunction_Methods *HumpMths = NULL;

void HumpEval(T *f, double *p, double *fp);
void HumpGrad(T *f, double *p, double *dfp);
void HumpHess(T *f, double *p, double *ddfp);

void HumpEval(T *f, double *p, double *fp)
  { nat_t pDim = f->d->fn.pDim;
    nat_t fDim = f->d->fn.fDim;
    int i;
    double h = 1.0;
    for (i = 0; i < pDim; i++) { double xi = p[i]; h *= 2.0*xi*(1-xi); }
    for (i = 0; i < fDim; i++) { fp[i] = h; }
  }
  
void HumpGrad(T *f, double *p, double *dfp)
  { affirm(FALSE, "not implemented yet"); } 

void HumpHess(T *f, double *p, double *ddfp)
  { affirm(FALSE, "not implemented yet"); }

// /*** sqry(p) **********************************/
// 
// static char *SqryDesc = "y^2";
// static char *SqryType = "SOF.Proc.sqry.";
// static SOProcFunction_Methods *SqryMths = NULL;
// 
// void SqryEval(T *f, double *p, double *fp);
// void SqryGrad(T *f, double *p, double *dfp);
// void SqryHess(T *f, double *p, double *ddfp);
// 
// void SqryEval(T *f, double *p, double *fp)
//   { double s = f->d->scale, y = p->c[1];
//     return s*y*y;
//   }
//   
// void SqryGrad(T *f, double *p, double *dfp)
//   { double s = f->d->scale, y = p->c[1];
//     return (r3_t){{ 0.0, 2.0*s*y, 0.0}};
//   }
// 
// void SqryHess(T *f, double *p, double *ddfp)
//   { double s = f->d->scale;
//     return (R3Hessian){{0.0, 0.0, s*2.0, 0.0, 0.0, 0.0}};
//   }  
// 
// /*** cubx(p) **********************************/
// 
// static char *CubxDesc = "x^3";
// static char *CubxType = "SOF.Proc.cubx.";
// static SOProcFunction_Methods *CubxMths = NULL;
// 
// void CubxEval(T *f, double *p, double *fp);
// void CubxGrad(T *f, double *p, double *dfp);
// void CubxHess(T *f, double *p, double *ddfp);
// 
// void CubxEval(T *f, double *p, double *fp)
//   { double s = f->d->scale, x = p->c[0];
//     return s*x*x*x;
//   }
//   
// void CubxGrad(T *f, double *p, double *dfp)
//   { double s = f->d->scale, x = p->c[0];
//     return (r3_t){{3.0*s*x*x, 0.0, 0.0}};
//   }
// 
// void CubxHess(T *f, double *p, double *ddfp)
//   { double s = f->d->scale, x = p->c[0];
//     return (R3Hessian){{s*6.0*x, 0.0, 0.0, 0.0, 0.0, 0.0}};
//   }
// 
// /*** qtrx(p) **********************************/
// 
// static char *QrtxDesc = "x^4";
// static char *QrtxType = "SOF.Proc.qrtx.";
// static SOProcFunction_Methods *QrtxMths = NULL;
// 
// void QrtxEval(T *f, double *p, double *fp);
// void QrtxGrad(T *f, double *p, double *dfp);
// void QrtxHess(T *f, double *p, double *ddfp);
// 
// void QrtxEval(T *f, double *p, double *fp)
//   { double s = f->d->scale, x = p->c[0];
//     double x2 = x*x;
//     return s*x2*x2;
//   }
//   
// void QrtxGrad(T *f, double *p, double *dfp)
//   { double s = f->d->scale, x = p->c[0];
//     return (r3_t){{4.0*s*x*x*x, 0.0, 0.0}};
//   }
// 
// void QrtxHess(T *f, double *p, double *ddfp)
//   { double s = f->d->scale, x = p->c[0];
//     return (R3Hessian){{s*12.0*x*x, 0.0, 0.0, 0.0, 0.0, 0.0}};
//   }
//   
// /*** quix(p) **********************************/
// 
// static char *QuixDesc = "x^5";
// static char *QuixType = "SOF.Proc.quix.";
// static SOProcFunction_Methods *QuixMths = NULL;
// 
// void QuixEval(T *f, double *p, double *fp);
// void QuixGrad(T *f, double *p, double *dfp);
// void QuixHess(T *f, double *p, double *ddfp);
// 
// void QuixEval(T *f, double *p, double *fp)
//   { double s = f->d->scale, x = p->c[0];
//     double x2 = x*x, x3 = x*x2;
//     return s*x2*x3;
//   }
//   
// void QuixGrad(T *f, double *p, double *dfp)
//   { double s = f->d->scale, x = p->c[0];
//     double x2 = x*x;
//     return (r3_t){{5.0*s*x2*x2, 0.0, 0.0}};
//   }
// 
// void QuixHess(T *f, double *p, double *ddfp)
//   { double s = f->d->scale, x = p->c[0];
//     double x2 = x*x, x3 = x*x2;
//     return (R3Hessian){{s*20.0*x3, 0.0, 0.0, 0.0, 0.0, 0.0}};
//   }
// 
// /*** sexx(p) **********************************/
// 
// static char *SexxDesc = "x^6";
// static char *SexxType = "SOF.Proc.sexx.";
// static SOProcFunction_Methods *SexxMths = NULL;
// 
// void SexxEval(T *f, double *p, double *fp);
// void SexxGrad(T *f, double *p, double *dfp);
// void SexxHess(T *f, double *p, double *ddfp);
// 
// void SexxEval(T *f, double *p, double *fp)
//   { double s = f->d->scale, x = p->c[0];
//     double x2 = x*x, x3 = x*x2;
//     return s*x3*x3;
//   }
//   
// void SexxGrad(T *f, double *p, double *dfp)
//   { double s = f->d->scale, x = p->c[0];
//     double x2 = x*x, x3 = x*x2;
//     return (r3_t){{6.0*s*x2*x3, 0.0, 0.0}};
//   }
// 
// void SexxHess(T *f, double *p, double *ddfp)
//   { double s = f->d->scale, x = p->c[0];
//     double x2 = x*x, x4 = x2*x2;
//     return (R3Hessian){{s*30.0*x4, 0.0, 0.0, 0.0, 0.0, 0.0}};
//   }
// 
// /*** sepx(p) **********************************/
// 
// static char *SepxDesc = "x^7";
// static char *SepxType = "SOF.Proc.sepx.";
// static SOProcFunction_Methods *SepxMths = NULL;
// 
// void SepxEval(T *f, double *p, double *fp);
// void SepxGrad(T *f, double *p, double *dfp);
// void SepxHess(T *f, double *p, double *ddfp);
// 
// void SepxEval(T *f, double *p, double *fp)
//   { double s = f->d->scale, x = p->c[0];
//     double x2 = x*x, x3 = x*x2, x4 = x2*x2;
//     return s*x3*x4;
//   }
//   
// void SepxGrad(T *f, double *p, double *dfp)
//   { double s = f->d->scale, x = p->c[0];
//     double  x2 = x*x, x3 = x*x2;
//     return (r3_t){{7.0*s*x3*x3, 0.0, 0.0}};
//   }
// 
// void SepxHess(T *f, double *p, double *ddfp)
//   { double s = f->d->scale, x = p->c[0];
//     double x2 = x*x, x3 = x*x2, x5 = x3*x2;
//     return (R3Hessian){{s*42.0*x5, 0.0, 0.0, 0.0, 0.0, 0.0}};
//   }
// 
// /*** octx(p) **********************************/
// 
// static char *OctxDesc = "x^8";
// static char *OctxType = "SOF.Proc.octx.";
// static SOProcFunction_Methods *OctxMths = NULL;
// 
// void OctxEval(T *f, double *p, double *fp);
// void OctxGrad(T *f, double *p, double *dfp);
// void OctxHess(T *f, double *p, double *ddfp);
// 
// void OctxEval(T *f, double *p, double *fp)
//   { double s = f->d->scale, x = p->c[0];
//     double x2 = x*x, x4 = x2*x2;
//     return s*x4*x4;
//   }
//   
// void OctxGrad(T *f, double *p, double *dfp)
//   { double s = f->d->scale, x = p->c[0];
//     double x2 = x*x, x3 = x*x2, x4 = x2*x2;
//     return (r3_t){{8.0*s*x3*x4, 0.0, 0.0}};
//   }
// 
// void OctxHess(T *f, double *p, double *ddfp)
//   { double s = f->d->scale, x = p->c[0];
//     double x2 = x*x, x3 = x*x2, x6 = x3*x3;
//     return (R3Hessian){{s*56.0*x6, 0.0, 0.0, 0.0, 0.0, 0.0}};
//   }
//   
// /*** band(p) **********************************/
// 
// static char *BandDesc = "1-x^2";
// static char *BandType = "SOF.Proc.band.";
// static SOProcFunction_Methods *BandMths = NULL;
// 
// void BandEval(T *f, double *p, double *fp);
// void BandGrad(T *f, double *p, double *dfp);
// void BandHess(T *f, double *p, double *ddfp);
// 
// void BandEval(T *f, double *p, double *fp)
//   { double s = f->d->scale, x = p->c[0];
//     return s*(1.0 - x*x);
//   }
//   
// void BandGrad(T *f, double *p, double *dfp)
//   { double s = f->d->scale, x = p->c[0];
//     return (r3_t){{-2.0*s*x, 0.0, 0.0}};
//   }
// 
// void BandHess(T *f, double *p, double *ddfp)
//   { double s = f->d->scale;
//     return (R3Hessian){{-s*2.0, 0.0, 0.0, 0.0, 0.0, 0.0}};
//   }
// 
// /*** expx(p) **********************************/
// 
// static char *ExpxDesc = "exp(x)";
// static char *ExpxType = "SOF.Proc.expx.";
// static SOProcFunction_Methods *ExpxMths = NULL;
// 
// void ExpxEval(T *f, double *p, double *fp);
// void ExpxGrad(T *f, double *p, double *dfp);
// void ExpxHess(T *f, double *p, double *ddfp);
// 
// void ExpxEval(T *f, double *p, double *fp)
//   { double s = f->d->scale, x = p->c[0];
//     return s*exp(x);
//   }
//   
// void ExpxGrad(T *f, double *p, double *dfp)
//   { double s = f->d->scale, x = p->c[0];
//     return (r3_t){{s*exp(x), 0.0, 0.0}};
//   }
// 
// void ExpxHess(T *f, double *p, double *ddfp)
//   { double s = f->d->scale, x = p->c[0];
//     return (R3Hessian){{s*exp(x), 0.0, 0.0, 0.0, 0.0, 0.0}};
//   }
//   
/*** sinx(p) **********************************/

/* A single cycle of a sinusoidal wave in the {X} direction. */

static char *SinxDesc = "sin(2*PI*x)";
static char *SinxType = "SOF.Proc.sinx.";
static SOProcFunction_Methods *SinxMths = NULL;

void SinxEval(T *f, double *p, double *fp);
void SinxGrad(T *f, double *p, double *dfp);
void SinxHess(T *f, double *p, double *ddfp);

void SinxEval(T *f, double *p, double *fp)
  { nat_t fDim = f->d->fn.fDim;
    int i;
    for (i = 0; i < fDim; i++){ fp[i] = sin(TWOPI*p[0]); }
  }
  
void SinxGrad(T *f, double *p, double *dfp)
  { affirm(FALSE, "not implemented yet"); }

void SinxHess(T *f, double *p, double *ddfp)
  { affirm(FALSE, "not implemented yet"); }

/*** sinx(p)*siny(p) **************************/

/* A 2×2 array of sinusoidal humps with alternating signs. */

static char *SxSyDesc = "sin(2*PI*x)*sin(2*PI*y)";
static char *SxSyType = "SOF.Proc.sxsy.";
static SOProcFunction_Methods *SxSyMths = NULL;

void SxSyEval(T *f, double *p, double *fp);
void SxSyGrad(T *f, double *p, double *dfp);
void SxSyHess(T *f, double *p, double *ddfp);

void SxSyEval(T *f, double *p, double *fp)
  { nat_t pDim = f->d->fn.pDim;
    nat_t fDim = f->d->fn.fDim;
    int i;
    affirm(pDim >= 2, "domain dimension too small");
    double x = p[0]; 
    double y = p[1];
    double v = sin(TWOPI*x)*sin(TWOPI*y);
    for (i = 0; i < fDim; i++){ fp[i] = v; }
  }
  
void SxSyGrad(T *f, double *p, double *dfp)
  { affirm(FALSE, "not implemented yet"); }

void SxSyHess(T *f, double *p, double *ddfp)
  { affirm(FALSE, "not implemented yet"); }

/*** sxcy(p) **************************/

/* Like {sxsy}, but shifted by a quarter-wave in {Y}. */

static char *SxCyDesc = "sin(2*PI*x)*cos(2*PI*y)";
static char *SxCyType = "SOF.Proc.sxcy.";
static SOProcFunction_Methods *SxCyMths = NULL;

void SxCyEval(T *f, double *p, double *fp);
void SxCyGrad(T *f, double *p, double *dfp);
void SxCyHess(T *f, double *p, double *ddfp);

void SxCyEval(T *f, double *p, double *fp)
  { nat_t pDim = f->d->fn.pDim;
    nat_t fDim = f->d->fn.fDim;
    affirm(pDim >= 2, "domain dimension too small");
    int i;
    double x = p[0]; 
    double y = p[1];
    double v = sin(TWOPI*x)*cos(TWOPI*y);
    for (i = 0; i < fDim; i++){ fp[i] = v; }
  }
  
void SxCyGrad(T *f, double *p, double *dfp)
  { affirm(FALSE, "not implemented yet"); }

void SxCyHess(T *f, double *p, double *ddfp)
  { affirm(FALSE, "not implemented yet"); }

// /*** cosx(p) **********************************/
// 
// static char *CosxDesc = "cos(x)";
// static char *CosxType = "SOF.Proc.cosx.";
// static SOProcFunction_Methods *CosxMths = NULL;
// 
// void CosxEval(T *f, double *p, double *fp);
// void CosxGrad(T *f, double *p, double *dfp);
// void CosxHess(T *f, double *p, double *ddfp);
// 
// void CosxEval(T *f, double *p, double *fp)
//   { double s = f->d->scale, x = p->c[0];
//     return s*cos(x);
//   }
//   
// void CosxGrad(T *f, double *p, double *dfp)
//   { double s = f->d->scale, x = p->c[0];
//     return (r3_t){{-s*sin(x), 0.0, 0.0}};
//   }
// 
// void CosxHess(T *f, double *p, double *ddfp)
//   { double s = f->d->scale, x = p->c[0];
//     return (R3Hessian){{-s*cos(x), 0.0, 0.0, 0.0, 0.0, 0.0}};
//   }  
// 
/*************************************************/ 

#define SOProcFunction_FileFormat "2003-04-28"

/* OVERRIDES FOR PARENT CLASS METHODS */
  
SOProcFunction *SOProcFunction_CopyMth(T *f)
  { SOProcFunction *fnew = SOProcFunction_New();
    fnew->d = SOProcFunction_Data_New();
    *(fnew->d) = *(f->d);
    fnew->type = f->type;
    fnew->m = f->m;
    return fnew;
  }

/* CLASS-SPECIFIC METHODS */
  
void SOProcFunction_WriteMth(T *f, FILE *wr)
  { filefmt_write_header(wr, "SOProcFunction", SOProcFunction_FileFormat);
    fprintf(wr, "type = %s\n",  f->type);
    fprintf(wr, "descr = %s\n",  f->d->descr);
    filefmt_write_footer(wr, "SOProcFunction");
    fflush(wr);
  }

/* OTHER PROCS */
  
SOProcFunction_Methods *SOProcFunction_Methods_New(void)
  { void *v = malloc(sizeof(SOProcFunction_Methods));
    return (SOProcFunction_Methods *)notnull(v, "out of mem for SOProcFunction_Methods");
  }

SOProcFunction_Data *SOProcFunction_Data_New(void)
  { void *v = malloc(sizeof(SOProcFunction_Data));
    return (SOProcFunction_Data *)notnull(v, "out of mem for SOProcFunction_Data");
  }

SOProcFunction *SOProcFunction_New(void)
  { void *v = malloc(sizeof(SOProcFunction));
    return (SOProcFunction *)notnull(v, "no mem for SOProcFunction");
  }


SOProcFunction *SOProcFunction_Make
  ( char *type,
    char *descr,
    SOProcFunction_Methods **mp,
    void (*eval)(T *f, double *p, double *fp),
    void (*grad)(T *f, double *p, double *dfp),
    void (*hess)(T *f, double *p, double *ddfp)
  )
{
  SOProcFunction *f = SOProcFunction_New();
  f->type = type;
  if ((*mp) == NULL)
    { SOProcFunction_Methods *m = SOProcFunction_Methods_New();
      /* Superclass methods: */
      m->fn.eval = (EvalMth *)eval;
      m->fn.grad = (GradMth *)grad;
      m->fn.hess = (HessMth *)hess;
      m->fn.write = (WriteMth *)&SOFunction_WriteMth; /* NOT SOProc! */
      m->fn.copy = (CopyMth *)&SOProcFunction_CopyMth;
      /* Class-specific methods */
      m->write = (WriteMth *)&SOProcFunction_WriteMth;
      (*mp) = m;
    }
  f->m = (*mp);
  f->d = SOProcFunction_Data_New();
  f->d->fn.pDim = 0; /* For now. */
  f->d->fn.fDim = 0; /* For now. */
  f->d->descr = descr;
  f->type = type;
  return f;
}

SOProcFunction *SOProcFunction_FromName(char *name, nat_t pDim, nat_t fDim)
  { SOProcFunction *f = SOProcFunction_New();
    #define MkFn SOProcFunction_Make
    if (strcmp(name, "unit") == 0)
      { f = MkFn(UnitType, UnitDesc, &UnitMths, UnitEval, UnitGrad, UnitHess); }
    else if (strcmp(name, "linx") == 0)
      { f = MkFn(LinxType, LinxDesc, &LinxMths, LinxEval, LinxGrad, LinxHess); }
    else if (strcmp(name, "liny") == 0)
      { f = MkFn(LinyType, LinyDesc, &LinyMths, LinyEval, LinyGrad, LinyHess); }
    else if (strcmp(name, "lind") == 0)
      { f = MkFn(LindType, LindDesc, &LindMths, LindEval, LindGrad, LindHess); }
    else if (strcmp(name, "cbsm") == 0)
      { f = MkFn(CbsmType, CbsmDesc, &CbsmMths, CbsmEval, CbsmGrad, CbsmHess); }
    else if (strcmp(name, "linz") == 0)
      { f = MkFn(LinzType, LinzDesc, &LinzMths, LinzEval, LinzGrad, LinzHess); }
    else if (strcmp(name, "sqrx") == 0)
      { f = MkFn(SqrxType, SqrxDesc, &SqrxMths, SqrxEval, SqrxGrad, SqrxHess); }
    else if (strcmp(name, "para") == 0)
      { f = MkFn(ParaType, ParaDesc, &ParaMths, ParaEval, ParaGrad, ParaHess); }
    else if (strcmp(name, "gaus") == 0)
      { f = MkFn(GausType, GausDesc, &GausMths, GausEval, GausGrad, GausHess); }
    else if (strcmp(name, "dgdx") == 0)
      { f = MkFn(DgDxType, DgDxDesc, &DgDxMths, DgDxEval, DgDxGrad, DgDxHess); }
    else if (strcmp(name, "mhat") == 0)
      { f = MkFn(MHatType, MHatDesc, &MHatMths, MHatEval, MHatGrad, MHatHess); }
    else if (strcmp(name, "hump") == 0)
      { f = MkFn(HumpType, HumpDesc, &HumpMths, HumpEval, HumpGrad, HumpHess); }
    else if (strcmp(name, "sinx") == 0)
      { f = MkFn(SinxType, SinxDesc, &SinxMths, SinxEval, SinxGrad, SinxHess); }
    else if (strcmp(name, "sxcy") == 0)
      { f = MkFn(SxCyType, SxCyDesc, &SxCyMths, SxCyEval, SxCyGrad, SxCyHess); }
    else if (strcmp(name, "sxsy") == 0)
      { f = MkFn(SxSyType, SxSyDesc, &SxSyMths, SxSyEval, SxSyGrad, SxSyHess); }
    else 
      { fprintf (stderr, "bad SOProcFunction function name = %s\n", name);
        affirm(FALSE, "aborted");
      };
    f->d->fn.pDim = pDim;
    f->d->fn.fDim = fDim;
    return f;
  }

SOProcFunction *SOProcFunction_Read(FILE *rd, nat_t pDim, nat_t fDim)
  { char *t, *descr;
    SOProcFunction *f;
    filefmt_read_header(rd, "SOProcFunction", SOProcFunction_FileFormat);
    t = nget_string(rd, "type"); fget_eol(rd);
    if (! isprefix("SOF.Proc.", t))
      { fprintf (stderr, "bad SOProcFunction type = %7s\n", t);
        affirm(FALSE, "expected \"SOF.Proc.\"");
      }
    /* Remove final {.} from subtype, and look it up: */
    { char *funcName = txtcat(t + 9, "#");
      affirm(strlen(funcName) > 0, "missing function name");
      affirm(funcName[strlen(funcName)-2] == '.', "missing \".\" in func type");
      funcName[strlen(funcName)-2] = '\000';
      f = SOProcFunction_FromName(funcName, pDim, fDim);
      free(funcName);
    }
    descr = nget_string(rd, "descr"); fget_eol(rd);
    if (strcmp(f->d->descr, descr) != 0)
      { fprintf (stderr, "description mismatch: \"%s\", \"%s\"\n", 
          f->d->descr, descr);
        affirm(FALSE, "aborted");
      }
    filefmt_read_footer(rd, "SOProcFunction");
    free(t);
    free(descr);
    return f;
  }

