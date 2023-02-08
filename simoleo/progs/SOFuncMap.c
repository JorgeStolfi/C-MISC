/* See SOFuncMap.h. */
/* Last edited on 2004-06-19 18:22:26 by stolfi */

#include <SOFuncMap.h>
#include <SOBasic.h>
#include <SOFunction.h>
#include <SOGrid.h>

#include <dg_grid.h>

#include <affirm.h>
#include <nat.h>
#include <string.h>
#include <math.h>

#define FM FuncMap
#define FD SOFuncMap_Data

dg_dim_t uDim;   /* Dimension of the function's value {u = f(p)}. */
dg_dim_t pDim;   /* Dimension of the argument points {p}. */
dg_dim_t vDim;   /* Dimension of the {FuncMap}'s result {v = FMap(f(p),p)}. */

#define FlinxDescr "3*x"
#define FlinxCoeff (0.5)
#define FlinxSol "linx"   /* u = x */ 
void Flinx(double *u, double *p, double *v);
void Flinx(double *u, double *p, double *v)
  { int i; for( i = 0; i < vDim; i++) v[i] = (double)3.0*p[X]; }

#define FsqrxDescr "6.5*x^2 - 2"
#define FsqrxCoeff (0.5)
#define FsqrxSol "sqrx"   /* u = x^2 */ 
void Fsqrx(double *u, double *p, double *v);
void Fsqrx(double *u, double *p, double *v)
  { double x = p[X];
    int i; for( i = 0; i < vDim; i++) v[i] = ((double)6.5*x*x - 2.0);
  }

#define FcubxDescr "12.5*x^3 - 6*x"
#define FcubxCoeff (0.5)
#define FcubxSol "cubx"   /* u = x^3 */ 
void Fcubx(double *u, double *p, double *v);
void Fcubx(double *u, double *p, double *v)
  { double x = p[X];
    int i; for( i = 0; i < vDim; i++) v[i] = ((double)12.5*x*x*x - (double)6.0*x);
  }

#define FqtrxDescr "20.5*x^4 - 12*x^2"
#define FqtrxCoeff (0.5)
#define FqtrxSol "qtrx"   /* u = x^4 */ 
void Fqtrx(double *u, double *p, double *v);
void Fqtrx(double *u, double *p, double *v)
  { double x = p[X];
    double x2 = (double)x*x, x4 = (double)x2*x2;
    int i; for( i = 0; i < vDim; i++) v[i] = ((double)20.5*x4 - (double)12.0*x2);
  }

#define FsepxDescr "x^7"
#define FsepxCoeff (0.5)
#define FsepxSol "sepx"  /* u = x^7 (???) */  
void Fsepx(double *u, double *p, double *v);
void Fsepx(double *u, double *p, double *v)
  { double x = p[X];
    double x2 = (double)x*x, x4 = (double)x2*x2;
    int i; for( i = 0; i < vDim; i++) v[i] = (double)x4*x2*x;
  }

#define FoctxDescr "x^8"
#define FoctxCoeff (0.5)
#define FoctxSol "octx"   /* u = x^8 */ 
void Foctx(double *u, double *p, double *v);
void Foctx(double *u, double *p, double *v)
  { double x = p[X];
    double x2 = (double)x*x, x4 = (double)x2*x2;
    int i; for( i = 0; i < vDim; i++) v[i] = (double)x4*x4;
  }

#define FbandDescr "2.5 - 6.5*x^2"
#define FbandCoeff (0.5)
#define FbandSol "band"   /* u = 1-x^2 */ 
void Fband(double *u, double *p, double *v);
void Fband(double *u, double *p, double *v)
  { double x = p[X];
    int i; for( i = 0; i < vDim; i++) v[i] = 2.5 - (double)6.5*x*x;
  }

#define FlinsumDescr "45 * x^3 + 45 * y^3 + 6 * x + 6 * y"
#define FlinsumCoeff (45.0)
#define FlinsumSol "cubsum"  
void Flinsum(double *u, double *p, double *v);
void Flinsum(double *u, double *p, double *v)
  { double x = p[X]; double y = p[Y];
    double x3 = (double)45.0 * x * x * x; double y3 = (double)45.0 * y * y * y;
    int i; for( i = 0; i < vDim; i++) v[i] = x3 + y3 + (double)6.0 * x + (double)6.0 * y;
  }

#define FsumprodDescr "-2 *( (x*(1-x))+(y*(1-y)) -8((x*(1-x))*(y*(1-y))) )"
#define FsumprodCoeff (80.0)
#define FsumprodSol "sumprod"  
void Fsumprod(double *u, double *p, double *v);
void Fsumprod(double *u, double *p, double *v)
  { double a = (double)p[0]*(1 - p[0]); double b = (double)p[1]*(1 - p[1]);
    int i; for( i = 0; i < vDim; i++) 
      v[i] = (double) -2 * (a + b - (double) 40 * (a * b));
  }

#define FlocaDescr "x*(186*x^4 -211*x^2 + 48)"
#define FlocaCoeff (0.5)
#define FlocaSol "loca"   /* u = x*(2*x^2-1)*(3*x^2 - 2) */ 
void Floca(double *u, double *p, double *v);
void Floca(double *u, double *p, double *v)
  { double x = p[X];
    double x2 = (double)x*x, x4 = (double)x2*x2;
    int i; for( i = 0; i < vDim; i++) 
      v[i] = (double)x*((double)-211.0*x2 + 48.0 + (double)186.0*x4);
  }

#define FsinxDescr "sin(x)"
#define FsinxCoeff (0.5)
#define FsinxSol "sinx"   /* u = sin(x) */ 
void Fsinx(double *u, double *p, double *v);
void Fsinx(double *u, double *p, double *v)
  { double x = p[X];
    int i; for( i = 0; i < vDim; i++) v[i] = sin(x);
  }

#define FcosxDescr "cos(x)"
#define FcosxCoeff (0.5)
#define FcosxSol "cosx"   /* u = cos(x) */ 
void Fcosx(double *u, double *p, double *v);
void Fcosx(double *u, double *p, double *v)
  { double x = p[X];
    int i; for( i = 0; i < vDim; i++) v[i] = cos(x);
  }

#define Fsinx_x_sinyDescr "-8*TWOPI^2*sin(x)*sin(y)+45*sin(x)*sin(y)"
#define Fsinx_x_sinyCoeff (45.0)
#define Fsinx_x_sinySol "sinx_x_siny" 
void Fsinx_x_siny(double *u, double *p, double *v);
void Fsinx_x_siny(double *u, double *p, double *v)
  { double x = p[X]; double y = p[Y];
    double srqPI = (double)TWOPI * TWOPI;
    double k = (double)-8 * srqPI + 45.0;  
    int i; 
    for( i = 0; i < vDim; i++) v[i] = (double)k * sin(TWOPI*x) * sin(TWOPI*y);
  }

#define FexpzDescr "exp(z)*(z^2 + 2*z - 0.5)"
#define FexpzCoeff (0.5)
#define FexpzSol "expz"    /* u = exp(z) */ 
void Fexpz(double *u, double *p, double *v);
void Fexpz(double *u, double *p, double *v)
  { double z = p[Z];
    int i; for( i = 0; i < vDim; i++) v[i] = exp(z)*((double)z*z + (double)2.0*z - 0.5);
  }

#define FexpxDescr "exp(x)*(x^2 + 2*x - 0.5)"
#define FexpxCoeff (0.5)
#define FexpxSol "expx"   /* u = exp(x) */ 
void Fexpx(double *u, double *p, double *v);
void Fexpx(double *u, double *p, double *v)
  { double x = p[X];
    int i; for( i = 0; i < vDim; i++) v[i] = exp(x)*((double)x*x + (double)2.0*x - 0.5);
  }

#define FmexpDescr "0.5*(exp(x) + u)*(x^2 + 2*x - 0.5)"
#define FmexpCoeff (0.5)
#define FmexpSol "expx"   /* u = exp(x) */ 
void Fmexp(double *u, double *p, double *v);
void Fmexp(double *u, double *p, double *v)
  { double x = p[X];
    int i; for( i = 0; i < vDim; i++) 
      v[i] = (double)0.5*(exp(x) + u[X])*((double)x*x + (double)2.0*x - 0.5);
  }

#define FdexpDescr "(0.75*exp(x) + 0.25*u)*(x^2 + 2*x - 0.5)"
#define FdexpCoeff (0.5)
#define FdexpSol "expx"   /* u = exp(x) */ 
void Fdexp(double *u, double *p, double *v);
void Fdexp(double *u, double *p, double *v)
  { double x = p[X];
    int i; for( i = 0; i < vDim; i++)
      v[i] = (double)((double)0.75*exp(x) + (double)0.25*u[X])*((double)x*x + (double)2.0*x - 0.5);
  }

#define FuexpDescr "(0.25*exp(x) + 0.75*u)*(x^2 + 2*x - 0.5)"
#define FuexpCoeff (0.5)
#define FuexpSol "expx"   /* u = exp(x) */ 
void Fuexp(double *u, double *p, double *v);
void Fuexp(double *u, double *p, double *v)
  { double x = p[X];
    int i; for( i = 0; i < vDim; i++)
      v[i] = (double)((double)0.25*exp(x) + (double)0.75*u[X])*((double)x*x + (double)2.0*x - 0.5);
  }

#define FoexpDescr "(0.25*exp(x) + 0.75*u)*(x^2 + 2*x + 7.0)"
#define FoexpCoeff (0.5)
#define FoexpSol "expx"   /* u = exp(x) */ 
void Foexp(double *u, double *p, double *v);
void Foexp(double *u, double *p, double *v)
  { double x = p[X];
    int i; for( i = 0; i < vDim; i++)
    v[i] = (double)((double)0.25*exp(x) + (double)0.75*u[X])*((double)x*x + (double)2.0*x + 7.0);
  }

#define FtexpDescr "(0.15*exp(x) + 0.85*u)*(x^2 + 2*x - 0.5)"
#define FtexpCoeff (0.5)
#define FtexpSol "expx"   /* u = exp(x) */ 
void Ftexp(double *u, double *p, double *v);
void Ftexp(double *u, double *p, double *v)
  { double x = p[X];
    int i; for( i = 0; i < vDim; i++)
    v[i] = (double)((double)0.15*exp(x) + (double)0.85*u[X])*((double)x*x + (double)2.0*x - 0.5);
  }

#define NoData (FD){(FM){NULL, uDim, pDim, vDim, TRUE}, NULL, NULL, 0.0, NULL}

SOFuncMap_Data SOFuncMap_FromName( char *name, dg_dim_t u, dg_dim_t p, dg_dim_t v )
{
  uDim = u; pDim = p; vDim = v;

  if (strcmp(name, "Flinx") == 0) 
    { return (FD){(FM){&Flinx, uDim, pDim, vDim, FALSE}, "Flinx", FlinxDescr, FlinxCoeff, FlinxSol}; }
  if (strcmp(name, "Fsqrx") == 0) 
    { return (FD){(FM){&Fsqrx, uDim, pDim, vDim, FALSE}, "Fsqrx", FsqrxDescr, FsqrxCoeff, FsqrxSol}; }
  if (strcmp(name, "Fcubx") == 0) 
    { return (FD){(FM){&Fcubx, uDim, pDim, vDim, FALSE}, "Fcubx", FcubxDescr, FcubxCoeff, FcubxSol}; }
  if (strcmp(name, "Fqtrx") == 0) 
    { return (FD){(FM){&Fqtrx, uDim, pDim, vDim, FALSE}, "Fqtrx", FqtrxDescr, FqtrxCoeff, FqtrxSol}; }
  if (strcmp(name, "Fsepx") == 0) 
    { return (FD){(FM){&Fsepx, uDim, pDim, vDim, FALSE}, "Fsepx", FsepxDescr, FsepxCoeff, FsepxSol}; }
  if (strcmp(name, "Foctx") == 0) 
    { return (FD){(FM){&Foctx, uDim, pDim, vDim, FALSE}, "Foctx", FoctxDescr, FoctxCoeff, FoctxSol}; }
  if (strcmp(name, "Fband") == 0) 
    { return (FD){(FM){&Fband, uDim, pDim, vDim, FALSE}, "Fband", FbandDescr, FbandCoeff, FbandSol}; }
  if (strcmp(name, "Flinsum") == 0) 
    { return (FD){(FM){&Flinsum, uDim, pDim, vDim, FALSE}, 
      "Flinsum", FlinsumDescr, FlinsumCoeff, FlinsumSol}; }
  if (strcmp(name, "Fsumprod") == 0) 
    { return (FD){(FM){&Fsumprod, uDim, pDim, vDim, FALSE}, 
      "Fsumprod", FsumprodDescr, FsumprodCoeff, FsumprodSol}; }
  if (strcmp(name, "Floca") == 0) 
    { return (FD){(FM){&Floca, uDim, pDim, vDim, FALSE}, "Floca", FlocaDescr, FlocaCoeff, FlocaSol}; }
  if (strcmp(name, "Fsinx") == 0) 
    { return (FD){(FM){&Fsinx, uDim, pDim, vDim, FALSE}, "Fsinx", FsinxDescr, FsinxCoeff, FsinxSol}; }
  if (strcmp(name, "Fcosx") == 0) 
    { return (FD){(FM){&Fcosx, uDim, pDim, vDim, FALSE}, "Fcosx", FcosxDescr, FcosxCoeff, FcosxSol}; }
  if (strcmp(name, "Fsinx_x_siny") == 0) 
    { return (FD){(FM){&Fsinx_x_siny, uDim, pDim, vDim, FALSE}, 
      "Fsinx_x_siny", Fsinx_x_sinyDescr, Fsinx_x_sinyCoeff, Fsinx_x_sinySol}; }
  if (strcmp(name, "Fexpz") == 0) 
    { return (FD){(FM){&Fexpz, uDim, pDim, vDim, FALSE}, "Fexpz", FexpzDescr, FexpzCoeff, FexpzSol}; }
  if (strcmp(name, "Fexpx") == 0) 
    { return (FD){(FM){&Fexpx, uDim, pDim, vDim, FALSE}, "Fexpx", FexpxDescr, FexpxCoeff, FexpxSol}; }
  if (strcmp(name, "Fmexp") == 0) 
    { return (FD){(FM){&Fmexp, uDim, pDim, vDim, FALSE}, "Fmexp", FmexpDescr, FmexpCoeff, FmexpSol}; }
  if (strcmp(name, "Fdexp") == 0) 
    { return (FD){(FM){&Fdexp, uDim, pDim, vDim, FALSE}, "Fdexp", FdexpDescr, FdexpCoeff, FdexpSol}; }
  if (strcmp(name, "Fuexp") == 0) 
    { return (FD){(FM){&Fuexp, uDim, pDim, vDim, FALSE}, "Fuexp", FuexpDescr, FuexpCoeff, FuexpSol}; }
  if (strcmp(name, "Foexp") == 0) 
    { return (FD){(FM){&Foexp, uDim, pDim, vDim, FALSE}, "Foexp", FoexpDescr, FoexpCoeff, FoexpSol}; }
  if (strcmp(name, "Ftexp") == 0) 
    { return (FD){(FM){&Ftexp, uDim, pDim, vDim, FALSE}, "Ftexp", FtexpDescr, FtexpCoeff, FtexpSol}; }
  affirm(FALSE, "unknown FuncMap name");
  return NoData;
}

