/* See SPFuncMap.h. */
/* Last edited on 2005-10-29 07:18:26 by stolfi */

#include <SPBasic.h>
#include <SPFuncMap.h>
#include <affirm.h>
#include <nat.h>
#include <bool.h>
#include <string.h>
#include <math.h>

#define FM FuncMap

/* RIGHT-HAND SIDES FOR HELMHOLTZ EQUATION PROBLEMS: */

#define FlinxDescr "2.5*x"
/* Solution: {u = x} for {c = 0.5}. */
double Flinx(double u, S2Point *p);
double Flinx(double u, S2Point *p)
  { double x = p->c[0];
    return 2.5*x;
  }
static FM FlinxMap;

#define FsqrxDescr "6.5*x^2 - 2"
/* Solution: {u = x^2} for {c = 0.5}. */
double Fsqrx(double u, S2Point *p);
double Fsqrx(double u, S2Point *p)
  { double x = p->c[0];
    return (6.5*x*x - 2.0);
  }
static FM FsqrxMap;

#define FcubxDescr "12.5*x^3 - 6*x"
/* Solution: {u = x^3} for {c = 0.5}. */
double Fcubx(double u, S2Point *p);
double Fcubx(double u, S2Point *p)
  { double x = p->c[0];
    return 12.5*x*x*x - 6.0*x;
  }
static FM FcubxMap;

#define FqrtxDescr "20.5*x^4 - 12*x^2"
/* Solution: {x^4} for {c = 0.5}. */
double Fqrtx(double u, S2Point *p);
double Fqrtx(double u, S2Point *p)
  { double x = p->c[0];
    double x2 = x*x, x4 = x2*x2;
    return 20.5*x4 - 12.0*x2;
  }
static FM FqrtxMap;

#define FsepxDescr "56.5*x^7 - 42*x^5"
/* Solution: {u = x^7} for {c = 0.5}. */
double Fsepx(double u, S2Point *p);
double Fsepx(double u, S2Point *p)
  { double x = p->c[0];
    double x2 = x*x, x4 = x2*x2;
    return x4*x*(56.5*x2 - 42);
  }
static FM FsepxMap;

#define FmsepDescr "56.5*(0.75*x^7 + 0.25*u) - 42*x^2*cbrt(x^2*u)"
/* Solution: {u = x^7} for {c = 0.5} (nonlinear). */
double Fmsep(double u, S2Point *p);
double Fmsep(double u, S2Point *p)
  { double x = p->c[0];
    double x2 = x*x, x4 = x2*x2, x5 = x4*x, x7 = x5*x2;
    return 56.5*(0.75*x7 + 0.25*u) - 42*x2*cbrt(x2*u);
  }
static FM FmsepMap;

#define FoctxDescr "72.5*x^8 - 56*x^6"
/* Solution: {u = x^8} for {c = 0.5}. */
double Foctx(double u, S2Point *p);
double Foctx(double u, S2Point *p)
  { double x = p->c[0];
    double x2 = x*x, x4 = x2*x2;
    return x4*x2*(72.5*x2 - 56);
  }
static FM FoctxMap;

#define FmoctDescr "x^6*(72.5*x^2 - 56) + x^2*(u - x^8)"
/* Solution: {u = x^8} for {c = 0.5} (nonlinear). */
double Fmoct(double u, S2Point *p);
double Fmoct(double u, S2Point *p)
  { double x = p->c[0];
    double x2 = x*x, x4 = x2*x2, x6 = x4*x2, x8 = x4*x4;
    return x6*(72.5*x2 - 56) + x2*(u - x8);
  }
static FM FmoctMap;

#define FbandDescr "2.5 - 6.5*x^2"
/* Solution: {u = 1-x^2} for {c = 0.5}. */
double Fband(double u, S2Point *p);
double Fband(double u, S2Point *p)
  { double x = p->c[0];
    return 2.5 - 6.5*x*x;
  }
static FM FbandMap;

#define FlocaDescr "x*(186*x^4 -211*x^2 + 48)"
/* Solution: ???. */
double Floca(double u, S2Point *p);
double Floca(double u, S2Point *p)
  { double x = p->c[0];
    double x2 = x*x, x4 = x2*x2;
    return x*(-211.0*x2 + 48.0 + 186.0*x4);
  }
static FM FlocaMap;

#define FsinxDescr "(1.5-x^2)*sin(x)+2*x*cos(x)"
/* Solution: {u = sin(x)} for {c = 0.5}. */
double Fsinx(double u, S2Point *p);
double Fsinx(double u, S2Point *p)
  { double x = p->c[0];
    return (1.5-x*x)*sin(x) + 2*x*cos(x);
  }
static FM FsinxMap;

#define FcosxDescr "(1.5-x^2)*cos(x)-2*x*sin(x)"
/* Solution: {u = cos(x)} for {c = 0.5}. */
double Fcosx(double u, S2Point *p);
double Fcosx(double u, S2Point *p)
  { double x = p->c[0];
    return (1.5-x*x)*cos(x) - 2*x*sin(x);
  }
static FM FcosxMap;

#define FmcosDescr "(1.5-x^2)*(0.75*cos(x)+0.25*u)-2*x*sin(x)"
/* Solution: {u = cos(x)} for {c = 0.5}. */
double Fmcos(double u, S2Point *p);
double Fmcos(double u, S2Point *p)
  { double x = p->c[0];
    return (1.5-x*x)*(0.75*cos(x) + 0.25*u) - 2*x*sin(x);
  }
static FM FmcosMap;

#define FexpzDescr "exp(z)*(z^2 + 2*z - 0.5)"
/* Solution: {u = exp(z)} for {c = 0.5}. */
double Fexpz(double u, S2Point *p);
double Fexpz(double u, S2Point *p)
  { double z = p->c[2];
    return exp(z)*(z*z + 2.0*z - 0.5);
  }
static FM FexpzMap;

#define FexpxDescr "exp(x)*(x^2 + 2*x - 0.5)"
/* Solution: {u = exp(x)} for {c = 0.5}. */
double Fexpx(double u, S2Point *p);
double Fexpx(double u, S2Point *p)
  { double x = p->c[0];
    return exp(x)*(x*x + 2.0*x - 0.5);
  }
static FM FexpxMap;

#define FmexpDescr "0.5*(exp(x) + u)*(x^2 + 2*x - 0.5)"
/* Solution: {u = exp(x)} for {c = 0.5} (iterative but linear). */
double Fmexp(double u, S2Point *p);
double Fmexp(double u, S2Point *p)
  { double x = p->c[0];
    return 0.5*(exp(x) + u)*(x*x + 2.0*x - 0.5);
  }
static FM FmexpMap;

#define FdexpDescr "(0.75*exp(x) + 0.25*u)*(x^2 + 2*x - 0.5)"
/* Solution: {u = exp(x)} for {c = 0.5} (iterative but linear). */
double Fdexp(double u, S2Point *p);
double Fdexp(double u, S2Point *p)
  { double x = p->c[0];
    return (0.75*exp(x) + 0.25*u)*(x*x + 2.0*x - 0.5);
  }
static FM FdexpMap;

#define FuexpDescr "(0.25*exp(x) + 0.75*u)*(x^2 + 2*x - 0.5)"
/* Solution: {u = exp(x)} for {c = 0.5} (iterative but linear). */
double Fuexp(double u, S2Point *p);
double Fuexp(double u, S2Point *p)
  { double x = p->c[0];
    return (0.25*exp(x) + 0.75*u)*(x*x + 2.0*x - 0.5);
  }
static FM FuexpMap;

#define FoexpDescr "(0.25*exp(x) + 0.75*u)*(x^2 + 2*x + 7.0)"
/* Solution: {u = exp(x)} for {c = 8.0} (iterative but linear). */
double Foexp(double u, S2Point *p);
double Foexp(double u, S2Point *p)
  { double x = p->c[0];
    return (0.25*exp(x) + 0.75*u)*(x*x + 2.0*x + 7.0);
  }
static FM FoexpMap;

#define FtexpDescr "(0.15*exp(x) + 0.85*u)*(x^2 + 2*x - 0.5)"
/* Solution: {u = exp(x)} for {c = 0.5} (iterative but linear). */
double Ftexp(double u, S2Point *p);
double Ftexp(double u, S2Point *p)
  { double x = p->c[0];
    return (0.15*exp(x) + 0.85*u)*(x*x + 2.0*x - 0.5);
  }
static FM FtexpMap;

/* TABLE INITIALIZATION */

static bool_t SPFuncMap_initialized = FALSE;

void SPFuncMap_Initialize(void);
void SPFuncMap_Initialize(void)
  {
    FlinxMap =  (FM){&Flinx, FALSE, FlinxDescr};
    FsqrxMap =  (FM){&Fsqrx, FALSE, FsqrxDescr};
    FcubxMap =  (FM){&Fcubx, FALSE, FcubxDescr};
    FqrtxMap =  (FM){&Fqrtx, FALSE, FqrtxDescr};
    FsepxMap =  (FM){&Fsepx, FALSE, FsepxDescr};
    FmsepMap =  (FM){&Fmsep, FALSE, FmsepDescr};
    FoctxMap =  (FM){&Foctx, FALSE, FoctxDescr};
    FmoctMap =  (FM){&Fmoct, FALSE, FmoctDescr};
    FbandMap =  (FM){&Fband, FALSE, FbandDescr};
    FlocaMap =  (FM){&Floca, FALSE, FlocaDescr};
    FsinxMap =  (FM){&Fsinx, FALSE, FsinxDescr};
    FcosxMap =  (FM){&Fcosx, FALSE, FcosxDescr};
    FmcosMap =  (FM){&Fmcos, FALSE, FmcosDescr};
    FexpzMap =  (FM){&Fexpz, FALSE, FexpzDescr};
    FexpxMap =  (FM){&Fexpx, FALSE, FexpxDescr};
    FmexpMap =  (FM){&Fmexp, FALSE, FmexpDescr};
    FdexpMap =  (FM){&Fdexp, FALSE, FdexpDescr};
    FuexpMap =  (FM){&Fuexp, FALSE, FuexpDescr};
    FoexpMap =  (FM){&Foexp, FALSE, FoexpDescr};
    FtexpMap =  (FM){&Ftexp, FALSE, FtexpDescr};
    SPFuncMap_initialized = TRUE;
  }

SPFuncMap SPFuncMap_FromName(char *name)
  {
    if (! SPFuncMap_initialized) { SPFuncMap_Initialize(); }
    if (strcmp(name, "Flinx") == 0) { return FlinxMap; }
    if (strcmp(name, "Fsqrx") == 0) { return FsqrxMap; }
    if (strcmp(name, "Fcubx") == 0) { return FcubxMap; }
    if (strcmp(name, "Fqrtx") == 0) { return FqrtxMap; }
    if (strcmp(name, "Fsepx") == 0) { return FsepxMap; }
    if (strcmp(name, "Fmsep") == 0) { return FmsepMap; }
    if (strcmp(name, "Foctx") == 0) { return FoctxMap; }
    if (strcmp(name, "Fmoct") == 0) { return FmoctMap; }
    if (strcmp(name, "Fband") == 0) { return FbandMap; }
    if (strcmp(name, "Floca") == 0) { return FlocaMap; }
    if (strcmp(name, "Fsinx") == 0) { return FsinxMap; }
    if (strcmp(name, "Fcosx") == 0) { return FcosxMap; }
    if (strcmp(name, "Fmcos") == 0) { return FmcosMap; }
    if (strcmp(name, "Fexpz") == 0) { return FexpzMap; }
    if (strcmp(name, "Fexpx") == 0) { return FexpxMap; }
    if (strcmp(name, "Fmexp") == 0) { return FmexpMap; }
    if (strcmp(name, "Fdexp") == 0) { return FdexpMap; }
    if (strcmp(name, "Fuexp") == 0) { return FuexpMap; }
    if (strcmp(name, "Foexp") == 0) { return FoexpMap; }
    if (strcmp(name, "Ftexp") == 0) { return FtexpMap; }
    affirm(FALSE, "unknown FuncMap name");
    return NoFMap;
  }

