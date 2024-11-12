

#define FLinxRHS "Flinx"
#define FlinxCoeff 0.5
#define FlinxSol "linx"   /* {u = x} */ 


#define FsqrxDescr "6.5*x^2 - 2"
#define FsqrxCoeff (0.5)
#define FsqrxSol "sqrx"   /* {u = x^2} */ 
double Fsqrx(double u, S2Point *p);
double Fsqrx(double u, S2Point *p)
  { double x = p->c[0];
    return (6.5*x*x - 2.0);
  }
static FD FsqrxData;

#define FcubxDescr "12.5*x^3 - 6*x"
#define FcubxCoeff (0.5)
#define FcubxSol "cubx"   /* {u = x^3} */ 
double Fcubx(double u, S2Point *p);
double Fcubx(double u, S2Point *p)
  { double x = p->c[0];
    return 12.5*x*x*x - 6.0*x;
  }
static FD FcubxData;

#define FqtrxDescr "20.5*x^4 - 12*x^2"
#define FqtrxCoeff (0.5)
#define FqtrxSol "qtrx"   /* {u = x^4} */ 
double Fqtrx(double u, S2Point *p);
double Fqtrx(double u, S2Point *p)
  { double x = p->c[0];
    double x2 = x*x, x4 = x2*x2;
    return 20.5*x4 - 12.0*x2;
  }
static FD FqtrxData;

#define FsepxDescr "x^7"
#define FsepxCoeff (0.5)
#define FsepxSol "sepx"  /* {u = x^7 (???)} */  
double Fsepx(double u, S2Point *p);
double Fsepx(double u, S2Point *p)
  { double x = p->c[0];
    double x2 = x*x, x4 = x2*x2;
    return x4*x2*x;
  }
static FD FsepxData;

#define FoctxDescr "x^8"
#define FoctxCoeff (0.5)
#define FoctxSol "octx"   /* {u = x^8} */ 
double Foctx(double u, S2Point *p);
double Foctx(double u, S2Point *p)
  { double x = p->c[0];
    double x2 = x*x, x4 = x2*x2;
    return x4*x4;
  }
static FD FoctxData;

#define FbandDescr "2.5 - 6.5*x^2"
#define FbandCoeff (0.5)
#define FbandSol "band"   /* {u = 1-x^2} */ 
double Fband(double u, S2Point *p);
double Fband(double u, S2Point *p)
  { double x = p->c[0];
    return 2.5 - 6.5*x*x;
  }
static FD FbandData;

#define FlocaDescr "x*(186*x^4 -211*x^2 + 48)"
#define FlocaCoeff (0.5)
#define FlocaSol "loca"   /* {u = x*(2*x^2-1)*(3*x^2 - 2)} */ 
double Floca(double u, S2Point *p);
double Floca(double u, S2Point *p)
  { double x = p->c[0];
    double x2 = x*x, x4 = x2*x2;
    return x*(-211.0*x2 + 48.0 + 186.0*x4);
  }
static FD FlocaData;

#define FsinxDescr "sin(x)"
#define FsinxCoeff (0.5)
#define FsinxSol "sinx"   /* {u = sin(x)} */ 
double Fsinx(double u, S2Point *p);
double Fsinx(double u, S2Point *p)
  { double x = p->c[0];
    return sin(x);
  }
static FD FsinxData;

#define FcosxDescr "cos(x)"
#define FcosxCoeff (0.5)
#define FcosxSol "cosx"   /* {u = cos(x)} */ 
double Fcosx(double u, S2Point *p);
double Fcosx(double u, S2Point *p)
  { double x = p->c[0];
    return cos(x);
  }
static FD FcosxData;

#define FexpzDescr "exp(z)*(z^2 + 2*z - 0.5)"
#define FexpzCoeff (0.5)
#define FexpzSol "expz"    /* {u = exp(z)} */ 
double Fexpz(double u, S2Point *p);
double Fexpz(double u, S2Point *p)
  { double z = p->c[2];
    return exp(z)*(z*z + 2.0*z - 0.5);
  }
static FD FexpzData;

#define FexpxDescr "exp(x)*(x^2 + 2*x - 0.5)"
#define FexpxCoeff (0.5)
#define FexpxSol "expx"   /* {u = exp(x)} */ 
double Fexpx(double u, S2Point *p);
double Fexpx(double u, S2Point *p)
  { double x = p->c[0];
    return exp(x)*(x*x + 2.0*x - 0.5);
  }
static FD FexpxData;

#define FmexpDescr "0.5*(exp(x) + u)*(x^2 + 2*x - 0.5)"
#define FmexpCoeff (0.5)
#define FmexpSol "expx"   /* {u = exp(x)} */ 
double Fmexp(double u, S2Point *p);
double Fmexp(double u, S2Point *p)
  { double x = p->c[0];
    return 0.5*(exp(x) + u)*(x*x + 2.0*x - 0.5);
  }
static FD FmexpData;

#define FdexpDescr "(0.75*exp(x) + 0.25*u)*(x^2 + 2*x - 0.5)"
#define FdexpCoeff (0.5)
#define FdexpSol "expx"   /* {u = exp(x)} */ 
double Fdexp(double u, S2Point *p);
double Fdexp(double u, S2Point *p)
  { double x = p->c[0];
    return (0.75*exp(x) + 0.25*u)*(x*x + 2.0*x - 0.5);
  }
static FD FdexpData;

#define FuexpDescr "(0.25*exp(x) + 0.75*u)*(x^2 + 2*x - 0.5)"
#define FuexpCoeff (0.5)
#define FuexpSol "expx"   /* {u = exp(x)} */ 
double Fuexp(double u, S2Point *p);
double Fuexp(double u, S2Point *p)
  { double x = p->c[0];
    return (0.25*exp(x) + 0.75*u)*(x*x + 2.0*x - 0.5);
  }
static FD FuexpData;

#define FoexpDescr "(0.25*exp(x) + 0.75*u)*(x^2 + 2*x + 7.0)"
#define FoexpCoeff (0.5)
#define FoexpSol "expx"   /* {u = exp(x)} */ 
double Foexp(double u, S2Point *p);
double Foexp(double u, S2Point *p)
  { double x = p->c[0];
    return (0.25*exp(x) + 0.75*u)*(x*x + 2.0*x + 7.0);
  }
static FD FoexpData;

#define FtexpDescr "(0.15*exp(x) + 0.85*u)*(x^2 + 2*x - 0.5)"
#define FtexpCoeff (0.5)
#define FtexpSol "expx"   /* {u = exp(x)} */
double Ftexp(double u, S2Point *p);
double Ftexp(double u, S2Point *p)
  { double x = p->c[0];
    return (0.15*exp(x) + 0.85*u)*(x*x + 2.0*x - 0.5);
  }
static FD FtexpData;

