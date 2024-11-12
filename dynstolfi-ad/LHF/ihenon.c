/*
* ia.c
* interval computation for henon map
* Luiz Henrique de Figueiredo
* 28 Jun 2004 09:03:56
*/

#include <stdio.h>

#define real double

typedef real Interval[2];
#define HI(x)	(x[0])
#define LO(x)	(x[1])

#if 0
#include <fenv.h>
#include <math.h>
#define fesetround(x) 	printf("fesetround: %d %d\n",x,(fesetround)(x))
#define ROUND_DOWN	fesetround(FE_DOWNWARD)
#define ROUND_UP	fesetround(FE_UPWARD)
#else
#include <fpu_control.h>
fpu_control_t cwup   =	_FPU_IEEE | _FPU_RC_UP;
fpu_control_t cwdown =	_FPU_IEEE | _FPU_RC_DOWN;
#define ROUND_DOWN	_FPU_SETCW(cwdown)
#define ROUND_UP	_FPU_SETCW(cwup)
#endif

static void Iprint(Interval x)
{
 printf("%p\t%.17g\t%.17g\t%s\n",(void*)x,LO(x),HI(x),(LO(x)==HI(x))?"true":"false");
}

static void Inew(Interval x, real lo, real hi)
{
 LO(x)=lo;
 HI(x)=hi;
}

static void Irational(Interval x, real p, real q)
{
 ROUND_DOWN;
 LO(x)=p/q;
#if 0
 HI(x)=nextafter(p/q,p/q+1);
 return;
#endif
 ROUND_UP;
 HI(x)=p/q;
}

void Isqr(Interval a, Interval b)
{
 real bmin=LO(b);
 real bmax=HI(b);
 if (bmin>=0)
 {
  ROUND_DOWN;
  LO(a)=bmin*bmin;
  ROUND_UP;
  HI(a)=bmax*bmax;
 }
 else if (bmax<=0)
 {
  ROUND_DOWN;
  LO(a)=bmax*bmax;
  ROUND_UP;
  HI(a)=bmin*bmin;
 }
 else
 {
  LO(a)=0.0;
  ROUND_UP;
  if (bmax<-bmin)
   HI(a)=bmin*bmin;
  else
   HI(a)=bmax*bmax;
 }
}

static void Ihenon(Interval x, Interval y, Interval a, Interval b)
{
 Interval x2,X;
 Isqr(x2,x);
 ROUND_DOWN;
 LO(X)=1.0-HI(a)*HI(x2)+LO(y);
 LO(y)=LO(b)*LO(x);
 ROUND_UP;
 HI(X)=1.0-LO(a)*LO(x2)+HI(y);
 HI(y)=HI(b)*HI(x);
 LO(x)=LO(X); HI(x)=HI(X);
}

static void doit(real xmin, real xmax, real ymin, real ymax)
{
 Interval a,b,x,y;
 Irational(a,14,10);	Iprint(a);
 Irational(b, 3,10);	Iprint(b);
 Inew(x,xmin,xmax); Iprint(x);
 Inew(y,ymin,ymax); Iprint(y);
 Ihenon(x,y,a,b); Iprint(x); Iprint(y);
 Irational(b, 1,8);	Iprint(b);
}

int main(void)
{
 doit(-1,1,-1,1);
 doit(0,0.1,0,0.1);
 return 0;
}
