/*
* xtrap.c
* find basin of oo
* Luiz Henrique de Figueiredo
* 20 May 2004 16:04:11
*/

#include <math.h>
#include <stdio.h>
#include <gp.h>

#define	a	(1.4)
#define	b	(0.3)

static void f(real* x, real* y)
{
 real X=*x;
 real Y=*y;
 *x=1.0-a*X*X+Y;
 *y=b*X;
}

#undef a
#undef b

static void image(real r)
{
 double t;
 double dt=2*3.14159265/32/r;
 double x0=0;
 double y0=0;
 gpmark(1,".");
 for (t=0; t<=6.3; t+=dt)
 {
  real x=x0+r*cos(t);
  real y=y0+r*sin(t);
  real xx=x;
  real yy=y;
  f(&x,&y);
  f(&x,&y);
  if (fabs(x)>=r || fabs(y)>=r)
  {
   continue;
   gpcolor(3); gpplot(xx,yy);
  }
  else
  {
   gpcolor(4); gpplot(xx,yy);
   gpcolor(xx>0 ? 2 : 5); gpplot(x,y);
  }
 }
 return;
 y0=y0-r;
 f(&x0,&y0);
 gpmark(3,"B");
 gpcolor(5); gpplot(x0,y0);
}

#define L 100

int main(void)
{
 real r;
 gpopen("henon attractor");
 gpwindow(-L,L,-L,L);
 gppalette(0,"black");
 gppalette(1,"white");
 gppalette(2,"grey");
 gppalette(3,"green");
 gppalette(4,"red");
 gppalette(5,"blue");
 for (r=0.1; r<=L; r++)
 {
  if(0)gpclear(0);
  image(r);
  gpwait(50);
 }
 gpclose(1);
 return 0;
}
