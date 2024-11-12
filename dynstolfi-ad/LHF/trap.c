/*
* henon.c
* henon attractor
* Luiz Henrique de Figueiredo
* 20 May 2004 16:04:11
*/

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

static void attractor(real x, real y)
{
 int n;
 for (n=0; n<3000; n++)
  f(&x,&y);
 for (n=0; n<3000; n++)
 {
  f(&x,&y);
  gpplot(x,y);
 }
}

static void sample(real px, real py, real qx, real qy, int n, int a, int b)
{
 int i;
 gpcolor(a);
 gpline(px,py,qx,qy);
 gpcolor(b);
 for (i=0; i<=n; i++)
 {
  real t=((real)i)/((real)n);
  real x=(1-t)*px+t*qx;
  real y=(1-t)*py+t*qy;
  f(&x,&y);
  gpplot(x,y);
 }
}

static void trap(real ax, real ay, real bx, real by, real cx, real cy, real dx, real dy, int n)
{
 int i;
 for (i=0; i<=n; i++)
 {
  real u=((real)i)/((real)n);
  real px=(1-u)*ax+u*bx;
  real py=(1-u)*ay+u*by;
  real qx=(1-u)*dx+u*cx;
  real qy=(1-u)*dy+u*cy;
  sample(px,py,qx,qy,n,0,2);
 }
 sample(ax,ay,bx,by,n,3,13);
 sample(bx,by,cx,cy,n,3,13);
 sample(cx,cy,dx,dy,n,3,13);
 sample(dx,dy,ax,ay,n,3,13);
 return;
 gpwait(-1);
 sample(ax,ay,bx,by,n,11,11);
 sample(bx,by,cx,cy,n,12,12);
 sample(cx,cy,dx,dy,n,13,13);
 sample(dx,dy,ax,ay,n,14,14);
}

int main(void)
{
 gpopen("henon attractor");
 gpwindow(-1.6,1.6,-0.5,0.5);
 gppalette(0,"black");
 gppalette(1,"white");
 gppalette(2,"red");
 gppalette(3,"green");
 gppalette(11,"red");
 gppalette(12,"green");
 gppalette(13,"blue");
 gppalette(14,"yellow");
#if 0
 gpcolor(2);
 attractor(0,0);
 gpclear(1);
#endif
 trap(-1.33,0.42,1.32,0.133,1.245,-0.14,-1.06,-0.5,1<<8);
 gpclear(1);
 trap(-1.33,0.42,1.32,0.133,1.245,-0.14,-1.06,-0.45,1<<8);
 gpclose(1);
 return 0;
}
