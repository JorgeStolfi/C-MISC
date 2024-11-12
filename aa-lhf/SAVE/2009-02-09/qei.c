/*
* qei.c
* enumeration of implicit curves using quadtrees and affine arithmetic
* Luiz Henrique de Figueiredo
* 18 Aug 93
*/

#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include "gp.h"
#include "arith.h"

void enumerate(real xmin, real xmax, real ymin, real ymax);
void explore(real xmin, real xmax, real ymin, real ymax);
void plot(real xmin, real xmax, real ymin, real ymax);
void show(int c, int t, real xmin, real xmax, real ymin, real ymax);
int haszero(real xmin, real xmax, real ymin, real ymax);
int toosmall(real xmin, real xmax, real ymin, real ymax);
void image(real xmin, real xmax, real ymin, real ymax, real* fmin, real* fmax);

#undef show
#define show(c,t,xmin,xmax,ymin,ymax)

real tol;

int main(int argc, char* argv[])
{
 if (argc>1)
  sscanf(argv[1],"%lg",&tol);
 else
  tol=0.01;
 enumerate(-0.5,1.5,-0.5,1.5);
 return 0;
}

void enumerate(real xmin, real xmax, real ymin, real ymax)
{
 gpopen("quadtree enumeration of implicit curves");
 gppalette(2,"red");
 gppalette(3,"grey");
 gpwindow(xmin,xmax,ymin,ymax);
 if (tol<0)				/* negative tol means in pixels */
 {
  real x0,x1,y;
  x0=0;    y=0; gpunview(&x0,&y);
  x1=-tol; y=0; gpunview(&x1,&y);
  tol=fabs(x1-x0);
 }
 explore(xmin,xmax,ymin,ymax);
 gpclose(1);
}

void explore(real xmin, real xmax, real ymin, real ymax)
{
 show(1,'p',xmin,xmax,ymin,ymax);
 if (haszero(xmin,xmax,ymin,ymax))
 {
  if (toosmall(xmin,xmax,ymin,ymax))
   plot(xmin,xmax,ymin,ymax);
  else
  {
   real xmid=(xmin+xmax)/2.0;
   real ymid=(ymin+ymax)/2.0;
   explore(xmin,xmid,ymid,ymax);
   explore(xmid,xmax,ymid,ymax);
   explore(xmin,xmid,ymin,ymid);
   explore(xmid,xmax,ymin,ymid);
  }
 }
 show(3,'p',xmin,xmax,ymin,ymax);
}

int haszero(real xmin, real xmax, real ymin, real ymax)
{
 real fmin,fmax;
 image(xmin,xmax,ymin,ymax,&fmin,&fmax);
 return (fmin<=0) && (fmax>=0);
}

int toosmall(real xmin, real xmax, real ymin, real ymax)
{
#if 0
 real fmin,fmax,tol=0.01;
 image(xmin,xmax,ymin,ymax,&fmin,&fmax);
 return ((fmax-fmin)<=tol);
#else
 return ((xmax-xmin)<=tol) && ((ymax-ymin)<=tol);
#endif
}

void image(real xmin, real xmax, real ymin, real ymax, real* fmin, real* fmax)
{
#if 0
#elif 1
 AAform x,y,f,x2,y2,a;
 aa_open();
 aa_interval(x,xmin,xmax);
 aa_interval(y,ymin,ymax);
 aa_sqr(x2,x); aa_sqr(y2,y); aa_add(a,x2,y2); aa_inv(f,a);
 aa_trans(x2,x,-1); aa_sqr(x2,x2); aa_add(a,x2,y2); aa_inv(a,a); aa_add(f,f,a);
 aa_trans(f, f, -7.8);
 aa_range(f,fmin,fmax);
#elif 0
 AAform x,y,f,x2,y2;
 aa_open();
 aa_interval(x,xmin,xmax);		/* x = [ xmin .. xmax ] */
 aa_interval(y,ymin,ymax);		/* y = [ ymin .. ymax ] */
 aa_sqr(x2,x);			        /* x2= x*x */
 aa_sqr(y2,y);			        /* y2= y*y */
 aa_add(f,x2,y2);			/* f=x2+y2 */
 aa_inv(f,f);			        /* f=1/f */
 aa_trans(f, f, -1);			/* f=x2+y2-1 */
 aa_range(f,fmin,fmax);			/* f= [ fmin .. fmax ] */
 aa_close();
#elif 0
 AAform x,y,f,x2,y2;
 aa_open();
 aa_interval(x,xmin,xmax);		/* x = [ xmin .. xmax ] */
 aa_interval(y,ymin,ymax);		/* y = [ ymin .. ymax ] */
 aa_mul(x2,x,x);			/* x2= x*x */
 aa_mul(y2,y,y);			/* y2= y*y */
 aa_add(f,x2,y2);			/* f=x2+y2 */
 aa_trans(f, f, -1);			/* f=x2+y2-1 */
 aa_range(f,fmin,fmax);			/* f= [ fmin .. fmax ] */
 aa_close();
#elif 0
 AAform x,y,f,x2,x3,y2;
 aa_open();
 aa_interval(x,xmin,xmax);		/* x = [ xmin .. xmax ] */
 aa_interval(y,ymin,ymax);		/* y = [ ymin .. ymax ] */
 aa_mul(y2,y,y);			/* y2= y*y */
 aa_mul(x2,x,x);			/* x2= x*x */
 aa_mul(x3,x,x2);			/* x3= x*x2 */
 aa_sub(f,y2,x3);			/* f=y2-x3 */
 aa_add(f,f,x);				/* f=f+x=y2-x3+x */
 aa_range(f,fmin,fmax);			/* f= [ fmin .. fmax ] */
 aa_close();
#elif 0
 AAform x,y,f,x2,y2,xy,xy2;
 aa_open();
 aa_interval(x,xmin,xmax);		/* x = [ xmin .. xmax ] */
 aa_interval(y,ymin,ymax);		/* y = [ ymin .. ymax ] */
 aa_mul(x2,x,x);			/* x2= x*x */
 aa_mul(y2,y,y);			/* y2= y*y */
 aa_mul(xy,x,y);			/* xy= x*y */
 aa_mul(xy2,xy,xy);			/* xy2= x*y*x*y */
 aa_scale(xy2,xy2,0.5);
 aa_add(f,x2,y2);			/* f=x2+y2 */
 aa_add(f,f,xy);			/* f=f+xy=x2+y2+xy */
 aa_sub(f,f,xy2);			/* f=f-xy2=x2+y2+xy-(xy)2/2 */
 aa_trans(f, f, -0.25);
 aa_range(f,fmin,fmax);			/* f= [ fmin .. fmax ] */
 aa_close();
#elif 1
 AAform x,y,f,x2,y2,a;
 aa_open();
 aa_interval(x,xmin,xmax);
 aa_interval(y,ymin,ymax);
 aa_sqr(x2,x); aa_sqr(y2,y); aa_add(a,x2,y2); aa_inv(f,a);
 aa_trans(x2,x,-1); aa_sqr(x2,x2); aa_add(a,x2,y2); aa_inv(a,a); aa_add(f,f,a);
 aa_sqr(x2,x);
 aa_trans(y2,y,-1); aa_sqr(y2,y2); aa_add(a,x2,y2); aa_inv(a,a); aa_add(f,f,a);
 aa_trans(f, f, -8.0);
 aa_range(f,fmin,fmax);
#if 0
aa_trace("x",x);
aa_trace("y",y);
aa_trace("f",f);
#endif
 aa_close();
#endif
}

void plot(real xmin, real xmax, real ymin, real ymax)
{
 real x=(xmin+xmax)/2.0;
 real y=(ymin+ymax)/2.0;
 show(3,'f',xmin,xmax,ymin,ymax);
 show(1,'p',xmin,xmax,ymin,ymax);
 gpcolor(2); gpplot(x,y);
#if 0
 gpwait(-1);
#endif
}

#undef show

void show(int c, int t, real xmin, real xmax, real ymin, real ymax)
{
 gpcolor(c);
 gpbegin(t);
  gppoint(xmin,ymin);
  gppoint(xmax,ymin);
  gppoint(xmax,ymax);
  gppoint(xmin,ymax);
 gpend();
#if 0
 gpwait(-1);
#endif
}
