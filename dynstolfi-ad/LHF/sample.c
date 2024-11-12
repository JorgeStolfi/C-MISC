/*
* henon.c
* henon attractor
* Luiz Henrique de Figueiredo
* 20 May 2004 16:04:11
*/

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

static void sample(real xmin, real xmax, real ymin, real ymax, int n)
{
 real dx=(xmax-xmin)/n;
 real dy=(ymax-ymin)/n;
 int i,j;
 for (i=0; i<=n; i++)
 {
  for (j=0; j<=n; j++)
  {
   real x=xmin+i*dx;
   real y=ymin+j*dy;
   gpcolor(2); gpplot(x,y);
   f(&x,&y);
   gpcolor(3); gpplot(x,y);
  }
 }
}
static void orbit(real x, real y, int n)
{
 int c=2;
 while (n--) { gpcolor(c); gpplot(x,y); f(&x,&y); c=5-c; }
}

int main(int argc, char* argv[])
{
 real r=0.1;
 int n=5;
 if (argc>1) sscanf(argv[1],"%lg",&r);
 if (argc>2) sscanf(argv[2],"%d",&n);
 gpopen("henon attractor");
 gpwindow(-1.5,1.5,-1.5,1.5);
 gpwindow(-100,100,-100,100);
 gppalette(0,"black");
 gppalette(1,"white");
 gppalette(2,"red");
 gppalette(3,"green");
 gppalette(4,"blue");
 gpmark(3,"B");
 sample(-1.5,1.6,-0.5,0.5,100);
 while (1)
 {
  real x,y;
  char* s=gpevent(1,&x,&y);
  if (s[0]=='k' && s[1]=='q') {gpclose(0); break;}
  if (s[0]=='k' && s[1]=='c') {gpclear(0); continue;}
  if (s[0]=='k' && s[1]=='w') {printf("%g\t%g\n",x,y); continue;}
  if (s[0]=='k' && s[1]=='o') {orbit(x,y,1000); continue;}
  if ((s[0]=='b' && s[2]=='-') || s[0]=='m')
  {
   orbit(x,y,2);
   continue;
  }
 }
 return 0;
}
