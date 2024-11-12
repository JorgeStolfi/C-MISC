/*
* ifs.c
* ifs interpreter
* Luiz Henrique de Figueiredo
* 06 Mar 93
*/

#include <stdlib.h>
#include <gp.h>

#define frand()		(rand()/((real)RAND_MAX))

typedef struct
{
 real a,b,c,d,e,f;			/* coefficients */
 real p;				/* probability */
} Transform;

Transform w[]=
{
#if 1
 {  0.85,  0.04, -0.04, 0.85, 0.0, 1.60, 0.85 },	/* fern */
 { -0.15,  0.28,  0.26, 0.24, 0.0, 0.44, 0.07 },
 {  0.20, -0.26,  0.23, 0.22, 0.0, 1.60, 0.07 },
 {  0.00,  0.00,  0.00, 0.16, 0.0, 0.00, 0.01 }
#elif 0
 {  0.5 ,  0.0 ,  0.0 , 0.5 , 0.5, 0.5 , 0.34 },
 {  0.5 ,  0.0 ,  0.0 , 0.5 , 1.0, 0.0 , 0.33 },
 {  0.5 ,  0.0 ,  0.0 , 0.5 , 0.0, 0.0 , 0.33 }
#elif 0
 {  0.42,  0.42, -0.42, 0.42, 0.0, 0.2 , 0.40 },
 {  0.42, -0.42,  0.42, 0.42, 1.0, 0.2 , 0.40 },
 {  0.1 ,  0.0 ,  0.0 , 0.1 , 0.0, 0.2 , 0.15 },
 {  0.0 ,  0.0 ,  0.0 , 0.5 , 0.0, 0.0 , 0.05 }
#endif
};

#define N	(sizeof(w)/sizeof(w[0]))

real p[N];				/* cummulative probabilities */

void main(void)
{
 int i,n;
 real t,x,y,x0,y0;
 real xmax,xmin,ymax,ymin;
 gpopen("iterated function system");
 gpwindow(-3,3,0,10);
 gppalette(0,"black");
 gppalette(1,"blue");
 gppalette(2,"red");
 gpcolor(1);
 for (t=0.0, i=0; i<N; i++)		/* compute cummulative probabilities */
  p[i]=(t+=w[i].p);
 x=y=0.0;
 xmax=ymax=-1e+38;
 xmin=ymin=+1e+38;
 for (n=0; ; n++)
 {
  x0=x; y0=y;
  t=frand();
  for (i=0; i<N; i++)
   if (t<=p[i]) break;
  t=w[i].a*x+w[i].b*y+w[i].e;
  y=w[i].c*x+w[i].d*y+w[i].f;
  x=t;
  if (x<xmin) xmin=x; else if (x>xmax) xmax=x;
  if (y<ymin) ymin=y; else if (y>ymax) ymax=y;
  gpcolor(1); gpplot(x,y);
  gpcolor(2); gpline(x0,y0,x,y);
  if ((n&1023)==0)
  {
   real x,y; char*s=gpevent(0,&x,&y);
   if (s!=NULL && *s=='k') break;
  }
 }
 gpclose(0);
 printf("n=%d\ngpwindow(%g,%g,%g,%g);\n",n,xmin,xmax,ymin,ymax);
}
