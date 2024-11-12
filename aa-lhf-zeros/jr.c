/*
* jr.c
* compute and plot joint range of two affine forms
* Luiz Henrique de Figueiredo -- bulk of code by Ken Clarkson.
* 13 Mar 2002 21:36:37
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <gp.h>
#include <jr.h>

#define NMAX (100)

static int ch2d(real** P, int n);

static real points[NMAX][2], *P[NMAX+1];	/* an extra position is used */
static int N=0;

static void addpoint(real x, real y)
{
 points[N][0]=x;
 points[N][1]=y;
 P[N]=points[N];
 ++N;
 return;
 printf("addpoint: N=%d\n",N);
}

static void show(int n, real* x, real* y)
{
 int i; 
 real xx=0.0;
 real yy=0.0;
 for (i=0; i<=n; i++)
 {
  xx+=x[i];
  yy+=y[i];
 }
 addpoint(xx,yy);
}

static void visit(int m, int n, real* x, real* y)
{
 if (m==0)
  show(n,x,y);
 else
 {
  visit(m-1,n,x,y);
  x[m]=-x[m]; y[m]=-y[m];
  visit(m-1,n,x,y);
  x[m]=-x[m]; y[m]=-y[m];
 }
}

real hull(int n, real* x, real* y)
{
 N=0;
 visit(n,n,x,y);
 N=ch2d(P,N);
 return 0;
}

void draw(int c, int out)
{
 int i;
 gpcolor(c);
 gpbegin('p');
 for (i=0; i<N; i++) gppoint(P[i][0],P[i][1]);
 gpend();
 if (!out) return;
if(0)return;
 printf("%% zonotope\n");
 for (i=0; i<N; i++) printf("%g %g %s\n",P[i][0],P[i][1],(i==0)?"M":"L");
 printf("E\n\n");
if(0)return;
 if (c!=0) gpcolor(9);
 for (i=0; i<N; i++) gpplot(P[i][0],P[i][1]);
}

/* Ken Clarkson [originally] wrote this.  Copyright (c) 1996 by AT&T. */

static int ccw(real** P, int i, int j, int k)
{
 real a = P[i][0] - P[j][0],
      b = P[i][1] - P[j][1],
      c = P[k][0] - P[j][0],
      d = P[k][1] - P[j][1];
 return (a*d-b*c) <= 0;		/* true if points i, j, k counterclockwise */
}

#define CMPM(c,A,B) \
	v = (*(real**)A)[c] - (*(real**)B)[c];\
	if (v>0) return 1;\
	if (v<0) return -1;

static int cmpl(const void *a, const void *b)
{
 real v;
 CMPM(0, a, b);
 CMPM(1, b, a);
 return 0;
}

static int cmph(const void *a, const void *b)
{
 return cmpl(b, a);
}

static int make_chain(real** V, int n, int (*cmp) (const void *, const void *))
{
 int i, j, s = 1;
 real *t;

 qsort(V, n, sizeof(real *), cmp);
 for (i = 2; i < n; i++)
   {
    for (j = s; j >= 1 && ccw(V, i, j, j - 1); j--);
    s = j + 1;
    t = V[s];
    V[s] = V[i];
    V[i] = t;
   }
 return s;
}

static int ch2d(real** P, int n)
{
 int u = make_chain(P, n, cmpl);/* make lower hull */
 if (!n)
  return 0;
 P[n] = P[0];
 return u + make_chain(P + u, n - u + 1, cmph);	/* make upper hull */
}
