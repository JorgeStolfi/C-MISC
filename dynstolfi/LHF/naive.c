/*
* henon.c
* henon attractor
* Luiz Henrique de Figueiredo
* 20 May 2004 16:04:11
*/

#include <gp.h>

#define	a	(1.4)
#define	b	(0.3)

void f(real* x, real* y)
{
 real X=*x;
 real Y=*y;
 *x=1.0-a*X*X+Y;
 *y=b*X;
}

void attractor(real x, real y)
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

int main(void)
{
 real x;
 gpopen("henon attractor");
 gpwindow(-1.6,1.6,-0.5,0.5);
 if(0)
 for (x=0.0; x<=1.0; x+=0.01)
  attractor(x,0.25);
 else
 attractor(0,0);
 gpclose(1);
 return 0;
}
