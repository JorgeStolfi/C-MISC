/* Example with two equations (circles) and square X range */
/* Last edited on 2009-02-11 15:50:05 by stolfi */

#include <f2circ.h>
#include <gp.h>

#define N          f2circ_N
#define M          f2circ_M
#define P          f2circ_P
#define func_args  f2circ_func_args
#define func_aa    f2circ_func_aa
#define func_show  f2circ_func_show

#define AXCTR -0.3
#define AYCTR -0.2
#define ARAD 1.2

#define BXCTR +0.3
#define BYCTR -0.3
#define BRAD 1.1

void f2circ_func_args(real xmin, real xmax, real ymin, real ymax, AAform X[], int S[])
  {
    aa_interval(X[0],xmin,xmax); S[0] = aa_last(); 
    aa_interval(X[1],ymin,ymax); S[1] = aa_last();
  }

void f2circ_func_aa(AAform X[], AAform Y[])
{
 AAform a,b;
 aa_clear(a); aa_clear(b); aa_clear(Y[0]);
 aa_trans(a,X[0],-(AXCTR)); aa_sqr(a,a);
 aa_trans(b,X[1],-(AYCTR)); aa_sqr(b,b);
 aa_add(Y[0],a,b);
 aa_trans(Y[0],Y[0],-(ARAD*ARAD));

 aa_clear(a); aa_clear(b); aa_clear(Y[1]);
 aa_trans(a,X[0],-(BXCTR)); aa_sqr(a,a);
 aa_trans(b,X[1],-(BYCTR)); aa_sqr(b,b);
 aa_add(Y[1],a,b);
 aa_trans(Y[1],Y[1],-(BRAD*BRAD));
}

void f2circ_func_show(void)
{
  gpcircle(AXCTR,AYCTR,ARAD);
  gpcircle(BXCTR,BYCTR,BRAD);
}
