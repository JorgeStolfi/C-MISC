/* Example with one equation (circle) and square X range */
/* Last edited on 2012-07-22 12:48:52 by stolfilocal */

#include <fcircle.h>
#include <gp.h>

#define N          fcircle_N
#define M          fcircle_M
#define P          fcircle_P
#define func_args  fcircle_func_args
#define func_aa    fcircle_func_aa
#define func_show  fcircle_func_show

#define XCTR -0.3
#define YCTR -0.2
#define RAD 1.2

void fcircle_func_args(real xmin, real xmax, real ymin, real ymax, AAform X[], int S[])
  {
    aa_interval(X[0],xmin,xmax); S[0] = aa_last(); 
    aa_interval(X[1],ymin,ymax); S[1] = aa_last();
  }

void fcircle_func_aa(AAform X[], AAform Y[])
  {
    AAform a,b;
    aa_clear(a); aa_clear(b); aa_clear(Y[0]);
    aa_trans(a,X[0],-(XCTR)); aa_sqr(a,a);
    aa_trans(b,X[1],-(YCTR)); aa_sqr(b,b);
    aa_add(Y[0],a,b);
    aa_trans(Y[0],Y[0],-(RAD*RAD));
  }

void fcircle_func_show(void)
  {
    gpcircle(XCTR,YCTR,RAD);
  }
