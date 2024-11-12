/* Example with one equation (circle) and non-square X range */
/* Last edited on 2012-07-22 12:49:20 by stolfilocal */

#include <fcirhex.h>
#include <gp.h>

#define N          fcirhex_N
#define M          fcirhex_M
#define P          fcirhex_P
#define func_args  fcirhex_func_args
#define func_aa    fcirhex_func_aa
#define func_show  fcirhex_func_show

#define XCTR -0.4
#define YCTR -0.1
#define RAD 1.2

void fcirhex_func_args(real xmin, real xmax, real ymin, real ymax, AAform X[], int S[])
  { /* Create an hexagon inside the box [xmin_xmax] x [ymin_ymax] */
    real xmid = (3*xmin + 2*xmax)/5;
    real ymid = (4*ymin + 3*ymax)/7;
    AAform a, b;
    aa_clear(a); 
    aa_interval(a, 0, 1); S[0] = aa_last();
    aa_interval(X[0],xmin,xmid); S[1] = aa_last();
    aa_interval(X[1],ymin,ymid); S[2] = aa_last();
    
    aa_clear(b); aa_scale(b, a, xmax-xmid);
    aa_add(X[0], X[0], b);
    aa_clear(b); aa_scale(b, a, ymax-ymid);
    aa_add(X[1], X[1], b);
  }

void fcirhex_func_aa(AAform X[], AAform Y[])
  {
    AAform a,b;
    aa_clear(a); aa_clear(b); aa_clear(Y[0]);
    aa_trans(a,X[0],-(XCTR)); aa_sqr(a,a);
    aa_trans(b,X[1],-(YCTR)); aa_sqr(b,b);
    aa_add(Y[0],a,b);
    aa_trans(Y[0],Y[0],-(RAD*RAD));
  }

void fcirhex_func_show(void)
  {
    gpcircle(XCTR,YCTR,RAD);
  }
