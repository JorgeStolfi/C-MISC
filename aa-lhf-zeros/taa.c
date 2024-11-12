/* Last edited on 2009-02-11 15:48:09 by stolfi */

#include <math.h>
/* #include <stddef.h> */
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <aa.h>
#include <gp.h>
#include <jr.h>
#include <aasolve.h>

#include <fcirhex.h>
#include <fcircle.h>
#include <f2circ.h>

/* Define N, M, P, func_args, func_aa, func_show: */
#define N          fcirhex_N
#define M          fcirhex_M
#define P          fcirhex_P
#define func_args  fcirhex_func_args
#define func_aa    fcirhex_func_aa
#define func_show  fcirhex_func_show

static void solve(real xmin, real xmax, real ymin, real ymax)
{
 int K, r, s;
 AAform X[N],Y[M];
 int S[P];
 aa_open();
 assert(N == 2);
 func_args(xmin, xmax, ymin, ymax, X, S);
 func_aa(X,Y);
 hull(aa_last(),X[0],X[1]); draw(4,0);
 gpwait(-1);
 for (r = 0; r < N; r++) { fprintf(stderr, "X[%d] = ", r); aa_trace("", X[r]); }
 for (s = 0; s < M; s++) { fprintf(stderr, "Y[%d] = ", s); aa_trace("", Y[s]); }
 K = aa_solve(N, X, M, Y, P, S);
 fprintf(stderr, "eliminated %d noises\n", K);
 for (r = 0; r < N; r++) { fprintf(stderr, "X[%d] = ", r); aa_trace("", X[r]); }
 for (s = 0; s < M; s++) { fprintf(stderr, "Y[%d] = ", s); aa_trace("", Y[s]); }
 hull(aa_last(),X[0],X[1]); draw(2,0);
 aa_close();
}

static void try(real xmin, real xmax, real ymin, real ymax)
{
 gpopen("solver");
 gppalette(0,"black");
 gppalette(1,"white");
 gppalette(2,"white");
 gppalette(3,"green");
 gppalette(4,"yellow");
 gpwindow(-2,2,-2,2);
 gpcolor(3);
 func_show();
 gpflush();
 solve(xmin,xmax,ymin,ymax);
 gpclose(1);
}

int main(int argc, char* argv[])
{
 try(0,1,0,1);
 return 0;
}
