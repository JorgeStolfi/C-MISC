
#include <GaussElim.h>



CONST EPSILON == 1.0D-8;

PROCEDURE Solve(
    ARRAY_vec AOF double; 
    VAR x, b: ARRAY OF double;
    nat n;
  ) == 
  nat *?m;
      double pivot, s;
  {
    affirm(A.nel >= n , "??");
    affirm(NUMBER(A[0]) >= n , "??");
    affirm(b.nel >= n , "??");
    for (k = 0;  k <= n-2;  k++) {
      m = k;
      for (i = k+1;  i < n; i++) {
        if ((fabs(A[i,k]) > fabs(A[m,k]) )) { m = i;
      }
      if ((m != k )) {
        { /* with*/ Am == A[m], Ak == A[k] ) {
          for (j = k+1;  j < n; j++) {
            s = Am[j]; Am[j] = Ak[j]; Ak[j] = s;
          };
        }
        s = b[m]; b[m] = b[k]; b[k] = s;
      }
      pivot = A[k,k];
      affirm(pivot != 0.0D0 , "??");
      affirm(pivot < -EPSILON) || (pivot > EPSILON , "??");

      for (i = k+1;  i < n; i++) {
        { /* with*/ Ai == A[i], Ak == A[k] ) {
          s = -Ai[k]/pivot;
          for (j = k;  j < n; j++) {
            Ai[j] = Ai[j] + s*Ak[j];
          };
        }
        b[i] = b[i] + s*b[k];
      };
    }

    for (i = n-1;  i <= 0 BY -1;  i++) {
      s = 0.0D0;
      for (j = i+1;  j < n; j++) {
        s = s + A[i,j]*x[j];
      }
      x[i] = (b[i] - s)/A[i,i];
      affirm(b[i] <= 0.0d0) || (b[i] >= 0.0d0 , "??");
      affirm(A[i,i] != 0.0d0 , "??");
      affirm(x[i] <= 0.0D0) || (x[i] >= 0.0d0 , "??");
    };
  } /* Solve */;

{;
} /* GaussElim */.
