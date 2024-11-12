/* See SPTimeBasisG3C1.h */
/* Last edited on 2005-08-27 21:18:53 by stolfi */ 

#include <SPTimeBasisG3C1.h>

#include <bool.h>
#include <affirm.h>

#include <stdio.h>

void SPTimeBasisG3C1_BasisCoeffs(int k, int u, double tStep, int diff, double *tc)
  { 
    double h = tStep;
    int tg = 3; /* Degree of {tau}. */
    switch (u)
      {
        case 0:
          switch(k)
            {
              case 0: 
                /* Coefficients of {tau[0,0](z-1)} */
                /* z^2(3 - 2z) = -2z^3 + 3z^2 */
                tc[0] = 00.0; tc[1] = 00.0; tc[2] = +3.0; tc[3] = -2.0;
                break;
              case 1:
                /* Coefficients of {tau[0,0](z)} */
                /* (z - 1)^2(2z + 1) = 2z^3 - 3z^2 + 1 */
                tc[0] = +1.0; tc[1] = 00.0; tc[2] = -3.0; tc[3] = +2.0;
                break;
              default:
                tc[0] = 00.0; tc[1] = 00.0; tc[2] = 00.0; tc[3] = 00.0;
            }
          break;
        case 1:
          switch(k)
            {
              case 0:
                /* Coefficients of {tau[0,1](z-1)}: */
                /* (z-1)z^2 = z^3 - z^2 */
                tc[0] = 00.0; tc[1] = 00.0; tc[2] = -h; tc[3] = +h;
                break;
              case 1:
                /* Coefficients of {tau[0,1](z)}: */
                /* z(1-z)^2 = z^3 - 2z^2 + z */
                tc[0] = 00.0; tc[1] = +h; tc[2] = -2*h; tc[3] = +h;
                break;
              default:
                tc[0] = 00.0; tc[1] = 00.0; tc[2] = 00.0; tc[3] = 00.0;
            }
          break;
        default:
          affirm(FALSE, "bad {u}");
      }

    /* Differentiate the time element with respect to time, {diff} times: */
    int tdg = tg; /* Degree of the differentiated polynomial. */
    int id, i;
    for (id = 0; id < diff; id++)
      { /* Diffrentiate the polynomial {tc[0..tdg]}.
          Note that the coefficients {tc} are for {z}, not {t},
          so we must multiply them by {dz/dt = 1/tStep}:
        */
        for (i = 1; i <= tdg; i++) { tc[i-1] = tc[i]*i/tStep; }
        if (tdg > 0) { tc[tdg] = 0.0; tdg--; }
      }
  }

double SPTimeBasisG3C1_BasisEval(int u, double t, double tj, double tStep, int diff)
  {
    /* Compute relative position {z} of {t} in interval {[tj-tStep __ tj]}: */
    double z = (t - tj)/tStep + 1;
    if ((z < 0) || (z > 2)) { return 0.0; }
    
    /* Find interval {[k-1 __ k]} and reduce {z} to that interval: */
    int k = 0;
    if (z > 1) { k = 1; z -= 1; }

    /* Get coeffs of {TD(tau[j,u],diff)} wrt {z} in that interval: */
    int tg = 3;
    int tdg = (diff > tg ? 0 : tg - diff);
    double tc[tg+1];
    SPTimeBasisG3C1_BasisCoeffs(k, u, tStep, diff, tc);
    /* Evaluate polynomial {tc[0..tdg]} at {z}: */
    int i;
    double tcVal = 0.0;
    for (i = tdg; i > 0; i--) { tcVal = (tcVal + tc[i])*z; }
    tcVal += tc[0];
    if (FALSE) 
      { int i;
        fprintf(stderr, "  u = %d  t = %g  tj = %g  diff = %d\n", u, t, tj, diff);
        fprintf(stderr, "  tc = (");
        for (i = tdg; i >= 0; i--) { fprintf(stderr, " %5g", tc[i]); }
        fprintf(stderr, " ) = %g\n", tcVal);
      }
    
    return tcVal;
  }

void SPTimeBasisG3C1_GaugeCoeffs(int v, double *pc)
  {
    switch(v)
      {
        case 0: 
          /* Coefficients of {pi[0,0](z-1)} = 1 */
          pc[0] = +1.0; pc[1] = 00.0;
          break;
        case 1:
          /* Coefficients of {pi[0,1](z-1)} = 2(z-1) + 1 = 2z - 1 */
          pc[0] = -1.0; pc[1] = +2.0;
          break;
        default:
          affirm(FALSE, "bad {v}");
      }
  }

double SPTimeBasisG3C1_GaugeEval(int v, double t, double tj, double tStep)
  {
    int pg = 1; /* Max degree of time gauge element. */
    /* Compute position of {t} in interval {[tj-tStep __ tj]}: */
    double z = (t - tj)/tStep + 1;
    if ((z < 0) || (z > 1)) { return 0.0; }
    /* Get coeffs of {pi[0,v]} wrt {z} in that interval: */
    double pc[pg+1];
    SPTimeBasisG3C1_GaugeCoeffs(v, pc);
    /* Evaluate polynomial {pc[0..pg]}: */
    return pc[0] + pc[1]*z;
  }

double SPTimeBasisG3C1_BasisGaugeDot(int k, int u, int v, double tStep, int diff, bool_t verbose)
  { 
    int tg = 3; /* Degree of {tau[0,u]}. */
    int tdg = (diff > tg ? 0 : tg - diff); /* Deg of {TD(tau[0,u],diff)}. */
    int pg = 1; /* Degree of {pi[0,v]}. */
    /* All coeffs below are for variable {z = t-1-k}: */
    /* Compute coeffs of {TD(tau[0,u],diff)} between {k-1} and {k}: */ 
    double tc[tg + pg + 1];
    SPTimeBasisG3C1_BasisCoeffs(k, u, tStep, diff, tc);
    
    /* Compute coeffs of {pi[0,v]} between {-1} and {0}: */ 
    double pc[pg + 1];
    SPTimeBasisG3C1_GaugeCoeffs(v, pc);
    
    /* Multiply {pc[0..pg]} into {tc[0..tdg]} yielding {tc[0..pg+tdg]}: */
    affirm(pg == 1, "bad pg");
    tc[tdg+1] = tc[tdg]*pc[1];
    int i;
    for (i = tdg; i > 0; i--) { tc[i] = tc[i-1]*pc[1] + tc[i]*pc[0]; }
    tc[0] = tc[0]*pc[0];
    tdg += pg;
    
    /* Integrate product: */ 
    double sum = 0;
    for (i = 0; i <= tdg; i++) { sum += tc[i]/(i + 1); }
    return tStep*sum;
  }

