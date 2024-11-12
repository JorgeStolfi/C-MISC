/* See SPTimeOnlyIntegral.h */
/* Last edited on 2005-08-21 14:02:04 by stolfi */

#include <SPTimeOnlyIntegral.h>
#include <SPBasic.h>
#include <vec.h>
#include <affirm.h>
#include <nat.h>
#include <math.h>
#include <stdlib.h>

void SPTimeOnlyIntegral_BySamples
  ( double func(double t),
    double_vec_t st,
    double_vec_t wt,
    double *sum,
    double *corr
  )
  { int j;
    for (j = 0; j < st.ne; j++)
      { double stj = st.e[j];
        double wtj = wt.e[j];
        double term = wtj*func(stj);
        /* Kahan's summation formula: */
        double tcorr = term - *corr;
        double newSum = *sum + tcorr;
        *corr = (newSum - *sum) - tcorr;
        *sum = newSum;
      }
  }
  
double SPTimeOnlyIntegral_OnInterval
  ( double func(double t),
    int smpOrder,
    double tMin,
    double tMax
  )
  { double sum = 0.0, corr = 0.0;
    int ns = 0;
    double stel[smpOrder]; double_vec_t st = {smpOrder, stel};
    double wtel[smpOrder]; double_vec_t wt = {smpOrder, wtel};
    SPTimeOnlyIntegral_GaussSampleInterval(tMin, tMax, smpOrder, &st, &wt, &ns);
    SPTimeOnlyIntegral_BySamples(func, st, wt, &sum, &corr);
    return sum;
  }

/* GENERATING SAMPLE POINTS AND WEIGHTS */

void SPTimeOnlyIntegral_GaussSampleInterval
  ( double tMin, double tMax, 
    int smpOrder, 
    double_vec_t *st, 
    double_vec_t *wt, 
    int *ns
  )
  { 
    /* ??? Find correct formulas! */
    int i;
    int m = *ns;
    double tWid = tMax - tMin;
    for (i = 0; i < smpOrder; i++)
      { /* Compute knots and weights relative to [-1 _ +1]: */
        st->e[m] = tMin + tWid*(i + 0.5)/smpOrder;
        wt->e[m] = tWid/smpOrder;
        m++;
      }
    (*ns) = m;
  }

