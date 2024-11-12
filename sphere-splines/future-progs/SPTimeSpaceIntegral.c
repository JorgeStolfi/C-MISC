/* See SPTimeSpaceIntegral.h */
/* Last edited on 2005-08-14 17:51:00 by stolfi */

#include <SPTimeSpaceIntegral.h>
#include <SPTimeOnlyIntegral.h>
#include <SPIntegral.h>
#include <SPTriang.h>
#include <SPBasic.h>
#include <vec.h>
#include <r3.h>
#include <r3x3.h>
#include <affirm.h>
#include <nat.h>
#include <math.h>
#include <stdlib.h>

void SPTimeSpaceIntegral_BySamples
  ( double func(S2Point *p, double t),
    S2Point_vec_t sp,
    double_vec_t wp,
    double_vec_t st,
    double_vec_t wt,
    double *sum, 
    double *corr
  )
  { int i,j;
    for (j = 0;  j < st.ne;  j++)
      { double stj = st.e[j];
        double wtj = wt.e[j];
        for (i = 0;  i < sp.ne;  i++)
          { S2Point *spi = &(sp.e[i]);
            double wpi = wp.e[i];
            double term = wpi*wtj*func(spi, stj);
            /* Kahan's summation formula: */
            double tcorr = term - *corr;
            double newSum = *sum + tcorr;
            *corr = (newSum - *sum) - tcorr;
            *sum = newSum;
          }
      }
  }
  
double SPTimeSpaceIntegral_OnSphereInterval
  ( double func(S2Point *p, double t),
    Triangulation *tri,
    int smpOrderTime,
    double tMin,
    double tMax
  )
  { int i;
    double sum = 0.0, corr = 0.0;
    int ns = 0;
    double stel[smpOrderTime]; double_vec_t st = {smpOrderTime, stel};
    double wtel[smpOrderTime]; double_vec_t wt = {smpOrderTime, wtel};
    if (tri == NULL) { tri = SPIntegral_GetDefaultTriangulation(); }
    SPTimeOnlyIntegral_GaussSampleInterval(tMin, tMax, smpOrderTime, &st, &wt, &ns);
    for (i = 0;  i < tri->side.ne;  i++)
      { Arc e = tri->side.e[i];
        Face *f = Left(e);
        SPTimeSpaceIntegral_BySamples(func, f->sp, f->wp, st, wt, &sum, &corr);
      }
    return sum;
  }
