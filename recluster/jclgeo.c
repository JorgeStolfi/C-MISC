/* Last edited on 2007-08-15 22:36:44 by stolfi */
/* See the authorship and copyright notice at the end of this file. */
/* Changes:

*/

#include <jclbasic.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <assert.h>

double L2ProbDist(double *a, double *b, long N)
  {
    double s = 0.;
    long i;
    for (i=0; i<N; i++) { double d = a[i] - b[i]; s += d*d; }
    return sqrt(s/2);
  }
      
double EarthMoverDist(double *a, double *b, long N)
  {
    double s = 0;
    double sa = 0;
    double sb = 0;
    long i;
    if (N==0) return (0);
    for (i=0; i<N; i++) 
      { double d;
        sa += a[i]; sb += b[i];
        /* "sa" is the earth mass in "a[0..i]", ditto for "sb" */
        d = fabs(sa - sb); 
        /* "d" is the mass that must be carried from "[i]" to "[i+1]" */
        s += d;
      }
    return sqrt(s/((double)(N-1)));
  }

void NormalizeProbs(double *a, long N, long R)
  { double s = 0;
    long i;
    for (i=0; i<N; i++) 
      { if (a[i] < 0) { ElementError("negative probability", R, i); }
        s += a[i];
      }
    if (s <= 0) { LineError("null distribution", R); }
    for (i=0; i<N; i++) a[i] /= s;
  }

void EstimateProbs(double *a, long N, long R)
  { double s = 0;
    long i;
    for (i=0; i<N; i++) 
      { if (a[i] < 0) { ElementError("negative count", R, i); }
        s += a[i];
      }
    s += (double)N;
    for (i=0; i<N; i++) a[i] = (a[i]+1)/s;
  }

