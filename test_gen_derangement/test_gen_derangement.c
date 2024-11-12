/* Last edited on 2013-01-12 15:36:46 by stolfilocal */
/* Checks the derangement generator used in the multiscale DNA matching package. */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#include <jsrandom.h>

void gen_biased_derangement(int n, int perm[]);
  /* Fills {perm[0..n-1]} with a derangement of {0..n-1} that is biased
    towards the identity, that is, {perm[i]} tends to be close to (but
    distinct from) {i}. Uses the {jsrandom.h} random generator, in its
    current state. */

int main(int argc, char **argv);

int main(int argc, char **argv)
  {
    /* Generate a random derangement: */
    int n = 100;
    int perm[n];
    srandom(46150001); /* Initialize the randomizer, for repeatability. */
    gen_biased_derangement(n, perm);
    
    /* Validate and compute the average relative displacement: */
    double sumr = 0;
    double sumd = 0;
    int i;
    for (i = 0; i < n; i++) 
      { fprintf(stderr, "%5d --> %5d\n", i, perm[i]);
        assert(perm[i] != i);
        double ri = ((double)perm[i])/((double)i+1);
        double di = (double)abs(i - perm[i]);
        sumr += ri;
        sumd += di;
      }
    fprintf(stderr, "\n");
    fprintf(stderr, "average relative position = %8.6f\n", sumr/n);
    fprintf(stderr, "average absolute displacement = %8.6f\n", sumd/n);
    return 0;
  }

void gen_biased_derangement(int n, int perm[])
  { 
    perm[0] = 0;
    int i;
    for(i = 1; i < n ; i++){
      /* Choose a random index {j} in {0..i-1}, make {{perm[i],perm[j] = {perm[j],i}}. */
      /* So bias {j} towards small {i-perm[j]} and {i-j}. */
      int j = i-1;
      while (j > 0) 
        { int di = i-perm[j];
          int dj = i-j;
          double d2 = (double)(di*di + dj*dj);
          double r = drandom();
          if (r*r*d2 < 6.0) { break; }
          j--;
        }
      assert((j >= 0) && (j < i));
      assert((perm[j] >= 0) && (perm[j] < i));
      /* Send {i} to {perm[j]}, {j} to {i}: */
      int t = perm[j];
      perm[j] = i;
      perm[i] = t;
    }
  }
