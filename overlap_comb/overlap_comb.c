/* Last edited on 2024-12-21 11:55:21 by stolfi */
/* Plots the smoothness of overlapped combinatorial pulses. */

/* 
  The binomial window weights {w[0..n-1]} are {w[i] = comb(n,i)/2^n}.
  The question is whether these weights are suitable for filtering 
  and downsampling; in particular, whether they are a partition
  of unit when used with steps larger than 1.
  
  Let's define a /binomial pulse/ the biinfinite signal {comb(n,k)},
  viewed as a function of {k\in\RZ} with {n} fixed. Let {X^n_d} be the
  sum of an infinite train of unscaled binomial pulses
  displaced by multiples of a step {d}. That is, {X^n_d[k] = SUM{
  comb(n,k-i*d) : i\in\RZ }}. The signal {X^n_d} is obviously periodic
  with period {d} and nowhere negative. The binomial weights are
  partition-of-unity iff {X^n_d[k]} is constant for all {k}.
  
  The convolution {X^n_d} is obviously constant for {d=1}, with total
  value {2^n}. On the other hand, for {d>=n} the binomial pulses do not
  overlap. For {d=n} the pulses are adjacent, and the convolution {x}
  ranges from 1 to {comb(n,n/2)}. For {d>n} the pulses are isolated some
  {x[j]} are zero.
  
  Somewhat surprisingly, the convolution {X^n_d} is constant also for {d=2}.
  That is because the sum of {comb(n,k)} for all even {k} is equal to
  the sum for all odd {k}, namely {2^[n-1}}. Thus the combinatorial
  weights are partition of unity also for downsampling a signal by 2:1.
  
  Even more surprisingly, for {d=3} the convolution is *almost* uniform.
  If it were uniform, all values should be equal to {(2^n)/d}. But {2^n}
  is not divisible by 3; indeed {2^n} is of the form {3*m+1} if {n} is
  even, and {3*m-1} if {n} is odd. In the first case the sequence
  {X^n_d} is {...,m+1,m,m,m+1,m,m,m+1,...} where the {m+1} values
  coincide with the centers of the binomial pulses ({j} congruent to
  {n/2} modulo 3). In the second case {X^n_d} is
  {...,m-1,m,m,m-1,m,m,...} where the {m} values lie on each side of the
  pulse centers ({j} congruent with {(n-1)/2} and {(n+1)/2} modulo 3).
  
  
  For {d} larger than 3 the convolution {x[j]} is 
  increasingly non-uniform, and the absolute variation increases
  with both {n} and {d}.  For {d=4} the difference is 
  {2^{\floor{n/2}}: 1 (for {n=1}),2,2,4,4,8,8,.... However the 
  relative variation decreases like {2/(1+2^{\floor{n/2})}.
  in general the maxima is at or adjacent to the center
  of the binomial pulses, and the minima halfway between them. */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <math.h>

#include <bool.h>
#include <jsfile.h>
#include <jsmath.h>

void gen_plot(FILE *wr, int n, int dmin, int dmax);
  /* Generates a plot file showing how well the sequence {x[k] = comb(n,k)}
    overlaps with itself under various displacements. */

int main(int argc, char **argv)
  {
    assert(argc == 5);
    int nmin = atoi(argv[1]); /* Min number of items. */
    int nmax = atoi(argv[2]); /* Max number of items. */
    int dmin = atoi(argv[3]); /* Min displacement btw overlaps. */
    int dmax = atoi(argv[4]); /* Max displacement. */

    int n;
    for (n = nmin; n <= nmax; n++)
      { char *fname = NULL;
        char *fname = jsprintf("out/comb-%04d.txt", n);
        FILE *wr = open_write(fname, TRUE);
        free(fname);
        gen_plot(wr, n, dmin, dmax);
        fclose(wr);
      }

    return 0;
  }

void gen_plot(FILE *wr, int n, int dmin, int dmax)
  {
    int d;
    for (d = dmin; d <= dmax; d++)
      { /* Output plot data with displacement {d}. */
        /* bool_t debug = ((d - dmin) <= 4) || ((dmax - d) <= 2); */
        bool_t debug = FALSE;
        int k;
        uint64_t smin = ~(uint64_t)0;
        uint64_t smax = 0;
        double hmin; /* Position of {smin} relative to pulse center (integer or half-integer). */
        double hmax; /* Position of {smax} relative to pulse center (integer or half-integer). */
        double kmed = ((double)n)/2; /* Position of pulse center rel to its first point. */
        for (k = 0; k < d; k++)
          { /* Compute sum {s} of {comb(n,k+r*d)} for all int {r}. */
            uint64_t s = 0;
            int rmin = -(k/d);
            int rmax = (n-k)/d;
            int r;
            for (r = rmin; r <= rmax; r++)
              { uint64_t c = comb(n, k+r*d);
                if (debug) { fprintf(stderr, "    n = %3d  k = %3d  comb(n,k) = %22llu\n", n,k+r*d,c); }
                s += c;
              }
            if (debug) { fprintf(stderr, "  n = %3d  d = %3d k = %3d  s = %22llu\n", n,d,k,s); }
            if (s < smin) { smin = s; hmin = k - kmed; }
            if (s > smax) { smax = s; hmax = k - kmed; }
          }
        /* Reduce {hmin,hmax} to the range {-d/2..+d/2} modulo d: */
        while (2*hmin < -d) { hmin += d; }
        while (2*hmin > +d) { hmin -= d; }
        while (2*hmax < -d) { hmax += d; }
        while (2*hmax > +d) { hmax -= d; }
        uint64_t sdif = smax - smin;
        double tot = exp2(n)/d;
        double fmin = ((double)smin)/tot;
        double fmax = ((double)smax)/tot;
        double f = ((double)(smax - smin))/((double)smax);
        fprintf(wr, 
          "%5d %5d  %22llu %22llu %22llu  %23.15e %23.15e  %23.15e %7.1f %7.1f\n", 
          n, d, smin,smax,sdif, fmin,fmax,f, hmin,hmax
        );
      }
  }
