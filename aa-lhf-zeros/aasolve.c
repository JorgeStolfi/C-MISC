/* zeros of function R^N -> R^M */
/* Last edited on 2024-12-21 11:57:11 by stolfi */

#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <values.h>

#include <aa.h>
#include <aasolve.h>

#define TINYVY 1.0e-6

int aa_solve(int N, AAform *X, int M, AAform *Y, int P, int *S)
  { int L = aa_last();
    int K = 0;
    int J[AA_N+1];
    
    /* The vector J[0..L] contains some permutation of the currently
      valid AA term indices 0..L; and entries J[0..P-1] are a
      permutation of the desired noises S[0..P-1]. At the end of the
      procedure, entries J[0..K] = S[0..K] are the indices of the
      noise symbols which were eliminated. Specifically, entry J[i],
      for i in 0..K-1, is the index of the noise symbol which was
      eliminated using equation Y[i]==0. Columns J[0],.. J[K-1] of
      Y[0],.. Y[K-1], in that order, form a KxK identity matrix; i.e.
      Y[i][J[j]] == (i==j) for any i,j in 0..K-1.

      Normally, the number K of discarded noises is equal to the number
      M of equations; but may be less, if the rows Y[0..M-1], excluded
      the column 0, are linearly dependent.  If this happens, then,
      after the recombination, the forms Y[K] through Y{M-1] will
      be all zero except for the constant term. */

    auto real select_next(void);
    
    real select_next(void)
      /* Permutes J[K..P-1] so that J[K] is the next noise symbol
        to be eliminated, if any. */
      { real vx[AA_N+1];  /* vx[jt] is norm_sqr of X[0..N-1][jt] */
        real vxtot = 0.0; /* SUM { vx[jt] : jt in J[K..P-1] } */
        real vy[AA_N+1];  /* vy[jt] is norm_sqr of Y[K..M-1][jt] */
        real vytot = 0.0; /* SUM { vy[jt] : jt in J[K..P-1] } */
        real maxeffect = -INFINITY;
        int t;
        /* We try to choose the noise which, if eliminated, would have
          maximally beneficial effect in the uncertainty of X[0..N-1],
          considering only columns J[K..P-1]. The following heuristic is
          just a guess. */
        
        fprintf(stderr, "selecting J[%d] among {", K);
        for (t = K; t < P; t++) { fprintf(stderr, " %d", J[t]); }
        fprintf(stderr, " }  norms: ");
        assert(K < M);
        assert(K < P);
        for (t = K; t < P; t++)
          { int jt = J[t]; 
            assert(jt != 0);
            { int r; real sum = 0.0;
              for (r = 0; r < N; r++) { real x = X[r][jt]; sum += x*x; }
              vx[jt] = sum; vxtot += sum;
            } 
            { int s; real sum = 0.0;
              for (s = K; s < M; s++) { real y = Y[s][jt]; sum += y*y; }
              vy[jt] = sum; vytot += sum;
            }
            fprintf(stderr, " (%g,%g)@%d", vx[jt], vy[jt], jt);
          }
        fprintf(stderr, "\n");
        fprintf(stderr, "est effects: ");
        for (t = K; t < P; t++)
          { int jt = J[t]; 
            if ((jt != 0) && (vy[jt] > TINYVY*vytot))
              { /* Guessed desirability of eliminating noise symbol jt. */
                real mk = (real)(M-K);
                real effect = vx[jt]*((mk + 1.0)/mk - vytot/vy[jt]);
                fprintf(stderr, " %g@%d", effect, jt);
                if (effect > maxeffect)
                  { J[t] = J[K]; J[K] = jt; maxeffect = effect; }
              }
          }
        fprintf(stderr, "\n");
        return vy[J[K]];
      }
      
    /* Initialize the column permutation vector J[0..L]: */
    { int i, k;
      /* Mark indices which are NOT to be considered (including 0): */
      for (i = 0; i <= L; i++) { J[i] = 1; }
      for (i = 0; i < P; i++)
        { assert((S[i] > 0) && (S[i] <= L) && J[S[i]]); J[S[i]] = 0; }
      /* Put those indices in J[P..L]: */ 
      k = L;
      for (i = L; i >= 0; i--) { if (J[i]) { J[k] = i; k--; } }
      assert(k == P-1);
      /* Copy the given indices into J[0..P-1]: */
      for (i = 0; i < P; i++) { J[i] = S[i]; }
    }

    K = 0;
    while(K < M)
      { /* At this point we have already eliminated each of the noise
          symbols J[0..K-1] from X[0..N-1] and from Y[0..M-1]; i.e.
          X[r][J[i]] and Y[s][J[i]] are zero for i in 0..K-1, r in 0..N-1
          and s in 0..M-1; except that Y[i] has unit dependency on noise
          symbol J[i] (i.e. Y[i][J[i]] = 1). */
          
        int jK, r, s;

        /* Permute J[K..P-1] so that J[K] is the next victim, if any. */
        if (select_next() <= 0.0) 
          { /* There are no more usable equations among Y[K..M-1]. */
            return K;
          }
        jK = J[K];
        assert((jK > 0) && (jK <= L)); 
        fprintf(stderr, "eliminating noise %d ", jK);
        /* Recombine equations Y[K..M-1] so that Y[K,jK] != 0,
          and Y[s,jK] == 0 for s in K+1..M-1. */
        for (s = K+1; s < M; s++)
          { /* Rotate rows Y[K] and Y[s]: */
            real cos =  Y[K][jK];
            real sin = -Y[s][jK];
            if (s!=0.0)
              { real w = sqrt(cos*cos + sin*sin);
                int t;
                cos /= w; sin /= w;
                /* Rotate un-eliminated columns (including col. 0) */
                for (t = K; t <= L; t++)
                  { int jt = J[t]; 
                    real a = cos*Y[K][jt] - sin*Y[s][jt]; 
                    real b = sin*Y[K][jt] + cos*Y[s][jt]; 
                    Y[K][jt]=a; Y[s][jt]=b;
                  }
                Y[s][jK]=0.0;
                Y[K][NORM]=NONORM; Y[s][NORM]=NONORM;
              }
          }

        /* Normalize equation Y[K] so that Y[K][jK] = 1: */
        assert(Y[K][jK] != 0.0);
        { real scale = 1.0/Y[K][jK]; int t;
          for (t = K; t <= L; t++) { Y[K][J[t]] *= scale; }
          Y[K][jK]=1.0;
          Y[K][NORM]=NONORM;
        }

        /* Now use equation Y[K] == 0 to eliminate noise symbol jK
          from X[0..N-1] and Y{0..M-1] (except from Y[K] itself): */
        fprintf(stderr, "using equation\nY[%d] = ", K);
        aa_trace("", Y[K]);
        for (r = 0; r < N; r++)
          { real scale = X[r][jK]; int t;
            for (t = K+1; t <= L; t++) 
              { int jt = J[t]; X[r][jt] -= scale * Y[K][jt]; }
            X[r][jK] = 0.0;
            X[r][NORM] = NONORM;
          }
        for (s = 0; s < M; s++)
          { if (s != K)
              { real scale = Y[s][jK]; int t;
                for (t = K+1; t <= L; t++)
                  { int jt = J[t]; Y[s][jt] -= scale * Y[K][jt]; }
                Y[s][jK] = 0.0;
                Y[s][NORM] = NONORM;
              }
          }

        /* We got rid of one more noise symbol, phew: */
        K++;
      }

    return K;
  }
