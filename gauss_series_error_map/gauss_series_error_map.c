#define PROG_NAME "gauss_series_error_map"
#define PROG_DESC "computes the error map for a regular Gaussian series"
#define PROG_VERS "1.0"

/* Last edited on 2024-11-08 18:52:11 by stolfi */
/* Created on 2007-07-11 by J. Stolfi, UNICAMP */

#define _GNU_SOURCE
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#include <affirm.h>
#include <bool.h>
#include <jsmath.h>
#include <jsfile.h>
#include <jsstring.h>
#include <argparser.h>
#include <gauss_distr.h>

#define PROG_HELP \
  "???"

#define PROG_INFO \
  "???"
  
/* !!! Rename {gauss_optimal_split} !!! */

/* PURPOSE

  We consider the problem of approximating functions with period 1
  by linear combinations of {N} Gaussians pulse trains with period
  1, where train {k} has a fixed deviation {d} and phase {k/N}. 
  
  APPROXIMATION SPACE
  
  Namely, the approximation space {\CA} consists of all functions {A}
  of the form
  
    { A(x) = SUM{ C[k] * B_{k,d}(x) : k \in 0..N-1 } }
    
  for any {x} in the period domain {D = [0_1]}, where
  the factors {C[0..N-1]} are arbitrary coefficients; and
  {B_{k,d}} is the infinite pulse train
  
    { B_{k,d}(x) = SUM{ G_{j + k/N,d}(x) : j \in \BZ } }
  
  where {G} is the Gaussian distribution with mean {a} and deviation
  {d}. That is, {G_{a,d}(x) = g((x-a)/d)/d} where {g} is the normal
  distribution, {g(z) = exp(-z^2/2)/sqrt(2*PI)}. Note that the integral
  of {B_{k,d}} in the period interval {D = [0 _ 1]} is 1.
  
  APPROXIMATION METHOD
  
  We use least squares as the approximation method.  Namely,
  for a given function {H} from some space {\CF} defined on {D}, we define the 
  approximation {H^a} to be the function in {\CA} that minimizes
  {<H - H^a|H - H^a>}, where {<|>} is an inner product on the
  space {\CA + \CF}.  It is known that {H^a} is the projection of {H}
  on the space {\CA} that is orthognal to {\CA} by that innner product;
  that is, satsfying {<H-H^a|A> = 0} for any {A\in \CA}.
  
  ERROR MAP
  
  We will use the error map technique of Gomide and Stolfi to evaluate 
  the quality of the least squares approximation at various points of the domain {D},
  for a whole space {\CF} of target functions to be approximated.
 
  Namely, for each {x} in {D}, we want to compute the root-mean-square value {E_{rms}(x)} of
  the error {F(x) - F^a(x)} as {F} ranges over {\CF_1}, the unit-energy subset of {\CF};
  where the energy of {F} is defined as {<<F|F>>} for some inner product {<<|>>} on {\CF}.  As shown
  by Gomide and Stolfi, {E_{rms}(x)} is related to the maximum error 
  
    { E_{max}(x) = \max\set{ |F(x) - F^a(x)| : F\in \CF_1 } }
    
  by a constant factor that does not depend on {x}.
  
  TARGET SPACE
  
  As the target space {\CF} we will use the space of Hartley series with period 1
  and {2*K+1} terms, where {2*K+1 <= N}.  Namely, a generic function {F} in {\CF} has the 
  form 
  
    { F(x) = SUM{ D[k]*H_{k}(x) : k \in -K..K }
  
  where {H_{k}} is a Hartley basis function of frequency {k}, 
  
    {H_{k}(x) = cos(2*PI*k*x + PI/4)}
  
  The norm {<<|>>} used to define the "energy" of a function {F} in {\CF}
  will be specified later.
  
  DISCRETIZATION
  
  We will work with discrete version of the functions, sampled
  at {M >> N} equally spaced points {xsmp[0..M-1]} over {D}.  
  That is, any function {F} defined over {D} will be represented
  by a column vector {F[0..M-1]} where {F[i] = F(xsmp[i])} for
  {i} in {0..M-1}.
  
  In particular, each basis element {B_{k,d}} of {\CA} can be viewed as column {k} of
  an {M\times N} matrix {B_{d}}, where {B_{d}[i,j] = B_{j,d}(xsmp[i]}}.
  
  Similarly, each Hartley basis element {H_{k}} is column {k} of an {M\times N} matrix {H}, 
  where {H[i,j] = H_{k}[xsmp[j]]}.
  
  The inner product {<|>} used in the approximation process is the 
  simple Cartesian inner product
  
    { <F|G> = F'*G = SUM{ F[i]*G[i] : i \in 0..M-1 }
    
  The inner product {<<|>>} used to define the energy of a target function
  will be a generic one, that can be expressed as {<<F|G>> = <L*F|L*G> = F'*(L'*L)*G} where 
  {L} is a non-singular {M\times M} matrix.  The matrix {L'*L} is symmetric and 
  positive definite; in fact, any symmetric positive definite marix can be
  factored in that form (Cholesky factorization), in many ways. 
  
  COMPUTING THE ERROR MAP
  
  To compute the error map, we need a basis {S} for {\CF} that is
  orthogonal under {<<|>>}.  If {<<|>>} is the ordinary inner product 
  {<|>}, we can use the sampled Hartley basis {H_{k}} for {k} in {-K..+K},
  as long as {M >= 2*K+1}.
  
  Then we need the projection matrix {R} of the approximation method,
  namely {B*(B'*B)^{-1}*B'} for the ordinary scalar product {<|>}. */
  
#define N_MAX (3)
  /* Max number of elements in the approximation basis {B}. */
  
#define M_MAX (10000)
  /* Max number of sampling args. */

#define D_MAX (0.800)
  /* Max deviation of the Gaussian pulses. */

typedef struct gto_options_t
  { char *outPrefix;                 /* Prefix for output file names. */
    int32_t nPulses;                 /* Dimension {N} of the approximation space. */
    int32_t maxFreq;                 /* Maximum absolute frequency {K} in the target space. */
    double dev;                      /* The deviation of each Gaussian pulse. */
    int32_t nSamples;                /* Number {M} of sample arguments. */
  } gto_options_t;
  /* The command line options.  */    

void gto_plot_funcs
  ( char *outPrefix,
    int32_t M, 
    double xsmp[], 
    double Hval[]
  );
  /* Writes to "{outPrefix}.txt" the sample arguments {xsmp[0..M-1]},
    one per line, each with the corresponding weights {wsmp[]}
    and the values of {H,Gm,Go,Gp,F} and {F-H}. */
    
gto_options_t *gto_parse_options(int argc, char **argv);
  /* Parses the command line optons. */

int main(int argc, char **argv);

/* IMPLEMENTATIONS */

int main(int argc, char **argv)
  {
    gto_options_t *o = gto_parse_options(argc, argv);

    /* Choose the sampling arguments and weights: */
    int32_t M = o->nSamples;
    double xsmp[M];  /* Sample arguments for inner products. */
    double wsmp[M];  /* Weights of those samples. */
    double xmax = o->maxA + gauss_distr_BIG_ARG; /* There should be nothing beyond this. */
    gto_compute_samples_and_weights(M, xmax, xsmp, wsmp);

    /* Sample the target function {H}: */
    double Hval[M];  /* Values of {H(xsmp[0..M-1])}. */
    gto_compute_H_values(M, xsmp, Hval);

    /* Find the optimium parameters: */
    gto_parms_t p;
    gto_find_optimum_parms(o, M, xsmp, wsmp, Hval, &p);
    
    /* Print the optimum parameters and mismatch: */
    double goal = gto_mean_sqr_mismatch(&p, M, xsmp, wsmp, Hval);
    gto_print_parms(stderr, 0, "optimum", &p, goal);
    { FILE *wr = open_write(txtcat(o->outPrefix, ".parms"), TRUE);
      gto_print_parms(wr, 0, "optimum", &p, goal);
      fclose(wr);
    }
    
    /* Write the plot data file: */
    gto_plot_funcs(o->outPrefix, &p, M, xsmp, wsmp, Hval);
    return 0;
  }
  
void gto_compute_samples_and_weights(int32_t M, double xmax, double xsmp[], double wsmp[])
  {
    double sum_w = 0.0;
    for (int32_t i = 0; i < M; i++)
      { xsmp[i] = (i*xmax)/(M-1);
        wsmp[i] = (i == 0 ? 1 : 2) * 1.0;
        assert(wsmp[i] >= 0.0);
        sum_w += wsmp[i];
      }
    /* Normalize the weights to unit sum, just for plot's sake: */
    int32_t Nsig = -1; /* Number of significative samples. */
    for (int32_t i = 0; i < M; i++) 
      { wsmp[i] /= sum_w;
        if ((wsmp[i] >= 1.0e-15) && (i >= Nsig)) { Nsig = i+1; }
      }
    fprintf(stderr, "M = %d xmax = %12.8f Nsig = %d\n", M, xmax, Nsig);
  }

void gto_compute_H_values(int32_t M, double xsmp[], double Hval[])
  {
    for (int32_t i = 0; i < M; i++) { Hval[i] = gauss_distr_PDF(xsmp[i], 0, 1); }
  }

void gto_find_optimum_parms
  ( gto_options_t *o, 
    int32_t M, 
    double xsmp[], 
    double wsmp[],
    double Hval[],
    gto_parms_t *p
  )
  {
    int32_t N = o->nTerms;
    bool_t debug_map = FALSE;
   
    /* Uses the Simplex Vertex-Edge (SVE) method for non-linear
      optimizaton to compute the parameters {A[1..N]} and {D[0..N]}.
      
      The SVE goal function is the weighted quadratic mismatch, for the
      optimum coefficients {C[0..N]} computed for the non-linear
      parameters by quadratic (least squares) optimization.

      The optimization parameters {sve_x[0..sve_n-1]} are
      {A[0],D[0],A[1],D[1],...,A[N],D[N]}, in that order, scaled and
      shifted from their user-specified ranges to {[0 _ 1]}; but
      skipping those parameters (such as {A[0]}) that have a trivial
      (single-value) range. */
    
    int32_t sve_n = gto_count_non_lin_parameters(N, o->rangeA, o->rangeD); 
    fprintf(stderr, "found %d non-linear  parameters to optimize\n", sve_n);
    double sve_x[sve_n];  /* The non-linear parameters, scaled and shifted. */
    
    auto void map_parms(gto_parms_t *pt, double xt[]);
    auto void unmap_parms(double xt[], gto_parms_t *pt);
      /* Map the approximation parameters {pt.D[0..N],pt.A[0..N]} to the 
        SVE minimization parameters {xt[0..sve_n-1]}, and vice-versa. */
      
    sign_t sve_dir = -1; /* Search for minimum.*/
    double *sve_ctr = NULL;    /* Center of search domain. */
    double sve_d1Max = 0.500;  /* Size of domain to search. */
    bool_t sve_dBox = TRUE;    /* Domain is a box of side {2*d1Max}. */
    double sve_rIni = 0.200;   /* Radius of initial simplex. */
    double sve_rMin = 0.002;   /* Min simplex radius. */
    double sve_rMax = 0.200;   /* Max simplex radius. */
    double sve_stop = 1.0e-10; /* Stop if displacement is less than this. */
    int32_t sve_maxIters = 50; /* Max ierations of the main loop. */
    bool_t sve_debug = FALSE;
    bool_t sve_debug_probes = FALSE;
    
    auto double sve_F(int ns, double xs[]);
      /* The goal function for {sve_minn_iterate}. */
      
    /* Initial guesses for the non-linear parameters: */
    gto_pick_initial_A_D_guesses(o, p);
    gto_compute_C_coeffs(p, o->noMiddle, M, xsmp, wsmp, Hval);
    map_parms(p, sve_x);
    double sve_Fx = sve_F(sve_n, sve_x); /* Value of the goal function for those parameters. */
    gto_print_parms(stderr, 2, "SVE - initial guess", p, sve_Fx);
   
    /* Call the optimizer: */
    sve_minn_iterate
      ( sve_n, sve_F, NULL, NULL, sve_x, &sve_Fx, sve_dir, 
        sve_ctr, sve_d1Max, sve_dBox, sve_rIni, sve_rMin, sve_rMax,
        sve_stop, sve_maxIters,
        sve_debug, sve_debug_probes
      );
    
    /* Final unmapping of parameters: */
    unmap_parms(sve_x, p);
    gto_compute_C_coeffs(p, o->noMiddle, M, xsmp, wsmp, Hval);
    gto_print_parms(stderr, 2, "SVE - result", p, sve_Fx);
      
    return;
    
    /* SVE PARAMETER MAPPING
    
      The SVE parameters {sve_x[0..sve_n-1]} range in {[0 _ 1]}. The
      mapping to {p.A[1..N]} and {p.D[0..N]} is such that the means
      {p.A[0..N]} always satisfy the range, order, and minimum
      separation constraints. Namely, if {x[i]} corresponds to the
      non-fixed parameter {p.A[k]}, then the range {[0 _ 1]} of the
      former is mapped to {[loAk _ hiAk]} where {loAk = p.A[k-1]+o.minSepA}
      and {hiAk = o.rangeA[k].end[1]}. The ranges
      {o.rangeA[0..N]} must be such that the second interval is always
      non-negligible. */
    
    /* Internal implementations: */
    void map_parms(gto_parms_t *pt, double xt[])
      { 
        int32_t kt = 0; /* Next SVE parameter to define. */
        if (debug_map) { fprintf(stderr, "map (N=%d): ", pt->N); }
        for (int32_t k = 0; k <= pt->N; k++)
          { /* Map the mean {A[k]} if not fixed: */
            interval_t *rAk = &(o->rangeA[k]);
            double loA = rAk->end[0];
            double hiA = rAk->end[1];
            if (loA < hiA)
              { assert(hiA - loA >= range_MIN);
                assert(kt < sve_n);
                xt[kt] = (pt->A[k] - loA)/(hiA - loA);
                if (debug_map) { fprintf(stderr, "(x[%d]=A[%d])", kt, k); }
                kt++;
              }
            else
              { assert(pt->A[k] == loA); }
            /* Map the deviation {D[k]} if not fixed: */
            interval_t *rDk = &(o->rangeD[k]);
            double loD = rDk->end[0];
            double hiD = rDk->end[1];
            if (loD < hiD)
              { assert(hiD - loD >= range_MIN);
                assert(kt < sve_n);
                xt[kt] = (pt->D[k] - loD)/(hiD - loD);
                if (debug_map) { fprintf(stderr, "(x[%d]=D[%d])", kt, k); }
                kt++;
              }
            else
              { assert(pt->D[k] == loD); }
          }
        if (debug_map) { fprintf(stderr, " kt=%d\n", kt); }
        assert(kt == sve_n);
      }
        
    void unmap_parms(double xt[], gto_parms_t *pt)
      { int32_t kt = 0; /* Next SVE parameter to unmap. */
        if (debug_map) { fprintf(stderr, "unmap (N=%d): ", pt->N); }
        for (int32_t k = 0; k <= pt->N; k++)
          { /* Map the mean {A[k]} if not fixed: */
            interval_t *rAk = &(o->rangeA[k]);
            double loA = rAk->end[0];
            double hiA = rAk->end[1];
            if (loA < hiA)
              { assert(kt < sve_n);
                pt->A[k] = loA + xt[kt]*(hiA - loA);
                if (debug_map) { fprintf(stderr, "(A[%d]=x[%d])", k, kt); }
                kt++;
              }
            else
              { pt->A[k] = loA; }
            /* Map the deviation {D[k]} if not fixed: */
            interval_t *rDk = &(o->rangeD[k]);
            double loD = rDk->end[0];
            double hiD = rDk->end[1];
            if (loD < hiD)
              { assert(kt < sve_n);
                pt->D[k] = loD + xt[kt]*(hiD - loD);
                if (debug_map) { fprintf(stderr, "(D[%d]=x[%d])", k, kt); }
                kt++;
              }
            else
              { pt->D[k] = loD; }
          }
        if (debug_map) { fprintf(stderr, " kt=%d\n", kt); }
        assert(kt == sve_n);
      }
    
    auto double sve_F(int ns, double xs[])
      { assert (ns == sve_n);
        unmap_parms(xs, p);
        double res = gto_mean_sqr_mismatch(p, M, xsmp, wsmp, Hval);
        gto_print_parms(stderr, 4, "sve_F", p, res);
        return res;
      }
  }

void gto_pick_initial_A_D_guesses(gto_options_t *o, gto_parms_t *p)
  { int32_t N = o->nTerms;
    p->N = N;
    for (int32_t k = 0; k <= N; k++)
      { /* Pick {A[k],D[k]} in the middle of their ranges: */
        p->D[k] = interval_mid(&(o->rangeD[k]));
        p->A[k] = interval_mid(&(o->rangeA[k]));
        p->C[k] = ((k == 0) && o->noMiddle ? 0.0 : NAN);
      }
    /* Paranoia*/
    gto_check_constraints(o, p);
  }

void gto_check_constraints(gto_options_t *o, gto_parms_t *p)
  { 
    int32_t N = o->nTerms;
    assert((N >= (o->noMiddle ? 1 : 0)) && (N <= nTerms_MAX));
    assert(p->N == N);
    
    /* Check the coefficients {p.C[0..N]}: */
    double Ctot = 0.0;
    for (int32_t k = 0; k <= N; k++)
      { double Ck = p->C[k];
        if ((k == 0) && o->noMiddle) { assert((! isnan(Ck)) && (Ck == 0.0)); }
        Ctot += (k == 0 ? 1.0 : 2.0) * Ck; /* Note: may become {NAN}. */
      }
    /* Check unit sum, if not {NAN}: */
    assert(isnan(Ctot) || (fabs(Ctot - 1.0) <= 1.0e-14));
    
    /* Check the deviations {p.D[0..N]}: */
    for (int32_t k = 0; k <= N; k++)
      { double Dk = p->D[k];
        assert(! isnan(Dk));
        interval_t rDk = o->rangeD[k];
        assert(rDk.end[0] >= o->minD);
        /* Check that {Dk} is in range: */
        assert((rDk.end[0] <= Dk) && (Dk <= rDk.end[1]));
      }
      
    /* Check the means {p.A[0..N]}: */
    assert(p->A[0] == 0.0); /* The mean {p.A[0]} must always be zero ... */
    assert((o->rangeA[0].end[0] == 0) && (o->rangeA[0].end[1] == 0)); /* ... and fixed there. */
    
    for (int32_t k = 1; k <= N; k++)
      { double Ak = p->A[k];
        assert(! isnan(Ak));
        /* Check the minimum separation: */
        assert(Ak > p->A[k-1] + o->minSepA);
        
        interval_t rAk = o->rangeA[k];
        assert(rAk.end[0] >= k*o->minSepA);
        assert(rAk.end[1] <= o->maxA - (N-k)*o->minSepA);
        
        /* Check that {Ak} is in range: */
        assert((rAk.end[0] <= Ak) && (Ak <= rAk.end[1]));
      }
  }

int32_t gto_count_non_lin_parameters(int32_t N, interval_t rangeA[], interval_t rangeD[])
  { int32_t sve_n = 0;
    for (int32_t k = 0; k <= N; k++)
      { interval_t *rAk = &(rangeA[k]);
        double loA = rAk->end[0];
        double hiA = rAk->end[1];
        if (loA < hiA)
          { assert(hiA - loA >= range_MIN);
            sve_n++;
          }
        interval_t *rDk = &(rangeD[k]);
        double loD = rDk->end[0];
        double hiD = rDk->end[1];
        if (loD < hiD)
          { assert(hiD - loD >= range_MIN);
            sve_n++;
          }
      }
    return sve_n;
  }

void gto_compute_C_coeffs
  ( gto_parms_t *p, 
    bool_t noMiddle,
    int32_t M,
    double xsmp[],
    double wsmp[],
    double Hval[]
  )
  { 
    int32_t N = p->N;
    assert(N >= (noMiddle ? 1 : 0));
    
    /* 
      Let {k0} be zero if {o->noMiddle} is false, and 1 if {o.noMiddle} is true.
      
      MISMATCH AS FUNCTION OF FREE COEFFICIENTS
    
      For given parameters {A[0..N]} and {D[0..N]}, the approximant {F}
      can be written as {F = SUM{ C[k]*G[k] : k \in k0..N} } where {G[0]
      = G_{A[0],D[0]}} and {G[k] = (G_{-A[k],D[k]} + G_{+A[k],D[k]}} for
      {k} in {1..N}.
      
      Let's rename {C[k0..N]} as {c[0..N0]}, where {N0 = N-k0}; that is,
      {c[r]} is the same as {C[r+k0]}, for all {r} in {0..N0}. We can
      then rewrite the approximant {F} as
      
        { F = SUM{ c[r]*g[r] : r \in 0..N0 } }
             
      where {g[r] = G[r+k0]}. The mismatch function {<F-H|F-H>} then too is a quadratic function 
      of the "free" coefficients {c[0..N0]} only, namely
      
        { Q(c[0..N0]) = 
            + <H|H> 
            - 2*SUM{ c[r]<g[r]|H> : r \in 0..N0  }
            + SUM{ c[s]*c[t]*<g[s]|g[t]> : s,t \in 0..N0 }
        }
        
      THE UNIT SUM CONSTRAINT
      
      Since the integral must be 1, the coefficients must satisfy 
      {SUM{V[k]*C[k] : k \in 0..N } = 1} where {V[0] = 1} and {V[1..N] = 2}.
      
      Since {C[0]} is zero if {k0 = 1}, we can rewrite this constraint 
      as {SUM{v[r]*r] : r \in 0..N0 } = 1} where {v[r] = V[r+k0]}.
      
      MINIMIZATION PROBLEM
      
      We then have the constrained least squares problem with {N1 = N0+1}
      coeeficients {c[0..N0]}, with moment matrix {A} ({N1×N1}) such that {A[i,j]} = 
      <g[i]|g[j]>} and right-hand side {B} ({N1×1}) such that 
      {B[i] = <g[i]|H>}; with the unit-sum constraint above.
      
      We then can solve the system to get {C[k0..N] = c[0..N0]},
      and set {C[0] = 0} if {noMiddle} is true. */
    
    /* Define the weights {V[0..N]} of the unit-sum constraint: */
    double V[N+1]; /* {V[k]} is the weight of {C[k]} in the unit-sum constraint. */
    for (int32_t k = 0; k <= N; k++) { V[k] = (k == 0 ? 1.0 : 2.0 ); }
      
    /* Determiner the free (non-fixed, non-zero) coeffs {c[0..N0]}: */
    int32_t k0 = (noMiddle ? 1 : 0); /* Free coeff {c[r]} is {C[k0+r]}. */
    int32_t N0 = N - k0; /* The free coeffs are indexed {0..N0}. */
    assert(N0 >= 0);

    double Gxi[N+1];  /* Work vector: sampled vals of the original basis functions. */
    double *g[N0+1];  /* Sampled basis elements of free coeffs. */

    /* Evaluate the basis functions {G[0..N]} at the sampling args: */ 
    for (int32_t r = 0; r <= N0; r++) { g[r] = notnull(malloc(M*sizeof(double)), "no mem"); }
    for (int32_t i = 0; i < M; i++)
      { /* Evaluate the original basis functions {G[0..N]} at the sampling points: */
        gto_eval_basis(xsmp[i], N, p->A, p->D, Gxi);
        /* Store them in {g[0..N0]}: */
        for (int32_t r = 0; r <= N0; r++) { g[r][i] = Gxi[r+k0]; }
      }
      
    /* Build the least squares system: */
    int32_t N1 = N0+1; /* Number of free coeffs. */
    double *lsq_A = notnull(malloc(N1*N1*sizeof(double)), "no mem");
    double *lsq_B = notnull(malloc(N1*sizeof(double)), "no mem");
    for (int32_t r = 0; r < N1; r++) 
      { for (int32_t t = r; t < N1; t++)
          { double Art = gto_inner_cum_prod(M, g[r], g[t], wsmp);
            lsq_A[r*N1 + t] = Art;
            if (r != t) { lsq_A[t*N1+ r] = Art; }
          }
        lsq_B[r] = gto_inner_cum_prod(M, g[r], Hval, wsmp);
      }

    /* Solve the system to obtain the free coeffs {c[0..N0]}: */
    double *lsq_U = notnull(malloc(N1*sizeof(double)), "no mem");
    double One = 1.0; /* RHS of the unit-sum constraint. */
    double *v = &(V[k0]); /* Coeffs of {c[0..R-1]} in the unit-sum constraint. */
    lsq_solve_system(N1, 1, lsq_A, lsq_B, 1,v, &One, lsq_U, NULL, FALSE);

    /* Copy {c[0..N0]} to {C[k0..N-1]}: */
    for (int32_t r = 0; r <= N0; r++) { p->C[r+k0] = lsq_U[r]; }

    /* Set {C[0] = 0} if so requested: */
    if (noMiddle) { p->C[0] = 0.0; }

    /* Release storage: */
    for (int32_t r = 0; r <= N0; r++) { free(g[r]); }
    free(lsq_A); free(lsq_B); free(lsq_U);
  }

double gto_mean_sqr_mismatch
  ( gto_parms_t *p, 
    int32_t M, 
    double xsmp[], 
    double wsmp[], 
    double Hval[]
  )
  { 
    /* Evaluate the approximation {F}: */
    double Fval[M]; /* Values of {F} at sample arguments. */
    gto_compute_F_values(p->N, p->A, p->D, p->C, M, xsmp, Fval);
    /* Compute the errors {F-H}: */
    double Ferr[M]; /* Values of {F-H} at sample arguments. */
    for (int32_t i = 0; i < M; i++) { Ferr[i] = Fval[i] - Hval[i]; }
    /* Average the squared error: */
    return gto_inner_cum_prod(M, Ferr, Ferr, wsmp);
  }
  
void gto_compute_F_values
  ( int32_t N,
    double A[], 
    double D[], 
    double C[], 
    int32_t M, 
    double xsmp[], 
    double Fval[]
  )
  {  
    double Gx[N+1];
    for (int32_t i = 0; i < M; i++)
      { gto_eval_basis(xsmp[i], N, A, D, Gx);
        double Fx = 0.0;
        for (int32_t k = 0; k <= N; k++) { Fx += C[k]*Gx[k]; }
        Fval[i] = Fx;
      }
  }

void gto_eval_basis
  ( double x, 
    int32_t N,
    double A[], 
    double D[], 
    double Gx[]
  )
  { double distr[2*N+1];
    gto_eval_distrs(x, N, A, D, distr);
    Gx[0] = distr[N];
    for (int32_t k = 1; k <= N; k++)
      { Gx[k] = distr[N-k] + distr[N+k]; }
  }

void gto_eval_distrs
  ( double x, 
    int32_t N,
    double A[], 
    double D[], 
    double distr[]
  )
  { assert (A[0] == 0.0);
    distr[N] = gauss_distr_PDF(x, 0.0, D[0]);
    for (int32_t k = 1; k <= N; k++)
      { distr[N-k] = gauss_distr_PDF(x, -A[k], D[k]);
        distr[N+k] = gauss_distr_PDF(x, +A[k], D[k]);
      }
  }

double gto_inner_cum_prod
  ( int32_t M, 
    double fval[], 
    double gval[],
    double wsmp[]
  )
  { double fsum = 0; /* Cumulant of {fval}. */
    double gsum = 0; /* Cumulant of {gval}. */
    double sum_fgw = 0.0; /* Sum of {fsum[i]*gsum[i]*wsmp[i]}. */
    double sum_w = 0.0;   /* Sum of {wsmp[i]}. */
    for (int32_t i = M-1; i >= 0; i--)
      { fsum += fval[i];
        gsum += gval[i];
        sum_fgw += wsmp[i]*fsum*gsum;
        sum_w += wsmp[i];
      }
    return sum_fgw/sum_w;
  }

void gto_print_parms(FILE *wr, int32_t ind, char *label, gto_parms_t *p, double goal)
  { fprintf(stderr, "%*s--- %s ------------------------------\n", ind, "", label);
    for (int32_t k = 0; k <= p->N; k++)
      { fprintf(wr, "%*s  C[%d] = %+20.16f", ind, "", k, p->C[k]);
        fprintf(wr, " A[%d] = %23.16f D[%d] = %23.16f\n", k, p->A[k], k, p->D[k]);
      }
    if (! isnan(goal)) { fprintf(wr, "%*s  goal = %20.16f\n", ind, "", goal); }
    char *dashes = txtrep("-", (uint32_t)(strlen(label) + 2));
    fprintf(stderr, "%*s---%s------------------------------\n", ind, "", dashes);
    free(dashes);
  }

void gto_plot_funcs
  ( char *outPrefix,
    gto_parms_t *p, 
    int32_t M, 
    double xsmp[], 
    double wsmp[], 
    double Hval[]
  )
  { 
    FILE *wr = open_write(txtcat(outPrefix, ".txt"), TRUE);
    int32_t N = p->N;
    /* Write negative samples too.  Note that all weights are folded, but one. */
    double *Fval = notnull(malloc(M*sizeof(double)), "no mem"); /* Sampled values of {F}. */
    double distr[2*N+1]; /* Temp: values of the distributions. */
    gto_compute_F_values(N, p->A, p->D, p->C, M, xsmp, Fval); 
    for (int32_t ii = -(M-1); ii <= (M-1); ii++)
      { int32_t i = (ii < 0 ? -ii : ii);
        double xi = (ii < 0 ? -1 : +1)*xsmp[i];
        double wi = (i == 0 ? 1.0 : 0.5)*wsmp[i]; /* Unfolded weight. */
        fprintf(wr, "%5d %+12.8f %16.12f", ii, xi, wi);
        fprintf(wr, "  %16.12f", Hval[i]);
        /* Compute the basis functions: */
        gto_eval_distrs(xi, N, p->A, p->D, distr);
        for (int32_t kk = -N; kk <= N; kk++)
           { int32_t k = (kk < 0 ? -kk : kk);
             fprintf(wr, "  %16.12f", p->C[k]*distr[N+kk]);
           }
        double Ferri = Fval[i] - Hval[i];
        fprintf(wr, "  %16.12f  %16.12f", Fval[i], Ferri);
        fprintf(wr, "\n");
      }
    free(Fval);
    fclose(wr);
  }
        
gto_options_t *gto_parse_options(int argc, char **argv)
  {
    /* Initialize argument parser: */
    argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, PROG_NAME " version " PROG_VERS ", usage:\n" PROG_HELP);
    argparser_set_info(pp, PROG_INFO);
    argparser_process_help_info_options(pp);

    /* Allocate the command line argument record: */
    gto_options_t *o = notnull(malloc(sizeof(gto_options_t)), "no mem");

    /* Parse keyword parameters: */
    argparser_get_keyword(pp, "-outPrefix");
    o->outPrefix = argparser_get_next_non_keyword(pp);
    
    argparser_get_keyword(pp, "-noMiddle");
    o->noMiddle = argparser_get_next_bool(pp);
    
    argparser_get_keyword(pp, "-nTerms");
    o->nTerms = (int32_t)argparser_get_next_int(pp, (o->noMiddle ? 1 : 0), nTerms_MAX);
    int32_t N = o->nTerms;
    
    argparser_get_keyword(pp, "-nSamples");
    o->nSamples = (int32_t)argparser_get_next_int(pp, 1, nSamples_MAX);
    
    argparser_get_keyword(pp, "-minD");
    double D_MIN = 20.0/(o->nSamples); /* Arbitrary limit for accurate sampling. */
    o->minD = argparser_get_next_double(pp, D_MIN, D_MAX);
    
    argparser_get_keyword(pp, "-maxA");
    o->maxA = argparser_get_next_double(pp, A_MIN, A_MAX);
    
    argparser_get_keyword(pp, "-minSepA");
    o->minSepA = argparser_get_next_double(pp, A_MIN, A_MAX/N);
    
    /* Parse the {A} parameter ranges: */
    gto_parse_range_options(pp, "-rangeA", 1, N, o->rangeA, A_MIN, A_MAX);
    o->rangeA[0] = (interval_t){{ 0, 0 }}; /* Term 0 is always at 0. */    
    /* Ensure order and separation constraints on {A[0..N]} can be met. */
    gto_tweak_A_ranges(pp, N, o->rangeA, o->minSepA); 

    /* Parse the {D} parameter ranges: */
    gto_parse_range_options(pp, "-rangeD", 0, N, o->rangeD, o->minD, D_MAX);
    
    /* Check for spurious arguments: */
    argparser_finish(pp);

    return o;
  }

void gto_parse_range_options
  ( argparser_t *pp, 
    char *key, 
    int32_t lok, 
    int32_t hik, 
    interval_t range[], 
    double vmin, 
    double vmax
  )
  { 
    assert(vmax - vmin >= 2.0*range_MIN);
    
    /* Set all to {NAN} to detect duplicates: */
    for (int32_t k = 0; k <= nTerms_MAX; k++) { range[k] = (interval_t){{ NAN, NAN }}; }
    
    /* Parse the {key} ranges: */
    while (argparser_keyword_present(pp, key))
      { int32_t k = (int32_t)argparser_get_next_int(pp, lok, hik);
        if (! isnan(range[k].end[0])) { argparser_error(pp, txtcat3("repeated \"",key,"\" spec")); }
        double vlo = argparser_get_next_double(pp, vmin, vmax);
        double vhi = argparser_get_next_double(pp, vlo, vmax);
        range[k] = (interval_t){{ vlo, vhi }};
      }
      
    /* Set unspecified ranges to their defaults: */
    for (int32_t k = lok; k <= hik; k++) 
      { if (isnan(range[k].end[0])) { range[k] = (interval_t){{ vmin, vmax }}; } }
  } 

void gto_tweak_A_ranges(argparser_t *pp, int32_t N, interval_t rangeA[], double minSepA)
  {
    /* The mean {A[0]} must be fixed at 0: */
    assert((rangeA[0].end[0] == 0) && (rangeA[0].end[1] == 0));
    
    interval_t newA[nTerms_MAX + 1]; /* Adjusted parameter ranges. */
    /* Ensure lower bounds are properly spaced: */
    newA[0].end[0] = 0.0; 
    for (int32_t k = 1; k <= N; k++)
      { double vlo = newA[k-1].end[0] + minSepA;
        newA[k].end[0] = fmax(vlo, rangeA[k].end[0]);
      }
    /* Ensure upper bounds are properly spaced: */
    newA[N].end[1] = rangeA[N].end[1];
    for (int32_t k = N-1; k >= 0; k--)
      { double vhi = newA[k+1].end[1] - minSepA;
        newA[k].end[1] = fmin(vhi, rangeA[k].end[1]);
      }
      
    /* Check changes and validity, and return result: */
    for (int32_t k = 0; k <= N; k++)
      { interval_t rAk = rangeA[k];
        assert(rAk.end[0] <= rAk.end[1]);
        bool_t fixed = (rAk.end[0] == rAk.end[0]); /* User specified that parameter is fixed. */
        interval_t nAk = newA[k];
        assert(nAk.end[0] >= rAk.end[0]);
        assert(nAk.end[1] <= rAk.end[1]);
        if ((nAk.end[0] != rAk.end[0]) || (nAk.end[0] != rAk.end[0]))
          { fprintf(stderr, "!! warning: range for A[%d] reduced from ", k);
            interval_gen_print(stderr, &(rAk), "%20.16f", "[", " _ ", "]");
            fprintf(stderr, " to ");
            interval_gen_print(stderr, &(newA[k]), "%20.16f", "[", " _ ", "]");
            fprintf(stderr, "\n");
            rangeA[k] = nAk; rAk = nAk;
          }
        assert(rAk.end[0] <= rAk.end[1]);
        
        if (rAk.end[0] > rAk.end[0])
          { /* Range was or became empty: */
            fprintf(stderr, "** error: empty range for A[%d] = ", k);
            interval_gen_print(stderr, &(rAk), "%20.16f", "[", " _ ", "]");
            fprintf(stderr, "\n");
            exit(1);
          }
        else if (fixed)
          { /* It was fixed, so it should have remained fixed: */
            assert(rAk.end[0] == rAk.end[0]);
          }
        else if ((rAk.end[1] - rAk.end[0]) < range_MIN)
          { /* It was not fixed, but become fixed or too small: */
            fprintf(stderr, "** error: range too small for A[%d] = ", k);
            interval_gen_print(stderr, &(rAk), "%20.16f", "[", " _ ", "]");
            fprintf(stderr, "\n");
            exit(1);
          }
      }
  
  /* Check that non-trivial ranges remain significant even if {A[k-1]} is maximum: */
  assert((rangeA[0].end[0] == 0) && (rangeA[0].end[1] == 0));
  for (int32_t k = 1; k <= N; k++)
    { interval_t rAk = rangeA[k];
      assert(rAk.end[0] <= rAk.end[1]);
      if (rAk.end[0] < rAk.end[1])
        { /* {A[k]} is not fixed. Check that it remains so, no matter what is {A[k-1]}: */
          double lo_sep = rangeA[k-1].end[1] + minSepA;
          double hi_sep = rAk.end[1];
          if (hi_sep - lo_sep < 1.0e-8)
            { fprintf(stderr, "** error: range for A[%d] may become empty or too small = ", k);
              interval_t sAk = (interval_t){{ lo_sep, hi_sep }};
              interval_gen_print(stderr, &(sAk), "%20.16f", "[", " _ ", "]");
              fprintf(stderr, "\n");
              exit(1);
            }
        }
      }
  }
