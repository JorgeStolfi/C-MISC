#define PROG_NAME "gauss_tripart_opt"
#define PROG_DESC "computes optimum approximation of a Gaussian by narrower Gaussians"
#define PROG_VERS "1.0"

/* Last edited on 2024-12-21 11:56:19 by stolfi */
/* Created on 2007-07-11 by J. Stolfi, UNICAMP */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#include <affirm.h>
#include <bool.h>
#include <jsmath.h>
#include <interval.h>
#include <interval_io.h>
#include <jsfile.h>
#include <jsstring.h>
#include <argparser.h>
#include <gauss_distr.h>
#include <sve_minn.h>
#include <lsq.h>

#define PROG_HELP \
  "???"

#define PROG_INFO \
  "???"
  
/* !!! Rename {gauss_optimal_split} !!! */

/* PURPOSE

  A Gaussian PDF is a function {G_{a,d}(x) = g((x-a)/d)/d}
  where {a} is the mean (average), {d} is the deviation,
  and {g(z) = exp(-z^2/2)/sqrt(2*PI)}.  Its integral over
  all {\BR} is 1.
  
  Our goal is to find the best approximation of the standard Gaussian {H
  = G_{0,1}} by a sum 
  
    {F = SUM{ C[k]*G_{A[k],D[k]}: k \in -N .. +N}
    
  where {N >= 0} and {C[k],A[k],D[k]} are parameters to be determined.
  
  CONSTRAINTS
  
  We want the solution to be symmetric, that is, {C[-k] = C[k]},
  {D[-k] = D[k]}, and {A[-k] = -A[k]}. Note that this condition implies
  {A[0] = 0}. Finally we want the integral to still be 1, which is the
  same as saying {SUM{ C[k] : k \in -N .. +N} = 1}.
  
  We are interested in solutions where the {A[k]} and {D[k]} are constrained to
  user-specified ranges.  Specifically, all the deviations {D[k]} 
  be much less than 1, possibly not exceeding {1/2}.  The optimum is then 
  expected to have {A[k]} bounded, with spacing less than 1.
  
  Without loss of generality, we will consider only solutions where {0 <
  A[1] < .. A[N]}.  In fact, we require {A[j] - A[i] >= (j-i)*minSepA} for
  any {0 <= i < j <= N}, where {minSepA} is a user-specified parameter.  
  
  We may want to consider special cases where some or
  all of the parameters {A[1..N]} and/or {D[0..N]} and/or {C[0..N]} are
  fixed. The lower bound on {D[k]} must always be large enough to ensure
  proper sampling. Note that {A[k]} and {D[k]} are irrelevant, and
  cannot be optimized, if {C[k]} is fixed as zero. However, it makes no
  sense to set {C[k]} to zero for any {k > 0}.
  
  OPTIMALITY CRITERION
  
  The optimality criterion is to minimize the weighted squared
  difference integral 
  
    {S(A,D,C) = \int_{-\oo}^{+\oo} w(x)*(MF(x) - MH(x))^2\, dx } 
    
  where {w(x)} is a weight function whose integral is 1. 
  
  The notations
  {MF} and {MH} above denote either {F} and {H},
  or their cumulative integrals {\int_{-\oo}^{x} F(z)\, dz} and
  {\int_{-\oo}^{x} H(z)\, dz}.  The latter choice means that the
  distributions {F} and {H} are compared with the mean-squared 
  Kormogorov distance.
  
  !!! Implement the mean square Kolmogorov distance !!!
  
  The integral {S(A,D,C)} can be written {<F-H||F-H>} where {<||>} is
  the function inner product
  
    { <f||g> = \int_{-\oo}^{+\oo} w(x)*f(x)*g(x)\, dx }
      
  This integral in turn will be approximated by {<F-H|F-H>} where 
  {<|>} is the weighted average DISCRETE inner product
  
    { <f|g>= (1/W) * SUM{ wsmp[i] * f(xsmp[i]) * g(xsmp[i]) : i \in 0..M-1 } } 
    
  for certain sampling arguments {xsmp[0..M-1]} and certain weights
  {wsmp[0..M-1]} whose sum is {W}. 
  
  DETERMINING THE COEFFICIENTS
  
  Since the goal function {S(A,D,C)} is quadratic on {C[0..N]}, for any
  given {A[0..N]} and {D[0..N]} the optimum values of these coefficients
  can be found by simple least squares, with the affine constraint {C[0]
  = (1-SUM{ C[1..N]}/2} (plus any other constraints that may arise from
  the specified {C} ranges). Therefore the problem is basically a
  nonlinear optimization on the parameters {A[1..N]} and
  {D[0..N]} (excluding those that are fixed or irrelevant). */
  
#define nTerms_MAX (3)
  /* Max number of terms on each side of the central one. */
  
#define nSamples_MAX (10000)
  /* Max number of sampling args. */

#define D_MAX (0.800)
  /* Max value of {D[0..N]}. */

#define A_MIN (0.00001)
  /* Min value of {A[1..N]}. */

#define A_MAX (3.0000)
  /* Max value of {A[1..N]}. */

#define range_MIN (1.0e-8)
  /* Minimum width of a parameter range. */

typedef struct gto_options_t
  { char *outPrefix;                 /* Prefix for output file names. */
    int32_t nTerms;                  /* Number of {F} terms on each side of the central term. */
    bool_t noMiddle;                 /* If true, the central term is omitted. */
    double minD;                     /* A general lower bound for {D[0..N]}. */
    double maxA;                     /* A general upper bound for {A[0..N]}. */
    interval_t rangeA[nTerms_MAX+1]; /* {rangeA[0..N]} is the allowed value range {A[1..N]}. */
    interval_t rangeD[nTerms_MAX+1]; /* {rangeD[0..N]} is the allowed value range of {D[0..N]}. */
    double minSepA;                  /* Minimum separation between {A} parameters. */
    int32_t nSamples;                /* Number of sample arguments. */
  } gto_options_t;
  /* The command line options.  The ranges are defined only for indices
    in {0..nTerms}.  The interval {rangeA[0]} is always {[0 _ 0]}. */    

typedef struct gto_parms_t
  { int32_t N;              /* Number of terms on each side of the central term. */
    double A[nTerms_MAX+1]; /* {A[0..N]} are the means of the terms. */
    double D[nTerms_MAX+1]; /* {D[0..N]} are the deviations. */
    double C[nTerms_MAX+1]; /* {C[0..N]} are the coefficients. */
  } gto_parms_t;
  /* The parameters of the approximation. */
     
void gto_compute_samples_and_weights(int32_t M, double xmax, double xsmp[], double wsmp[]);
  /* Computes the sample arguments {xsmp[0..M-1]} and the
    corresponding inner product weights {wsmp[0..M-1]}.
    
    In the current implementation, The sampling arguments will be evenly
    spaced in {[-xMax _ +xmax]},inclusive the ends. All the weights will be
    equal to some positive constant {K}. */

void gto_compute_H_values(int32_t M, double xsmp[], double Hval[]);
  /* Computes {Hval[i] = H(xsmp[i])} for {i} in {0..M-1}. */
    
void gto_tweak_A_ranges(argparser_t *pp, int32_t N, interval_t rangeA[], double minSepA);
  /* Tweaks the {A} parameter ranges {rangeA[0..N]} so that the 
    ordering and separation constraints can be satisfied 
    by the mapping and unmapping formulas. */ 

void gto_find_optimum_parms
  ( gto_options_t *o, 
    int32_t M, 
    double xsmp[], 
    double wsmp[],
    double Hval[],
    gto_parms_t *p
  );
  /* Determines the parameters {p.A[1..N]} and {p.D[0..N]}
    and the coefficients {p.C[0..N]} of the optimum approximation {F},
    where {N = p.N}.  
    
    If any computed coefficient {p.C[k]} is zero, then the corresponding
    computed parameters {p.A[k]} and {p.D[k]} will be chosen arbitrarily.

    The parameters are restricted to the ranges specified in
    {o.rangeA[0..N]} and {o.rangeD[0..N]}. Moreover, the means
    {p.A[0..N]} should satisfy the minimum separation {minSepA}
    specified by the user; see {gto_check_constraints}.
    
    The optimality criterion is to minimize the quadratic mismatch
    {<F-H|F-H>}, where {<|>} is the inner product computed by 
    {gto_inner_cum_prod}, with arguments {M,wsmp[0..M-1]} and functions
    sampled at the points {xsmp[0..M=1]}. Assumes
    that {Hval[0..M-1]} are the sampled values of {H}. */
   
void gto_pick_initial_A_D_guesses(gto_options_t *o, gto_parms_t *p);
  /* Chooses initial guesses for the parameters {p->A[0..N]}
    and {p->D[0..N]}, where {N = o->nTerms}, consistent with the
    ranges {o->rangeA[0..N]} and {o->rangeD[0..N}. 
    
    Also sets {p->N} to {N}, and {p->C[0..N]} to {NAN}. */
    
void gto_check_constraints(gto_options_t *o, gto_parms_t *p);
  /* Check that the parameters {p.A[0..N]}, {p.D[0..N]}, and {p.C[0..N]}
    are in the respective ranges {o.rangeA[0..N]} and {o.rangeD[0..N]}; 
    where {N} is {p.N}. Aborts with error if these constraints are 
    not satisfied.
    
    Also checks the consistency of the parameter ranges with {o.minD}
    and {o.maxA}.  
    
    The coefficients {p.C[0..N]} may be {NAN}. If they are not {NAN},
    the sum {p.C[0] + 2 * SUM{ p.C[k] : k \in 1..N }} must be 1.
    
    Also checks the ordering and minimum separation constraints of the
    parameters {p.A[0..N]}. Specifically, {A[0]} must be fixed at 0, and
    the parameters {A[1..N]} must be increasing and satisfy {A[i] -
    A[i-1] >= minSepA}. 
    
    Moreover, if a range {o.rangeA[k]} is not 
    trivial, it must be such that {p.A[k]} still has a non-trivial
    actual range (accounting for minimum separation) 
    even if {p.A[k-1]} is maximum and {p.A[k+1]} is minimum.*/
    
int32_t gto_count_non_lin_parameters(int32_t N, interval_t rangeA[], interval_t rangeD[]);
   /* Returns the number of parameters {A[0..N]} and {D[0..N]} that are not fixed
     (that is, that have non-trivial ranges). */

void gto_compute_C_coeffs
  ( gto_parms_t *p, 
    bool_t noMiddle,
    int32_t M,
    double xsmp[],
    double wsmp[],
    double Hval[]
  );
  /* Computes the parameters {p.C[0..N]} of {F} that minimize the 
    mismatch {<F-H|F-H>} between {H} and {F}, for the current values
    of {p.A[0..N]} and {p.D[0..N]}; where {N = p.N}.
    
    The computed coefficients will satisfy the unit-sum constraint {C[0]
    + 2*SUM{C[1..N]} = 1}. Moreover, if {noMiddle} is {TRUE}, the
    coefficient {C[0]} is fixed at 0.
    
    The optimality criterion is to minimize the quadratic mismatch
    {<F-H|F-H>}, where {<|>} is the inner product computed by
    {gto_inner_cum_prod} with arguments {M,wsmp[0..M-1]} and functions sampled
    at the points {xsmp[0..M-1]}. Assumes that {Hval[0..M-1]} are the
    sampled values of {H}.  */

double gto_mean_sqr_mismatch
  ( gto_parms_t *p, 
    int32_t M, 
    double xsmp[], 
    double wsmp[], 
    double Hval[]
  );
  /* Computes the weighted mean quadratic mismatch between the cumulants {SH} and {SF}
    of {H} and {F}, namely {<F-H|F-H>}, where {<|>} is the inner product computed by 
    {gto_inner_cum_prod} with arguments {M,wsmp[0..M-1]} and functions sampled
    at the points {xsmp[0..M-1]}. Assumes that {Hval[0..M-1]} are the
    sampled values of {H}, and {p} are the parameters of {F}. */ 

void gto_compute_F_values
  ( int32_t N,
    double A[], 
    double D[], 
    double C[], 
    int32_t M, 
    double xsmp[], 
    double Fval[]
  );
  /* Computes the values {Fval[0..M-1]} of the approximant {F} at {xsmp[0..M-1]},
    for the given approximation parameters {A[0..N]} and {D[0..N]}
    and coefficients {C[0..N]}. Ignores {A[k]} and {D[k]} if {C[k]} is zero. */

void gto_print_parms(FILE *wr, int32_t ind, char *label, gto_parms_t *p, double goal);
  /* Writes the given parameters to {wr}, indented by {ind} columns.  
    Also writes the goal function {goal} if it is not {NAN}. */
  
void gto_eval_basis
  ( double x, 
    int32_t N,
    double A[], 
    double D[], 
    double Gx[]
  );
  /* Evaluates the basis distributions of {F} at {x}. Specifically, 
    sets {Gx[k] = G_{A[k],D[k]}(x)} for {k} in {0..N}.
    
    Namely, sets {G[0] = G_{A[0],D[0]}(x)
    } and {G[k] = (G_{-A[k],D[k]} + G_{+A[k],D[k]}} 
    for {k} in {1..N}. */

void gto_eval_distrs
  ( double x, 
    int32_t N,
    double A[], 
    double D[], 
    double distr[]
  );
  /* Evaluates the distributions {distr[0..2*N]} used in the approximant
    at the argument {x}.  Note that {distr} must  have {2*N+1} elements.
    
    Namely, sets {dist[N]} to G_{A[0],D[0]}(x), and {distr[N±k]}
    to G_{±A[k],D[k]} for {k} in {1..N}. */

double gto_inner_cum_prod
  ( int32_t M, 
    double fval[], 
    double gval[],
    double wsmp[]
  );
  /* Computes the cumulative inner product {<f|g>} of two functons, given their
    sampled values {fval[0..M-1]} and {gval[0..M-1]}.  
    
    The inner product is defined as the weighted average of {fsum[i]*gsum[i]},
    that is, {SUM{wsmp[i] *fsum[i] * gsum[i]} / SUM{wsmp[i]}}; where
    {fsum[i] = SUM{ fval[j] : j \in 0..i } and ditto for {gsum[i]}. */
    
void gto_plot_funcs
  ( char *outPrefix,
    gto_parms_t *p, 
    int32_t M, 
    double xsmp[], 
    double wsmp[], 
    double Hval[]
  );
  /* Writes to "{outPrefix}.txt" the sample arguments {xsmp[0..M-1]},
    one per line, each with the corresponding weights {wsmp[]}
    and the values of {H,Gm,Go,Gp,F} and {F-H}. */
    
gto_options_t *gto_parse_options(int argc, char **argv);
  /* Parses the command line optons. */

void gto_parse_range_options
  ( argparser_t *pp, 
    char *key, 
    int32_t lok, 
    int32_t hik, 
    interval_t range[], 
    double vmin, 
    double vmax
  );
  /* Parses zero or more of the command line options "{key} {k} {vlo} {vhi}" 
    for various {k} in {lok..hik}, and stores each interval {[vlo _ vhi]} into {range[k]}.
    Requires {vmin <= vlo <= vhi <= vmax} and fails if the same {k} 
    appears more than once.  Intervals {range[lok .. hik]} that are not specified 
    are set to {[vmin _ vmax]}. All other intervals in {range[0 .. nTerms_MAX]} 
    are set of {[NAN _ NAN]}. */

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
    for (uint32_t i = 0;  i < M; i++)
      { xsmp[i] = xmax*((2.0*i)/(M-1) - 1.0);
        wsmp[i] = 1.0;
        assert(wsmp[i] >= 0.0);
        sum_w += wsmp[i];
      }
    /* Normalize the weights to unit sum, just for plot's sake: */
    int32_t Nsig = -1; /* Number of significative samples. */
    for (uint32_t i = 0;  i < M; i++) 
      { wsmp[i] /= sum_w;
        if ((wsmp[i] >= 1.0e-15) && (i >= Nsig)) { Nsig = i+1; }
      }
    fprintf(stderr, "M = %d xmax = %12.8f Nsig = %d\n", M, xmax, Nsig);
  }

void gto_compute_H_values(int32_t M, double xsmp[], double Hval[])
  {
    for (uint32_t i = 0;  i < M; i++) { Hval[i] = gauss_distr_PDF(xsmp[i], 0, 1); }
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
    double *ctr = NULL;        /* Center of search domain. */
    double sve_d1Max = 0.500;  /* Size of domain to search. */
    bool_t sve_dBox = TRUE;    /* Domain is a box of side {2*d1Max}. */
    double sve_rIni = 0.200;   /* Radius of initial simplex. */
    double sve_rMin = 0.002;   /* Min simplex radius. */
    double sve_rMax = 0.200;   /* Max simplex radius. */
    double sve_minStep = 1.0e-10; /* Stop if displacement is less than this. */
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
        ctr, sve_d1Max, sve_dBox, sve_rIni, sve_rMin, sve_rMax,
        sve_minStep, sve_maxIters,
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
        for (uint32_t k = 0;  k <= pt->N; k++)
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
        for (uint32_t k = 0;  k <= pt->N; k++)
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
    for (uint32_t k = 0;  k <= N; k++)
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
    for (uint32_t k = 0;  k <= N; k++)
      { double Ck = p->C[k];
        if ((k == 0) && o->noMiddle) { assert((! isnan(Ck)) && (Ck == 0.0)); }
        Ctot += (k == 0 ? 1.0 : 2.0) * Ck; /* Note: may become {NAN}. */
      }
    /* Check unit sum, if not {NAN}: */
    assert(isnan(Ctot) || (fabs(Ctot - 1.0) <= 1.0e-14));
    
    /* Check the deviations {p.D[0..N]}: */
    for (uint32_t k = 0;  k <= N; k++)
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
    
    for (uint32_t k = 1;  k <= N; k++)
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
    for (uint32_t k = 0;  k <= N; k++)
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
    for (uint32_t k = 0;  k <= N; k++) { V[k] = (k == 0 ? 1.0 : 2.0 ); }
      
    /* Determiner the free (non-fixed, non-zero) coeffs {c[0..N0]}: */
    int32_t k0 = (noMiddle ? 1 : 0); /* Free coeff {c[r]} is {C[k0+r]}. */
    int32_t N0 = N - k0; /* The free coeffs are indexed {0..N0}. */
    assert(N0 >= 0);

    double Gxi[N+1];  /* Work vector: sampled vals of the original basis functions. */
    double *g[N0+1];  /* Sampled basis elements of free coeffs. */

    /* Evaluate the basis functions {G[0..N]} at the sampling args: */ 
    for (uint32_t r = 0;  r <= N0; r++) { g[r] = notnull(malloc(M*sizeof(double)), "no mem"); }
    for (uint32_t i = 0;  i < M; i++)
      { /* Evaluate the original basis functions {G[0..N]} at the sampling points: */
        gto_eval_basis(xsmp[i], N, p->A, p->D, Gxi);
        /* Store them in {g[0..N0]}: */
        for (uint32_t r = 0;  r <= N0; r++) { g[r][i] = Gxi[r+k0]; }
      }
      
    /* Build the least squares system: */
    int32_t N1 = N0+1; /* Number of free coeffs. */
    double *lsq_A = notnull(malloc(N1*N1*sizeof(double)), "no mem");
    double *lsq_B = notnull(malloc(N1*sizeof(double)), "no mem");
    for (uint32_t r = 0;  r < N1; r++) 
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
    for (uint32_t r = 0;  r <= N0; r++) { p->C[r+k0] = lsq_U[r]; }

    /* Set {C[0] = 0} if so requested: */
    if (noMiddle) { p->C[0] = 0.0; }

    /* Release storage: */
    for (uint32_t r = 0;  r <= N0; r++) { free(g[r]); }
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
    for (uint32_t i = 0;  i < M; i++) { Ferr[i] = Fval[i] - Hval[i]; }
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
    for (uint32_t i = 0;  i < M; i++)
      { gto_eval_basis(xsmp[i], N, A, D, Gx);
        double Fx = 0.0;
        for (uint32_t k = 0;  k <= N; k++) { Fx += C[k]*Gx[k]; }
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
    for (uint32_t k = 1;  k <= N; k++)
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
    for (uint32_t k = 1;  k <= N; k++)
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
    for (uint32_t k = 0;  k <= p->N; k++)
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
    double *Fval = notnull(malloc(M*sizeof(double)), "no mem"); /* Sampled values of {F}. */
    double distr[2*N+1]; /* Temp: values of the distributions. */
    gto_compute_F_values(N, p->A, p->D, p->C, M, xsmp, Fval); 
    for (uint32_t i = 0;  i < M; i++)
      { double xi = xsmp[i];
        double wi = wsmp[i];
        fprintf(wr, "%5d %+12.8f %16.12f", i, xi, wi);
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
    for (uint32_t k = 0;  k <= nTerms_MAX; k++) { range[k] = (interval_t){{ NAN, NAN }}; }
    
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
    for (uint32_t k = 1;  k <= N; k++)
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
    for (uint32_t k = 0;  k <= N; k++)
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
  for (uint32_t k = 1;  k <= N; k++)
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
