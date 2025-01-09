#define PROG_NAME "compute_gauss_resampling_weights"
#define PROG_DESC "computes optimum weights for Gaussian resampling of a Gaussian-sampled signal"
#define PROG_VERS "1.0"

/* Last edited on 2024-12-01 00:26:48 by stolfi */
/* Created on 2007-07-11 by J. Stolfi, UNICAMP */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <complex.h>

#include <affirm.h>
#include <bool.h>
#include <jsmath.h>
#include <wt_table.h>
#include <wt_table_binomial.h>
#include <rn.h>
#include <rmxn.h>
#include <gausol_solve.h>

/* GAUSSIAN KERNEL

  A /Gaussian kernel/ is any Gaussian distribution with zero mean, 
  
  | { g_\sigma(z) = exp(-(z/\sigma)^2/2)/(\sigma \sqrt{2\pi})}.
  
  The parameter {\sigma} is the standard deviation of the 
  distribution, and is the distance of its inflection points.
  We can also write
  
  | { g_\sigma(z) =  \lambda^{z^2/\sigma^2}/\sigma }.
  
  where {\lambda = 1/\sqrt{2e\pi}}

  GAUSSIAN-SAMPLED SIGNALS

  We assume given a digital signal {s[*]}, nominally bi-infinite,
  obtained from some continuous signal {S(*)} by sampling with 
  a Gaussian kernel with standard deviation {\sigma}. That is
  
  |  { s[i] = \integral{S(x) g(x-i,\sigma)}{x} } 
  
  where {g} is the 
  
  

  Our goal is to find the best matching discrete bounded-suport filters
  for 2-downsampling with {nd} taps and 2-upsampling with {nu} taps, for
  various values of {nd,nu}. These filters are used to convert between a
  /fine/ sequence {x} and a /coarse/ sequence {y} whose samples are
  spaced twice as far apart as those of {x}. The filters are denoted
  {D,U}, so that {Y = D(x)} and {x = U(y)}.

  The filters are defined by two sets of weights {wd[0..nd-1]} and
  {wu[0..nu-1]}. The downsampling formula is
  {y[i]=SUM{wd[k]*x[2*i+k-rd]} where {rd=floor(nd/2)}. The upsampling
  formula is {x[i] = SUM{wu[k]*y[(i+k-ru)/2]} where {ru=floor(nu/2)} and
  only the terms with {i+k-ru} even are considered. Alternately, the
  upsampling formula can be viewed as splatting each {y[j]} multiplied
  by {wu[k]} onto {x[2*j+k-ru]} for {k} in {0..nu-1}.

  The two sets of weights are assumed to be reversal-symmetric, so there
  are {md=ceil(nd/2)} unknowns in {wd}, and {mu=ceil(nu/2)} unknowns in
  {wu}. E.g. for {nd=5,nu=7} the weights would be {wd=[a,b,c,b,a]}
  ({mu=3}) and {wu=[e,f,g,h,g,f,e]} ({mu=4}).

  The parameters {nd,nu} must have the same parity. If {nd,nu} are odd,
  the coarse sample {y[i]} will be aligned with sample {x[2*i]}. If
  {nd,nu} are even, the coarse sample {y[i]} will be located between
  samples {x[2*i]} and {x[2*i+1]}.

  The downsampling filter weights must add to 1, so that a constant
  unit signal produces a constant unit signal on the output. For the
  same reason, the upsampling filter weights must add to 2.
  
  The quality of the downsampling filter {D} is measured by its ability
  to suppress aliasing of frequencies in the range {[1/4 _ 1/2)], and
  preserve the rest, in some frequency-sensitive metric. The quality of
  an upsampling filter {U} is measured by the absence of frequencies in
  that range from its output. We can also evaluate both together by
  comparing {y' = D(U(y))} for some unit-norm sinusoid {y} with
  frequency in {[0,1/2)}, and {x' = U(D(x))} for some unit-norm sinusoid
  {x} with frequency in {[0,1/4]}. */

typedef struct opt_params_t
  { /* Client-specified parameters: */
    int nd;       /* Width of the downsampling filter window {wd}. */
    int nu;       /* Width of the upsampling filter window {wu}. */

    /* Derived parameters: */
    int md; /* Number of independent coefs in {wd} ({ceil(nd/2)}). */
    int mu; /* Number of independent coefs in {wu} ({ceil(nu/2)}. */

    int rd; /* Index shift for downsampling ({floor(nd/2)}). */
    int ru; /* Index shift for upsampling ({floor(nu/2)}. */

    /* Client-specified optimization options: */
    bool_t dopt;  /* Optimize the frequency response of {D}. */
    bool_t uopt;  /* Optimize the frequency response of {U}.*/

    /* Paramters of the constrained least squares system: */
    int nw; /* Total number of indep filter coefs being optimized. */
    int nc; /* Number of indep affine constraints over {wd,wu}. */
    int ne; /* Number of unknowns and equations in linear system ({nw+nc}). */
    
    int id; /* If {dopt}, the indep coefs of {wd} are variables {id..id+md-1}. */ 
    int iu; /* If {uopt}, the indep coefs of {wu} are variables {iu..iu+mu-1}. */ 
    int ic; /* The Lagrange multipliers are variables {ic..ic+nc-1} */ 
  } opt_params_t;


void optimize_many_cases(int min_nd, int max_nd);
  /* Tries to find good weightsfor several values of {nd,nu}. */

void test_specific_weights(void);
  /* Analyzes some specific filters. */

void find_specific_upsampling(void);
  /* Finds optimum upsampler {U} for a given downsampler {D}
    in the sense of {D(U(y))} being closer to {y}. */

void find_best_weights(opt_params_t *o);
    
opt_params_t *opt_parameters(int nd, bool_t dopt, int nu, bool_t uopt);
  /* Creates an {opt_params_t} for downsampling filter
    of {nd} taps and upsampling filters with minimum {nu}
    taps, with optimization choices {dopt,uopt}.  Computes the 
    parameters {.md,.mu,.mv,.nc,.ne}. */

void build_system(opt_params_t *o, double M[], double b[], double s[]);
  /* Builds the constrained least squares system for the given filter configuration and goal. */
  
void add_D_goal_term(opt_params_t *o, double M[], double b[], double f);
void add_U_goal_term(opt_params_t *o, double M[], double b[], double f);
  /* Adds to {M,b} the terms of the unconstrained least squares system 
     corresponding to frequency {f} (which must be in the range {[0 _ 1/2)}. */
     
void append_D_unit_sum_constraints(opt_params_t *o, double M[], double b[]);
void append_U_unit_sum_constraints(opt_params_t *o, double M[], double b[]);
  /* Inserts in {M,b} the equations for the unit-sum constraints. */

double spectrum_profile(double f);
  /* Returns the frequency response of the sampling process, namely 
    the amplitude {M(f)} of a discrete Fourier helix with frequency {f} 
    that is assumed to have unit importance in a discrete signal.
    
    The frequency {f} must be in the range {[-1/2 _ 1/2]}
    since other frequencies are aliased to that interval. 
    The response {M(f)} may be assumed to be real and symmetric
    around frequency 0, with {M(0)=1}. */

double complex compute_D_transfer(int nd, double wd[], double f);
  /* Returns the complex gain {GD(f)} for frequency {f} of the
    downsampling filter {D} with weights {wd[0..nd-1]}. Namely, if {x}
    is a Fourier helix with frequency {f} and amplitude 1, 
    {x[i]=cexp(2*PI*I*f*i)}, the downsampled signal {y = D(x)} will be a
    Fourier helix with frequency {2*f} and amplitude {GD(f)}. */

void rmxn_add_to_block(int m, int n, double *M, int ib, int jb, int mb, int nb, double s, double *B);
  /* Adds the matrix {B} scaled by {s} to the submatrix of {M} with size
    {mb x nb} size that starts on row {ib} and column {jb}. Assumes that
    {M} has size {m x n}, that {B} has size {mb x nb}, that {M} and {B}
    are disjoint, and that the submatrix lies entirely inside {M}. */

void extract_wd_from_solution(opt_params_t *o, double s[], double wd[]);
void extract_wu_from_solution(opt_params_t *o, double s[], double wu[]);
  /* Extracts a specific set of weights ({w}, {wu} or {wv}) from the 
    solution vector {s}, taking into account the symmetries. */

void print_system(FILE *wr, int m, int n, double M[], double b[]);
  /* Prints the system {M*z = b} to {wr}, assumed to have
    {m} equations and {n} unknowns.  Assumes {M} mas {m*n} elements,
    linearized by rows, and {b} has {m} elements. */

void print_solution(FILE *wr, opt_params_t *o, double s[]);
  /* Prints the solution {s} of the optimization prolem. */

void print_filter_weights(FILE *wr, char *name, int n, double w[]);
  /* Prints to {wr} the set of filter weights {w[k]}
    labeled "{name}[{k}]" for {k} in {0..n-1}. */

int main(int argc, char **argv);

int main(int argc, char **argv)
  {
    find_specific_upsampling();
    test_specific_weights();
    if (TRUE) 
      { int min_nd =  1; /* Min taps in downsampling filter. */
        int max_nd = 10; /* Max taps in downsampling filter. */
        optimize_many_cases(min_nd, max_nd);
      }
    return 0;
  }
  
void test_specific_weights(void)
  { 
    /* Test combinatorial downsampling filters: */
    int nd;
    for (nd = 13; nd <= 13; nd++)
      { fprintf(stderr, "--- binomial downsampling filter with nd = %d ---\n", nd);
        double wd[nd];
        wt_table_binomial_fill(nd, wd, NULL);
        wt_table_normalize_sum(nd, wd);
        int nf = 128;
        int kf;
        for (kf = 0; kf <= nf; kf++)
          { double f = 0.5*((double)kf)/((double)nf);
            double complex DTf = compute_D_transfer(nd, wd, f);
            double Mf = spectrum_profile(f);
            double Mg = spectrum_profile(2*f);
            fprintf
              ( stderr, 
                "  %4d %16.14f   %+17.14f %+17.14f  %18.14f  %18.14f  %18.14f\n", 
                kf, f, creal(DTf), cimag(DTf), cabs(DTf), Mf, cabs(DTf)*Mf/Mg
              );
          }
      }
  }
  
void find_specific_upsampling(void)
  { 
    int nd = 5;
    int nu = 5;
    fprintf(stderr, "--- trying to find the best U(%d) for binomial D(%d) ---\n", nu,nd);
    
    opt_params_t *o = opt_parameters(nd, FALSE, nu, TRUE);

    /* Get the binomial weights {wd[0..nd-1]}: */
    double wd[nd];
    wt_table_binomial_fill(nd, wd, NULL);
    wt_table_normalize_sum(nd, wd);
    
    /* Index shifs: */
    int rd = o->rd;
    int ru = o->ru;
    
    int mu = o->mu;          /* Unknowns in {wu}. */
    char* unames = "abcdefgij"; /* Names of those unknowns. */
    
    /* Combined filter {D(U())}: */
    int nz = nd + nu - 1; /* Total width (always odd). */
    int rz = nz/2;        /* Index shift. */

    /* Build a system {A*w=b} that says {y'==y} if {y' = D(U(y))} and {y = [...,0,0,0,1,0,0,0,...]} */
    bool_t req_normu = TRUE;       /* TRUE to add the even/odd normaization on {wu}. */
    int nc = (req_normu ? 1 : 0);  /* Number of extra constraints. */
    int ne = nz + nc;              /* Add one equation for sum constraint on {wu}. */
    int nw = mu;                   /* Number of unknowns. */
    double scale = pow(2,nd-1);    /* Scale all eqs by this for clarity. */
    double A[ne*nw];               /* Matrix of goal system. */
    double h[ne];                  /* RHS of goal system: */
    rmxn_zero(ne,nw,A);
    rn_zero(ne,h);
    int iy,kd,ku;
    /* Compute the right-hand sides {h[0..nz-1]}. */
    for (iy = 0; iy < nz; iy++)
      { /* Equation {iy} says that the contrib of {y[iy]} to {y'[rz]} should be {(iy==rz)}: */
        h[iy] = scale*(iy == rz);
      }
    /* Enumerate all paths through the {D(U())} filter to obtain the matrix {A}: */
    for (kd = 0; kd < nd; kd++)
      { /* Find which element {ix} of {x} is used by tap {kd} of {D} when producing {y[rz]}}: */
        int ix = 2*rz+kd-rd;
        fprintf(stderr, "  y[%d] += wd[%d]*x[%d]\n", rz, kd, ix);
        /* Scan the taps of the {wu} filter: */
        for (ku = 0; ku < nu; ku++)
          { /* Find which element {iy} of {y} is used by tap {ku} of {U} when producing {x[ix]}}. */
            /* Note that only half of them exist: */
            int iy2 = ix+ku-ru;
            if ((iy2 % 2) == 0) 
              { int iy = iy2/2;
                fprintf(stderr, "    x[%d] += wu[%d]*y[%d]", ix, ku, iy);
                /* We are looking at the path from {y[iy]} to {y'[rz]} through {wd[kd]} and {wu[ku]}. */
                assert((iy >= 0) && (iy < nz));
                /* Decide which unknown {wu[k]} is: */
                int tu = (ku <= ru ? ku : nu-1-ku);
                fprintf(stderr, "  (wu[%d] = %c)\n", ku, unames[tu]); 
                int iw = 0 + tu;
                A[iy*nw + iw] += scale*wd[kd];
              }
          } 
      }
    int ie = nz; /* Next available row in {A,h}. */
    if (req_normu)
      { /* Append a constraint that says that even and odd {wu}s must add to the same: */
        for (ku = 0; ku < nu; ku++) 
          { /* Find which unknown {iw} is {wu[ku]}: */
            int iw =  0 + (ku <= ru ? ku : nu-1-ku);
            /* The sum of even {ku}s must equal the sum of the odd {ku}s: */
            double sg = ((ku % 2) == 0 ? +1 : -1);
            A[ie*nw + iw] += scale*sg;
          }
        h[ie] = 0;
      }
    /* Print the overconstrained system: */
    print_system(stderr, ne,nw,A,h);
      
    /* Build the least squares system {M*x=b}: */
    double M[nw*nw];
    double b[nw];
    rmxn_tr_mul(ne,nw,nw,A,A,M);
    rmxn_tr_mul(ne,nw,1,A,h,b);
    print_system(stderr, nw,nw,M,b);
    fprintf(stderr, "determinant = %15.8f\n", rmxn_det(nw,M));
    
    double w[nw];
    uint32_t r1;
    gausol_solve(nw, nw, M, 1, b, w, TRUE,TRUE, 0.0, NULL, &r1);
    affirm(r1 == nw, "system is indeterminate/impossible");
    fprintf(stderr, "\n");
    print_filter_weights(stderr, "w", nw, w);
    fprintf(stderr, "\n");
    
    /* Extract weights: */
    double wu[nu];
    for (ku = 0; ku < nu; ku++) { wu[ku] = (ku <= ru ? w[0 + ku] : w[0 + nu-1-ku]); }
    print_filter_weights(stderr, "wu", nu, wu);
    fprintf(stderr, "\n");
  }

void optimize_many_cases(int min_nd, int max_nd)
  {
    int nd;
    for (nd = min_nd; nd <= max_nd; nd++)
      { int min_nu = (int)imax(nd % 2, nd - 2); /* Min taps in upsampling filter. */
        int max_nu = nd + 4;               /* Max taps in upsampling filter. */
        int nu;
        for (nu = min_nu; nu <= max_nu; nu += 2)
          { 
            bool_t dopt =  TRUE;  /* Optimize {D}? */
            bool_t uopt = FALSE;  /* Optimize {U}? */
            opt_params_t *o = opt_parameters(nd, dopt, nu, uopt);
            find_best_weights(o);
          }
      }
  }

void find_best_weights(opt_params_t *o)
  { 
    int ne = o->ne; /* Number of unknowns, incuding Lagrange multipliers. */
    double M[ne*ne]; /* System's matrix. */
    double b[ne]; /* System's RHS. */
    double s[ne]; /* Solution (weights and Lagrange multipliers). */
    build_system(o, M,b,s);
    uint32_t r2;
    gausol_solve(ne, ne, M, 1, b, s, TRUE,TRUE, 0.0, NULL, &r2);
    affirm(r2 == ne, "system is indeterminate/impossible");
    print_solution(stdout,o,s);
  }

opt_params_t *opt_parameters(int nd, bool_t dopt, int nu, bool_t uopt)
  {
    opt_params_t *o = notnull(malloc(sizeof(opt_params_t)), "no mem");
    fprintf(stderr, "building optimization parameters");
    fprintf(stderr, " nd = %d  dopt = %c", nd, "FT"[dopt]);
    fprintf(stderr, " nu = %d  uopt = %c", nu, "FT"[uopt]);
    fprintf(stderr, "\n");
    
    demand((nd%2) == (nu%2), "{nd,nu} mut have the same parity");
    
    /* Given filter parameters: */
    o->nd = nd;
    o->nu = nu;
    
    /* Independent filter coeffs: */
    o->md = (nd+1)/2;
    o->mu = (nu+1)/2;
    
    /* Index shifts: */
    o->rd = nd/2;
    o->ru = nu/2;

    /* Given optimization options: */
    o->dopt = dopt;
    o->uopt = uopt;

    /* Linear system parameters: */
    o->nw = 0;
    o->nc = 0;
    o->id = 0;
    if (dopt) 
      { o->iu = o->id + o->md;
        o->nw += o->md;
        /* Normalization may follow from freq response goals. */
      }
    else
      { o->iu = 0; }
    if (uopt)
      { o->ic = o->iu + o->mu;
        o->nw += o->mu;
        /* Normalization may follow from freq response goals. */
      }
    else
      { o->ic = o->iu; }
    /* There does not seem to be a need for noralization constraints. */
    o->ne = o->nw + o->nc;
    return o;
  }

void build_system(opt_params_t *o, double M[], double b[], double s[])
  {
    /* 
      The system has the form {[[G,C'];[C,0]] [t;1] = [w;z]} where 
      
      * {H} is the Hessian of the quadratic penalty function to be minimized
      * {t} is the linear part of that quadratic function, negated
      * {C} is the matrix of unit-sum constraints
      * {w} is the column of unknown weights
      * {z} is the column of Lagrange multipliers. 
      
      The quadratic penalty function to be minimized is the sum of 
      many terms, one for each sample frequency {f} {[0 _ 1/4]}.
      Therefore we add the Hessians of these partial penalties to get {H}. 
    */
    int ne = o->ne;
    rmxn_zero(ne, ne, M);
    rn_zero(ne,b);
    int nf = 128; /* Should be enough. */
    int kf;
    for (kf = 0; kf < 128; kf++)
      { double f = 0.5*((double)kf)/((double)nf);
        if (o->dopt) { add_D_goal_term(o,M,b,f); }
        if (o->uopt) { add_U_goal_term(o,M,b,f); }
      }
    if (o->dopt) { append_D_unit_sum_constraints(o,M,b); }
    if (o->uopt) { append_U_unit_sum_constraints(o,M,b); }
  }
  
void add_D_goal_term(opt_params_t *o, double M[], double b[], double f)
  { 
    /* Compute the matrix {G} that maps the independent {D} filter coeffs to 
       the complex gain of the {D} filter for Fourier helix of frequency {f}: */
       
    int nd = o->nd; /* Number taps in {D} filter. */
    int md = o->md; /* Number of independent {D} filter coeffs. */
    double G[2*md]; /* A {2 × md} matrix; rows are real part, imag part. */
    int i;
    double wd[nd]; /* Filter weights. */
    for (i = 0; i < md; i++)
      { /* Create a weight vector {wd} where only the unknown coef {i} is 1: */
        rn_zero(nd,wd);
        wd[i] = 1;
        wd[nd-1-i] = 1; /* Note that it may be the same as {wd[i]}. */
        /* Compute the filter's gain: */
        double complex DTf = compute_D_transfer(nd, wd, f);
        G[0*md + i] = creal(DTf);
        G[1*md + i] = cimag(DTf);
      }
      
    /* Compute the Hessian part {H} of the goal: */
    double H[md*md];
    rmxn_tr_mul(2,md,md,G,G,H);
      
    /* Compute the frequency {g} in {y} created by frequency {f} in {x}: */
    double g = (f >= 0.25 ? 2*f - 0.5 : 2*f);
    /* Compute the reference amplitudes for frequency {f} and {g}: */
    double Mf = spectrum_profile(f);
    double Mg = spectrum_profile(g);
    double S = Mf/Mg;
    rmxn_add_to_block(o->ne,o->ne,M,o->id,o->id,md,md,2*S*S,H);
    if (f < 0.25)
      { /* The goal is to preserve the frequency profile, namely {(Mf/Mg)^2*|G*wd|^2 ~ 1}: */
        /* Add the independent term: */
        for (i = 0; i < md; i++) { b[o->id + i] += 2*S*G[0*md + i]; }
      }
    else
      { /* The goal is to avoid aliasing, namely {Mf^2*|G*wd|^2/Mg^2 ~ 0}: */
        /* The independent term is zero.*/
      }
  }
  
void add_U_goal_term(opt_params_t *o, double M[], double b[], double f)
  { 
    affirm(FALSE, "not implemented");
  }

void append_D_unit_sum_constraints(opt_params_t *o, double M[], double b[])
  { 
    fprintf(stderr, "no D unit sum constraints used\n");
  }

void append_U_unit_sum_constraints(opt_params_t *o, double M[], double b[])
  { 
    fprintf(stderr, "no U unit sum constraints used\n");
  }

void rmxn_add_to_block(int m, int n, double *M, int ib, int jb, int mb, int nb, double s, double *B)  
  { 
    demand(m >= 0, "invalid {m}");
    demand(n >= 0, "invalid {n}");
    demand(mb >= 0, "invalid {mb}");
    demand(nb >= 0, "invalid {nb}");
    demand((ib >= 0) && (ib + mb <= m), "invalid row subrange");
    demand((jb >= 0) && (jb + nb <= n), "invalid column subrange");
    int i,j;
    for (i = 0; i < mb; i++)
      { int iM = ib + i;
        for (j = 0; j < nb; j++) 
          { int jM = jb + i;
            M[iM*n+jM] += s*B[i*nb+j];
          }
      }
  }

double complex compute_D_transfer(int nd, double wd[], double f)
  {
    int rd = nd/2;
    double complex DT = 0;
    int k;
    for (k = 0; k < nd; k++)
      { double complex bk = cexp(2*M_PI*I*f*((double)k-rd));
        if (FALSE) { fprintf(stderr, "    %4d  %18.14f  %+17.14f %+17.14f\n", k, wd[k], creal(bk), cimag(bk)); }
        DT += wd[k]*bk;
      }
    return DT;
  }

// double complex compute_U_transfer(int nu, double wu[], int nv, double wv[], double f)
//   {
//     int rd = nd/2;
//     double complex UT = 0;
//     int k;
//     for (k = 0; k < nu; k++)
//       { double complex bk = cexp(2*M_PI*I*f*((double)k-rd) ???);
//         UT += wd[k]*bk;
//       }
//     for (k = 0; k < nv; k++)
//       { double complex bk = cexp(2*M_PI*I*f*((double)k-rd) ???);
//         UT += wd[k]*bk;
//       }
//     return UT;
//   }
 
void print_system(FILE *wr, int m, int n, double M[], double b[])
  {
    fprintf(wr, "\n");
    fprintf(wr, "system is %d x %d\n", m,n);
    int i, j;
    for (i = 0; i < m; i++)
      { double *Mi = &(M[i*n]);
        double *bi = &(b[i]);
        for (j = 0; j < n; j++)
          { if (Mi[j] != 0)
              { fprintf(wr, " %+8.3f*w%d", Mi[j], j); }
            else
              { fprintf(wr, " %11s", ""); }
          }
        fprintf(wr, " = %+8.3f", *bi);
        fprintf(wr, "\n\n");
      }
     fprintf(wr, "\n");
  }
  
double spectrum_profile(double f)
  { 
    /* !!! TO THINK !!! */
    
    double sigma = 0.128;
    
    auto double G(double z);
    
    double G(double z)
      { return exp(-(z/sigma)*(z/sigma)/2); }
      
    if (TRUE)
      { /* Folded-over Gaussian: */
        int k = 0;
        double sum0 = 1;
        double sumf = G(f);
        bool_t converged = FALSE;
        for (k = 1; !converged; k++)
          { double g0 = G(0-k) + G(0+k);
            double gf = G(f-k) + G(f+k);
            double prev0 = sum0;
            double prevf = sumf;
            sum0 += g0;
            sumf += gf;
            converged = ((sum0 == prev0) && (sumf == prevf)); 
          }
        return sumf/sum0;
      }
    else
      { /* A raised cosine window: */
        double Mmin = 0.0;
        return (1-Mmin)*(1 + cos(2*M_PI*f))/2 + Mmin;
      }
  }

void print_solution(FILE *wr, opt_params_t *o, double s[])
  {
    if (o->dopt)
      { int nd = o->nd;
        double wd[nd];
        extract_wd_from_solution(o,s,wd);
        print_filter_weights(wr, "wd", nd, wd);
        fprintf(wr, "\n");
      }
    if (o->uopt)
      { int nu = o->nu;
        double wu[nu];
        extract_wu_from_solution(o,s,wu);
        print_filter_weights(wr, "wu", nu, wu);
        fprintf(wr, "\n");
      }
  }
  
void extract_wd_from_solution(opt_params_t *o, double s[], double wd[])
  {
    int nd = o->nd;
    int id = o->id;
    int md = o->md;
    int kd;
    for (kd = 0; kd < nd; kd++) { wd[kd] = ( kd < md ? s[id + kd] : s[id + (nd-1-kd)]); }
  }

void extract_wu_from_solution(opt_params_t *o, double s[], double wu[])
  {
    int nu = o->nu;
    int iu = o->iu;
    
    int mu = o->mu;
    
    int ku;
    for (ku = 0; ku < nu; ku++) { wu[ku] = ( ku < mu ? s[iu + ku] : s[iu + (nu-1-ku)]); }
  }

void print_filter_weights(FILE *wr, char *name, int n, double w[])
  { 
    int i;
    for (i = 0; i < n; i++)
      { fprintf(wr, "  %s[%d] = %+17.14f\n", name, i, w[i]); }
  }

