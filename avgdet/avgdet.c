/* Last edited on 2024-12-21 11:57:06 by stolfi */ 

/* Determines the RMS value of the determinant of a near-singular matrix. */

#include <stdint.h>
#include <stdio.h>
#include <math.h>

#include <jsrandom.h>
#include <affirm.h>
#include <rn.h>
#include <rmxn.h>
#include <rmxn_test_tools.h>

#define N_MAX 8
#define N_TRIALS 100000
#define EPS_MIN 1.0e-11
#define EPS_MAX 1.0e-8

void avgdet_find_avg_det
  ( int32_t n,
    double eps,
    double *avg_P,
    double *gva_P,
    bool_t verbose
  );
  /* Computes and prints the RMS average values {avg} of {det(A)}
    and {gva} of {det(B)} when {A} is a normal matrix at distance {eps} from a
    singular {n×n} matrix, and {B} is the inverse of {A}
    normalized. Returns {avg} and {gva} in {*avg_P} and {*gva_P}.
    
    More precisely, {A} is derived from a singular {n×n} matrix {S},
    with unit RMS element value and {n>=2}, by adding to each element of
    {S} a random perturbation in the range {[-eps _ +eps]} and then
    re-normalizing the result to unit RMS elem value. {B} is then the
    inverse of {A}, scaled to unit RMS elem value. */

void avgdet_perturb(int32_t n, double *A, double eps, bool_t verbose);
  /* Adds to each element of the  {n×n} array {A} a random
    perturbation in the range {[-eps _ +eps]}. */

void advget_normal_invert(int32_t n, double *A, double *B, bool_t verbose);
  /* Stores into {B} the inverse of the {n×n} matrix {A}, and then
    scales {B} so that it has unit RMS elem value. */
    
void avgdet_print_array(int32_t n, double *A, char *name);
  /* Prints {A} on {stderr}. */

int32_t main(int32_t argn, char **argv)
  {
    srand(1993);
    srandom(1993);
    
    double dr[N_MAX+1]; /* Value of {avgDet/eps} for each {n}. */
    dr[0] = 0;
    dr[1] = 0;
    fprintf(stderr, "  %2s %12s %12s %12s %12s %12s\n", "n", "eps", "avg", "gva", "avg/eps", "gva/sep");
    for (int32_t n = 2; n <= N_MAX; n++)
      { double eps = 0.0;
        while (eps <= 1.001*EPS_MAX)
          { bool_t verbose = FALSE;
            double avg, gva;
            avgdet_find_avg_det(n, eps, &avg, &gva, verbose);
            
            fprintf(stderr, "  %2d %12.4e %12.4e %12.4e", n, eps, avg, gva);
            if (eps != 0)
              { double ref = eps;
                double fer = pow(eps, n-1);
                fprintf(stderr, " %12.4f %12.4f", avg/ref, gva/fer);
              }
            fprintf(stderr, "\n");

            /* Save last entry of this {n} for finite difference analysis: */
            dr[n] = avg/eps;

            eps = (eps == 0 ? EPS_MIN : eps*10);
          }
        fprintf(stderr, "\n");
      }
    /* Print finite differences: */
    fprintf(stderr, "finite differences of {avgDet/eps}:\n");
    for (int32_t kd = N_MAX; kd >= 1; kd--)
      { /* At this stage {dr[1..kd]} are the differences. */
        for (int32_t j = 2; j <= kd; j++)
          { fprintf(stderr, " %+10.6f", dr[j]);
            dr[j-1] = dr[j] - dr[j-1];
          }
        fprintf(stderr, "\n");
      }
    return 0;
  }

void avgdet_find_avg_det
  ( int32_t n,
    double eps,
    double *avg_P,
    double *gva_P,
    bool_t verbose
  )
  { 
    if (verbose) { fprintf(stderr, "  computing avg det for eps = %12.4e\n", eps); }
    double *A = talloc(n*n, double); 
    double *B = talloc(n*n, double); /* Inverse of {A}. */
    double sum_dA2 = 0, sum_dB2 = 0;
    for (uint32_t k = 0;  k < N_TRIALS; k++)
      { bool_t sub_verbose = (verbose && (k <= 2));
        if (sub_verbose) { fprintf(stderr, "    trial %d ...\n", k); }

        avgdet_throw_singular(n, A, sub_verbose);
        if (sub_verbose) { avgdet_print_array(n, A, "normalized singular {A}"); }
        rmxn_test_tools_check_all_different(n, n, A, "bug!");

        avgdet_perturb(n, A, eps, sub_verbose);
        if (sub_verbose) { avgdet_print_array(n, A, "normalized perturbed {A}"); }
        
        /* Compute determinants: */
        double dA = rmxn_det(n, A);
        if (sub_verbose) { fprintf(stderr, "    det(A) = %24.16e\n", dA); }
        sum_dA2 += dA*dA;
        
        if (eps > 0)
          { advget_normal_invert(n, A, B, sub_verbose);
            if (sub_verbose) { avgdet_print_array(n, B, "normalized inverse {B}"); }

            double dB = rmxn_det(n, B);
            if (sub_verbose) { fprintf(stderr, "    det(B) = %24.16e\n", dB); }
            sum_dB2 += dB*dB;
          }
      }
    double avg = sqrt(sum_dA2/N_TRIALS);
    double gva = sqrt(sum_dB2/N_TRIALS);
    
    free(A);
    free(B);
    
    (*avg_P) = avg;
    (*gva_P) = gva;
  }

void avgdet_throw_singular(int32_t n, double *A, bool_t verbose)
  { rmxn_throw_singular(n, A);
    if (verbose) { avgdet_print_array(n, A, "raw singular {A}"); }
    double Aenorm = rmxn_norm(n, n, A)/n; /* RMS of each element. */
    if (verbose) { fprintf(stderr, "    RMS of element = %24.16e\n", Aenorm); }
    rmxn_scale(n, n, 1.0/Aenorm, A, A);
    Aenorm = rmxn_norm(n, n, A)/n;
    if (verbose) { fprintf(stderr, "    after scaling =  %24.16e\n", Aenorm); }
  }    

void avgdet_perturb(int32_t n, double *A, double eps, bool_t verbose)
  { 
    for (uint32_t i = 0;  i < n; i++)
      { for (uint32_t j = 0;  j < n; j++)
         { A[i*n + j] += (2*drandom() - 1)*eps; }
      }
  }

void advget_normal_invert(int32_t n, double *A, double *B, bool_t verbose)
  { 
    rmxn_inv(n, A, B);
    if (verbose) { avgdet_print_array(n, B, "raw inverse {B}"); }
    double Benorm = rmxn_norm(n, n, B)/n;
    if (verbose) { fprintf(stderr, "    RMS of element = %24.16e\n", Benorm); }
    if ((! isfinite(Benorm)) || (Benorm < 1.0e-200))
      { avgdet_print_array(n, A, "    perturbed normalized {A}");
        avgdet_print_array(n, B, "    raw inverse {B}");
        fprintf(stderr, "    RMS of element = %24.16e\n", Benorm);
      }
    demand(isfinite(Benorm) && (Benorm >= 1.0e-200), "inverse is bad or too small");
    rmxn_scale(n, n, 1.0/Benorm, B, B);
    Benorm = rmxn_norm(n, n, B)/n;
    if (verbose) { fprintf(stderr, "    after scaling =  %24.16e\n", Benorm); }
  }

void avgdet_print_array(int32_t n, double *A, char *name)
  { fprintf(stderr, "    %s:\n", name);
    rmxn_gen_print(stderr, n, n, A, "%+10.7f", "", "", "", "       [ ", " ", " ]\n");
    fputc('\n', stderr);
  }
