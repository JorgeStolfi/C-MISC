/* Last edited on 2024-11-12 08:33:23 by stolfi */ 

/* Determines the RMS value of the determinant of a near-singular matrix. */

#define _GNU_SOURCE
#include <stdint.h>
#include <stdio.h>
#include <math.h>

#include <jsrandom.h>
#include <affirm.h>
#include <rn.h>
#include <rmxn.h>

#define N_MAX 8
#define N_TRIALS 10000
#define TOL_MIN 1.0e-12
#define TOL_MAX 1.0e-8

double avgdet_find_avg_det(int32_t n, double mag, double tol, bool_t verbose);
  /* Computes and prints the RMS average value {avg} of {det(B)}
    when {B} is at relative distance {tol} from a singular {n×n} matrix.  
    Returns {avg}.
    
    More precisely, {B} is derived from a singular {n×n} matrix {A}, with RMS
    element value {mag} and {n>=2}, by adding to each element of {A} a random
    perturbation in {[-tol _ +tol]} times {mag}. */
    
void avgdet_throw_singular(int32_t n, double *A, double mag, bool_t verbose);
  /* Fills the {n×n} array {A} with a matrix that is singular 
    apart from roundoff errors, whose elements have RMS {mag}. */

void avgdet_perturb(int32_t n, double *A, double eps, bool_t verbose);
  /* Adds to each element of the  {n×n} array {A} a random
    perturbation in the range {[-eps _ +eps]}. */

void avgdet_print_array(int32_t n, double *A, char *name);
  /* Prints {A} on {stderr}. */

int32_t main(int32_t argn, char **argv)
  {
    double dr[N_MAX+1]; /* Value of {avgDet/tol} for each {n}. */
    dr[0] = 0;
    dr[1] = 0;
    fprintf(stderr, "  %2s %12s %12s %12s %12s %12s %12s\n", "n", "mag", "tol", "f", "avg", "avg/f", "avg/f/tol");
    for (int32_t n = 2; n <= N_MAX; n++)
      { double tol = 0.0;
        double mag = 1.0;
        while (tol <= 1.001*TOL_MAX)
          { bool_t verbose = FALSE;
            double avg = avgdet_find_avg_det(n, mag, tol, verbose);

            /* Predicted RMS value of {det/tol}: */
            double f = pow(mag,n);

            fprintf(stderr, "  %2d %12.4f %12.4e %12.4e %12.4e %12.4e", n, mag, tol, f, avg, avg/f);
            if (tol != 0) { fprintf(stderr, " %12.4f", avg/f/tol); }
            fprintf(stderr, "\n");

            /* Save last entry of this {n} for finite difference analysis: */
            dr[n] = avg/f/tol;

            tol = (tol == 0 ? TOL_MIN : tol*10);
            mag = 2*mag;
          }
        fprintf(stderr, "\n");
      }
    /* Print finite differences: */
    fprintf(stderr, "finite differences of {avgDet/tol}:\n");
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

double avgdet_find_avg_det(int32_t n, double mag, double tol, bool_t verbose)
  { 
    if (verbose) { fprintf(stderr, "  computing avg det for tol = %12.4e\n", tol); }
    double *A = talloc(n*n, double);
    double sum_d2 = 0;
    for (int32_t k = 0; k < N_TRIALS; k++)
      { bool_t sub_verbose = (verbose && (k <= 2));
        if (sub_verbose) { fprintf(stderr, "    trial %d ...\n", k); }

        avgdet_throw_singular(n, A, mag, sub_verbose);
        if (verbose) { avgdet_print_array(n, A, "original {A}"); }

        avgdet_perturb(n, A, tol*mag, sub_verbose);
        if (verbose) { avgdet_print_array(n, A, "perturbed {A}"); }
        
        /* RMS value of each element: */
        double d = rmxn_det(n, A);
        if (sub_verbose) { fprintf(stderr, "    det = %24.16e\n", d); }
        sum_d2 += d*d;
      }
    double avg = sqrt(sum_d2/N_TRIALS);
    return avg;
  }

void avgdet_throw_singular(int32_t n, double *A, double mag, bool_t verbose)
  { demand(n >= 2, "invalid size");
    /* Clear out row {n-1}: */
    double *Alast = &(A[(n-1)*n]);
    for (int32_t j = 0; j < n; j++) { Alast[j] = 0.0; }
    /* Fill rows {0..n-2} with random elems, mix into row {n-1}: */
    for (int32_t i = 0; i < n-1; i++) 
      { double ri = 2*drand48() - 1;
        for (int32_t j = 0; j < n; j++)
          { A[i*n + j] = (2*drand48() - 1); 
            Alast[j] += ri * A[i*n + j];
          }
      }
    double s = rmxn_norm(n, n, A)/n; /* RMS of each element. */
    if (verbose) { fprintf(stderr, "    RMS of element = %24.16e\n", s); }
    rmxn_scale(n, n, mag/s, A, A);
    s = rmxn_norm(n, n, A)/n;
    if (verbose) { fprintf(stderr, "    after scaling =  %24.16e\n", s); }
  }

void avgdet_perturb(int32_t n, double *A, double eps, bool_t verbose)
  { 
    for (int32_t i = 0; i < n; i++)
      { for (int32_t j = 0; j < n; j++)
         { A[i*n + j] += (2*drand48() - 1)*eps; }
      }
    if (verbose) { avgdet_print_array(n, A, "perturbed"); }
  }

void avgdet_print_array(int32_t n, double *A, char *name)
  { fprintf(stderr, "    %s:\n", name);
    rmxn_gen_print(stderr, n, n, A, "%+10.7f", "", "", "", "       [ ", " ", " ]\n");
    fputc('\n', stderr);
  }
