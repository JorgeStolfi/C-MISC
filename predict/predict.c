#define PROG_NAME "predict"
#define PROG_DESC "extrapoate time-series data"
#define PROG_VERS "1.1"

/* Copyright © 2008 by the State University of Campinas (UNICAMP). */
/* Last edited on 2023-02-25 15:39:02 by stolfi */


#include <bool.h>
#include <affirm.h>
#include <rmxn.h>
#include <gauss_elim.h>
#include <rn.h>
#include <jsfile.h>
#include <vec.h>

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <values.h>

/* GENERAL PARAMETERS */

typedef struct options_t 
  { int m;  /* Number of coefficients in the formula. */
  } options_t;

/* INTERNAL PROTOTYPES */

options_t *parse_options(int argc, char **argv);
  /* Parses the command line options, returns an {options_t}. */

double_vec_t read_data(FILE *rd);
  /* Reads from {rd} a sequence of numeric data values from {stdin}, 
    returns them as a {double_vec_t} value. */

int main (int argc, char **argv);

void analyze(int n, double X[], int m, double C[]);
  /* Finds the formula {Y[i] = SUM{C[k]*X[i+k] : k \in 0..m-1} that 
    best predicts the value of {X[i+m]} for {i} in {0..n-m-1},
    in the least-squares sense. */
  
void predict(int n, double X[], int m, double C[], double Y[]);
  /* Computes the predicted values {Y[0..n-m-1]} given the values
    {X[0..n-1]} and the predictor coefficients {C[0..m-1]}. */

/* IMPLEMENTATIONS */

int main (int argc, char **argv)
  { options_t *o = parse_options(argc, argv);
    int m = o->m;
    double_vec_t Xv = read_data(stdin);
    int n = Xv.ne;
    double *X = Xv.e;
    demand(n >= 2*m, "insufficient data");
    double C[m];
    analyze(n, X, m, C);
    double Y[n-m];
    predict(n, X, m, C, Y);
    write_output(n, X, m, C, Y);
    fclose(stderr);
    fclose(stdout);
    return (0);
  }
  
void analyze(int n, double X[], int m, double C[])
{
  /* Build the system {A*D = 0} */
  double *A = rmxn_alloc(m+1,m+1);
  int i, j, k;
  for (i = 0; i <= m; i++)
    { for (j = 0; j <= i; j++)
        { double sum = 0.0;
          for (k = 0; k < n-m; k++) { sum += X[i+k]*X[j+k]; }
          A[i*m+j] = A[j*m+i] = sum;
        }
    }
}
  
void find_compact_pulse(int ns, int nt, double X[])
  { /* Working storage for the goal function: */
    fftw_complex *in = fftw_malloc(sizeof(fftw_complex)*nt); 
    fftw_complex *ot = fftw_malloc(sizeof(fftw_complex)*nt);
    fftw_plan plan = fftw_plan_dft_1d(nt, in, ot, FFTW_FORWARD, FFTW_ESTIMATE);

    /* Compute the max freq allowed {fMax}: */
    demand(ns < nt, "{ns} must be less than {nt}");
    demand(ns % 2 == 0, "{ns} must be even");
    int fMax = (nt - ns)/2;
    
    /* Compute the number {nk} of constraints (freqs to kill): */
    int nk = ns - 1;
    assert(nk == ns-1); /* We should have just one freedom left. */
    
    /* Build the linear system that says that hi-freq terms are zero. */
    int rows = nk + 1;    /* {nk} freq-kill eqs plus unit-sum eq. */
    int cols = ns;        /* {ns} unknowns (samples). */
    double A[rows*cols];  /* The constraint coefficient matrix. */
    double B[rows];       /* The independent-term vector. */
    double R[rows];       /* Residuals of constraints. */
    double Y[cols];       /* Estimated error in {X} computed from residuals. */
    int iu = nk;          /* Row of unit-sum constraint in {M} */
    int i, j;
    /* Frequency-kill constraints: */
    for (i = 0; i < nk; i++)
      { for (j = 0; j < ns; j++)
          { int f = fMax + 1 + i;
            A[i*cols + j] = sin(M_PI/4 + (2*M_PI*f*j)/nt);
          }
        B[i] = 0;
      }
    /* Unit-sum constraint: */
    for (j = 0; j < ns; j++) { A[iu*cols + j] = 1; }
    B[iu] = 1;
    
    /* Solve the system: */
    int r1 = gsel_solve(rows, cols, A, 1, B, X, 0.0);
    demand(r1 == rows, "indeterminate system");
    
    int iter;
    for (iter = 0; iter < 4; iter++)
      { /* Print the current solution: */
        rn_gen_print(stderr, ns, X, "%20.16f", "X =\n  ", "\n  ", "\n");

        /* Compute the residual error: */
        gsel_residual(rows, cols, A, 1, B, X, R);

        /* Print the residual: */
        rn_gen_print(stderr, rows, R, "%20.16f", "R =\n  ", "\n  ", "\n");

        /* Solve the residual system: */
        int r2 = gsel_solve(rows, cols, A, 1, R, Y, 0.0);
        demand(r2 == rows, "indeterminate system");

        /* Print correction: */
        rn_gen_print(stderr, ns, Y, "%20.16f", "Y =\n  ", "\n  ", "\n");

        /* Subtract correction: */
        rn_sub(cols, X, Y, X);
      }
      
    /* Print final solution: */
    rn_gen_print(stderr, ns, X, "%20.16f", "X =\n  ", "\n  ", "\n");
    
    /* Write the pulse for plotting: */
    { FILE *wr = open_write("pulse-ot.dat", TRUE);
      int i;
      double Lipp = INFINITY; /* {log(X[i-2])} */
      double Lip = INFINITY;  /* {log(X[i-1])} */
      for (i = -2; i < nt+2; i++)
        { int k = ((i % nt) + nt) % nt;
          double Xi = (k < ns ? X[k] : 0);
          double Yi = (k < ns ? Y[k] : 0);
          /* Logs and derivatives of logs: */
          double Li = (Xi <= 0 ? log(1e-300) : log(Xi));
          double DLi = (Lip == INFINITY ? 0 : Li - Lip);
          double DDLi = (Lipp == INFINITY ? 0 : (Li - 2*Lip + Lipp)/2);
          fprintf
            ( wr, "%5d  %20.16f  %+12.5f %+12.5f %+12.5f  %20.16f\n", 
              i, Xi, Li, DLi, DDLi, Yi
            );
          Lipp = Lip; Lip = Li;
        }
      fclose(wr);
    }
    
    /* Write the power spectrum for plotting: */
    { double P[fMax+1];
      pulse_spectrum(ns, nt, X, in, ot, plan, P);
      FILE *wr = open_write("spectrum-ot.dat", TRUE);
      int f;
      double Lfpp = INFINITY; /* {log(P[f-2])} */
      double Lfp = INFINITY;  /* {log(P[f-1])} */
      for (f = 0; f <= nt/2; f++) 
        { /* Use the bilateral spectrum, it plots nicer than the folded one: */
          double Pf = (f == 0 ? P[f] : P[f]/2); 
          /* Use the bilateral spectrum, it plots nicer than the folded one: */
          /* Logs and derivatives of logs: */
          double Lf = (Pf <= 0 ? log(1e-300) : log(Pf));
          double DLf = (Lfp == INFINITY ? 0 : Lf - Lfp);
          double DDLf = (Lfpp == INFINITY ? 0 : (Lf - 2*Lfp + Lfpp)/2);
          fprintf
            ( wr, "%5d  %20.16f  %+12.5f %+12.5f %+12.5f\n", 
              f, Pf, Lf, DLf, DDLf
            );
          Lfpp = Lfp; Lfp = Lf;
        }
      fclose(wr);
    }
  }

void pulse_spectrum
  ( int ns, 
    int nt,
    double X[], 
    fftw_complex in[], 
    fftw_complex ot[], 
    fftw_plan plan, 
    double P[]
  )
  { /* Compute the Fourier transform: */
    int i;
    for (i = 0; i < nt; i++) 
      { in[i][0] = (i < ns ? X[i] : 0);
        in[i][1] = 0;
      }
    fftw_execute(plan);

    /* Extract the power spectrum from FFT. */
    fourier_to_spectrum(nt, ot, P);
  }

void fourier_to_spectrum(int nt, fftw_complex ot[], double P[])
  {  int fMax = nt/2;  /* Max frequency */
     auto double cmod2(fftw_complex c); /* Complex modulus of {c}, squared. */
     double cmod2(fftw_complex c) { return c[0]*c[0] + c[1]*c[1]; }
     /* We must combine Fourier terms with same freq. */
     /* They are {ot[i]} and {ot[n-i]}, modulo {n}. */
     /* Frequency 0 has a single term: */
     P[0] = cmod2(ot[0])/nt;
     /* Other freqs have two terms, except possibly for {fMax}: */
     int f;
     for (f = 1; f <= fMax; f++) 
       { double Pi = cmod2(ot[f]); 
         int j = nt - f;
         if (j != f) { Pi += cmod2(ot[j]); }
         P[f] = Pi/nt;
       }
  }
