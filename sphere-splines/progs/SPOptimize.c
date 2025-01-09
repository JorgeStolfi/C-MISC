/* See SPOptimize.h */
/* Last edited on 2024-12-21 11:30:20 by stolfi */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <vec.h>
#include <rn.h>
#include <affirm.h>
#include <nat.h>
#include <bool.h>
#include <jsrandom.h>

#include <SPOptimize.h>
#include <SPBasic.h>

void SPOptimize_Method1
  ( double func(double_vec_t x),
    double_vec_t x,
    double *fx,
    int niter, 
    double minStep,
    double maxStep, 
    char *plotName,
    bool_t verbose
  )
  {
    int iter;
    FILE *optPlot = NULL;
    
    #define SCALEFACT 2.0
    #define DECAY 0.90
    #define EXTRAPOLATE PHI
    #define INTERPOLATE (1.0/PHI)
    #define NBEESTEPS 100
    
    int NX = x.ne; /* Number of variables. */
    
    int NS = (int)(ceil(log(maxStep/minStep)/log(SCALEFACT))); /* Num of scales. */
    double_vec_t succRate = double_vec_new(NS); /* Success rate at each scale. */
    nat_vec_t sorted = nat_vec_new(NS); /* Scales sorted by success rate. */
    
    auto void ChooseNextPoint(void);
    auto void HandleFailure(double step);
    auto void HandleSuccess(double step);
      /* The steps of the optimization method. */
    
    auto int scale_index(double d);
      /* Clips the given {d} to the range {[minStep..maxStep]}
        and converts it into an integer logarithmic scale index. */
        
    auto double nominal_step(int scale);
      /* The nominal step size for the given scale. */
        
    auto void update_ratio(double d, double succ);
      /* Updates the success ratio for steps of the same scale as
        {d}, to account for a new success ({succ == 1.0}) or
        failure ({succ == 0.0}). */
        
    auto void reset_ratios(double d);
      /* Resets the success ratios to account for the fact
        that the current optimum moved by {d}. */
    
    auto void bubble_scale(int scale);
      /* Re-sorts the {sorted} list after a change in {succRate[scale]}. */
    
    int scale_index(double d)
      { d = fabs(d);
        if (d <= minStep)
          { return 0; }
        else if (d >= maxStep)
          { return NS-1; }
        else
          { double rscale = log(d/minStep)/log(SCALEFACT);
            int scale = (int)floor(rscale);
            return (scale < 0 ? 0 : (scale >= NS ? NS-1 : scale));
          }
      }
      
    double nominal_step(int scale)
      { double step = 0.99999 * minStep * pow(SCALEFACT, scale + 1);
        affirm(scale_index(step) == scale, "inconsistent scale_index");
        return step; 
      }
      
    void update_ratio(double d, double succ)
      { int scale = scale_index(d);
        /* Update success ratio: */
        succRate.e[scale] *= DECAY;
        succRate.e[scale] += (1.0 - DECAY)*succ;
        bubble_scale(scale);
      }
    
    void reset_ratios(double d)
      { int scale;
        for (scale = 0; scale < NS; scale++)
          { double step = nominal_step(scale);
            double decay = step/(step + d);
            succRate.e[scale] *= decay;
            succRate.e[scale] += (1.0 - decay)*0.5;
            bubble_scale(scale);
          }
      }
    
    void bubble_scale(int scale)
      { int i;
        /* Fix order of {scale} in the success ranking of scales: */
        for (i = 0; i < NS; i++) { if (sorted.e[i] == scale) { break; } }
        affirm(i < NS, "buggy scale ranking");
        while ((i > 0) && (succRate.e[scale] > succRate.e[sorted.e[i-1]]))
          { sorted.e[i] = sorted.e[i-1]; i--; }
        sorted.e[i] =  scale;
        while ((i < NS-1) && (succRate.e[scale] < succRate.e[sorted.e[i+1]]))
          { sorted.e[i] = sorted.e[i+1]; i++; }
        sorted.e[i] =  scale;
      }

    double_vec_t xini = double_vec_new(NX); /* Initial point: */
    double_vec_t xold = double_vec_new(NX); /* Previous best {x}. */
    double_vec_t xnew = double_vec_new(NX); /* Next trial: */
    double fxnew; /* {func(xnew)} */

    double xstep = 0.0;  /* Current extrapolation/interpolation step. */
    bool_t xinter;       /* Interpolation mode flag. */
    
    void ChooseNextPoint(void)
      { int k;
        /* Start at current optimum: */
        for (k = 0; k < NX; k++) { xnew.e[k] = x.e[k]; }
        if (fabs(xstep) > 0.0)
          { /* Move further along the previous jump: */
            double d = rn_L_inf_dist(NX, xold.e, x.e);
            if (verbose) 
              { fprintf(stderr, "  xstep = %12.8f", xstep); }
            else
              { int xscale = scale_index(xstep);
                fprintf(stderr, "(X%d)", xscale);
              }
            for (k = 0; k < NX; k++) 
              { xnew.e[k] += xstep/d * (x.e[k] - xold.e[k]); }
          }
        else
          { /* Pick the most promising step size: */
            int rscale = sorted.e[0];
            double rstep = nominal_step(rscale);
            /* Randomly perturb the coefficients by about {step}: */
            if (verbose) 
              { fprintf(stderr, "  rstep = %12.8f", rstep); }
            else
              { fprintf(stderr, "(R%d)", rscale); }
            for (k = 0; k < NX; k++) 
              { xnew.e[k] += rstep*(2*drandom() - 1.0); }
          }
      }

    void HandleSuccess(double step)
      { /* Last step reduced the current minimum: */
        int k;
        /* Reset success ratios of smaller step sizes: */
        reset_ratios(step);
        /* Update the success ratio for this step size: */
        update_ratio(step, 1.0);
        /* Give a chance also to the next largest step size: */
        update_ratio(step*SCALEFACT, 0.75);
        /* Retry this direction again next time: */
        xstep = hypot(EXTRAPOLATE*step, minStep);
        if (xstep > maxStep) { xstep = maxStep; }
        if (verbose) { fprintf(stderr, "xstep = %12.8f\n", xstep); }
        xinter = FALSE;
        /* Save last optimum and move to current one: */ 
        for (k = 0; k < NX; k++)
          { xold.e[k] = x.e[k]; x.e[k] = xnew.e[k]; }
        *fx = fxnew; 
      }
    
    void HandleFailure(double step)
      { /* Last step failed. */
        /* Update the success ratio for this step size: */
        update_ratio(step, 0.0); 
        /* Adjust extrapolation/interpolation step: */
        if (fabs(xstep) > 0.0)
          { if (! xinter)
            { /* Start interpolating, forward first: */
              xstep = - fabs(xstep)/EXTRAPOLATE;
              xinter = TRUE;
            }
            /* Shrink and reverse the interpolation step: */
            xstep = - xstep * INTERPOLATE;
            if (fabs(xstep) < minStep) { xstep = 0.0; }
          }
      }

    double start; /* CPU clock at start of iterations */
    int k, i;

    /* Initialize random generators and scale ranking: */
    srandom(1947); srand(1947);
    for (i = 0; i < NS; i++) { succRate.e[i] = 0.5; sorted.e[i] = i; }
    /* Save initial solution: */
    for (k = 0; k < NX; k++)
      { xold.e[k] = xini.e[k] = x.e[k]; }
    /* Open plot progress file, if any: */
    if ((plotName != NULL) && (plotName[0] != '\000'))
      { optPlot = open_write(txtcat(plotName, "-opt.plt"), TRUE); 
        fprintf(optPlot, "# %5s %12s  %12s %12s  %16s %16s\n",
          "iter", "time", "xstep", "step", "fbest", "ftest"
        );
      }
    start = user_cpu_time_usec();
    for (iter = 0; iter < niter; iter++)
      { if (verbose) { fprintf(stderr, "\niter %4d", iter); }
        ChooseNextPoint();
        fxnew = func(xnew);
        { double d = rn_L_inf_dist(NX, x.e, xnew.e);
          if (fxnew < *fx)
            { fprintf(stderr, "- \n"); HandleSuccess(d); }
          else
            { fprintf(stderr, "+ "); HandleFailure(d); }
          /* Report progress: */
          if (verbose) 
            { fprintf(stderr, "\niter %4d", iter); 
              fprintf(stderr, "  fxnew = %16.8e  fx = %16.8e", fxnew, *fx);
            }

          /* Plot progress: */
          if (optPlot != NULL) 
            { double time = user_cpu_time_usec()-start;
              fprintf(optPlot, "  %5d %12.6f  %12.8f %12.8f  %16.8e %16.8e\n", 
                iter, time/1000000.0,  xstep, d,  *fx, fxnew
              ); 
            }
        }
      }
    if (optPlot != NULL) { fclose(optPlot); }

    /* Beeline initial-final plot */
    if ((plotName != NULL) && (*plotName != '\000'))
      { int istep;
        double d = rn_dist(NX, xini.e, x.e);
        FILE *beePlot = open_write(txtcat(plotName, "-sec.plt"), TRUE); 
        nat_t over = NBEESTEPS/5;
        fprintf(beePlot, "#  %8s %10s %16s\n", "frac", "dist", "f");
        for (istep = -over; istep <= NBEESTEPS+over; istep++)
          { double t = ((double)istep)/((double)NBEESTEPS), s = 1-t;
            rn_mix(NX, s, xini.e, t, x.e, xnew.e);
            fxnew = func(xnew);
            fprintf(beePlot, "  %8.5f %10.8f %16.8e\n", t, t*d, fxnew);
          }
        fclose(beePlot);
      }
    fprintf(stderr, "\n");
    free(xnew.e); free(xini.e); free(xold.e);
    free(succRate.e); free(sorted.e);
  }
