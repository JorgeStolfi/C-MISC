/* SPOptimize.h -- General-purpose optimization heuristics. */
/* Last edited on 2005-06-06 11:15:51 by stolfi */

#ifndef SPOptimize_H
#define SPOptimize_H

#include <vec.h>
#include <bool.h>

void SPOptimize_Method1
  ( double func(double_vec_t x),
    double_vec_t x,
    double *fx,
    int niter, 
    double minStep,
    double maxStep, 
    char *plotName,
    bool_t verbose
  );
  /* Tries to adjust the parameter vector {x} so as to minimize {func(x)},
    by an ad-hoc heuristic.
    
    Upon entry, {x} must contain the initial guess, and {fx} the value
    of {func(x)}. Upon exit {x} will be the best solution found, and
    {fx} the corresponding function value.
    
    Performs up to {niter} evaluations of {func}. Probes are at least
    {minStep} and at most {maxStep} from current best solution.
    
    If {plotName} is not "", also writes a plot of the progress to
    file {plotName} plus extension ".plt". */

#endif
