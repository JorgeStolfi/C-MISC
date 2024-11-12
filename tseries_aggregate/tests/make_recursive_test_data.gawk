#! /usr/bin/gawk -f
# Last edited on 2014-01-14 20:12:18 by stolfilocal

BEGIN \
  {
    abort = -1;
    srand(4615);

    # Generates a series of data samples {X[0..NX-1]} by a linear recurrence formula.
    # The formula is {X[i] = coeff[1] + SUM{ coef[r+1]*X[i-r] : r \in 1..NH } + D * RND(i)}
    # where {coeff[1..NH+1]} is a user-given vector,
    # and {RND(i)} is a random normal variable with zero mean and unit deviation.
    # The first {NH} elements {X[0..NH-1]} are user-specified.

    if (NX == "") { arg_error("must define {NX}"); } # Number of data samples to generate
    if (C == "") { arg_error("must define {C}"); } # Recurrence coeffs, comma-separated.
    if (X == "") { arg_error("must define {X}"); } # Initial values, comma-separated.
    if (D == "") { arg_error("must define {D}"); } # Deviation of random component.

    NX = NX + 0;

    # Split {C} into model coefficients {coef[1..NC]}: 
    # Coefficient {C[1]} is the independent term: 
    NC = split(C, coef, /[ ]*[,][ ]*/);
    if (NC < 2) { arg_error("must give at elast two coeffs in {C}"); }
    if (NX < NC) { arg_error("too few data points"); }
    NH = NC-1; # Number of historical values of {X} in recurrence formula. 

    # Get initial {X} values {init[1..NI]}: 
    NI = split(X, init, /[ ]*[,][ ]*/);
    if (NI != NH) { arg_error("wrong number of initial values"); }

    # History values are {Z[1..NC]}.
    # {Z[0]} is always 1.
    # The last {NH=NC-1} values are {Z[r] = X[i-r]} for {r} in {1..NH}.
    split("", Z);
    Z[0] = 1;
    for (i = 0; i < NX; i++)
      { if (i < NH) 
          { Xi = init[i+1]; }
        else
          { Xi = 0;
            for (k = 0; k <= NH; k++) { Xi = Xi + coef[k+1]*Z[k]; }
          }
          Xi = Xi + D * nrandom();
          printf "%5d %16.12f\n", i, Xi;
        # Shift {Z[1..NH]} and store {Xi}:
        for (k = NH; k > 1; k--) { Z[k] = Z[k-1]; }
        Z[1] = Xi;
      }
     abort = 0;
     exit (abort);
  }
     
(abort >= 0) { exit(abort); }

END \
  { if (abort >= 0) { exit(abort); } }

function nrandom(   n,i,sum,avg,dev) 
  {
    # Normal random variable with zero mean and unit deviation.
    n = 20; # Number of uniform variables to add. 
    sum = 0;
    for (i = 0; i < n; i++) { sum = sum + rand(); }
    avg = 0.5*n; # Expcted value of {sum}.
    dev = sqrt(n/12.0);
    return (sum - avg)/dev;
  }

function arg_error(msg)
  { printf "** %s\n", msg > "/dev/stderr"; 
    abort = 1;
    exit(abort);
  } 
   
   
