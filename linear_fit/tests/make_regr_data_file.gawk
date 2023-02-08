#! /usr/bin/gawk -f
# Last edited on 2023-02-02 11:13:21 by stolfi

# Creates a test data file for {linfit}. */
# Each line will have "{ID} {Z[i]} {W[i]} {

BEGIN {
  abort = -1;
  srand(4615*417);
  
  # Get the file name.  Must be "rangen_NX${NX}_wt${weighted}_ui${unit}"
  if (fname == "") { arg_error("must define {fname}"); }
  
  # Obtain the file paramaters from the file name: 
  nf = split(fname, fpart, /[_.\/]/);
  if (nf != 4) { arg_error("invalid {fname}"); }
  type = fpart[1];
  if (type != "rangen") { arg_error("invalid data type"); }
  NX = peel("nx", fpart[2],1,99);
  weighted = peel("wt", fpart[3],0,1);
  unit = peel("ui", fpart[4],0,1);
  
  # Write header:
  printf "# Generated by {make_regr_data_file.gawk} \"%s\"\n", fname 
  printf "# Data values {Z[i]} created by this formula\n" 
  split("", C)
  NXR = (unit ? NX-1 : NX); # Number of random {X} terms. */
  for (kx = 0; kx < NX; kx++) {
    C[kx] = 20*rand() - 10;
    printf "#   %+16.8f", C[kx];
    if (kx < NXR) { printf " * X[%2d]", kx; }
    printf "\n"
  }
  noise = 0.1 + 0.9*rand();
  printf "#   %+16.8f * RND()\n", noise;
  printf "#\n";

  # Generate and write the data:
  NZ = 1000; # Number of cases.
  split("", X);
  for (kt = 0; kt < NZ; kt++) {
    # Choose the independent variables {X[0..NX-1]} and compute the value {Z} to be fitted: */
    Z = 0.0;
    for (kx = 0; kx < NX; kx++) {
      if (kx < NXR) 
        { X[kx] = 2*rand() - 1.0; }
      else
        { X[kx] = 1; }
      Z += C[kx]*X[kx];
    }
    Z += noise*grand();
    printf "%5d %+16.8f", kt, Z;
    if (weighted) { 
      W = 0.1 + 0.9*rand();
      printf " %12.8f", W;
    }
    printf "  ";
    for (kx = 0; kx < NXR; kx++) {
      printf " %+16.8f", X[kx];
    }
    if (unit) { printf " 1"; }
    printf "\n";
  }
}

function peel(tag,x,xmin,xmax,   n) {
  if (substr(x,1,length(tag)) != tag) { arg_error(("bad field in file name \"" x "\"")); }
  n = substr(x,3) + 0;
  if ((n < xmin) || (n > xmax)) { arg_error(("value out of range \"" x "\"")); }
  return n
}

function grand(   kg,ng,s) {
  # Approximation of a random value with Gaussian (normal) distribution.
  ng = 20;
  s = 0;
  for (kg = 0; kg < ng; kg++) {
    s += 2*rand() - 1;
  }
  return s*sqrt(3/ng);
}
          
function arg_error(msg)
  { printf "** %s\n", msg > "/dev/stderr"; 
    abort = 1;
    exit(abort);
  } 

