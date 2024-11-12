#! /usr/bin/gawk -f
# Last edited on 2013-12-17 02:35:14 by stolfilocal

BEGIN{
  abort = -1;
  
   # Generates a sequence of values {X[0..NX-1]}  by adding multiple sinusoids
   # {X[i] = SUM{ ampl[r]*sin(2*PI*(i*freq[r]/NX + phase[r]) : r \in 1..NF }}
   # where {ampl[1..NF],freq[1..NF],phase[1..NF]} are user-specified.
   
   # !!! Add random noise term !!!

   if (NX == "") { arg_error("must define {NX}"); }
   if (A == "") { arg_error("must define {A}"); }
   if (F == "") { arg_error("must define {F}"); }
   if (P == "") { arg_error("must define {P}"); }
   
   NX = NX + 0; 
   
   # Get model coefficients {ampl[1..NF],freq[1..NF],phase[1..NF]}: 
   NF = split(A, ampl, /[ ]*[,][ ]*/);
   NF2 = split(F, freq, /[ ]*[,][ ]*/);
   NF3 = split(P, phase, /[ ]*[,][ ]*/);
   if (NF2 != NF) { arg_error("inconsistent elem count {A,F}"); }
   if (NF3 != NF) { arg_error("inconsistent elem count {A,P}"); }
   for (k = 1; k <= NF; k++) 
     { printf "%3d %12.8f %12.8f %12.8f\n", k, ampl[k], freq[k], phase[k] > "/dev/stderr"; }
   
   pi = 3.14159265358979323844;
   
   for (i = 0; i < NX; i++)
     { Xi = 0;
       for (k = 1; k <= NF; k++) { Xi = Xi + ampl[k]*sin(2*pi*(i*freq[k]/NX + phase[k])); }
       printf "%5d %16.10f\n", i, Xi;
     }
    abort = 0;
    exit (abort);
  }
     
(abort >= 0) { exit(abort); }

END { if (abort >= 0) { exit(abort); } }

function arg_error(msg)
  { printf "** %s\n", msg > "/dev/stderr"; 
    abort = 1;
    exit(abort);
  } 
   
   
