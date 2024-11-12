#! /usr/bin/gawk -f
# Last edited on 2012-12-23 12:10:14 by stolfilocal

# Given a list of decimal numbers, tries to find matching fractions with common divisor.

BEGIN {
  split("", frac);
  nf = 0;
}

# Remove #-comments
//{ gsub(/[\011]/," ", $0); gsub(/[ ]*[\#].*$/, "", $0); }

# Ignore blank lines
/^[ ]*$/ { next; }

# Save numbers in list
/^[ ]*[-+]?[0-9]*([0-9]|[.][0-9]+)(|[eE][-+]?[0-9]+)[ ]*$/ { 
  v = $1 + 0;
  if (v < 0) { frac[nf] = -v; } else { frac[nf] = v; }
  nf++;
}

# Process them:
END {
    printf "found %d fractions\n", nf > "/dev/stderr";
  if (nf == 0) { printf "no fractions?\n" > "/dev/stderr"; exit(1); }
  
  split("",den);
  split("",err);
  nden = 0;
  for (d = 1; d < 1000; d++)
    { 
      terr2 = 0;
      for (i = 0; i < nf; i++)
        { ni = int(frac[i]*d + 0.5);
          ei = frac[i] - ni/d;
          terr2 = terr2 + ei*ei;
        }
      avgerr = sqrt(terr2/nf);
      relerr = avgerr*d;
      if (relerr < 0.05)
        { den[nden] = d; err[nden] = relerr; nden++; }
    }
  /* Sort {den,err} by increasing error: */
  for(i = 0; i < nden; i++)
    { tden = den[i]; terr = err[i];
      j = i; 
      while((j > 0) && (err[j-1] > terr))
        { err[j] = err[j-1]; den[j] = den[j-1]; j--; }
      den[j] = tden; err[j] = terr;
    }
  for (i = 0; i < nden; i++)
    { printf  "%6d %8.6f ", den[i],err[i] > "/dev/stderr";
      for (j = 0; j < nf; j++) { printf " %13.6f", frac[j]*den[i] > "/dev/stderr"; }
      printf "\n" > "/dev/stderr";
    }
}
