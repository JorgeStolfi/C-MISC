#! /usr/bin/gawk -f
# Last edited on 2004-12-23 01:58:16 by stolfi

BEGIN{
  Expt=2.5; Coef=3; Zmin=7;
  for (i=1; i<=10; i++) 
    { Zi=exp(log(i/Coef)/Expt)+Zmin; 
      printf "%24.16e\n", Zi;
    }
}

