#! /bin/gawk -f
# Last edited on 2009-11-10 00:44:15 by stolfi

# Condenses the output of {tmaze_stats.c} into power-of-two size bins.

BEGIN{ lo = 0; hi = 1; clr();  }

/^[ ]*[0-9]/ {
  sz = $1; act = $2; dct = $3; tsz = $4; act_cum = $5; tsz_cum = $6; ect = $7;
  if (sz == 0) { next; }
  if (sz > hi) { prt(); clr(); lo = hi+1; hi = 2*hi; }
  bin_act += act;
  bin_tsz += tsz;
  next;
}

//{ print; next; }

END{ prt(); }

function clr(){ bin_act = 0; bin_tsz = 0; }

function prt(){ printf "  %9d %9d  %11.1f %11.1f\n", lo, hi, bin_act, bin_tsz; }


