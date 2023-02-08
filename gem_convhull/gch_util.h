#ifndef util_H
#define util_H
/* Last edited on 2014-07-22 02:21:48 by stolfilocal */

#define MAX(a,b) (((a) > (b)) ? (a) : (b))
#define MIN(a,b) (((a) < (b)) ? (a) : (b))

void gch_int_array_sort(int v[], int n, int ixv[]);
  /* Stores into {ixv[0..n-1]} a permutation of {0..n-1} such that
    the integers {v[ixv[0..n-1]]} are sorted in non-decreasing order. */

int gch_parse_int_vals(char *str, int nmax, int vals[]);
  /* Parses the string {str} into zero or more integers, stores them in {val[0..n-1]},
    returns the count {n} of integers parsed.  Fails if {n > nmax}. */

#endif
