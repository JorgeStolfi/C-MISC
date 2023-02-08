/* See {gch_util.h} */
/* Last edited on 2014-07-22 02:22:25 by stolfilocal */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include <gch_util.h>

void gch_int_array_sort(int *v, int n, int *ixv)
  {
    auto int cmp(const void *a, const void *b);
      /* Given the addresses of two elements {ixv[ia]}, {ixv[ib]},
        returns {-1,0,+1} depending on order of {v[ixv[ia]]} and {v[ixv[ib]]}. */

    int cmp(const void *a, const void *b)
      { int ixa = *((int*)a);
        int ixb = *((int*)b);
        if (v[ixa] < v[ixb])
          { return -1; }
        else if (v[ixa] > v[ixb])
          { return +1; }
        else
          { return 0; }
      }

    int i;
    for (i = 0; i < n; i++) { ixv[i] = i; }

    qsort(ixv, n, sizeof(int), cmp);
  }


int gch_parse_int_vals(char *str, int nmax, int vals[])
  { int n = 0;
    char *p = str;
    while (((*p) != '\0') && ((*p) != '\n'))
      { /* Skip blanks: */
        while ((*p) == ' ') { p++; }
        if (((*p) != '\0') && ((*p) != '\n'))
          { /* Parse the number: */
            if (n >= nmax) { fprintf(stderr, "** too many vertices\n"); assert(0); }
            int nr = sscanf(p, "%d", &(vals[n]));
            if (nr != 1) { fprintf(stderr, "** format error\n"); assert(0); }
            n++;
            /* Skip the number: */
            while (((*p) != ' ') && ((*p) != '\0') && ((*p) != '\n')) { p++; }
          }
      }
    return n;
  }
