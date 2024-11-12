/* Checks "arith" routines */

#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <values.h>
#include "arith.h"
#include "fn.h"

/* PROTOTYPES */

void eval_whole(
  real xmin, real xmax, 
  real ymin, real ymax, 
  real zmin, real zmax, 
  real *fmin, real *fmax
);

void eval_fine(
  int n,
  real xmin, real xmax, 
  real ymin, real ymax, 
  real zmin, real zmax, 
  real *fmin, real *fmax
);

int nice_numer(int n);

void split_range(real rmin, real rmax, int n, int i, real *xmin, real *xmax);

/* PROCEDURES */

int nice_number(int n)
{
  if (n<0) n = -n;
  while ((n>0) && ((n % 10) == 0)) n /= 10;
  return ((n < 10) || ((n < 50) && ((n % 5) == 0)));
}

void split_range(real rmin, real rmax, int n, int i, real *xmin, real *xmax)
{
  *xmin = rmin + (rmax - rmin) * (real)i/(real)n;
  *xmax = rmin + (rmax - rmin) * (real)(i+1)/(real)n;
}

#define NC 9

int main(int argc, char* argv[])
{
  int level, n;
  if (argc>1 && sscanf(argv[1],"%d", &level)==1) /* ok */; else level = 4;
  n=1<<level;
  {
    int x,y,z;
    int ntests = 0;
    int bug, bad, nbug = 0;
    real accbad = 0.0;
    real totbits = 0.0;
    real avgacc;
    
    fprintf(stderr, "function = \"%s\"\n", fnname());
    fndescr(stderr);
    
    for (x=0; x<NC; x++)
      for (y=0; y<NC; y++)
        for (z=0; z<NC; z++)
          { real xmin, xmax, ymin, ymax, zmin, zmax, wmin, wmax, fmin, fmax;
            real wwid, fwid, acc;
            split_range(-3.0, 3.0, NC, x, &xmin, &xmax);
            split_range(-3.0, 3.0, NC, y, &ymin, &ymax);
            split_range(-3.0, 3.0, NC, z, &zmin, &zmax);
            eval_whole(xmin, xmax, ymin, ymax, zmin, zmax, &wmin, &wmax);
            eval_fine (n, xmin, xmax, ymin, ymax, zmin, zmax, &fmin, &fmax);
            wwid = wmax - wmin;
            fwid = fmax - fmin;
            bug = (((fmin - wmin) < - wwid * 1.e-5 ) || ((fmax - wmax) > wwid * 1.e-5));
            /* if (bug) { printf("** error **\n"); nbug++; } */
            acc = log((wmax - wmin)/((fmax - fmin) + 1.0e-30))/log(2.0);
            bad = acc > accbad;
            if (bad) { accbad = acc; }
            if (bug | bad | nice_number(ntests))
              { printf("%3d %3d %3d", x, y, z);
                printf(" [%12.6f __ %12.6f] = %12.6f", wmin, wmax, wmax - wmin);
                printf(" [%12.6f __ %12.6f] = %12.6f", fmin, fmax, fmax - fmin);
                printf(" acc = %8.3f\n", acc);
              }
            if (acc > 0.0) totbits += acc;
            ntests++;
            if (bug && (nbug>100)) exit(1); 
          }
    avgacc = totbits/(real)(ntests);
    printf("\n");
    printf("average bit loss = %9.3f  max  %9.3f\n", avgacc, accbad);
    return (0);
  }
}

void eval_fine(
  int n,
  real xmin, real xmax, 
  real ymin, real ymax, 
  real zmin, real zmax, 
  real *fmin, real *fmax
)
{
  int x,y,z;
  *fmin = MAXREAL; 
  *fmax = -MAXREAL;

  for (x=0; x<n; x++)
    for (y=0; y<n; y++)
      for (z=0; z<n; z++)
        { real xlo, xhi, ylo, yhi, zlo, zhi, flo, fhi;
          split_range(xmin, xmax, n, x, &xlo, &xhi);
          split_range(ymin, ymax, n, y, &ylo, &yhi);
          split_range(zmin, zmax, n, z, &zlo, &zhi);
          eval_whole(xlo, xhi, ylo, yhi, zlo, zhi, &flo, &fhi);
          if (flo < *fmin) *fmin = flo;
          if (fhi > *fmax) *fmax = fhi;
        }
}

void eval_whole(
  real xmin, real xmax, 
  real ymin, real ymax, 
  real zmin, real zmax, 
  real *fmin, real *fmax
)
{AAform x,y,z,f;
 aa_open();
 aa_interval(x,xmin,xmax); aa_debug(x); /* x = [ xmin .. xmax ] */
 aa_interval(y,ymin,ymax); aa_debug(y); /* y = [ ymin .. ymax ] */
 aa_interval(z,zmin,zmax); aa_debug(z); /* z = [ zmin .. zmax ] */
 fneval(x, y, z, f); aa_debug(f);
 aa_range(f, fmin, fmax);               /* f= [ fmin .. fmax ] */
 aa_close();
}

