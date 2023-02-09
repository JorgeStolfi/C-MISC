/* Last edited on 2008-07-14 22:21:39 by stolfi */
/* See the authorship and copyright notice at the end of this file. */
/* Changes:

*/

#include <jclimage.h>
#include <jclbasic.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <jclerror.h>

/* INTERNAL PROTOTYPES */

#define Gamma  (2.0)
#define LogGamma (0.69314718055994530942)
  
void SelectColors(pixel *map, long NC)
  { 
    long ip;
    for (ip=0; ip<NC; ip++)
      { 
        pixel pi;
        double y, R, G, B, tR, tB;
        double d = ((double)ip)/((double)(NC-1));
        d = (d<1?(d>0?d:1.0):0.0);
        if (Gamma != 1) y = ( d==0 ? 1 : (exp(d*LogGamma)-Gamma)/(1-Gamma) );

        tR = y*y*(3-2*y);
        tB = y;

        R = 0.03*y + 0.97*(1-cos(3*M_PI*tR))/2;
        B = (1-cos(7*M_PI*tB))/2;
        G = (y - 0.299*R - 0.114*B)/0.588;
        if ((G < 0) || (G > 1)) fprintf(stderr, "y = %f G = %f\n", y, G);
        pi.R = (int)(R*255.0 + 0.5);
        pi.G = (int)(G*255.0 + 0.5);
        pi.B = (int)(B*255.0 + 0.5);
        map[ip] = pi;
      }
  }
    
void WritePseudoColorImage(FILE* f, byte *g, long M, long N, pixel *map, long NC)
  {
    long k = 0;
    long nBytes = M*N;
    /* Write PPM header: */
    fprintf(f, "P6\n%ld %ld\n255\n", N, M);
    /* Write distances as color values: */
    for(k = 0; k < nBytes; k++)
      { byte b = g[k];
        pixel pi = map[(b < NC ? b : NC-1)];
        putc(pi.R, f); putc(pi.G, f), putc(pi.B, f);
      }
    fflush(f);
  }

byte *MatrixToImage(
    double **d, long *row, long M, long *col, long N, 
    double dMin, double dMax, long NC
  )
  {
    long i, j;
    byte *g = (byte *)Alloc(M*N*sizeof(byte));
    long k = 0;
    double hMin = log(dMin);
    double hMax = log(dMax);
    hScale = ((double)(NC-1.000001))/(hMax - hMin);
    for(i=0;i<M;i++)
      { for (j=0;j<N;j++) 
          { 
            double dd = d[row[i]][col[j]];
            if (dd <= 0)
              { g[k] = 0; }
            else
              { double hh = log(dd);
                long ip = (long)((hh-hMin)*hScale + 1);
                g[k] = (byte)(ip < NC ? (ip > 0 ? ip : 0) : NC-1);
              }
            k++;
          }
      }
    return g;
  }
  
