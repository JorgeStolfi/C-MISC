/* Last edited on 2007-01-04 04:31:31 by stolfi */ 

#include <jsrandom.h>

#include <math.h>
#include <stdio.h>

#define NN 20
#define NT 50000

int main(int argn, char **argv)
  { int n, k, i;
    for (n = 0; n < NN; n++)
      { double sum = 0;
        for (k = 0; k < NT; k++)
          { double d2 = 0;
            for (i = 0; i < n; i++)
              { double xi = frandom();
                d2 += xi * xi;
              }
            sum += sqrt(d2);
          }
        printf("n = %3d  avg = %.2f  rms = %.2f\n", 
          n, sum/((double)NT), sqrt(((double) n)/3)
        );
      }
    return 0;
  }
                
   
