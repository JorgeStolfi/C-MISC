/* Last edited on 2023-01-14 12:21:47 by stolfi */

#include <memory.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <arith.h>
#include <fn.h>

typedef struct{
  float x, y, z, q;
} Charge;

/* PROTOTYPES */

/* IMPLEMENTATIONS */

#define NQ 8

static Charge ch[NQ] =
  {
/*
    (Charge){+0.230, +0.270, -0.160, -0.335},
    (Charge){-0.606, +0.690, +0.360, +0.515},
    (Charge){+1.460, +1.490, +1.440, -0.385},
*/
    (Charge){+0.540, +1.400, +1.050, +1.000},
    (Charge){+0.200, +0.950, -0.080, +1.000},
    (Charge){-1.080, 0+.330, -0.780, +1.000},
    (Charge){+0.960, +1.410, -0.160, -1.000},
    (Charge){-1.130, -0.957, -1.370, -0.500},
    (Charge){-0.040, +0.090, -0.020, -0.500},
    (Charge){+0.930, +0.060, +0.770, -0.500},
    (Charge){-1.250, -1.210, -1.290, -0.500}
  };

char *fnname(void)
{ 
  return("fncharges"); 
} 
  
void fndescr(FILE *f)
{
  fprintf(f, "GIVEN (x, y, z) DO\n");
  fprintf(f, "  f = sum{ Q_i / sqrt((x-xi)^2 + (y-yi)^2 + (z-zi)^2);\n");
  fprintf(f, "RETURN (f)\n");
}

void fneval(AAform x, AAform y, AAform z, AAform f)
{
  AAform dx, dy, dz, s;
  int k;
  
  aa_interval(f, 0.0, 0.0);
  
  for (k = 0; k< NQ; k++)
    { aa_trans(dx, x, -ch[k].x); aa_sqr(s, dx);
      aa_trans(dy, y, -ch[k].y); aa_sqr(dy, dy); aa_add(s, s, dy);
      aa_trans(dz, z, -ch[k].z); aa_sqr(dz, dz); aa_add(s, s, dz);
      aa_sqrt(s, s);
      aa_inv(s, s);
      aa_scale(s, s, ch[k].q);
      aa_add(f, f, s);
    }
}

