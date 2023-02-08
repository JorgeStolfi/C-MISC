/* See {pov_utils.h} */
/* Last edited on 2005-08-26 15:24:03 by stolfi */

#include <pov_utils.h>

#include <r3.h>
#include <affirm.h>
#include <bool.h>
#include <sign.h>

#include <values.h>
#include <stdio.h>
#include <math.h>

void pov_sphere(FILE *fpov, r3_t *P, double radius, char *txname)
  { 
    fprintf(fpov, "sphere { ");
    pov_point(fpov, P);
    fprintf(fpov, ", %.3f", radius);
    fprintf(fpov, " texture { %s } }\n", txname);
  }

void pov_cylinder(FILE *fpov, r3_t *P, r3_t *Q, double radius, char *txname)
  { 
    fprintf(fpov, "cylinder { ");
    pov_point(fpov, P);
    fprintf(fpov, ", ");
    pov_point(fpov, Q);
    fprintf(fpov, ", %.3f", radius);
    fprintf(fpov, " texture { %s } }\n", txname);
  }

void pov_triangle(FILE *fpov, r3_t *P, r3_t *Q, r3_t *R, char *txname)
  { 
    fprintf(fpov, "triangle { ");
    pov_point(fpov, P);
    fprintf(fpov, ", ");
    pov_point(fpov, Q);
    fprintf(fpov, ", ");
    pov_point(fpov, R);
    fprintf(fpov, " texture { %s } }\n", txname);
  }
    
void pov_point(FILE *fpov, r3_t *P)
  {
    fprintf(fpov, "<");
    fprintf(fpov, " %.2f", P->c[0]);
    fprintf(fpov, ", %.2f", P->c[1]);
    fprintf(fpov, ", %.2f", P->c[2]);
    fprintf(fpov, " >");
  }
