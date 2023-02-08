/* Basic operations for generating POV-Ray models. */
/* Last edited on 2005-08-26 15:21:50 by stolfi */

#ifndef pov_utils_H
#define pov_utils_H

#include <r3.h>
#include <argparser.h>
#include <stdio.h>

void pov_sphere(FILE *fpov, r3_t *P, double radius, char *txname);
  /* Outputs to {fpov} a call to the POV-Ray "sphere" primitive,
    given the Cartesian coordinates of the center and the
    radius. The sphere will be painted with the POV-Ray texture
    called {txname}. */

void pov_cylinder(FILE *fpov, r3_t *P, r3_t *Q, double radius, char *txname);
  /* Outputs to {fpov} a call to the POV-Ray "cylinder" primitive,
    given the Cartesian coordinates of the basis centers and the
    radius. The cylinder will be painted with the POV-Ray texture
    called {txname}. */

void pov_triangle(FILE *fpov, r3_t *P, r3_t *Q, r3_t *R, char *txname);
  /* Outputs to {fpov} a call to the POV-Ray "triangle" primitive,
    given the Cartesian coordinates of the corner points. The triangle
    will be painted with the POV-Ray texture called {txname}. */

void pov_point(FILE *fpov, r3_t *P);
  /* Writes to {fpov} the coordinates of the point {P}, in the POV-Ray format. */

#endif
