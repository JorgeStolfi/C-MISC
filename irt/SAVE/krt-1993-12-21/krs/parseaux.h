#ifndef PARSEAUX_H
#define PARSEAUX_H

#input "shape.h"

/* Tese routines return nonzero result if request is invalid. */

int setdefaults (void);

int makeobject ( int nsurf, t_solid solid );

int makelight (
  double bright,
  double x, double y, double z
);

int makesurface (
  int nsurf,
  double ar, double ag, double ab,
  double dr, double dg, double db,
  double sr, double sg, double sb,
  double coef, 
  double refl,
  double transp
);

#endif
