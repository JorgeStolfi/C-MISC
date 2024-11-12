#ifndef SHADE_H
#define SHADE_H

#include "raymath.h"
#include "light.h"
#include "object.h"

typedef struct
  {
    double  r, g, b;                /* red, green, blue.    */
  } 
  t_color;

typedef struct
  {
    double  ar, ag, ab;             /* ambient r, g, b      */
    double  dr, dg, db;             /* diffuse r, g, b      */
    double  sr, sg, sb;             /* specular r, g, b     */
    double  coef;                   /* specular coef        */
    double  refl;                   /* reflection 0-1       */
    double  transp;                 /* transparency 0-1     */
  } 
  t_surface;

void shade(
    t_3d *pos, t_3d *ray, t_3d *nrm,
    t_object *obj,
    t_color *color
)

#endif
