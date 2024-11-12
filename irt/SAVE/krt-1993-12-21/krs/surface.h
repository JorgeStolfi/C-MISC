#ifndef SURFACE_H
#define SURFACE_H

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

#endif
