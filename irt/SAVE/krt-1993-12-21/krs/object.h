#ifndef OBJECT_H
#define OBJECT_H

#include "shape.h"
#include "surface.h"

typedef struct
  {
    int       id;     /* object number         */
    t_solid   solid;  /* geometry of object */
    t_surface *surf;  /* surface properties */
  } 
  t_object;

#endif
