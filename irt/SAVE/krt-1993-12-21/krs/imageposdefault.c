#include <math.h>
#include "imageposdefault.h"
#include "imagepos.h"
#include "raymath.h"

void imageposdefault(t_image *image)
{
  t_3d eyep, lookp, zenith;
  int nx, ny;
  double hfov, vfov;
  
  eyep.x = 100.0;
  eyep.y = 0.0;
  eyep.z = 0.0;
  lookp.x = 0.0;
  lookp.y = 0.0;
  lookp.z = 0.0;
  zenith.x = 0.0;
  zenith.y = 1.0;
  zenith.z = 0.0;
  sizex = 512;
  sizey = 512;
  hfov = 50.0;
  vfov = 50.0;

  imagepos(
    &eyep, &lookp, &zenith,
    hfov, vfov,
    nx, ny,
    image
  );
}

    
