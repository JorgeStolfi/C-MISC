#include <math.h>
#include "krtimagepos.h"
#include <r3.h>

#define DEGREETORADIAN (3.1415926/180.)

void imagepos( 
  r3_t obsp,                 /* Nominal World coordinates of observer */
  r3_t sctr,                 /* World coordinates of image's center */
  r3_t up,                   /* Observer's approximate "up" direction */
  double hfov, double vfov,  /* Camera's field of view, in degrees */
  int nx, int ny,            /* Image dimensions, in pixels */
  t_image *image             /* Out: image position and size */
)
{
    int i,j;
    r3_t sx, sy, sz;
    double dist, pixsize;

    sz[0] = obsp[0] - sctr[0];
    sz[1] = obsp[1] - sctr[1];
    sz[2] = obsp[2] - sctr[2];
    dist = r3_normalize(sz);
    image->sz = sz;

    /* sx = up cross sz */
    r3_cross(sx, up, sz);
    r3_normalize(sx);
    image->sx = sx;
    
    /* sy = sx cross sz */
    r3_cross(sy, sx, sz);
    r3_normalize(sy);
    image->sy = sy;

    /* scale sx, sy so that pixels are 1x1: */
    image->nx = nx;
    image->dx = 2.0 * dist * tan( 0.5 * hfov * DEGREETORADIAN ) / nx;
    image->firstx = (0.5 - (((float)nx)/2.0))*image->dx;

    image->ny = ny;
    image->dx = 2.0 * dist * tan( 0.5 * vfov * DEGREETORADIAN ) / ny;
    image->firsty = (0.5 - (((float)ny)/2.0))*image->dy;

}

    
