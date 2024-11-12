#ifndef KRTIMAGEPOS_H
#define KRTIMAGEPOS_H

#include <r3.h>

typedef struct
  {
    r3_t sctr;     /* Center of image */
    r3_t sx;       /* Screen X axis vector */
    r3_t sy;       /* Screen Y axis vector */
    r3_t sz;       /* Screen Z axis vector, pointing TOWARDS observer */
    
    int nx;        /* Number of pixels in the X direction */
    int ny;        /* Number of pixels in the Y direction */
    
    double dx;     /* Pixel size in X direction */
    double dy;     /* Pixel size in Y direction */
    
    double firstx; /* Screen X of center of leftmost pixel column */
    double firsty; /* Screen Y of center of lowest pixel row */
  } 
  t_image;

void imagepos( 
  r3_t obsp,                 /* Nominal World coordinates of observer */
  r3_t sctr,                 /* World coordinates of image's center */
  r3_t up,                   /* Observer's approximate "up" direction */
  double hfov, double vfov,  /* Camera's field of view, in degrees */
  int nx, int ny,            /* Image dimensions, in pixels */
  t_image *image             /* Out: image position and size */
);

#endif
