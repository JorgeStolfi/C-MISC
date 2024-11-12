/* h3.h --- points and planes of oriented projective 3-space */

#ifndef H3_H
#define H3_H

#include <r3.h>
#include <stdio.h>

/*
  Note that an h3_point_t or h3_plane_t has the same 
  size and shape as an r4_t. */

typedef struct { double c[4]; } h3_point_t;

typedef struct { double c[4]; } h3_plane_t;

void h3_from_cart (r3_t *c, h3_point_t *p);
  /*
    Stores in $*p$ the point whose Cartesian coordinates are $*c$ */

void h3_to_cart (h3_point_t *p, r3_t *c);
  /*
    Stores in $*c$ the Cartesian coordinates of point $*p$ */

void h3_infty (r3_t *dir, h3_point_t *p);
  /*
    Stores in $*p$ the point at infinity whose (hither) direction is $*dir$. */

void h3_mix (double at, h3_point_t *a, double bt, h3_point_t *b, h3_point_t *p);
  /*
    Stores in $p$ the linear combination $at*a + bt*b$ */

void h3_inf_reduce(h3_point_t *a);
  /*
    Scales the homogeneous coordinates of $a$ so that 
    the largest one in absolute value is $\pm 1$. */

void h3_print_point (FILE *f, h3_point_t *a);
  /* 
    Prints $a$ to file $f$ as [w x y z]. */

void h3_print_plane (FILE *f, h3_point_t *a);
  /* 
    Prints $a$ to file $f$ as <W X Y Z>. */

#endif
