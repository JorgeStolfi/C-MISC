/* rgb.h --- RGB color triples (scaled to [0__1]) */

#ifndef RGB_H
#define RGB_H

/*
  Note that an rgb_t has the same size and shape as an r3_t. */

typedef struct { double c[3]; } rgb_t;

void rgb_black (rgb_t *r);
  /*
    Sets $*r$ to the black color (0,0,0) */

void rgb_white (rgb_t *r);
  /*
    Sets $*r$ to the white color (1,1,1) */

void rgb_scale (double s, rgb_t *a, rgb_t *r);
  /*
    Sets $r$ to color $a$ scaled by $s$. */

void rgb_mix_in (double s, rgb_t *a, rgb_t *r);
  /*
    Adds color $a$ scaled by $s$ to color $r$. */

void rgb_weigh (rgb_t *a, rgb_t *b, rgb_t *r);
  /*
    Sets $*r$ to the componentwise product of $*a$ and $*b$ */

void rgb_weigh_in (rgb_t *a, rgb_t *b, rgb_t *r);
  /*
    Adds the componentwise product of $*a$ and $*b$ to $r$. */

int rgb_is_black (rgb_t *r);
  /*
    TRUE iff all components of $*r$ are zero. */

int rgb_is_almost_black (rgb_t *r);
  /*
    TRUE iff all components of $*r$ are less than 0.001. */

#endif
