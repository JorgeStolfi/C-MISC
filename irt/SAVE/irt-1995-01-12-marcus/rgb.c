/* See rgb.h */

#include "rgb.h"

void rgb_black (rgb_t *r)
  { r->c[0] = 0.0; r->c[1] = 0.0; r->c[2] = 0.0; }

void rgb_white (rgb_t *r)
  { r->c[0] = 1.0; r->c[1] = 1.0; r->c[2] = 1.0; }

int rgb_is_black (rgb_t *r)
  { return((r->c[0] == 0.0) && (r->c[1] == 0.0) && (r->c[2] == 0.0)); }

int rgb_is_almost_black (rgb_t *r)
  { return((r->c[0] < 0.001) && (r->c[1] < 0.001) && (r->c[2] < 0.001)); }

void rgb_weigh (rgb_t *a, rgb_t *b, rgb_t *r)
  { r->c[0] = a->c[0] * b->c[0];
    r->c[1] = a->c[1] * b->c[1];
    r->c[2] = a->c[2] * b->c[2];
  }

void rgb_weigh_in (rgb_t *a, rgb_t *b, rgb_t *r)
  { r->c[0] += a->c[0] * b->c[0];
    r->c[1] += a->c[1] * b->c[1];
    r->c[2] += a->c[2] * b->c[2];
  }

void rgb_scale (double s, rgb_t *a, rgb_t *r)
  { r->c[0] = a->c[0] * s;
    r->c[1] = a->c[1] * s;
    r->c[2] = a->c[2] * s;
  }

void rgb_mix_in (double s, rgb_t *a, rgb_t *r)
  { r->c[0] += a->c[0] * s;
    r->c[1] += a->c[1] * s;
    r->c[2] += a->c[2] * s;
  }

