/* See h3.h */

#include "h3.h"
#include <r3.h>
#include <math.h>
#include <ioprotos.h>
#include <stdio.h>

void h3_from_cart (r3_t *c, h3_point_t *p)
  { p->c[0] = 1.0;
    p->c[1] = c->c[0];
    p->c[2] = c->c[1];
    p->c[3] = c->c[2];
  }

void h3_to_cart (h3_point_t *p, r3_t *c)
  { double w = p->c[0];
    c->c[0] = p->c[1] / w;
    c->c[1] = p->c[2] / w;
    c->c[2] = p->c[3] / w;
  }

void h3_infty (r3_t *dir, h3_point_t *p)
  { p->c[0] = 0.0;
    p->c[1] = dir->c[0];
    p->c[2] = dir->c[1];
    p->c[3] = dir->c[2];
  }

void h3_mix (double at, h3_point_t *a, double bt, h3_point_t *b, h3_point_t *p)
  { p->c[0] = at * a->c[0] + bt * b->c[0];
    p->c[1] = at * a->c[1] + bt * b->c[1];
    p->c[2] = at * a->c[2] + bt * b->c[2];
    p->c[3] = at * a->c[3] + bt * b->c[3];
  }

void h3_inf_reduce(h3_point_t *a)
  { double m = 0.0;
    double ai;
    ai = fabs(a->c[0]); if (ai > m) m = ai;
    ai = fabs(a->c[1]); if (ai > m) m = ai;
    ai = fabs(a->c[2]); if (ai > m) m = ai;
    ai = fabs(a->c[3]); if (ai > m) m = ai;
    if (m != 1.0)
      {
	a->c[0] /= m;
	a->c[1] /= m;
	a->c[2] /= m;
	a->c[3] /= m;
      }
  }

void h3_print_point (FILE *f, h3_point_t *a)
  { fprintf(f, 
      "[%16.8e %16.8e %16.8e %16.8e]", 
      a->c[0], a->c[1], a->c[2], a->c[3]
    );
  }

void h3_print_plane (FILE *f, h3_point_t *a)
  { fprintf(f, 
      "<%16.8e %16.8e %16.8e %16.8e>", 
      a->c[0], a->c[1], a->c[2], a->c[3]
    );
  }
