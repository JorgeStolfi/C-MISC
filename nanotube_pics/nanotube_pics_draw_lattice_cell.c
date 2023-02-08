// Last edited on 2019-05-06 12:46:18 by jstolfi

#include <nanotube_pics_defs.h>

void nanotube_pics_draw_lattice_cell
  ( PSStream *ps,
    nanotube_pics_style_t *sty,
    double x, 
    double y, 
    double a1x,
    double a1y,
    double wx,
    double wy,
    int32_t stage,
    double dim
  )
  {
    double dx = sty->dx;
    double dy = sty->dy;
    
    nanotube_pics_draw_atom(ps,sty, x +   dx/4, y + dy/6, +dx/2, +dy/2, a1x,a1y, wx,wy, stage, 0, dim);
    nanotube_pics_draw_atom(ps,sty, x + 3*dx/4, y - dy/6, +dx/2, -dy/2, a1x,a1y, wx,wy, stage, 1, dim);
  }
