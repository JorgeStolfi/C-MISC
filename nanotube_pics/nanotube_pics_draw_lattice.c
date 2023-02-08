// Last edited on 2019-05-05 05:55:34 by jstolfi
 
#include <nanotube_pics_defs.h>
 
void nanotube_pics_draw_lattice
  ( PSStream *ps,
    nanotube_pics_style_t *sty,
    double a1x,
    double a1y,
    double wx,
    double wy,
    int32_t stage,
    double dim
  )
  { 
    double xtot = sty->xtot;
    double ytot = sty->ytot;

    double dx = sty->dx;
    double dy = sty->dy;
    
    double b0x = sty->b0x;
    double b0y = sty->b0y;

    // Approx number of graphene grid cells in the picture:
    int32_t NX = (int32_t)floor(xtot/dx)+1;
    int32_t NY = (int32_t)floor(ytot/dy)+1;

    for (int32_t iy = -3; iy < NY + 3; iy++) { 
      double y = b0y + dy*iy;
      for (int32_t ix = -3; ix < NX + 3; ix++) { 
        double x = b0x + dx*(ix - 0.5*(iy % 2));
        nanotube_pics_draw_lattice_cell(ps,sty, x,y, a1x,a1y, wx,wy, stage,dim);
      }
    }
  }
