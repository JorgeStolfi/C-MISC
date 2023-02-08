// Last edited on 2019-05-06 18:34:34 by jstolfi

#include <nanotube_pics_defs.h>

void nanotube_pics_draw_ray
  ( PSStream *ps,
    nanotube_pics_style_t *sty,
    double cx,
    double cy,
    double wx,
    double wy,
    bool_t dash, 
    int32_t stage,
    double dim,
    nanotube_pics_elem_style_t *esty
  )
  {
    double R = hypot(sty->xtot, sty->ytot);

    double wL = sqrt(wx*wx + wy*wy);
    double wnx = wx/wL; double wny = wy/wL;

    double p1x = cx + R*wnx; 
    double p1y = cy + R*wny;
    
    bool_t arrow = FALSE;
    nanotube_pics_draw_segment(ps,sty, cx,cy, p1x,p1y, dash,arrow,stage,dim,esty);
  }  
