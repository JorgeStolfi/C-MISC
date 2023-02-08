// Last edited on 2019-05-06 10:47:42 by jstolfi

#include <nanotube_pics_defs.h>

void nanotube_pics_draw_perp_line
  ( PSStream *ps,
    nanotube_pics_style_t *sty,
    double cx,
    double cy,
    double tx,
    double ty,
    bool_t dash, 
    int32_t stage,
    double dim,
    nanotube_pics_elem_style_t *esty
  )
  {
    
    double R = hypot(sty->xtot, sty->ytot);

    double tL = sqrt(tx*tx + ty*ty);
    double tnx = tx/tL; double tny = ty/tL;

    double p0x = cx + R*tny; 
    double p0y = cy - R*tnx;
    
    double p1x = cx - R*tny; 
    double p1y = cy + R*tnx;
    
    bool_t arrow = FALSE;
    nanotube_pics_draw_segment(ps,sty, p0x,p0y, p1x,p1y, dash,arrow,stage,dim,esty);
  }  
