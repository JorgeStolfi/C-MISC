// Last edited on 2019-05-06 11:25:38 by jstolfi
 
#include <nanotube_pics_defs.h>
 
void nanotube_pics_draw_dot 
  ( PSStream *ps,
    nanotube_pics_style_t *sty, 
    double cx,
    double cy,
    int32_t stage,
    double dim,
    nanotube_pics_elem_style_t *esty
  )
  {
    nanotube_pics_elem_style_t fsty;
    double ex,ey;
    nanotube_pics_stage_fudge(sty, stage,dim,esty, &ex,&ey,&fsty);
    
    pswr_set_pen(ps, fsty.dcol.c[0], fsty.dcol.c[1], fsty.dcol.c[2], fsty.lw, 0,0);
    pswr_set_fill_color(ps, fsty.fcol.c[0], fsty.fcol.c[1], fsty.fcol.c[2]);
    pswr_circle(ps, cx+ex,cy+ey, fsty.rdot, TRUE, TRUE);
  }  
