// Last edited on 2019-05-06 18:49:37 by jstolfi

#include <nanotube_pics_defs.h>

void nanotube_pics_draw_segment
  ( PSStream *ps,
    nanotube_pics_style_t *sty,
    double ppx,
    double ppy,
    double qqx,
    double qqy,
    bool_t dash,
    bool_t arrow,
    int32_t stage, 
    double dim,
    nanotube_pics_elem_style_t *esty
  )
  { 
    nanotube_pics_elem_style_t fsty;
    double ex,ey;
    nanotube_pics_stage_fudge(sty, stage,dim,esty, &ex,&ey,&fsty);
    double dashLength = (dash ? 12.0*px : 0.0);
    double dashSpace = (dash? 8.0*px : 0.0);
    
    pswr_set_pen(ps, fsty.dcol.c[0], fsty.dcol.c[1], fsty.dcol.c[2], fsty.lw, dashLength, dashSpace);
    pswr_segment(ps, ppx+ex,ppy+ey, qqx+ex,qqy+ey);
    if (arrow > 0)
      { 
        // To get the shadow and halo in the proper place we need to fill and draw:
        double arlen = fsty.arlen;     // Length in mm (with no extra).
        double arwid = fsty.arwid;     // Width in mm (with no extra).
        
        pswr_set_pen(ps, fsty.dcol.c[0], fsty.dcol.c[1], fsty.dcol.c[2], fsty.lw, 0.0,0.0);
        pswr_set_fill_color(ps, fsty.fcol.c[0], fsty.fcol.c[1], fsty.fcol.c[2]);
        pswr_arrowhead(ps, ppx+ex,ppy+ey, qqx+ex,qqy+ey, arwid,arlen,1.0, TRUE, TRUE);
      }
  }

