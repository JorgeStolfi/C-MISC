// Last edited on 2019-05-06 19:53:47 by jstolfi

#include <nanotube_pics_defs.h>

void nanotube_pics_draw_strip_edges
  ( PSStream *ps,
    nanotube_pics_style_t *sty,
    double cx,
    double cy,
    double tx,
    double ty,
    int32_t stage
  )
  {
    nanotube_pics_elem_style_t *edgsty = &(sty->edgsty); 
    bool_t dash = TRUE;
    double dim = 0.0;
    
    nanotube_pics_draw_perp_line(ps,sty, cx,cy,       +tx,+ty, dash,stage,dim,edgsty);
    nanotube_pics_draw_perp_line(ps,sty, cx+tx,cy+ty, -tx,-ty, dash,stage,dim,edgsty);
  }  

