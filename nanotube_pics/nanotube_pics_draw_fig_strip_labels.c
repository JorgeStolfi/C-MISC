// Last edited on 2019-05-06 22:48:37 by jstolfi

#include <nanotube_pics_defs.h>

void nanotube_pics_draw_fig_strip_labels
  ( PSStream *ps,
    nanotube_pics_style_t *sty, 
    double a1x,
    double a1y,
    double a2x,
    double a2y,
    int32_t stage
  ) 
  {
    nanotube_pics_elem_style_t *labsty = &(sty->labsty);
    nanotube_pics_elem_style_t *rfvsty = &(sty->rfvsty);
    
    double wx = a2x-a1x;
    double wy = a2y-a1y;

    // Plot dots at the reference atoms:
    double dim = 0.0;
    nanotube_pics_draw_dot(ps,sty, a1x,a1y, stage,dim,rfvsty);
    nanotube_pics_draw_dot(ps,sty, a2x,a2y, stage,dim,rfvsty);
    
    // Plot the circumference arrow:
    bool_t dash = FALSE;
    bool_t arrow = TRUE;
    nanotube_pics_draw_segment(ps,sty, a1x,a1y, a2x,a2y, dash,arrow,stage,dim,rfvsty);

    // Write labels:
    (void)nanotube_pics_draw_label(ps,sty, a1x,a1y, -wx,-wy, 130*px, "A1", 0.5,0.5, stage,dim,labsty);
    (void)nanotube_pics_draw_label(ps,sty, a2x,a2y, +wx,+wy,  90*px, "A2", 0.5,0.5, stage,dim,labsty);

   }

