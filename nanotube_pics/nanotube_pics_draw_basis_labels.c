// Last edited on 2019-05-06 19:13:06 by jstolfi

#include <nanotube_pics_defs.h>

void nanotube_pics_draw_basis_labels
  ( PSStream *ps,
    nanotube_pics_style_t *sty, 
    double a0x,
    double a0y,
    double ux,
    double uy,
    double vx,
    double vy,
    int32_t stage
  ) 
  {
    nanotube_pics_elem_style_t *labsty = &(sty->labsty);
    nanotube_pics_elem_style_t *bassty = &(sty->bassty);
    
    double aux = a0x+ux;
    double auy = a0y+uy;

    double avx = a0x+vx;
    double avy = a0y+vy;
 
    // Plot dots at the reference atom:
    double dim = 0.0;
    nanotube_pics_draw_dot(ps,sty, a0x,a0y, stage,dim,bassty);
    nanotube_pics_draw_dot(ps,sty, aux,auy, stage,dim,bassty);
    nanotube_pics_draw_dot(ps,sty, avx,avy, stage,dim,bassty);
    
    // Plot the basis arrows:
    bool_t dash = FALSE;
    bool_t arrow = TRUE;
    nanotube_pics_draw_segment(ps,sty, a0x,a0y, aux,auy, dash,arrow,stage,dim,bassty);
    nanotube_pics_draw_segment(ps,sty, a0x,a0y, avx,avy, dash,arrow,stage,dim,bassty);

    // Write labels:
    double mux = a0x+0.5*ux;
    double muy = a0y+0.5*uy;
    (void)nanotube_pics_draw_label(ps,sty, mux,muy, +uy,-ux, 40*px, "u", 0.5,1.0, stage,dim,labsty);
    
    double mvx = a0x+0.5*vx;
    double mvy = a0y+0.5*vy;
    (void)nanotube_pics_draw_label(ps,sty, mvx,mvy, -vy,+vx, 30*px, "v", 1.0,0.0, stage,dim,labsty);
  }

