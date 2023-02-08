// Last edited on 2019-05-06 21:09:25 by jstolfi

#include <nanotube_pics_defs.h>

void nanotube_pics_draw_fig_master_labels
  ( PSStream *ps,
    nanotube_pics_style_t *sty, 
    double a1x,
    double a1y,
    double ux,
    double uy,
    double vx,
    double vy,
    int32_t stage
  ) 
  {
    fprintf(stderr, "!! hi\n");
    
    nanotube_pics_elem_style_t *vatsty = &(sty->vatsty);
    nanotube_pics_elem_style_t *vbtsty = &(sty->vbtsty);
    nanotube_pics_elem_style_t *vntsty = &(sty->vntsty);
    nanotube_pics_elem_style_t *edgsty = &(sty->edgsty);
    double dim = 0.0;
    
    // Draw the sector boundary lines:
    bool_t rdash = TRUE;
    nanotube_pics_draw_ray(ps,sty, a1x,a1y, ux,uy, rdash, stage,dim,edgsty);
    nanotube_pics_draw_ray(ps,sty, a1x,a1y, vx,vy, rdash, stage,dim,edgsty);
    nanotube_pics_draw_ray(ps,sty, a1x,a1y, ux+vx,uy+vy, rdash, stage,dim,edgsty);

    // Plot dot at the reference atom:
    nanotube_pics_draw_dot(ps,sty, a1x,a1y, stage,dim,vatsty);
    
    // Enumerate the valid pairs
    int32_t N = 10;
    for (int32_t iu = 1; iu < N; iu++) {
      for (int32_t iv = 0; iv < N; iv++) {
        nanotube_pics_elem_style_t *esty = (iv <= iu ? vatsty : vbtsty);
        double wx = iu*ux + iv*vx; double a2x = a1x + wx;
        double wy = iu*uy + iv*vy; double a2y = a1y + wy;
        
        // Draw a dot on the atom: 
        nanotube_pics_draw_dot(ps,sty, a2x,a2y, stage,dim,esty);

        if 
          ( ((iu == 3) && (iv == 3)) || 
            ((iu == 3) && (iv == 1)) || 
            ((iu == 1) && (iv == 3)) || 
            ((iu == 4) && (iv == 0))
          )
          {
            // Plot the arrow:
            bool_t vdash = FALSE;
            bool_t varrow = TRUE;
            nanotube_pics_draw_segment(ps,sty, a1x,a1y, a2x,a2y, vdash,varrow,stage,dim,esty);
        
            // Draw the label:
            char *xpair = NULL;
            asprintf(&xpair, "(%d,%d)", iu, iv);
            (void)nanotube_pics_draw_label(ps,sty, a2x,a2y, 0,0,0, xpair, 0.5,-0.7, stage,dim,vntsty);
            free(xpair);
          }
      }
    }
  }

