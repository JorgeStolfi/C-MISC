// Last edited on 2019-05-06 11:36:03 by jstolfi

#include <nanotube_pics_defs.h>

void nanotube_pics_draw_atom
  ( PSStream *ps,
    nanotube_pics_style_t *sty,
    double cx,
    double cy,
    double sx,
    double sy,
    double a1x,
    double a1y,
    double wx,
    double wy,
    int32_t stage,
    int32_t class,
    double dim 
  )
  {
    double dotw = wx*wx + wy*wy;
    bool_t vis;
    if (dotw == 0) {
      vis = TRUE;
    } else {
      double dotp = (cx - a1x)*wx + (cy - a1y)*wy;
      vis = ((dotp > -0.001*dotw) && (dotp < 0.999*dotw));
    }

    if (vis)
      { 
        // Bonds:
        double mbd = 1.2;      //  Extra length of bond.
        double hbx = mbd*sx/2; double hby = -mbd*sy/3;     //  Vector towards "horiz" neighbor.
        double vbx = 0;        double vby = mbd*2*sy/3;    //  Vector towards "vert" neighbor.
        
        nanotube_pics_elem_style_t *bdsty = &(sty->bdsty);
        bool_t dash = FALSE;
        bool_t arrow = FALSE;
        nanotube_pics_draw_segment(ps,sty, cx, cy, cx+hbx, cy+hby, dash,arrow,stage,dim,bdsty);
        nanotube_pics_draw_segment(ps,sty, cx, cy, cx-hbx, cy+hby, dash,arrow,stage,dim,bdsty);
        nanotube_pics_draw_segment(ps,sty, cx, cy, cx+vbx, cy+vby, dash,arrow,stage,dim,bdsty);

        nanotube_pics_elem_style_t *atsty = &(sty->atsty[class]);
        nanotube_pics_draw_dot(ps,sty,  cx,cy, stage,dim,atsty);
      }
  }

