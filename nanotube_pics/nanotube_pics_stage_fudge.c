// Last edited on 2019-05-06 11:59:10 by jstolfi
 
#include <nanotube_pics_defs.h>
 
void nanotube_pics_stage_fudge
  ( nanotube_pics_style_t *sty,
    int32_t stage,
    double dim,
    nanotube_pics_elem_style_t *esty,
    double *exP,
    double *eyP,
    nanotube_pics_elem_style_t *fsty
  )
  {  
    (*fsty) = (*esty); // Most fields don't change.
    
    if (stage == 0)
      { fsty->lw = esty->lw + 2*esty->rwh;
        fsty->dcol = sty->whcol;
        fsty->fcol = sty->whcol;
        (*exP) = 0.0;  (*eyP) = 0.0;
      }
    else if (stage == 1)
      { fsty->lw = esty->lw + 2*sty->rsh;
        fsty->dcol = sty->shcol;
        fsty->fcol = sty->shcol;
        (*exP) = sty->shx; (*eyP) = sty->shy;
      }
    else if (stage == 2)
      { fsty->lw = esty->lw;
        fsty->dcol = frgb_mix(1-dim, &(esty->dcol), dim, &(sty->bgcol));
        fsty->fcol = frgb_mix(1-dim, &(esty->fcol), dim, &(sty->bgcol));
        (*exP) = 0.0; (*eyP) = 0.0;
      }
  }
