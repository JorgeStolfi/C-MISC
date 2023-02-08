// Last edited on 2019-05-06 23:14:22 by jstolfi
 
#include <nanotube_pics_defs.h>

double nanotube_pics_draw_label
  ( PSStream *ps,
    nanotube_pics_style_t *sty, 
    double qx,
    double qy, 
    double tx,
    double ty, 
    double dst,
    char *lab,
    double xalg, 
    double yalg,
    int32_t stage,
    double dim,
    nanotube_pics_elem_style_t *esty
  )
  {   
    // Normalize the displacement vector:
    double tL = sqrt(tx*tx + ty*ty);
    if (tL > 0) { tx = tx/tL; ty = ty/tL; }
    
    // Compute the stroke width for shadows and halos:
    nanotube_pics_elem_style_t fsty;
    double ex,ey;
    nanotube_pics_stage_fudge(sty, stage,dim,esty, &ex,&ey,&fsty);
    
    pswr_set_pen(ps, fsty.dcol.c[0], fsty.dcol.c[1], fsty.dcol.c[2], fsty.lw, 0.0, 0.0);
    pswr_set_fill_color(ps, fsty.fcol.c[0], fsty.fcol.c[1], fsty.fcol.c[2]);
    pswr_set_label_font(ps, fsty.font, fsty.fontsize);
    bool_t fill = TRUE;
    bool_t draw = TRUE; // Needed for halos and shadows. Set {lw} to 0 in {esty} for naked chars.
    pswr_fill_draw_label(ps, lab, qx+dst*tx+ex, qy+dst*ty+ey, 0.0, xalg, yalg, fill,draw);
    // Hack: Estimate width of label assuming Arial font:
    double wd = 0.0;
    double dch = 2.0*fsty.fontsize*px;  // Estimated width of 1 character.
    char *p = lab;
    while ((*p) != 0)
      { // Check for characters that matter:
        if (((*p) == 'i') || ((*p) == 'l') || ((*p) == 'j'))
          { wd += 0.4*dch; }
        else if (((*p) == ' ') || ((*p) == '.') || ((*p) == ',') || ((*p) == 't') || ((*p) == 'f'))
          { wd += 0.5*dch; }
        else if ((*p) == '-')
          { wd += 0.6*dch; }
        else if ((*p) == '\260')
          { wd += 0.7*dch; }
        else if ((*p) == 's')
          { wd += 0.9*dch; }
        else if ((*p) == 'w')
          { wd += 1.2*dch; }
        else if ((*p) == 'm')
          { wd += 1.4*dch; }
        else  
          { wd += dch; }
        p++;
      }
    return wd;
  }

