// Last edited on 2019-05-13 02:44:55 by jstolfi

#include <nanotube_pics_defs.h>

void nanotube_pics_draw_fig_key
  ( PSStream *ps,
    nanotube_pics_style_t *sty, 
    int32_t n,
    int32_t m,
    double kxmin,
    double kxmax,
    double kymin,
    double kymax,
    int32_t stage
  ) 
  {
    nanotube_pics_elem_style_t *kyvsty = &(sty->kyvsty);
    nanotube_pics_elem_style_t *kynsty = &(sty->kynsty);
    nanotube_pics_elem_style_t *kyssty = &(sty->kyssty);
    
    // Font spacing constants:
    double dly = 4.5*kynsty->fontsize*px; // Interline distance.
    double dsx = 2.5*kynsty->fontsize*px; // Width of one blank space.
    
    // Reference point on baseline of first line, right margin:
    double kcy = kymax - 1.0*dly;
    double kcx = kxmax - 1.0*dsx;
    
    // Background rectangle:
    pswr_set_fill_color(ps, sty->bgcol.c[0], sty->bgcol.c[1], sty->bgcol.c[2]);
    pswr_rectangle(ps, kxmin, kxmax, kymin, kymax, TRUE, FALSE);
       
    // Current point at right margin of key.
    double cy = kcy;
    double cx = kcx;
    double dcx; // Estimated width of stuff just written.
    double duy; // Adjustment for descender/ascender of unit.
    
    // Nanotube type {n,m} (name in Latin italic font, no unit):
    cx = kcx;
    
    dcx = nanotube_pics_draw_key_line(ps,sty, cx,cy, "m", (double)m, 0, "",0, stage,kyvsty,kynsty,kyvsty);
    cx = cx - dcx - 1.5*dsx;

    dcx = nanotube_pics_draw_key_line(ps,sty, cx,cy, "n", (double)n, 0, "",0, stage,kyvsty,kynsty,kyvsty); 
    cx = cx - dcx;

    cx = kcx; cy = cy - dly;

    // Diameter (name in Latin italic font, unit " pm" in upright font, like digits):
    double diam = 78.3*sqrt(n*n + m*m + n*m);
    duy = -0.7*kynsty->fontsize*px;  
    dcx = nanotube_pics_draw_key_line(ps,sty, cx,cy, "d", diam, 1, " pm",duy, stage,kyvsty,kynsty,kynsty); 
    cx = cx - dcx;
    
    cx = kcx; cy = cy - dly;

    // Angle  (name and unit in Symbol font):
    double wx = n + 0.5*m;
    double wy = m*0.5*sqrt(3);
    double alpha = atan2(wy,wx)*180.0/M_PI;
    duy = +2.0*kynsty->fontsize*px; 
    dcx = nanotube_pics_draw_key_line(ps,sty, cx,cy, "\141", alpha, 1, "\260",duy, stage,kyssty,kynsty,kyssty); 
    cx = cx - dcx;
    
    cx = kcx; cy = cy - dly;

  }

