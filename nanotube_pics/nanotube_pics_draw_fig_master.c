// Last edited on 2019-05-13 02:29:30 by jstolfi

#include <nanotube_pics_defs.h>

void nanotube_pics_draw_fig_master(void)
  {
    fprintf(stderr,  "!! Generating master picture\n");

    double xtot = 120*mm;
    double ytot = 120*mm;
    nanotube_pics_style_t *sty = nanotube_pics_def_style(xtot,ytot);
    
    char *fname = "out/nanotube_strip_master.eps";
    PSStream *ps = nanotube_pics_fig_open(sty, fname); 
    
    double orgx = sty->orgx;
    double orgy = sty->orgy;

    double dx = sty->dx;
    double dy = sty->dy;

    // Fudamental vectors of graphene lattice:
    double ux = dx;     double uy = 0;
    double vx = 0.5*dx; double vy = dy;

    // Reference atom of graphene to show sector:
    double a1x = orgx + 1.0*dx;
    double a1y = orgy + 2*dy;

    // Origin to show the fundamental vectors:
    double a0x = orgx + 1*dx;  
    double a0y = orgy + 8*dy;

    // ----------------------------------------------------------------------
    // DRAWING

    // Draw dimmed graphene lattice:
    double dim = 0.80;
    nanotube_pics_draw_lattice(ps,sty,     0.0,0.0, 0.0,0.0, 2, dim);

    // Plot label/arrow with halos:
    nanotube_pics_draw_fig_master_labels(ps,sty, a1x,a1y, ux,uy, vx,vy, 0);
    nanotube_pics_draw_basis_labels(ps,sty, a0x,a0y, ux,uy, vx,vy, 0);

    nanotube_pics_draw_fig_master_labels(ps,sty, a1x,a1y, ux,uy, vx,vy, 2);
    nanotube_pics_draw_basis_labels(ps,sty, a0x,a0y, ux,uy, vx,vy, 2);

    nanotube_pics_fig_close(ps); 
  }
