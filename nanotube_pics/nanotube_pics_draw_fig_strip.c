// Last edited on 2019-05-14 22:01:27 by jstolfi

#include <nanotube_pics_defs.h>

void nanotube_pics_draw_fig_strip(int32_t n, int32_t m)
  {
    fprintf(stderr,  "!! Generating strip picture for n = %d m = %d\n", n, m);

    double xtot = ((n <= 3) || (n+m <= 3) ? 120 : (n == 4 ? 160 : 200))*mm;
    double ytot = 120*mm; 
    nanotube_pics_style_t *sty = nanotube_pics_def_style(xtot,ytot);
    
    char *fname = NULL;
    if (m == 0)
      { asprintf(&fname, "out/nanotube_strip_%+03d_%03d.eps", n, m); }
    else
      { asprintf(&fname, "out/nanotube_strip_%+03d_%+03d.eps", n, m); }
    PSStream *ps = nanotube_pics_fig_open(sty, fname); 
    
    // Dimensions of lattice cell (2 atoms):
    double dx = sty->dx;
    double dy = sty->dy;

    // Fudamental vectors of graphene lattice:
    double ux = dx;     double uy = 0;
    double vx = 0.5*dx; double vy = dy;

    // Atom where to show the basis:
    double a0x = sty->orgx + 0.5*dx;
    double a0y = sty->orgy + (m < 0 ? 7.0 : 1.0)*dy;

    // Unrolling vector:
    double wx = n*ux + m*vx; 
    double wy = n*uy + m*vy;
    // double wL = hypot(wx, wy);

    // Key position:
    nanotube_pics_elem_style_t *kynsty = &(sty->kynsty);
    double dkly = 4.5*kynsty->fontsize*px; // Estimated line spacing of key.
    double dksx = 2.5*kynsty->fontsize*px; // Estimated width of one blank space.
    double kxmax = xtot - 4*px;
    double kxmin = kxmax - 10*dksx;
    double kymax = (m < 0 ? 3.6*dkly + 4*px : ytot - 4*px);
    double kymin = kymax - 3.6*dkly;

    // Strip reference atom (class 0) on graphene lattice:
    double a1x, a1y;
    double krx = kxmin;
    double kry = (m < 0 ? kymin : kymax);
    nanotube_pics_fig_strip_choose_a1(sty, a0x,a0y, ux,uy, vx,vy, krx,kry, wx,wy, &a1x,&a1y);
    
    // Strip edge reference:
    double edx = a1x;
    double edy = a1y;
    
    // Congruent atom (class 0) on graphene lattice:
    double a2x = a1x + wx;
    double a2y = a1y + wy;

    // ----------------------------------------------------------------------
    // DRAWING

    // Draw dimmed graphene lattice:
    double dim = 0.80;
    nanotube_pics_draw_lattice(ps,sty,     0.0,0.0, 0.0,0.0, 2, dim);

    // Draw the halo of the strip atoms and bonds, and of the strip edge lines:
    nanotube_pics_draw_strip_edges(ps,sty, edx,edy, wx,wy,   0);
    nanotube_pics_draw_lattice(ps,sty,     a1x,a1y, wx,wy,   0, 0.0);

    // Draw shadows of the strip atoms and bonds:
    nanotube_pics_draw_lattice(ps,sty,     a1x,a1y, wx,wy,   1, 0.0);

    // Draw the strip lattice proper:
    nanotube_pics_draw_strip_edges(ps,sty, edx,edy, wx,wy,   2);
    nanotube_pics_draw_lattice(ps,sty,     a1x,a1y, wx,wy,   2, 0.0);

    // Plot label/arrow with halos:
    nanotube_pics_draw_fig_strip_labels(ps,sty, a1x,a1y,a2x,a2y, 0);
    nanotube_pics_draw_basis_labels(ps,sty, a0x,a0y, ux,uy, vx,vy, 0);

    nanotube_pics_draw_fig_strip_labels(ps,sty, a1x,a1y,a2x,a2y, 2);
    nanotube_pics_draw_basis_labels(ps,sty, a0x,a0y, ux,uy, vx,vy, 2);

    // Plot key with halos:
    nanotube_pics_draw_fig_key(ps,sty, n,m, kxmin,kxmax, kymin,kymax, 0);
    nanotube_pics_draw_fig_key(ps,sty, n,m, kxmin,kxmax, kymin,kymax, 2);
                                            
    // (void)nanotube_pics_draw_label(ps,sty, 10*px,ytot-2*dkly, 0,0,0, "|ssssssssss|", 0,0, 2,0.0,&(kynsty));
    // (void)nanotube_pics_draw_label(ps,sty, 10*px,ytot-3*dkly, 0,0,0, "|          |", 0,0, 2,0.0,&(kynsty));
    // (void)nanotube_pics_draw_label(ps,sty, 10*px,ytot-4*dkly, 0,0,0, "|..........|", 0,0, 2,0.0,&(kynsty));
    // (void)nanotube_pics_draw_label(ps,sty, 10*px,ytot-5*dkly, 0,0,0, "|----------|", 0,0, 2,0.0,&(kynsty));
    // (void)nanotube_pics_draw_label(ps,sty, 10*px,ytot-6*dkly, 0,0,0, "|2222222222|", 0,0, 2,0.0,&(kynsty));

    free(fname);
    nanotube_pics_fig_close(ps); 
  }
