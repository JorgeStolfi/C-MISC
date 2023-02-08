

#include <nanotube_pics_defs.h>

PSStream *nanotube_pics_fig_open(nanotube_pics_style_t *sty, char *fname)
  {
    double xtot = sty->xtot;
    double ytot = sty->ytot;
    
    char *prefix = "nanotube_diag";
    FILE *wr = open_write(fname, TRUE);
    bool_t eps = TRUE;
    PSStream *ps = pswr_new_stream(prefix, wr, eps, NULL, NULL, FALSE, xtot, ytot); 
    pswr_new_picture(ps, 0,xtot, 0, ytot);
    pswr_set_canvas_window(ps, 0,xtot, 0,ytot);
    pswr_set_client_window(ps, 0,xtot, 0,ytot);

    // Paint background:
    double tad = 2*px;
    pswr_set_fill_color(ps, sty->bgcol.c[0], sty->bgcol.c[1], sty->bgcol.c[2]);
    pswr_rectangle(ps, 0-tad, xtot+tad, 0-tad, ytot+tad, TRUE, FALSE);
    
    return ps;
  }
