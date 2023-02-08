// Last edited on 2019-05-14 21:30:38 by jstolfi
// Generates Postscript pictures for the "carbon nanotube" Wikipedia article.

#include <nanotube_pics_defs.h>

int main(int argc, char**argv)
  { 
    // nanotube_pics_draw_fig_master();
    nanotube_pics_draw_fig_strip(+1,00);
    nanotube_pics_draw_fig_strip(+1,+1);
    nanotube_pics_draw_fig_strip(+2,00);
    nanotube_pics_draw_fig_strip(+2,+1);
    nanotube_pics_draw_fig_strip(+2,+2);
    nanotube_pics_draw_fig_strip(+3,00);
    nanotube_pics_draw_fig_strip(+3,+1);
    //    nanotube_pics_draw_fig_strip(+3,+2);
    //    nanotube_pics_draw_fig_strip(+3,+3);
    //    nanotube_pics_draw_fig_strip(+4,00);
    nanotube_pics_draw_fig_strip(+4,+1);
    //    nanotube_pics_draw_fig_strip(+4,+2);
    //    nanotube_pics_draw_fig_strip(+4,+3);
    //    nanotube_pics_draw_fig_strip(+4,+4);
    //    nanotube_pics_draw_fig_strip(+5,+0);
    //    nanotube_pics_draw_fig_strip(+5,+1);
    //    nanotube_pics_draw_fig_strip(+5,+2);
    //    nanotube_pics_draw_fig_strip(+5,+3);
    //    nanotube_pics_draw_fig_strip(+5,+4);
    //    nanotube_pics_draw_fig_strip(+5,+5);
    nanotube_pics_draw_fig_strip(+4,-1);
    nanotube_pics_draw_fig_strip(+1,+3);
    nanotube_pics_draw_fig_strip(-1,+4);
    return 0;
  }
