/* See {mesh_plot.h} */
/* Last edited on 2013-10-02 03:15:00 by stolfilocal */

#include <basic.h>
#include <triangulate.h>
#include <mesh.h>
#include <mesh_plot.h>
#include <plot_utils.h>
#include <argparser_geo.h>

#include <r2.h>
#include <r3.h>
#include <r4.h>
#include <r3x3.h>
#include <r4x4.h>
#include <rmxn.h>
#include <hr3.h>
#include <pswr.h>
#include <affirm.h>
#include <argparser.h>
#include <bool.h>
#include <sign.h>
#include <frgb.h>

#include <values.h>
#include <stdio.h>
#include <math.h>

void plot_mesh(
    PSStream *fps,
    hr3_pmap_t *map,         /* Perspective projection matrix. */
    mesh_t *tri,             /* Reference triangulation, of NULL */
    int N,                   /* Mesh subdivision parameter. */
    frgb_t *color,           /* Color to use for {fMin}. */
    r3_t *dLight,            /* Direction towards main light source. */
    double shadow            /* Amount of darkening by shadow. */
  )
  {
    auto void TrianglePaint(r3_t *P, r3_t *Q, r3_t *R);
      /* Paints color-mapped values of an affine scalar field in a 
        flat triangle, given the corner points and their 
        function values. */
    
    /* auto void PlotEdge(qarc_t e); */
      /* Draws the edge {e}. */

    int i;

    qarc_vec_t side = tri->side;
    pswr_comment(fps, "painting the faces");
    fprintf(stderr, "painting %d faces\n", side.ne/2);
    enum_trianglets(tri, N, 0.0001, &TrianglePaint);

    /* Draw the edges: */
    
    qarc_vec_t arc = tri->arc;
    pswr_comment(fps, "drawing the edges");
    fprintf(stderr, "drawing %d edges\n", arc.ne/2);
    pswr_set_pen(fps, 0.0,0.0,0.0,  0.1, 0.0,0.0);
    for (i = 0;  i < arc.ne; i += 2)
      { qarc_t e = arc.e[i];
        if (! quad_arc_is_null(e)) 
          { r3_t *P = &(ORGP(e)); 
            r3_t *Q = &(DESTP(e)); 
            draw_edge(fps, P, Q, map);
          }
      }

    qarc_vec_t out = tri->out;
    pswr_comment(fps, "drawing the vertices");
    fprintf(stderr, "drawing %d vertices\n", out.ne/2);
     for (i = 0;  i < out.ne; i ++)
      { qarc_t e = out.e[i]; 
        if (! quad_arc_is_null(e))
          { r3_t *P = &(ORGP(e)); 
            draw_vertex(fps, P, map);
          }
      }

    void TrianglePaint(r3_t *P, r3_t *Q, r3_t *R)
      { paint_triangle(fps, P, Q, R, map, color, dLight, shadow); }
  }
