/* See {mesh_pov.h} */
/* Last edited on 2013-10-02 03:14:34 by stolfilocal */

#include <basic.h>
#include <triangulate.h>
#include <mesh.h>
#include <mesh_pov.h>
#include <pov_utils.h>

#include <r3.h>
#include <affirm.h>
#include <bool.h>
#include <sign.h>
#include <frgb.h>

#include <values.h>
#include <stdio.h>
#include <math.h>

/* INTERNAL PROTOTYPES */

void pov_mesh(FILE *fpov, mesh_t *tri, int N, double radius)
  {
    auto void do_pov_triangle(r3_t *P, r3_t *Q, r3_t *R);
      /* Paints color-mapped values of an affine scalar field in a 
        flat triangle, given the corner points and their 
        function values. */

    fprintf(fpov, "\n// Faces\n\n");
    fprintf(stderr, "painting %d faces\n", tri->side.ne);
    enum_trianglets(tri, N, 0.0, &do_pov_triangle);

    void do_pov_triangle(r3_t *P, r3_t *Q, r3_t *R)
      { pov_triangle(fpov, P, Q, R, "tx_mesh_f"); }

    if (radius > 0.0)
      { int NE = tri->arc.ne/2;
        fprintf(fpov, "\n// Edges\n\n");
        fprintf(stderr, "drawing %d edges\n", NE);
        int i;
        for (i = 0; i < NE; i ++) 
          { qarc_t a = tri->arc.e[2*i];
            r3_t *P = &(ORGP(a));
            r3_t *Q = &(DESTP(a));
            pov_cylinder(fpov, P, Q, radius, "tx_mesh_e");
          }
      }

    if (radius > 0.0)
      { int NV = tri->out.ne;
        fprintf(fpov, "\n// Vertices\n\n");
        fprintf(stderr, "plotting %d vertices\n", NV);
        int i;
        for (i = 0; i < NV; i ++) 
          { qarc_t a = tri->out.e[i];
            r3_t *P = &(ORGP(a));
            pov_sphere(fpov, P, radius, "tx_mesh_v");
          }
      }
  }

