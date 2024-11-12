/* See {mesh_plot.h} */
/* Last edited on 2005-08-20 15:09:04 by stolfi */

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

/* INTERNAL PROTOTYPES */

typedef void triangle_visit_proc_t(r3_t *P, r3_t *Q, r3_t *R);
  /* A procedure that can process a triangular sphere fragment. */

void generic_plot
  ( PSStream *fps,
    hr3_pmap_t *map,             /* Perspective projection matrix. */
    mesh_t *tri,                 /* Reference triangulation. */
    int N,                       /* Mesh refinement parameter. */
    triangle_visit_proc_t *proc  /* Called for each sphere fragment. */
  );
  /* Generates the plotting mesh and calls {proc} on each piece. The
    corners of each piece are in the orientation specified by 
    the initial edge {.
    
    If {tri} is not null, each face of {tri} is covered by a separate
    mesh that lies slightly inside the face; the meshes of adjacent faces
    do NOT share any points.  */
    
void paint_faces(
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
    
    void TrianglePaint(r3_t *P, r3_t *Q, r3_t *R)
      { paint_triangle(fps, P, Q, R, map, color, dLight, shadow); }
  
    int i;

    quad_arc_vec_t side = tri->side;
    pswr_comment(fps, "painting faces");
    fprintf(stderr, "painting %d faces\n", side.nel/2);
    generic_plot(fps, map, tri, N, &TrianglePaint);

    /* Draw the edges: */
    auto void PlotEdge(quad_arc_t e);
    
    void PlotEdge(quad_arc_t e)
      { r3_t *P = &(quad_org(e)->curr); 
        r3_t *Q = &(quad_dst(e)->curr); 
        draw_edge(fps, P, Q, map);
      }

    quad_arc_vec_t arc = tri->arc;
    pswr_comment(fps, "drawing the mesh");
    fprintf(stderr, "drawing %d edges\n", arc.nel/2);
    pswr_set_pen(fps, 0.0,0.0,0.0,  0.1, 0.0,0.0);
    for (i = 0;  i < arc.nel; i += 2)
      { quad_arc_t e = arc.el[i];
        if (e != NULL_REF)
          { PlotEdge(e); }
      }
  }
  
void generic_plot
  ( PSStream *fps,
    hr3_pmap_t *map,             /* Perspective projection matrix. */
    mesh_t *tri,                 /* Mesh to plot. */
    int N,                       /* Mesh refinement parameter. */
    triangle_visit_proc_t *proc  /* Called for each triangle fragment. */
  )
  { 
    auto void processTriangle(r3_t *P, r3_t *Q, r3_t *R);
      /* Subdivides the given triangle into {N^2} trianglets and calls {proc}
        on them. */

    r3_t S[N + 1]; /* Saved points. */
    
    void processTriangle(r3_t *P, r3_t *Q, r3_t *R)
      { double fN = (double)N;
        int i;
        for (i = 0; i <= N; i++)
          { int j;
            for (j = 0; j <= N-i; j++) 
              { int k = N - i - j;
                r3_t T = (r3_t)
                  {{(P->c[0]*i + Q->c[0]*j + R->c[0]*k)/fN,
                    (P->c[1]*i + Q->c[1]*j + R->c[1]*k)/fN,
                    (P->c[2]*i + Q->c[2]*j + R->c[2]*k)/fN
                  }};
                r3_dir(&T, &T); 
                if (i > 0)
                  { /* Plot triangles using points from previous row: */
                    /* Be sure to preserve orientation rel. to {P,Q,R}: */
                    r3_t *V = &(S[j]);  
                    r3_t *W = &(S[j+1]);
                    if (j > 0)
                      { r3_t *U = &(S[j-1]);
                        proc(V, U, &T);
                      }
                    proc(W, V, &T);
                  }
                /* Save point: */
                S[j] = T;
              }
          }
      }

    /* Process each part of the given triangulation */
    int NT = tri->side.nel;
    double epsilon = 0.0001; /* Relative perturbation to avoid edges */
    int i, j;
    /* Plot each triangle: */
    for (i = 0; i < NT; i++)
      { quad_arc_t e = tri->side.el[i];
        r3_t P = quad_dst(quad_onext(e))->curr;
        r3_t Q = quad_org(e)->curr;
        r3_t R = quad_dst(e)->curr;

        /* Perturb points so that they lie slightly inside the triangle: */
        for (j = 0; j < 3; j++)
          { double Bj = (P.c[j] + Q.c[j] + R.c[j])/3.0;
            P.c[j] = (1.0-epsilon)*P.c[j] + epsilon*Bj;
            Q.c[j] = (1.0-epsilon)*Q.c[j] + epsilon*Bj;
            R.c[j] = (1.0-epsilon)*R.c[j] + epsilon*Bj;
          }
        r3_dir(&P, &P);
        r3_dir(&Q, &Q); 
        r3_dir(&R, &R);

        processTriangle(&P, &Q, &R);
      }
  }

