/* Plotting a mesh of triangles. */
/* Last edited on 2005-08-20 14:42:45 by stolfi */

#ifndef mesh_plot_H
#define mesh_plot_H

#include <basic.h>
#include <triangulate.h>
#include <mesh.h>
#include <plot_utils.h>

#include <bool.h>
#include <sign.h>
#include <pswr.h>
#include <frgb.h>
#include <r3.h>
#include <r4.h>
#include <r4x4.h>
#include <hr3.h>
#include <argparser.h>

void paint_faces(
    PSStream *fps,
    hr3_pmap_t *map,        /* Perspective projection matrix. */
    mesh_t *tri,             /* Mesh to plot. */
    int N,                   /* Mesh subdivision parameter. */
    frgb_t *color,           /* Base color of mesh. */
    r3_t *dLight,            /* Direction towards main light source. */
    double shadow            /* Amount of darkening by shadow. */
  );  
  /* Paints the faces of the mesh {tri} with color {color}.
  
    Each triangle is split into {N^2} trianglets.
    
    Parts of the sphere away from the direction {dLight} are darkened
    slightly, and parts facing towards it are lightened. The {shadow}
    parameter specifies the maximum relative change in the color's
    intensity, in either sense. */
 
#endif
