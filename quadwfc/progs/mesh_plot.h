/* Plotting a mesh of triangles. */
/* Last edited on 2023-02-12 11:45:00 by stolfi */

#ifndef mesh_plot_H
#define mesh_plot_H

#include <basic.h>
#include <triangulate.h>
#include <mesh.h>
#include <plot_utils.h>

#include <bool.h>
#include <sign.h>
#include <epswr.h>
#include <frgb.h>
#include <r3.h>
#include <r4.h>
#include <r4x4.h>
#include <hr3.h>
#include <argparser.h>

void plot_mesh(
    epswr_figure_t *eps,
    hr3_pmap_t *map,         /* Perspective projection matrix. */
    mesh_t *tri,             /* Mesh to plot. */
    int N,                   /* Mesh subdivision parameter. */
    frgb_t *color,           /* Base color of mesh. */
    r3_t *dLight,            /* Direction towards main light source. */
    double shadow            /* Amount of darkening by shadow. */
  );  
  /* Paints the mesh {tri}.
  
    The faces of the mesh {tri} are painted with color {color}.
    Each triangle is split into {N^2} trianglets.
    Parts of the sphere away from the direction {dLight} are darkened
    slightly, and parts facing towards it are lightened. The {shadow}
    parameter specifies the maximum relative change in the color's
    intensity, in either sense.
    
    The edges of the mesh are drawn in black.
    
    The vertices of the mesh are plotted as dots, filled with 
    the color {dcolor} and outlined in black. */

#endif
