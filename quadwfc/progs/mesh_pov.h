/* Generating POV-ray description of a mesh of triangles. */
/* Last edited on 2005-08-26 15:45:12 by stolfi */

#ifndef mesh_pov_H
#define mesh_pov_H

#include <basic.h>
#include <triangulate.h>
#include <mesh.h>
#include <pov_utils.h>

#include <bool.h>
#include <sign.h>
#include <r3.h>
#include <r4.h>
#include <r4x4.h>
#include <hr3.h>
#include <argparser.h>

void pov_mesh(FILE *fpov, mesh_t *tri, int N, double radius);  
  /* Writes to {fpov} a POV-ray description of the mesh {tri}.
    Each triangle is split into {N^2} trianglets.
    Edges and vertices are rendered as cylinders and balls
    with the given {radius}. */

#endif
