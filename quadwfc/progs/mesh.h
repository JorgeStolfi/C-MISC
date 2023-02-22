/* Triangular meshes. */
/* Last edited on 2023-02-12 23:50:24 by stolfi */

#ifndef mesh_H
#define mesh_H

#include <basic.h>
#include <quad.h>
#include <wavefront.h>

#include <stdio.h>

typedef struct mesh_t 
  {
    quad_arc_vec_t out;  /* Element {out[i]} is one arc out of site number {i}. */
    quad_arc_vec_t arc;  /* Elems {arc[2*i]} and {arc[2*i+1]} are the two arcs of edge {i}. */
    quad_arc_vec_t side; /* Element {side[i]} is one arc whose left face is triangle {i}. */
  } mesh_t;

mesh_t *mesh_from_topology(quad_arc_t e);
  /* Enumerates all edges, vertices, and faces of the polyhedron
    hanging from {e}, assigns consecutive numbers to them, then packs
    the tables of those elements into a {mesh_t} record. */

void write_mesh(FILE *wr, mesh_t *tri);
  /* Writes the mesh {tri} to file {wr}, in a format that 
    can be read back by {read_mesh}. */

mesh_t *read_mesh(FILE *rd);
  /* Reads a mesh description from file {rd}. */
  
typedef void triangle_visit_proc_t(r3_t *P, r3_t *Q, r3_t *R);
  /* A procedure that can process a triangle. */

void enum_trianglets(mesh_t *tri, int N, double epsilon, triangle_visit_proc_t *proc);
  /* Subdivides each non-hole face of {tri} into {N^2} trianglets, and
    calls {proc} on each piece. If {epsilon>0}, the corners of each
    face are pulled inwards by {epsilon} (relative), so that the
    meshes of adjacent faces do NOT share any points. */

#endif
