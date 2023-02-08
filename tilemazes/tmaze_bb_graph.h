#ifndef tmaze_bb_graph_H
#define tmaze_bb_graph_H

/* Graph representation of blip-blop mazes. */
/* Last edited on 2009-11-09 23:37:09 by stolfi */

#include <stdint.h>
#include <limits.h>

#include <bool.h>
#include <dgraph.h>
#include <tmaze_bb.h>

dgraph_t tmaze_bb_graph_make(tmaze_t *M);
  /* Creates the undirected graph for the Blip-Blop maze {M}. 
    
    For each side of each cell there is a /joint vertex/, which is
    shared by the two cells separated by that side. The vertex is said
    to have type {H} or {V} depending on whether the the side is
    palallel to the vertical {X} or {Y} axis, respectively. The roads
    in each tile are represented by two undirected edges, each
    connecting a distinct {V} joint of the tile to a distinct {H}
    joint. */

dgraph_vertex_index_t tmaze_bb_graph_vertex(tmaze_t *M, int x, int y, tmaze_dir_t dir);
  /* Returns the index of the joint vertex of the graph of maze {M}
    associated to the {dir} side of the cell at column {x} and row {y}. */

void tmaze_bb_graph_get_cell_and_direction_from_vertex
  ( tmaze_t *M, 
    dgraph_vertex_index_t v, 
    int *xP, 
    int *yP,
    tmaze_dir_t *dirP
  );
  /* Given a cell joint vertex {v} of of the graph of maze {M},
    returns the column {*xP} and row {*yP} of the corresponding cell
    in the maze, and the direction {*dirP} of the joint.
    
    The latter may be {tmaze_dir_W} (in which case {*xP,*yP} is the
    cell just East of the joint, ranging in {0..M.nxWE-1,0..M.ny-1})
    or {tmaze_dir_S} (in which case {*xP,*yP} is the cell just North
    of the joint, ranging in {0..M.nx-1,0..M.nySN-1}). */

void tmaze_bb_graph_get_cell_and_dirs_from_edge
  ( tmaze_t *M, 
    dgraph_vertex_index_t v1, 
    dgraph_vertex_index_t v2,
    int *xP, 
    int *yP, 
    tmaze_dir_t *dir1P,
    tmaze_dir_t *dir2P
  );
  /* Given the two end vertices {v1,v2} of an edge of the graph of
    {M}, returns the column {*xP} and row {*yP} of the corresponding
    cell in the maze, and the directions {*dir1P} and {*dir2P} of the
    sides of that cell that contain the vertices {v1} and
    {v2},respectively. The indices {*xP,*yP} will range in
    {0..nx-1,0..ny-1}, irrespective of {M.torus}. The procedure fails
    when {M} is toroidal and either {nx} or {ny} is 1, since in that
    case the cell cannot be determined from {v1} and {v2}. */

/* COMPONENT SIZES */

void tmaze_bb_graph_count_components_by_size
  ( dgraph_t *G, 
    tmaze_t *M, 
    tmaze_size_t *msP,
    tmaze_comp_count_t **ctP
  );
  /* Takes the graph {G} of maze {M}. Returns in {*msP} the maximum
    size (that is, vertex count) {ms} of any connected component of
    {G}, and in {*ctP} the address of a newly allocated vector
    {ct[0..ms]} such that {ct[k]} is the number of connected
    components of {G} with exactly {k} vertices, for each {k}.
    
    Note that each component of {G} is either a closed simple cycle or
    a simple path, and has at most {2*nt} edges. Note also that the
    vector {ct} will have {ms+1} elements. */
    
double tmaze_bb_graph_predicted_comp_count(tmaze_t *M, tmaze_size_t size);
  /* Statistically expected count of components with {size} vertices
    in a maze with {M}'s size and topology but with random tiles.
    
    Currently tries to estimate only /plain/ components. In the open
    topology, the plain components are whole closed cycles (not open
    paths). In toroidal topology, they are closed cycles that do not
    wind around the torus. It is accurate only for large mazes, and
    returns NAN for most larger {size}s. */

void tmaze_bb_graph_get_edge_components
  ( dgraph_t *G, 
    tmaze_t *M, 
    dgraph_edge_count_t *neP,
    tmaze_comp_count_t *ncP,
    tmaze_size_t **szP,
    dgraph_vertex_index_t **rtP
  );
  /* Takes the graph {G} of maze {M}. Stores into {*neP} the number of
    directed edges {ne} of {G} (which is {2*M.nt}), into {*ncP} the
    number of connected components of {G}, and into {*szP} and {*rtP} 
    the addresses of two newly allocated vectors {sz[0..nv-1]} and
    {rt[0..ne-1]}, where {nv} is the number of vertices in {G}. 
    
    Then fills {sz} so that {sz[v]} is the size of
    the component that contains vetex {v}, for each {v} in {0..nv-1}.
    
    Also fills {rt} so that {rt[e]} is the root vertex of the
    component of {G} that contains edge {e} of {G}, for each {e} in
    {0..ne-1}. The indexing of edges in {rt} is NOT the same as in
    {G}: instead, the data for edges in the cell of {M} with index {k}
    are in {rt[2*k]} for the edge incident to the South side, and
    {rt[2*k+1} for the edge incident to the North side. */
    
#endif
