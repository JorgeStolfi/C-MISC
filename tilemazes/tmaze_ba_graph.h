#ifndef tmaze_ba_graph_H
#define tmaze_ba_graph_H

/* Graph representation of Brasilia Airport mazes. */
/* Last edited on 2009-11-09 23:35:24 by stolfi */

#include <stdint.h>
#include <limits.h>

#include <bool.h>
#include <dgraph.h>
#include <tmaze_ba.h>

dgraph_t tmaze_ba_graph_make(tmaze_t *M);
  /* Creates the undirected graph for the Brasilia Airport maze {M}. 
    
    For each cell, the graph has a /center vertex/. For each side of
    each cell there is a /joint vertex/, except that adjacent cells
    share the correponding joint vertices. The T-path of the tile in
    each cell is represented by three graph edges connecting the
    center vertex to the three appropriate joint vertices. */

dgraph_vertex_index_t tmaze_ba_graph_center_vertex(tmaze_t *M, int x, int y);
  /* Returns the number of the vertex of the graph of maze {M}
    that reprsents the center of the cell at column {x} and row {y}. */

bool_t tmaze_ba_graph_is_center_vertex(tmaze_t *M, dgraph_vertex_index_t v);
  /* TRUE iff {v} is a cell center vertex of the graph of maze {M}. */

dgraph_vertex_index_t tmaze_ba_graph_joint_vertex(tmaze_t *M, int x, int y, tmaze_dir_t dir);
  /* Returns the index of the joint vertex of the graph of maze {M}
    associated to the {dir} side of the cell at column {x} and row {y}. */

bool_t tmaze_ba_graph_is_joint_vertex(tmaze_t *M, dgraph_vertex_index_t v);
  /* TRUE iff {v} is a cell joint vertex of the graph of maze {M}. */

void tmaze_ba_graph_sort_edge_endpoints
  ( tmaze_t *M, 
    dgraph_vertex_index_t v1, 
    dgraph_vertex_index_t v2, 
    dgraph_vertex_index_t *vCP, 
    dgraph_vertex_index_t *vJP
  );
  /* Given the indices {v1,v2} of the end vertices of an edge of the
    maze {M}, in any order, returns in {*vCP} the vertex which is a
    cell-center, and in {*vJP} the other vertex, which must be a joint
    vertex of the same cell. */

void tmaze_ba_graph_cell_from_center_vertex
  ( tmaze_t *M, 
    dgraph_vertex_index_t v, 
    int *xP, 
    int *yP
  );
  /* Given a cell center vertex {v} of of the graph of maze {M},
    returns the column {*xP} and row {*yP} of the corresponding
    cell. */

void tmaze_ba_graph_get_cell_and_direction_from_joint_vertex
  ( tmaze_t *M, 
    dgraph_vertex_index_t v, 
    int *xP, 
    int *yP,
    tmaze_dir_t *dirP
  );
  /* Given a cell joint vertex {v} of of the graph of maze {M},
    returns the column {*xP} and row {*yP} of the corresponding cell
    in the maze, and the direction {*dirP} of the joint.
    
    The latter may be {tmaze_dir_W} (in which case {*xP,*yP} is the cell
    just East of the joint, ranging in {0..M.nxWE-1,0..M.ny-1}) or
    {tmaze_dir_S} (in which case {*xP,*yP} is the cell just North of the
    joint, ranging in {0..M.nx-1,0..M.nySN-1}). */

void tmaze_ba_graph_get_cell_and_direction_from_edge
  ( tmaze_t *M, 
    dgraph_vertex_index_t v1, 
    dgraph_vertex_index_t v2,
    int *xP, 
    int *yP, 
    tmaze_dir_t *dirP
  );
  /* Given the two end vertices {v1,v2} of an edge of the graph of
    {M}, returns the column {*xP} and row {*yP} of the corresponding
    cell in the maze, and the direction {*dirP} of the edge. The
    latter may be {tmaze_dir_N}, {tmaze_dir_W}, {tmaze_dir_S},or
    {tmaze_dir_E}. In any case {*xP,*yP} will range in
    {0..nx-1,0..ny-1}. */

/* COMPONENT SIZES */

void tmaze_ba_graph_count_components_by_size
  ( dgraph_t *G, 
    tmaze_t *M, 
    tmaze_size_t *msP,
    tmaze_comp_count_t **ctP
  );
  /* Takes the graph {G} of maze {M}. Returns in {*msP} the maximum
    size (that is, tile count) {ms} of any connected component of
    {G}, and in {*ctP} the address of a newly allocated vector
    {ct[0..ms]} such that {ct[k]} is the number of connected
    components of {G} that cover exactly {k} tiles (i.e. that have {k}
    cell-center vertices), for each {k} in {0..ms}. In particular
    {ct[0]} is the number of isolated joint vertices. Note that {ct}
    has {ms+1} elements, not {ms}. */
    
double tmaze_ba_graph_predicted_comp_count(tmaze_t *M, tmaze_size_t size);
  /* Statistically expected count of components with {size} vertices
    in a maze with {M}'s size and topology but with random tiles.
    
    Currently tries to estimate only components that do not
    wind around the torus. It is accurate only for large mazes, and
    returns NAN for most larger {size}s. */

#endif
