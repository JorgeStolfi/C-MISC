/* See {tmaze_ba.h} */
/* Last edited on 2009-11-09 23:39:07 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <math.h>

#include <bool.h>
#include <affirm.h>
#include <dgraph.h>

#include <tmaze_ba.h>
#include <tmaze_ba_graph.h>

dgraph_vertex_index_t tmaze_ba_graph_center_vertex(tmaze_t *M, int x, int y)
  {
    tmaze_cell_index_t k = tmaze_cell_index(x, y, M->nx, M->ny);
    /* The center of tile {k} is vertex {t} (for both topologies): */
    return (dgraph_vertex_index_t)k;
  }

bool_t tmaze_ba_graph_is_center_vertex(tmaze_t *M, dgraph_vertex_index_t v)
  { dgraph_vertex_count_t nvC = M->nx * M->ny;  /* Number of center vertices. */
    return (v < nvC);
  }

dgraph_vertex_index_t tmaze_ba_graph_joint_vertex(tmaze_t *M, int x, int y, tmaze_dir_t dir)
  { dgraph_vertex_count_t nvC = M->nt;  /* Number of center vertices. */
    dgraph_vertex_count_t nvWE = M->nxWE * M->ny; /* Number of W-joints. */
    dgraph_vertex_count_t nvSN = M->nx * M->nySN; /* Number of S-joints. */
    switch (dir)
      { 
      case tmaze_dir_N:
        { /* The North joint of a tile in row {y} is the South joint
            of the tile on the row {y+1} of the same column, wrapping
            around if the topology is toroidal: */
          y = (y + 1) % M->nySN;
          goto dir_S;
        }
      case tmaze_dir_E:
        { /* The East joint of a tile in column {x} is the West joint
            of the tile in column {x+1} of the same row, wrapping
            around if the topology is toroidal: */
          x = (x + 1) % M->nxWE;
          goto dir_W;
        }
      case tmaze_dir_W:
      dir_W:
        { /* The West joint of a tile is numbered {nvC + k},
            where {k} is the tile index computed as if the maze
            had {nxWE} columns and {ny} rows. */
          tmaze_cell_index_t k = tmaze_cell_index(x, y, M->nxWE, M->ny);
          assert(k < nvWE);
          return nvC + k;
          break;
        }
      case tmaze_dir_S:
      dir_S:
        { /* The South joint of a tile is numbered {nvC + nvW + k},
            where {k} is the tile index computed as if the maze
            had {nx} columns and {nySN} rows. */
          tmaze_cell_index_t k = tmaze_cell_index(x, y, M->nx, M->nySN);
          assert(k < nvSN);
          return nvC + nvWE + k;
          break;
        }
      default: assert(FALSE);
      }
  }

bool_t tmaze_ba_graph_is_joint_vertex(tmaze_t *M, dgraph_vertex_index_t v)
  { dgraph_vertex_count_t nvC = M->nt;  /* Number of center vertices. */
    dgraph_vertex_count_t nvWE = M->nxWE * M->ny; /* Number of W-joints. */
    dgraph_vertex_count_t nvSN = M->nx * M->nySN; /* Number of S-joints. */
    return ((v >= nvC) && (v < nvC + nvWE + nvSN));
  }
  
void tmaze_ba_graph_sort_edge_endpoints
  ( tmaze_t *M, 
    dgraph_vertex_index_t v1, 
    dgraph_vertex_index_t v2, 
    dgraph_vertex_index_t *vCP, 
    dgraph_vertex_index_t *vJP
  )
{ if (tmaze_ba_graph_is_center_vertex(M, v1) && tmaze_ba_graph_is_joint_vertex(M,v2))
      { (*vCP) = v1; (*vJP) = v2; }
  else if (tmaze_ba_graph_is_center_vertex(M,v2) && tmaze_ba_graph_is_joint_vertex(M,v1))
      { (*vCP) = v2; (*vJP) = v1; }
    else
      { /* Invalid edge: */ assert(FALSE); }
  }

void tmaze_ba_graph_cell_from_center_vertex
  ( tmaze_t *M, 
    dgraph_vertex_index_t v, 
    int *xP, 
    int *yP
  )
  { dgraph_vertex_count_t nvC = M->nt;  /* Number of center vertices. */
    if(v < nvC)
      { /* Get the tile index from the center vertex index: */
        tmaze_cell_index_t k = (tmaze_cell_index_t)v;  
        /* Get the coords from the tile index: */
        tmaze_cell_position(k, M->nx, M->ny, xP, yP);
      }
    else
      { demand(FALSE, "vertex is not a tile center"); }
  }

void tmaze_ba_graph_get_cell_and_direction_from_joint_vertex
  ( tmaze_t *M, 
    dgraph_vertex_index_t v, 
    int *xP, 
    int *yP,
    tmaze_dir_t *dirP
  )
  {
    dgraph_vertex_count_t nvC = M->nt;  /* Number of center vertices. */
    dgraph_vertex_count_t nvWE = M->nxWE * M->ny; /* Number of W-joints. */
    dgraph_vertex_count_t nvSN = M->nx * M->nySN; /* Number of S-joints. */

    if (v < nvC)
      { demand(FALSE, "vertex is not a joint"); }
    else if (v < nvC + nvWE)
      { /* Vertex is a W-joint: */
        /* Get the tile index from the joint vertex: */
        tmaze_cell_index_t k = (tmaze_cell_index_t)(v - nvC);  
        /* Get the coords from the tile index: */
        tmaze_cell_position(k, M->nxWE, M->ny, xP, yP);
        (*dirP) = tmaze_dir_W;
      }
    else if (v < nvC + nvWE + nvSN)
      { /* Vertex is an S-joint: */
        /* Get the tile index from the joint vertex: */
        tmaze_cell_index_t k = (tmaze_cell_index_t)(v - nvC - nvWE);  
        /* Get the coords from the tile index: */
        tmaze_cell_position(k, M->nx, M->nySN, xP, yP);
        (*dirP) = tmaze_dir_S;
      }
    else
      { demand(FALSE, "invalid joint vertex index"); }
  }

void tmaze_ba_graph_get_cell_and_direction_from_edge
  ( tmaze_t *M, 
    dgraph_vertex_index_t v1, 
    dgraph_vertex_index_t v2,
    int *xP, 
    int *yP, 
    tmaze_dir_t *dirP
  )
  {
    /* Sort the edge's endpoints: */
    dgraph_vertex_index_t vC, vJ;
    tmaze_ba_graph_sort_edge_endpoints(M, v1, v2, &vC, &vJ);
    
    /* Get the tile column and row {xC,yC} from the center vertex: */
    int xC, yC;
    tmaze_ba_graph_cell_from_center_vertex(M, vC, &xC, &yC);
    assert(xC < M->nx);
    assert(yC < M->ny);
    
    /* Get the tile column and row {xJ,yJ} from the joint vertex: */
    int xJ, yJ;
    tmaze_dir_t dirJ;
    tmaze_ba_graph_get_cell_and_direction_from_joint_vertex(M, vJ, &xJ, &yJ, &dirJ);
    
    /* Get the edge's direction: */
    tmaze_dir_t dir;
    if ((xJ == xC) && (yJ == yC))
      { /* West or South joint: */
        dir = dirJ;
      }
    else if (xJ == xC)
      { /* Must be a North joint: */
        assert(dirJ == tmaze_dir_S);
        assert(yJ == ((yC + 1) % M->nySN));
        dir = tmaze_dir_N;
      }
    else if (yJ == yC)
      { /* Must be an East joint: */
        assert(dirJ == tmaze_dir_W);
        assert(xJ == ((xC + 1) % M->nxWE));
        dir = tmaze_dir_E;
      }
      
    /* Consistency checking: */
    assert(vC == tmaze_ba_graph_center_vertex(M, xC, yC));
    assert(vJ == tmaze_ba_graph_joint_vertex(M, xC, yC, dir));
    
    /* Return results: */
    (*xP) = xC; 
    (*yP) = yC;
    (*dirP) = dir;
  }

dgraph_t tmaze_ba_graph_make(tmaze_t *M)
  {
    demand(M->nx >= 2, "maze width must be at least 2");
    demand(M->ny >= 2, "maze height must be at least 2");

    demand(M->nySN <= dgraph_MAX_VERTEX_COUNT/M->nxWE/3, "too many vertices");
    dgraph_vertex_count_t nvC = M->nt;    /* Number of center vertices. */
    dgraph_vertex_count_t nvWE = M->nxWE * M->ny;  /* Number of W-joints. */
    dgraph_vertex_count_t nvSN = M->nx * M->nySN;  /* Number of S-joints. */
    dgraph_vertex_count_t nv = nvC + nvWE + nvSN;  /* Number of vertices in graph. */

    demand(nvC <= dgraph_MAX_EDGE_COUNT/6, "too many edges");
    dgraph_edge_count_t ne = 6*nvC;  /* Number of edges in graph. */

    /* Create an undirected graph with {nv} vertices {ne} undefined edges: */
    dgraph_t G = dgraph_new(nv, nv, ne);

    /* Define the edges: */
    dgraph_edge_count_t me = 0; /* Number of edges added so far. */
    int x; /* Tile x index (column), in {0..nx-1}. */
    int y; /* Tile y index (row), in {0..ny-1}. */
    for (y = 0; y < M->ny; y++)
      { for (x = 0; x < M->nx; x++)
          { /* Get the tile's index: */
            tmaze_cell_index_t k = tmaze_cell_index(x, y, M->nx, M->ny);
            /* Vertex numbers for this tile: */
            dgraph_vertex_index_t vC = tmaze_ba_graph_center_vertex(M, x, y);
            dgraph_vertex_index_t vN = tmaze_ba_graph_joint_vertex(M, x, y, tmaze_dir_N);
            dgraph_vertex_index_t vS = tmaze_ba_graph_joint_vertex(M, x, y, tmaze_dir_S);
            dgraph_vertex_index_t vW = tmaze_ba_graph_joint_vertex(M, x, y, tmaze_dir_W);
            dgraph_vertex_index_t vE = tmaze_ba_graph_joint_vertex(M, x, y, tmaze_dir_E);
            /* Get the orientation of the tile's T-road: */
            tmaze_ba_tile_t tile = M->tile[k];
            /* Add the edges that are not excluded: */
            if (tile != tmaze_ba_tile_N) { me = dgraph_add_undirected_edge(&G, vC, vN, me); }
            if (tile != tmaze_ba_tile_S) { me = dgraph_add_undirected_edge(&G, vC, vS, me); }
            if (tile != tmaze_ba_tile_W) { me = dgraph_add_undirected_edge(&G, vC, vW, me); }
            if (tile != tmaze_ba_tile_E) { me = dgraph_add_undirected_edge(&G, vC, vE, me); }
          }
      }
    assert(me == ne);
    dgraph_trim(&G, ne);
    return G;
  }
 
void tmaze_ba_graph_count_components_by_size
  ( dgraph_t *G, 
    tmaze_t *M, 
    tmaze_size_t *msP,
    tmaze_comp_count_t **ctP
  )
  {
    tmaze_cell_count_t nt = M->nt;  /* Number of tiles in maze. */
    
    assert(G->cols == G->rows);
    dgraph_vertex_count_t nv = G->cols;  /* Total number of vertices. */
    
    /* Find a maximally connected spanning forest of {G}: */
    dgraph_vertex_index_t *parent = dgraph_find_spanning_forest(G);

    /* For each root {v}, compute the number {tsize[v]} of tile centers in its tree: */
    dgraph_vertex_count_t *tsize = notnull(malloc(nv*sizeof(dgraph_vertex_count_t)), "no mem");
    dgraph_vertex_index_t v;
    for (v = 0; v < nv; v++) { tsize[v] = 0; }
    for (v = 0; v < nv; v++) 
      { dgraph_vertex_index_t r = dgraph_find_root(parent, v); 
        assert(r < nv);
        if (tmaze_ba_graph_is_center_vertex(M, v)) { tsize[r]++; }
      }

    /* Compute the max possible number of tile centers {ms} in any component: */
    tmaze_size_t ms = 0;
    for (v = 0; v < nv; v++) { if (tsize[v] > ms) { ms = tsize[v]; } }
    assert(ms <= nt); /* Since there is one tile center vertex per tile. */
    
    /* Count components by number of tile centers in tree: */
    tmaze_comp_count_t *ct = notnull(malloc((ms+1)*sizeof(dgraph_vertex_count_t)), "no mem");
    tmaze_size_t sz;
    for (sz = 0; sz <= ms; sz++) { ct[sz] = 0; }
    for (v = 0; v < nv; v++) 
      { if (parent[v] == v)
          { /* Vertex {v} is a root; get the number {sz} of tile centers in its tree: */
            sz = tsize[v];
            assert(sz <= ms);
            ct[sz]++;
          }
      }

    /* Cleanup and return: */
    free(tsize);
    free(parent);
    (*msP) = ms;
    (*ctP) = ct;
  }

double tmaze_ba_graph_predicted_comp_count(tmaze_t *M, tmaze_size_t size)
  {
    int nx = M->nx;
    int ny = M->ny;
    int nt = M->nt; /* Number of cells in maze. */
    
    /* Compute vertex counts for the graph: */
    dgraph_vertex_count_t nvC = nt;            /* Num of cell vertices. */
    dgraph_vertex_count_t nvWE = M->nxWE * ny; /* Num of vert joints. */
    dgraph_vertex_count_t nvSN = nx * M->nySN; /* Num of horz joints. */
    dgraph_vertex_count_t nvJ = nvWE + nvSN;   /* Num of joint vertices. */
    
    tmaze_size_t ms = nt;  /* Max size of any component. */
    
    double exp_ct;  /* Statistically expected count of plain components in torus. */

    if (size == 0) 
      { exp_ct = (1/16.0)*nvJ; /* Expected num of isolated joint vertices. */
      }
    else if (size == 1)
      { exp_ct = (1/64.0)*nvC; /* Expected num of 1-cell components. */
      }
    else if (size == 2)
      { exp_ct = (9/2048.0)*nvC; /* Expected num of 2-cell components. */
      }
    else if (size == 3)
      { exp_ct = (19/16384.0)*nvC; /* Expected num of 3-cell components. */
      }
    else if (size > ms)
      { exp_ct = 0; }
    else
      { exp_ct = NAN; }

    return exp_ct;
  }
