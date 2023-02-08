/* See {tmaze_bb.h} */
/* Last edited on 2009-11-09 23:39:59 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <math.h>

#include <bool.h>
#include <affirm.h>
#include <dgraph.h>

#include <tmaze_bb.h>
#include <tmaze_bb_graph.h>

dgraph_vertex_index_t tmaze_bb_graph_vertex(tmaze_t *M, int x, int y, tmaze_dir_t dir)
  { dgraph_vertex_count_t nvWE = M->nxWE * M->ny; /* Number of W-joints. */
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
        { /* The West joint of a tile is numbered {k},
            where {k} is the tile index computed as if the maze
            had {nxWE} columns and {ny} rows. */
          tmaze_cell_index_t k = tmaze_cell_index(x, y, M->nxWE, M->ny);
          assert(k < nvWE);
          return k;
          break;
        }
      case tmaze_dir_S:
      dir_S:
        { /* The South joint of a tile is numbered {nvW + k},
            where {k} is the tile index computed as if the maze
            had {nx} columns and {nySN} rows. */
          tmaze_cell_index_t k = tmaze_cell_index(x, y, M->nx, M->nySN);
          assert(k < nvSN);
          return nvWE + k;
          break;
        }
      default: assert(FALSE);
      }
  }

void tmaze_bb_graph_get_cell_and_direction_from_vertex
  ( tmaze_t *M, 
    dgraph_vertex_index_t v, 
    int *xP, 
    int *yP,
    tmaze_dir_t *dirP
  )
  {
    dgraph_vertex_count_t nvWE = M->nxWE * M->ny; /* Number of W-joints. */
    dgraph_vertex_count_t nvSN = M->nx * M->nySN; /* Number of S-joints. */

    if (v < nvWE)
      { /* Vertex is a W-joint: */
        /* Get the tile index from the joint vertex: */
        tmaze_cell_index_t k = (tmaze_cell_index_t)v;  
        /* Get the coords from the tile index: */
        tmaze_cell_position(k, M->nxWE, M->ny, xP, yP);
        (*dirP) = tmaze_dir_W;
      }
    else if (v < nvWE + nvSN)
      { /* Vertex is an S-joint: */
        /* Get the tile index from the joint vertex: */
        tmaze_cell_index_t k = (tmaze_cell_index_t)(v - nvWE);  
        /* Get the coords from the tile index: */
        tmaze_cell_position(k, M->nx, M->nySN, xP, yP);
        (*dirP) = tmaze_dir_S;
      }
    else
      { demand(FALSE, "invalid joint vertex index"); }
  }

void tmaze_bb_graph_get_cell_and_dirs_from_edge
  ( tmaze_t *M, 
    dgraph_vertex_index_t v1, 
    dgraph_vertex_index_t v2,
    int *xP, 
    int *yP, 
    tmaze_dir_t *dir1P,
    tmaze_dir_t *dir2P
  )
  {
    /* Check whether the problem is soluble: */
    demand((! M->torus) || ((M->nx >= 2) && (M->ny >= 2)), "ambiguous");
    
    /* Get the tile column, row, and side {x1,y1,dir1} from {v1}: */
    int x1, y1;
    tmaze_dir_t dir1;
    tmaze_bb_graph_get_cell_and_direction_from_vertex(M, v1, &x1, &y1, &dir1);
    assert(x1 < M->nx); assert(y1 < M->ny);
    assert((dir1 == tmaze_dir_W) || (dir1 == tmaze_dir_S));
    
    /* Get the tile column, row, and side {x2,y2,dir2} from {v2}: */
    int x2, y2;
    tmaze_dir_t dir2;
    tmaze_bb_graph_get_cell_and_direction_from_vertex(M, v2, &x2, &y2, &dir2);
    assert(x2 < M->nx); assert(y2 < M->ny);
    assert((dir2 == tmaze_dir_W) || (dir2 == tmaze_dir_S));
    
    /* Determine the actual indices {x,y} and fix directions {dir1,dir2}: */
    int x, y;
    if (x2 == tmaze_cell_row_col_inc(x1, 1, M->nx, M->torus))
      { /* Column {x2} is East of {x1}: */
        assert(dir2 == tmaze_dir_W);
        x = x1; dir2 = tmaze_dir_E;
      }
    else if (x2 == tmaze_cell_row_col_inc(x1, -1, M->nx, M->torus))
      { /* Column {x2} is West of {x1}: */
        assert(dir1 == tmaze_dir_W);
        x = x2; dir1 = tmaze_dir_E;
      }
    else
      { assert(x1 == x2);
        x = x1;
      }
    if (y2 == tmaze_cell_row_col_inc(y1, 1, M->ny, M->torus))
      { /* Column {y2} is North of {y1}: */
        assert(dir2 == tmaze_dir_S);
        y = y1; dir2 = tmaze_dir_N;
      }
    else if (y2 == tmaze_cell_row_col_inc(y1, -1, M->ny, M->torus))
      { /* Column {y2} is South of {y1}: */
        assert(dir1 == tmaze_dir_S);
        y = y2; dir1 = tmaze_dir_N;
      }
    else
      { assert(y1 == y2);
        y = y1;
      }
      
    /* Consistency checking: */
    assert(v1 == tmaze_bb_graph_vertex(M, x, y, dir1));
    assert(v2 == tmaze_bb_graph_vertex(M, x, y, dir2));
      
    /* Return results: */
    (*xP) = x; 
    (*yP) = y;
    (*dir1P) = dir1;
    (*dir2P) = dir2;
  }

dgraph_t tmaze_bb_graph_make(tmaze_t *M)
  {
    demand(M->nx >= 2, "maze width must be at least 2");
    demand(M->ny >= 2, "maze height must be at least 2");

    /* There is one vertex on each tile side. Better safe than sorry: */
    demand(M->nySN <= dgraph_MAX_VERTEX_COUNT/2/M->nxWE, "too many vertices");
    dgraph_vertex_count_t nvWE = M->nxWE * M->ny;  /* Number of W joint vertices. */
    dgraph_vertex_count_t nvSN = M->nx * M->nySN;  /* Number of S joint vertices. */
    dgraph_vertex_count_t nv = nvWE + nvSN;  /* Number of vertices in graph. */

    /* There are two edges in each tile: */
    demand(M->nx <= dgraph_MAX_EDGE_COUNT/2/M->ny, "too many edges");
    dgraph_edge_count_t ne = 4*M->nx*M->ny;  /* Number of directed edges in graph. */

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
            dgraph_vertex_index_t vN = tmaze_bb_graph_vertex(M, x, y, tmaze_dir_N);
            dgraph_vertex_index_t vS = tmaze_bb_graph_vertex(M, x, y, tmaze_dir_S);
            dgraph_vertex_index_t vW = tmaze_bb_graph_vertex(M, x, y, tmaze_dir_W);
            dgraph_vertex_index_t vE = tmaze_bb_graph_vertex(M, x, y, tmaze_dir_E);
            /* Get the orientation of the tile: */
            tmaze_bb_tile_t tile = M->tile[k];
            /* Add the diagonal roads: */
            if (tile == tmaze_bb_tile_BLIP)
              { me = dgraph_add_undirected_edge(&G, vN, vE, me);
                me = dgraph_add_undirected_edge(&G, vS, vW, me);
              }
            else
              { me = dgraph_add_undirected_edge(&G, vN, vW, me);
                me = dgraph_add_undirected_edge(&G, vS, vE, me);
              }
          }
      }
    assert(me == ne);
    dgraph_trim(&G, ne);
    return G;
  }

void tmaze_bb_graph_count_components_by_size
  ( dgraph_t *G, 
    tmaze_t *M, 
    tmaze_size_t *msP,
    tmaze_comp_count_t **ctP
  )
  {
    tmaze_cell_count_t nt = M->nt;  /* Number of tiles in maze. */
    
    /* Each tile should have 4 directed edges: */
    assert(G->ents == 4*nt);
    
    /* The matrix must be square: */
    assert(G->cols == G->rows);
    dgraph_vertex_count_t nv = G->cols;  /* Total number of vertices. */
    
    bool_t debug = FALSE;
    
    /* Find a maximally connected spanning forest of {G}: */
    dgraph_vertex_index_t *parent = dgraph_find_spanning_forest(G);

    /* For each root {v}, compute the number {tsize[v]} of vertices in its tree: */
    dgraph_vertex_count_t *tsize = notnull(malloc(nv*sizeof(dgraph_vertex_count_t)), "no mem");
    dgraph_vertex_index_t v;
    for (v = 0; v < nv; v++) { tsize[v] = 0; }
    for (v = 0; v < nv; v++) 
      { dgraph_vertex_index_t r = dgraph_find_root(parent, v); 
        if(debug) { fprintf(stderr, " %3d --> %3d\n", v, r); }
        assert(r < nv);
        tsize[r]++;
      }

    /* Compute the max possible number of vertices {ms} in any component: */
    tmaze_size_t ms = 0;
    for (v = 0; v < nv; v++) { if (tsize[v] > ms) { ms = tsize[v]; } }
    
    assert(ms <= 2*nt);
      /* Since the vertices have undirected degree 2 (interior sides)
        or 1 (boundary sides), each component is either a simple cycle
        or a simple path. If it is a cycle, it has at most {2*M.nt}
        edges and therefore at most {2*nt} vertices. If the component
        is a open path, it means that {M.torus} is false, and
        therefore at least 4 of its edges end on the maze boundary.
        Therefore the path may have at most {2*nt-2} edges, that is,
        at most {2*nt-1} vertices. */
    
    /* Count components by number of vertices: */
    tmaze_comp_count_t *ct = notnull(malloc((ms + 1)*sizeof(dgraph_vertex_count_t)), "no mem");
    tmaze_size_t sz;
    for (sz = 0; sz <= ms; sz++) { ct[sz] = 0; }
    for (v = 0; v < nv; v++) 
      { if (parent[v] == v)
          { /* Vertex {v} is a root; get the number {sz} of vertices in its tree: */
            sz = tsize[v];
            assert(sz <= ms);
            ct[sz]++;
          }
      }
      
    /* Consistency checks: */
    for (sz = 0; sz <= ms; sz++)
      { if (sz <= 1)
          { /* Blip-Blop maze graphs have no 0-vertex or 1-vertex components. */
            assert(ct[sz] == 0);
          }
        else if ((sz % 2) == 1)
          { /* Odd-sz components are expected only along edges or wrapping the torus: */
            if (M->torus && ((M->nx % 2) == 0) && ((M->ny % 2) == 0)) { assert(ct[sz] == 0); }
          }
        else if (sz == 2)
          { /* There are no 2-vertex components in the torus: */
            if (M->torus) { assert(ct[sz] == 0); }
          }
        else if (sz == 6)
          { /* There should be no plain components with 6 vertices: */
            if (M->torus && (M->nx > 3) && (M->ny > 3)) { assert(ct[sz] == 0); }
          }
        else if (sz == 10)
          { /* There should be no plain components with 10 vertices: */
            if (M->torus && (M->nx > 5) && (M->ny > 5)) { assert(ct[sz] == 0); }
          }
      }

    /* Cleanup and return: */
    free(tsize);
    free(parent);
    (*msP) = ms;
    (*ctP) = ct;
  }
 
double tmaze_bb_graph_predicted_comp_count(tmaze_t *M, tmaze_size_t size)
  {
    int nx = M->nx;
    int ny = M->ny;
    int nt = M->nt; /* Number of cells in maze. */
    
    /* Compute vertex counts for the graph: */
    dgraph_vertex_count_t nvC = nt;            /* Num of cell vertices. */
    dgraph_vertex_count_t nvWE = M->nxWE * ny; /* Num of vert joints. */
    dgraph_vertex_count_t nvSN = nx * M->nySN; /* Num of horz joints. */
    dgraph_vertex_count_t nvJ = nvWE + nvSN;   /* Num of joint vertices. */
    dgraph_vertex_count_t nv = nvC + nvJ;      /* Tot num of vertices. */

    tmaze_size_t ms = nv;  /* Max num of vertices in any component. */
    
    double exp_ct;  /* Statistically expected count of plain components in torus. */

    if (size <= 1)
      { /* The Blip-Blop maze has no 0-vertex or 1-vertex components. */
        exp_ct = 0;
      }
    else if ((size % 2) == 1)
      { /* Odd-size components are expected only along edges or wrapping the torus: */
        exp_ct = 0; /* Ignore edge effects. */
      }
    else if (size == 2)
      { /* There are no 2-vertex components in the torus: */
        exp_ct = 0;
      }
    else if (size == 4)
      { /* A plain 4-vertex component is 4 tiles ((/,\),(\,/)): {1/2^4}: */
        exp_ct = nt*1.0/16.0;
      }
    else if (size == 6)
      { /* There should be no plain components with 6 vertices: */
        exp_ct = 0;
      }
    else if (size == 8)
      { /* A plain 8-vertex component is a Dumbbell: 
          Dumbbell ({3×3}): ((*,/,\),(/,/,/),(\,/,*)) or its transposal -- {2/2^7}: */
        exp_ct = nt*2/128.0;
      }
    else if (size == 10)
      { /* There should be no plain components with 10 vertices: */
        exp_ct = 0; 
      }
    else if (size == 12)
      { /* A whole nonwrapping 12-vertex component is a bottle or fly or small cross: */
        /* Bottle ({4×4}):  ((*,*,/,\),(*,/,/,/),(/,/,/,*),(\,/,*,*)) or its transp: {2/2^{10}}. */
        /* Fly ({3×4}):     ((*,/,\,*),(/,/,\,\),(\,/,\,/)) or its vflip or transp: {4/2^{10}}. */
        /* SmCross ({4×4}): ((*,/,\,*),(/,/,\,\),(\,\,/,/),(*,\,/,*)): {1/2^{12}}. */
        exp_ct = nt*(8+16+1)/4096.0;
      }
    else if (size > ms)
      { exp_ct = 0; }
    else
      { exp_ct = NAN; }

    return exp_ct;
  }


void tmaze_bb_graph_get_edge_components
  ( dgraph_t *G, 
    tmaze_t *M, 
    dgraph_edge_count_t *neP,
    tmaze_comp_count_t *ncP,
    tmaze_size_t **szP,
    dgraph_vertex_index_t **rtP
  )
  {
    tmaze_cell_count_t nt = M->nt;    /* Number of tiles in maze. */
    
    /* Each tile should have 4 directed edges: */
    assert(G->ents == 4*nt);
    
    /* The matrix must be square: */
    assert(G->cols == G->rows);
    dgraph_vertex_count_t nv = G->cols;  /* Total number of vertices. */
    
    /* Find a maximally connected spanning forest of {G}: */
    dgraph_vertex_index_t *parent = dgraph_find_spanning_forest(G);
    
    /* For each root {v}, compute the number {sz[v]} of vertices in its tree: */
    tmaze_size_t *sz = notnull(malloc(nv*sizeof(dgraph_vertex_count_t)), "no mem");
    dgraph_vertex_index_t v;
    for (v = 0; v < nv; v++) { sz[v] = 0; }
    for (v = 0; v < nv; v++) 
      { dgraph_vertex_index_t r = dgraph_find_root(parent, v); 
        assert(r < nv);
        sz[r]++;
      }
    /* Copy {sz[r]} of each root {r} to all vertices in its component, and count components: */
    tmaze_comp_count_t nc = 0;
    for (v = 0; v < nv; v++) 
      { dgraph_vertex_index_t r = dgraph_find_root(parent, v); 
        assert(r < nv);
        if (r == v) 
          { /* {v} is root: */ nc++; }
        else
          { /* {v} is non-root: */ sz[v] = sz[r]; }
      }
      
    /* Find the root vertex {rt[e]} of the component that contains each edge {e} of {G}: */
    dgraph_edge_count_t ne = 2*nt; /* Number of undirected edges in graph. */
    dgraph_vertex_index_t *rt = notnull(malloc(ne*sizeof(dgraph_vertex_index_t)), "no mem");
    /* Fill {rt} with invalid values, for consistency checking: */
    int e;
    for (e = 0; e < ne; e++) { rt[e] = nv; }
    /* Scan the edges of {G} and fill {rt}: */
    dgraph_edge_index_t p;
    for (p = 0; p < G->ents; p++)
      { /* We should have stored only TRUE entries: */
        assert(G->e[p].val); 
        /* Get the origin and destination vertices of edge {p}: */
        dgraph_vertex_index_t v1 = G->e[p].row; assert(v1 < nv);
        dgraph_vertex_index_t v2 = G->e[p].col; assert(v2 < nv);
        /* Get the root of the component that contains {v1}: */
        dgraph_vertex_index_t r1 = dgraph_find_root(parent, v1); assert(r1 < nv);
        dgraph_vertex_index_t r2 = dgraph_find_root(parent, v2); assert(r2 < nv);
        assert(r1 == r2);
        /* Find the edge {(v1,v2)} in the maze: */
        int x, y;
        tmaze_dir_t dir1, dir2;
        tmaze_bb_graph_get_cell_and_dirs_from_edge(M, v1, v2, &x, &y, &dir1, &dir2);
        /* Compute the index {e} of the edge in {rt}: */
        tmaze_cell_index_t k = tmaze_cell_index(x, y, M->nx, M->ny);
        if ((dir1 == tmaze_dir_S) || (dir2 == tmaze_dir_S))
          { e = 2*k; }
        else if ((dir1 == tmaze_dir_N) || (dir2 == tmaze_dir_N))
          { e = 2*k + 1; }
        else
          { assert(FALSE); }
        /* Save the root of {e}'s component: */
        rt[e] = r1;
      }
      
    /* Consistency check: */
    for (e = 0; e < ne; e++) { assert(rt[e] < nv); }
    
    /* Cleanup and return: */
    free(parent);
    (*neP) = ne;
    (*ncP) = nc;
    (*szP) = sz;
    (*rtP) = rt;
  }
