/* See {tmaze_ba.h} */
/* Last edited on 2023-02-03 23:41:35 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <assert.h>

#include <bool.h>
#include <affirm.h>
#include <epswr.h>

#include <tmaze_ba.h>

/* INTERNAL PARAMETERS */

#define tmaze_ba_ROAD_RADIUS (0.5/3.0)
  /* Half-width of the T-road in Brasilia Airport mazes. */

#define tmaze_ba_CROSS_RADIUS (0.5/6.0)
  /* Half-width of the arms of the cross used to indicate 
    empty tileset in Brasilia Airport maze patterns. */

/* IMPLEMENTATIONS */

tmaze_t tmaze_ba_make(int32_t nx, int32_t ny, bool_t torus)
  {
    tmaze_t M = tmaze_make(nx, ny, torus);
    return M;
  }

void tmaze_ba_random(tmaze_t *M, int32_t seed)
  {
    tmaze_random(M, seed, tmaze_ba_tileset_FOUR);
  }

tmaze_t tmaze_ba_from_string(char *string, bool_t torus)
  { 
    /* Get the maze size from the string: */
    int32_t nx, ny;
    tmaze_size_from_string(string, &nx, &ny);

    /* Now allocate the maze: */
    tmaze_t M = tmaze_make(nx, ny, torus);
    
    /* Fill {M.tile} from {string}, from TOP TO BOTTOM: */
    int32_t x; /* Cell x index (column), in {0..nx-1}. */
    int32_t y; /* Cell y index (row), in {0..ny-1}. */
    char *p = string;
    for (y = ny-1; y >= 0; y--)
      { for (x = 0; x < nx; x++)
          { int32_t k = tmaze_cell_index(x, y,nx, ny);
            char c = (*p); p++;
            demand(c != 0, "string is too short");
            tmaze_tile_t t;
            switch (c)
              {
                case 'S': t = tmaze_ba_tile_S; break;
                case 'N': t = tmaze_ba_tile_N; break;
                case 'W': t = tmaze_ba_tile_W; break;
                case 'E': t = tmaze_ba_tile_E; break;
                default: 
                  fprintf(stderr, "** invalid char in string = '%c'\n", c);
                  exit(1);
              }
            M.tile[k] = t;
          }
        char s = (*p); p++;
        if (s != ';')
          { fprintf(stderr, "** expecting ';', found '%c'\n", s);
            exit(1);
          }
      }
    if ((*p) != 0)
      { fprintf(stderr, "** excess char in string = '%c'\n", (*p));
        exit(1);
      }
    return M;
  }

tmaze_pattern_t tmaze_ba_pattern_from_string(char *string, bool_t torus, bool_t sub, int32_t val)
  { 
    /* Get the pattern size from the string: */
    int32_t nx, ny;
    tmaze_size_from_string(string, &nx, &ny);
          
    /* Now allocate the pattern: */
    tmaze_pattern_t M = tmaze_make_pattern(nx, ny, torus, sub, val);
    
    /* Fill {M.tset} from {string}, from TOP TO BOTTOM: */
    int32_t x; /* Cell x index (column), in {0..nx-1}. */
    int32_t y; /* Cell y index (row), in {0..ny-1}. */
    char *p = string;
    for (y = ny-1; y >= 0; y--)
      { for (x = 0; x < nx; x++)
          { int32_t k = tmaze_cell_index(x, y, nx, ny);
            char c = (*p); p++;
            demand(c != 0, "string is too short");
            tmaze_tileset_t S;
            switch (c)
              {
                case 'S': S = (1 << tmaze_ba_tile_S); break;
                case 'N': S = (1 << tmaze_ba_tile_N); break;
                case 'W': S = (1 << tmaze_ba_tile_W); break;
                case 'E': S = (1 << tmaze_ba_tile_E); break;
                case '*': S = tmaze_ba_tileset_FOUR; break;
                case 'x': S = tmaze_tileset_NONE; break;
                default: 
                  fprintf(stderr, "** invalid char in string = '%c'\n", c);
                  exit(1);
              }
            M.tset[k] = S;
          }
        char s = (*p); p++;
        if (s != ';')
          { fprintf(stderr, "** expecting ';', found '%c'\n", s);
            exit(1);
          }
      }
    if ((*p) != 0)
      { fprintf(stderr, "** excess char in string = '%c'\n", (*p));
        exit(1);
      }
    return M;
  }

tmaze_size_t tmaze_ba_tot_size(tmaze_t *M)
  { 
    return M->nt;
  }

tmaze_size_t tmaze_ba_max_comp_size(tmaze_t *M)
  { 
    return M->nt;
  }

void tmaze_ba_print(FILE *wr, tmaze_t *M)
  {
    int32_t x; /* Cell x index (column), in {0..nx-1}. */
    int32_t y; /* Cell y index (row), in {0..ny-1}. */
    if (M->torus)
      { fprintf(wr, "%4s ", "");
        for (x = 0; x < M->nx; x++) { fputc('|', wr); }
        fputc('\n', wr);
      }
    /* Print rows in the correct order: */
    for (y = M->ny - 1; y >= 0; y--)
      { fprintf(wr, "%4d ", y);
        for (x = 0; x < M->nx; x++)
          { /* Compute the cell's index in the array: */
            tmaze_cell_index_t k = tmaze_cell_index(x, y, M->nx, M->ny);
            /* Print the tile's orientation: */
            fputc(tmaze_ba_tile_chars[M->tile[k]], wr);
          }
        if (M->torus) { fputc('-', wr); } 
        fputc('\n', wr);
      }
  }

void tmaze_ba_plot_tile
  ( epswr_figure_t *eps,
    tmaze_ba_tile_t tile,
    int32_t x,
    int32_t y,
    bool_t open[],
    bool_t fill,
    bool_t draw
  )
  { int32_t nc = 8; /* Number of corners in the T-road's outline. */
    double xp[nc], yp[nc]; /* X and Y coordinates of corners of T-road. */
    tmaze_ba_path_tee_make(xp, yp, nc, tmaze_ba_ROAD_RADIUS); 

    /* Rotate path from the canonical {tmaze_tile_N} orientation to {tile}'s: */
    tmaze_dir_t dir = (tmaze_dir_t)tile;
    int32_t k = dir - tmaze_dir_N;
    tmaze_path_rotate(xp, yp, nc, k);
    
    /* Translate the path from the origin to the cell's center: */
    tmaze_path_translate(xp, yp, nc, x + 0.5, y + 0.5);
    
    if (fill) { epswr_polygon(eps, TRUE, xp, yp, nc, TRUE, FALSE, TRUE); }
    
    if (draw)
      { /* Determine which legs of the T-road are open: */
        bool_t cL = open[(dir + 1) % 4]; /* Left arm of 'T'. */
        bool_t cM = open[(dir + 2) % 4]; /* Base of 'T'. */
        bool_t cR = open[(dir + 3) % 4]; /* Right arm of 'T'. */
        /* Draw the outlines of the T-road inside the cell: */
        epswr_segment(eps, xp[0], yp[0], xp[1], yp[1]);
        epswr_segment(eps, xp[1], yp[1], xp[2], yp[2]);
        if (! cR) { epswr_segment(eps, xp[2], yp[2], xp[3], yp[3]); }
        epswr_segment(eps, xp[3], yp[3], xp[4], yp[4]);
        if (! cL) { epswr_segment(eps, xp[4], yp[4], xp[5], yp[5]); }
        epswr_segment(eps, xp[5], yp[5], xp[6], yp[6]);
        epswr_segment(eps, xp[6], yp[6], xp[7], yp[7]);
        if (! cM) { epswr_segment(eps, xp[7], yp[7], xp[0], yp[0]); }
      }
  }
  
void tmaze_ba_path_tee_make
  ( double xp[], 
    double yp[], 
    int32_t nc, 
    double r
  )
  { /* This procedure generates 8 corners, so: */
    assert(nc == 8); 
    int32_t kc = 0; /* Counts corners appended to the T-road's outline. */
    xp[kc] = r;    yp[kc] = -0.5; kc++;
    xp[kc] = r;    yp[kc] = -r;   kc++;
    xp[kc] = +0.5; yp[kc] = -r;   kc++;
    xp[kc] = +0.5; yp[kc] = +r;   kc++;
    xp[kc] = -0.5; yp[kc] = +r;   kc++;
    xp[kc] = -0.5; yp[kc] = -r;   kc++;
    xp[kc] = -r;   yp[kc] = -r;   kc++;
    xp[kc] = -r;   yp[kc] = -0.5; kc++;
    assert(kc == nc);
  }
  
void tmaze_ba_plot_maze
  ( char *dir,
    char *name,
    tmaze_t *M,
    bool_t cells,  /* TRUE paints the background of each tile. */
    bool_t grids,  /* TRUE paints the cell gridlines. */
    double tsize   /* Cell size (mm) */
  )
  {
    int32_t val = 0;         /* Number of ports to show on each side. */
    double margin = 2.0; /* Figure margin in mm. */
    
    auto void plot_tile(epswr_figure_t *eps, tmaze_t *M, int32_t x, int32_t y, bool_t fill, bool_t draw);
      /* Plots the tileset in column {x} and  row {y} of maze pattern {M}, with given style. */

    void plot_tile(epswr_figure_t *eps, tmaze_t *M, int32_t x, int32_t y, bool_t fill, bool_t draw)
      { tmaze_ba_plot_maze_tile(eps, M, x, y, fill, draw); }

    tmaze_plot_maze
      ( dir, name,
        M,
        cells,
        grids,
        &plot_tile,
        val,
        tsize,
        margin
      );
  }
  
void tmaze_ba_plot_maze_tile
  ( epswr_figure_t *eps,
    tmaze_t *M, 
    int32_t x, 
    int32_t y,
    bool_t fill, 
    bool_t draw
  )
  {
    /* Get the cell index: */
    tmaze_cell_index_t k = tmaze_cell_index(x, y, M->nx, M->ny);
    /* Determine which sides of the tile should be left open: */
    bool_t open[4];
    tmaze_ba_get_open_ends_of_maze_tile(M, k, open);
    /* Plot the tile: */
    tmaze_ba_plot_tile(eps, M->tile[k], x, y, open, fill, draw);
  }
 
void tmaze_ba_plot_pattern_tileset
  ( epswr_figure_t *eps,
    tmaze_pattern_t *M, 
    int32_t x, 
    int32_t y,
    bool_t fill, 
    bool_t draw
  )
  {
    /* Get the cell's index: */
    tmaze_cell_index_t k = tmaze_cell_index(x, y, M->nx, M->ny);
    /* Get the set of allowed tiles {S} for this cell: */
    tmaze_tileset_t S = M->tset[k];

    if (tmaze_tileset_contains(S, tmaze_ba_tileset_FOUR))
      { /* Free cell -- leave it blank: */ }
    else if (S == tmaze_tileset_NONE)
      { /* Invalid cell -- paint a square dot: */
        tmaze_ba_plot_cross(eps, x, y, fill, draw);
      }
    else 
      { /* Must be a singleton set, get its direction: */
        tmaze_tile_t tile = tmaze_min_tile(S);
        tmaze_tileset_t S1 = (1 << tile);
        if (S != S1)
          { fprintf(stderr, "** S = %02x\n", S);
            fprintf(stderr, "** tile = %d  (1<<tile) = %02x\n", tile, S1);
            demand(S == (1 << tile), "can plot only trivial, ALL, or NONE tilesets");
          }
        /* Determine which sides of the cell should be left open: */
        bool_t open[4];
        tmaze_ba_get_open_ends_of_pattern_tileset(M, k, open);
        /* Plot the tile: */
        tmaze_ba_plot_tile(eps, tile, x, y, open, fill, draw);
      }
  }
  
void tmaze_ba_get_open_ends_of_maze_tile
  ( tmaze_t *M, 
    tmaze_cell_index_t k, 
    bool_t open[]
  )
  { /* Enumerate the four neighbors of cell {k}: */
    int32_t id;
    for (id = 0; id < 4; id++)
      { /* Get the direction {dir} towards neighbor: */
        tmaze_dir_t dir = (tmaze_dir_t)id; 
        /* Get the index {kN} of the {dir}-neighbor of cell {t}: */
        tmaze_cell_index_t kN = 
          tmaze_cell_neighbor_index(k, dir, M->nx, M->ny, M->torus);
        /* Decide whether any end of tile {k} on side {dir} should be left open: */
        if (kN == tmaze_cell_index_NONE)
          { /* There is no neighbor in that direction: */
            open[id] = FALSE;
          }
        else
          { /* Get the orientation {dirN} of the neighboring tile: */
            tmaze_ba_tile_t tileN = M->tile[kN];
            tmaze_dir_t dirN = (tmaze_dir_t)tileN;
            /* Leave end open unless neighbor is oriented opposite to {id}: */
            tmaze_dir_t opp = (tmaze_dir_t)((id + 2) % 4); /* opposite of {dir}. */
            open[id] = (dirN != opp);
          }
      }
  }
    
void tmaze_ba_get_open_ends_of_pattern_tileset
  ( tmaze_pattern_t *M, 
    tmaze_cell_index_t k, 
    bool_t open[]
  )
  {
    /* Enumerate the four neighbors of cell {k}: */
    int32_t id;
    for (id = 0; id < 4; id++)
      { /* Get the direction {dir} towards the neighbor: */
        tmaze_dir_t dir = (tmaze_dir_t)id;
        /* Get the index {kN} of the {dir}-neighbor of cell {t}: */
        tmaze_cell_index_t kN = 
          tmaze_cell_neighbor_index(k, dir, M->nx, M->ny, M->torus);
        /* Decide whether any end of tile {k} on side {dir} should be left open: */
        if (kN == tmaze_cell_index_NONE)
          { /* There is no neighbor in that direction. */
            /* Road is possibly open if sub-board, closed if rectangle: */
            open[id] = M->sub;
          }
        else
          { /* Get the tileset {SN} of the neighboring cell: */
            tmaze_tileset_t SN = M->tset[kN];
            if (SN == tmaze_tileset_NONE)
              { /* Neighbor is invalid, may as well close the road: */
                open[id] = FALSE; 
              }
            else
              { /* Get one tile {tileN} from {SN}: */
                tmaze_tile_t tileN = tmaze_min_tile(SN);
                demand(tileN < 4, "invalid tile in tileset");
                /* Check whether {SN == {tileN}}: */
                if (SN == (1 << tileN))
                  { /* Yes; get the orientation {dirN} of {tileN}: */
                    tmaze_dir_t dirN = (tmaze_dir_t)tileN;
                    /* Leave end open unless neighbor is oriented opposite to {id}: */
                    tmaze_dir_t opp = (tmaze_dir_t)((id + 2) % 4); /* opp. of {dir}. */
                    open[id] = (dirN != opp);
                  }
                else
                  { /* Neighbor's tileset is neither empty nor singleton. */
                    /* Since it *may* connect, better leave the {dir}-leg open: */
                    open[id] = TRUE;
                  }
              }
          }
      }
  }   

void tmaze_ba_plot_cross
  ( epswr_figure_t *eps,
    int32_t x, 
    int32_t y,
    bool_t fill, 
    bool_t draw
  )
  {
    tmaze_plot_cross(eps, x, y, tmaze_ba_CROSS_RADIUS, fill, draw);
  }
