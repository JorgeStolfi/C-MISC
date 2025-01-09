/* See {tmaze_bb.h} */
/* Last edited on 2024-12-21 11:31:17 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <assert.h>

#include <bool.h>
#include <affirm.h>
#include <frgb.h>
#include <jsrandom.h>

#include <tmaze_bb.h>

/* INTERNAL PARAMETERS */

/* #define tmaze_bb_ROAD_RADIUS (0.0517767) */
#define tmaze_bb_ROAD_RADIUS ((M_SQRT2 - 1)/4)
  /* Half-width of the roads in Blip Blop mazes. */

/* #define tmaze_bb_CROSS_RADIUS (0.0517767) */
#define tmaze_bb_CROSS_RADIUS ((M_SQRT2 - 1)/4)
  /* Half-width of the arms of the cross used to indicate 
    empty tileset in Blip Blop maze patterns. */

/* IMPLEMENTATIONS */

tmaze_bb_tile_t tmaze_bb_get_tile_from_dirs(tmaze_dir_t dir1, tmaze_dir_t dir2)
  {
    /* Swap {dir1,dir2} if necessary so that {dir1} is West or East: */
    if ((dir1 == tmaze_dir_N) || (dir1 == tmaze_dir_S))
      { tmaze_dir_t dirt = dir1; dir1 = dir2; dir2 = dirt; }
    
    /* Now figure out the tile from the directions: */
    tmaze_bb_tile_t tile;
    if (dir1 == tmaze_dir_W)
      { if (dir2 == tmaze_dir_N)
          { tile = tmaze_bb_tile_BLIP; }
        else if (dir2 == tmaze_dir_S)
          { tile = tmaze_bb_tile_BLOP; }
        else 
          { assert(FALSE); }
      }
    else if (dir1 == tmaze_dir_E)
      { if (dir2 == tmaze_dir_N)
          { tile = tmaze_bb_tile_BLOP; }
        else if (dir2 == tmaze_dir_S)
          { tile = tmaze_bb_tile_BLIP; }
        else 
          { assert(FALSE); }
      }
    else
      { assert(FALSE); }
      
    return tile;
  }

tmaze_t tmaze_bb_make(int nx, int ny, bool_t torus)
  {
    tmaze_t M = tmaze_make(nx, ny, torus);
    return M;
  }

void tmaze_bb_random(tmaze_t *M, int seed)
  {
    tmaze_random(M, seed, tmaze_bb_tileset_BOTH);
  }

tmaze_t tmaze_bb_from_string(char *string, bool_t torus)
  { 
    /* Get the maze size from the string: */
    int nx, ny;
    tmaze_size_from_string(string, &nx, &ny);
          
    /* Now allocate the maze: */
    tmaze_t M = tmaze_make(nx, ny, torus);
    
    /* Fill {M.tile} from {string}, from TOP TO BOTTOM: */
    int x; /* Cell x index (column), in {0..nx-1}. */
    int y; /* Cell y index (row), in {0..ny-1}. */
    char *p = string;
    for (y = ny-1; y >= 0; y--)
      { for (x = 0; x < nx; x++)
          { int k = tmaze_cell_index(x, y, M.nx, M.ny);
            char c = (*p); p++;
            demand(c != 0, "string is too short");
            tmaze_tile_t t;
            switch (c)
              {
                case '>': t = tmaze_bb_tile_BLIP; break;
                case '<': t = tmaze_bb_tile_BLOP; break;
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

tmaze_pattern_t tmaze_bb_pattern_from_string
  ( char *string, 
    bool_t torus,
    bool_t sub,
    int val
  )
  { 
    /* Get the pattern size from the string: */
    int nx, ny;
    tmaze_size_from_string(string, &nx, &ny);
          
    /* Now allocate the pattern: */
    tmaze_pattern_t M = tmaze_make_pattern(nx, ny, torus, sub, val);
    
    /* Fill {M.tset} from {string}, from TOP TO BOTTOM: */
    int x; /* Cell x index (column), in {0..nx-1}. */
    int y; /* Cell y index (row), in {0..ny-1}. */
    char *p = string;
    for (y = ny-1; y >= 0; y--)
      { for (x = 0; x < nx; x++)
          { int k = tmaze_cell_index(x, y, M.nx, M.ny);
            char c = (*p); p++;
            demand(c != 0, "string is too short");
            tmaze_tileset_t S;
            switch (c)
              {
                case '>': S = (1 << tmaze_bb_tile_BLIP); break;
                case '<': S = (1 << tmaze_bb_tile_BLOP); break;
                case '*': S = tmaze_bb_tileset_BOTH; break;
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
  
tmaze_size_t tmaze_bb_tot_size(tmaze_t *M)
  { 
    /* The size of a component is the count of its graph vertices. So: */
    demand(M->nySN <= tmaze_MAX_SIZE/2/M->nxWE, "too many vertices in graph");
    tmaze_size_t nvWE = M->nxWE * M->ny;  /* Number of W joint vertices. */
    tmaze_size_t nvSN = M->nx * M->nySN;  /* Number of S joint vertices. */
    tmaze_size_t nv = nvWE + nvSN;        /* Number of vertices in graph. */
    return nv;
  }

tmaze_size_t tmaze_bb_max_comp_size(tmaze_t *M)
  { 
    return tmaze_bb_tot_size(M);
  }

void tmaze_bb_print(FILE *wr, tmaze_t *M)
  {
    int x; /* Cell x index (column), in {0..nx-1}. */
    int y; /* Cell y index (row), in {0..ny-1}. */
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
            fputc(tmaze_bb_dir_chars[M->tile[k]], wr);
          }
        if (M->torus) { fputc('-', wr); } 
        fputc('\n', wr);
      }
  }

void tmaze_bb_plot_tile
  ( epswr_figure_t *eps,
    tmaze_bb_tile_t tile,
    int x,
    int y,
    bool_t curved,
    frgb_t *rgb_S,  /* Fill color to use for rod that ends on South side. */
    frgb_t *rgb_N,  /* Fill color to use for rod that ends on North side. */
    bool_t open[],
    bool_t fill,
    bool_t draw
  )
  { int nc = 8; /* Number of points in the outline of one elbow road. */
    double xp[nc], yp[nc]; /* X and Y coordinates of points of elbow road. */
    
    int p; /* Which road (0 ends at North, 1 at South). */
    for (p = 0; p < 2; p++)
      { 
        if (curved)
          { tmaze_bb_path_elbow_bezier_make(xp, yp, nc, tmaze_bb_ROAD_RADIUS); }
        else
          { tmaze_bb_path_elbow_poly_make(xp, yp, nc, tmaze_bb_ROAD_RADIUS); }

        /* Rotate path from the North-West to the needed orientation: */
        int k = 2*p + (tile == tmaze_bb_tile_BLOP ? 1 : 0);
        tmaze_path_rotate(xp, yp, nc, k);

        /* Translate the path from the origin to the cell's center: */
        tmaze_path_translate(xp, yp, nc, x + 0.5, y + 0.5);

        if (fill)
          { /* Set the appropriate fill color, if so requested: */
            if ((rgb_S != NULL) || (rgb_N != NULL))
              { frgb_t *rgb = (p == 0 ? rgb_N : rgb_S);
                epswr_set_fill_color(eps, rgb->c[0], rgb->c[1], rgb->c[2]);
              }
            /* Fill the road: */
            if (curved)
              { epswr_bezier_polygon(eps, TRUE, xp, yp, nc/4, TRUE, FALSE, TRUE); }
            else
              { epswr_polygon(eps, TRUE, xp, yp, nc, TRUE, FALSE, TRUE); }
          }

        if (draw)
          { /* Draw the outlines of the road inside the cell: */
            if (curved)
              { epswr_curve(eps, xp[0], yp[0], xp[1], yp[1], xp[2], yp[2], xp[3], yp[3]); }
            else
              { epswr_segment(eps, xp[0], yp[0], xp[1], yp[1]);
                epswr_segment(eps, xp[1], yp[1], xp[2], yp[2]);
                epswr_segment(eps, xp[2], yp[2], xp[3], yp[3]);
              }
            /* Get {dir34}, direction of side that contains segm. {[3]--[4]}: */ 
            tmaze_dir_t dir34 = (tmaze_dir_E + k) % 4; 
            if (! open[dir34]) { epswr_segment(eps, xp[3], yp[3], xp[4], yp[4]); }
            if (curved)
              { epswr_curve(eps, xp[4], yp[4], xp[5], yp[5], xp[6], yp[6], xp[7], yp[7]); }
            else
              { epswr_segment(eps, xp[4], yp[4], xp[5], yp[5]);
                epswr_segment(eps, xp[5], yp[5], xp[6], yp[6]);
                epswr_segment(eps, xp[6], yp[6], xp[7], yp[7]);
              }
            /* Get {dir70}, direction of side that contains segm. {[7]--[0]}: */ 
            tmaze_dir_t dir70 = (tmaze_dir_N + k) % 4; 
            if (! open[dir70]) { epswr_segment(eps, xp[7], yp[7], xp[0], yp[0]); }
          }
      }
  }

void tmaze_bb_path_elbow_poly_make(double xp[], double yp[], int nc, double r)
  { /* This procedure generates 8 corners, so: */
    assert(nc == 8); 
    double h = 0.5;
    tmaze_path_quarter_circle_poly_make(&(xp[0]), &(yp[0]), 4);
    tmaze_path_scale(&(xp[0]), &(yp[0]), 4, -(h+r), -(h+r));
    tmaze_path_translate(&(xp[0]), &(yp[0]), 4, h, h);
    tmaze_path_quarter_circle_poly_make(&(xp[4]), &(yp[4]), 4);
    tmaze_path_scale(&(xp[4]), &(yp[4]), 4, -(h-r), -(h-r));
    tmaze_path_translate(&(xp[4]), &(yp[4]), 4, h, h);
    tmaze_path_reverse(&(xp[4]), &(yp[4]), 4);
  }

void tmaze_bb_path_elbow_bezier_make(double xp[], double yp[], int nc, double r)
  { /* This procedure generates 2 Bezier arcs, so: */
    assert(nc == 8); 
    double h = 0.5;
    tmaze_path_quarter_circle_bezier_make(&(xp[0]), &(yp[0]), 4);
    tmaze_path_scale(&(xp[0]), &(yp[0]), 4, -(h+r), -(h+r));
    tmaze_path_translate(&(xp[0]), &(yp[0]), 4, h, h);
    tmaze_path_quarter_circle_bezier_make(&(xp[4]), &(yp[4]), 4);
    tmaze_path_scale(&(xp[4]), &(yp[4]), 4, -(h-r), -(h-r));
    tmaze_path_translate(&(xp[4]), &(yp[4]), 4, h, h);
    tmaze_path_reverse(&(xp[4]), &(yp[4]), 4);
  }
  
void tmaze_bb_plot_maze
  ( char *dir,
    char *name,
    tmaze_t *M,
    bool_t cells,  /* TRUE paints the background of each tile. */
    bool_t grids,  /* TRUE paints the cell grid. */
    double tsize,  /* Cell size (mm) */
    bool_t curved  /* If TRUE uses curved roads, FALSE uses polygonal ones. */
  )
  {
    int val = 0;         /* Number of ports to show on each side. */
    double margin = 2.0; /* Figure margin in mm. */

    auto void plot_tile(epswr_figure_t *eps, tmaze_t *M, int x, int y, bool_t fill, bool_t draw);
      /* Plots the tileset in column {x} and  row {y} of maze pattern {M}, with given style. */

    void plot_tile(epswr_figure_t *eps, tmaze_t *M, int x, int y, bool_t fill, bool_t draw)
      { tmaze_bb_plot_maze_tile(eps, M, x, y, curved, NULL, NULL, fill, draw); }

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
  
void tmaze_bb_plot_maze_tile
  ( epswr_figure_t *eps,
    tmaze_t *M, 
    int x, 
    int y,
    bool_t curved,
    frgb_t *rgb_S,  /* Fill color to use for rod that ends on South side. */
    frgb_t *rgb_N,  /* Fill color to use for rod that ends on North side. */
    bool_t fill, 
    bool_t draw
  )
  {
    /* Get the cell index: */
    tmaze_cell_index_t k = tmaze_cell_index(x, y, M->nx, M->ny);
    /* Determine which sides of the tile should be left open: */
    bool_t open[4];
    tmaze_bb_get_open_ends_of_maze_tile(M, k, open);
    /* Plot the tile: */
    tmaze_bb_plot_tile(eps, M->tile[k], x, y, curved, rgb_S, rgb_N, open, fill, draw);
  }
 
void tmaze_bb_plot_pattern_tileset
  ( epswr_figure_t *eps,
    tmaze_pattern_t *M, 
    int x, 
    int y,
    bool_t curved,
    bool_t fill, 
    bool_t draw
  )
  {
    /* Get the cell's index: */
    tmaze_cell_index_t k = tmaze_cell_index(x, y, M->nx, M->ny);
    /* Get the set of allowed tiles {S} for this cell: */
    tmaze_tileset_t S = M->tset[k];

    if (tmaze_tileset_contains(S, tmaze_bb_tileset_BOTH))
      { /* Free cell -- leave it blank: */ }
    else if (S == tmaze_tileset_NONE)
      { /* Invalid cell -- paint a square dot: */
        tmaze_bb_plot_cross(eps, x, y, fill, draw);
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
        tmaze_bb_get_open_ends_of_pattern_tileset(M, k, open);
        /* Plot the tile: */
        tmaze_bb_plot_tile(eps, tile, x, y, curved, NULL, NULL, open, fill, draw);
      }
  }
  
void tmaze_bb_get_open_ends_of_maze_tile
  ( tmaze_t *M, 
    tmaze_cell_index_t k, 
    bool_t open[]
  )
  { /* Enumerate the four neighbors of cell {k}: */
    int id;
    for (id = 0; id < 4; id++)
      { /* Get the direction {dir} towards neighbor: */
        tmaze_dir_t dir = (tmaze_dir_t)id; 
        /* Get the index {kN} of the {dir}-neighbor of cell {k}: */
        tmaze_cell_index_t kN = 
          tmaze_cell_neighbor_index(k, dir, M->nx, M->ny, M->torus);
        /* Decide whether any end of tile {k} on side {dir} should be left open: */
        open[id] = (kN != tmaze_cell_index_NONE);
      }
  }
    
void tmaze_bb_get_open_ends_of_pattern_tileset
  ( tmaze_pattern_t *M, 
    tmaze_cell_index_t k, 
    bool_t open[]
  )
  {
    /* Enumerate the four neighbors of cell {k}: */
    int id;
    for (id = 0; id < 4; id++)
      { /* Get the direction {dir} towards the neighbor: */
        tmaze_dir_t dir = (tmaze_dir_t)id;
        /* Get the index {kN} of the {dir}-neighbor of cell {k}: */
        tmaze_cell_index_t kN = 
          tmaze_cell_neighbor_index(k, dir, M->nx, M->ny, M->torus);
        /* Decide whether any end of tile {k} on side {dir} should be left open: */
        open[id] = (kN != tmaze_cell_index_NONE) || M->sub;
      }
  }   

void tmaze_bb_plot_cross
  ( epswr_figure_t *eps,
    int x, 
    int y,
    bool_t fill, 
    bool_t draw
  )
  {
    tmaze_plot_cross(eps, x, y, tmaze_bb_CROSS_RADIUS, fill, draw);
  }
