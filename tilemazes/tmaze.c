/* See {tmaze.h} */
/* Last edited on 2023-10-01 19:37:31 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <sys/stat.h>
#include <sys/types.h>

#include <bool.h>
#include <epswr.h>
#include <affirm.h>
#include <jsrandom.h>

#include <tmaze.h>
#include <tmaze_gen.h>

bool_t tmaze_tile_is_in_set(tmaze_tile_t t, tmaze_tileset_t S)
  { return ((S & (1 << t)) != 0); }
  
bool_t tmaze_tileset_contains(tmaze_tileset_t S, tmaze_tileset_t P)
  { return (P & (~ S)) == 0; }

tmaze_tile_t tmaze_min_tile(tmaze_tileset_t S)
  { 
    S = (S & tmaze_tileset_ALL);
    demand(S != 0, "the empty set has no minimum"); 
    int iel = 0;
    while ((S & 1) == 0) { S = (S >> 1); iel++; }
    assert(iel <= tmaze_MAX_TILE_TYPES);
    return (tmaze_tile_t)iel;
  }

tmaze_tile_t tmaze_pick_tile(tmaze_tileset_t S)
  { /* Collect the elements of {S} in {el[0..nel-1]}: */
    int nel = 0;
    tmaze_tile_t el[tmaze_MAX_TILE_TYPES];
    int t;
    for (t = 0; t < tmaze_MAX_TILE_TYPES; t++)
      { if (tmaze_tile_is_in_set((tmaze_tile_t)t, S))
          { el[nel] = (tmaze_tile_t)t; nel++; }
      }
    /* Check for empty input: */
    demand(nel > 0, "cannot pick from an empty set");
    /* Now pick one element at random: */
    int k = int32_abrandom(0, nel-1);
    return el[k];
  }

tmaze_cell_index_t tmaze_cell_index(int x, int y, int nx, int ny)
  {
    return x + y * nx;
  }

void tmaze_cell_position(tmaze_cell_index_t k, int nx, int ny, int *xP, int *yP)
  {
    (*xP) = k % nx;
    (*yP) = k / nx;
  }

int tmaze_cell_row_col_inc(int z, int d, int n, bool_t torus)
  {
    demand((z >= 0) && (z < n), "invalid row/col index");
    z += d;
    if (torus)
      { z = z%n; if (z < 0) { z += n; } }
    else
      { if ((z < 0) || (z >= n)) { z = -1; } }
    return z;
  }

tmaze_cell_index_t tmaze_cell_neighbor_index
  ( tmaze_cell_index_t k, 
    tmaze_dir_t dir,
    int nx,
    int ny,
    bool_t torus
  )
  {
    int x, y;
    tmaze_cell_position(k, nx, ny, &x, &y);
    switch (dir)
      { 
      case tmaze_dir_N:
        { y = tmaze_cell_row_col_inc(y, +1, ny, torus); break; }
      case tmaze_dir_W:
        { x = tmaze_cell_row_col_inc(x, -1, nx, torus); break; }
      case tmaze_dir_S:
        { y = tmaze_cell_row_col_inc(y, -1, ny, torus); break; }
      case tmaze_dir_E:
        { x = tmaze_cell_row_col_inc(x, +1, nx, torus); break; }
      }
    if ((x < 0) || (x >= nx) || (y < 0) | (y >= ny))
      { return tmaze_cell_index_NONE; }
    else
      { return tmaze_cell_index(x, y, nx, ny); }
  }

tmaze_t tmaze_make(int nx, int ny, bool_t torus)
  { 
    demand(nx > 0, "maze width must be positive");
    demand(ny > 0, "maze height must be positive");

    /* Alocate and fill the header record: */
    tmaze_t M;
    M.nx = nx;
    M.ny = ny;
    M.torus = torus;
    
    
    /* Compute number of cells and of vertical and horizontal cell joints: */
    demand(ny <= tmaze_MAX_CELL_COUNT/nx, "too many cells");
    M.nt = nx*ny;                    /* Number of cells in maze. */
    M.nxWE = (torus ? nx : nx + 1);  /* Number of vert cell joints per row. */
    M.nySN = (torus ? ny : ny + 1);  /* Number of horiz cell joints per col. */
    
    /* Allocate the cell array and fill it with tile 0: */
    M.tile = notnull(malloc(M.nt*sizeof(tmaze_tile_t)), "no mem");
    tmaze_fill(&M, tmaze_tile_0);
    return M;
  }

void tmaze_fill(tmaze_t *M, tmaze_tile_t t)
  {
    int x; /* Cell x index (column), in {0..nx-1}. */
    int y; /* Cell y index (row), in {0..ny-1}. */
    for (y = 0; y < M->ny; y++)
      { for (x = 0; x < M->nx; x++)
          { int k = tmaze_cell_index(x, y, M->nx, M->ny);
            M->tile[k] = t;
          }
      }
  }

void tmaze_random(tmaze_t *M, int seed, tmaze_tileset_t S)
  { srandom(seed);
    int x; /* Cell x index (column), in {0..nx-1}. */
    int y; /* Cell y index (row), in {0..ny-1}. */
    for (y = 0; y < M->ny; y++)
      { for (x = 0; x < M->nx; x++)
          { int k = tmaze_cell_index(x, y, M->nx, M->ny);
            M->tile[k] = tmaze_pick_tile(S);
          }
      }
  }

void tmaze_size_from_string(char *string, int *nxP, int *nyP)
  {
    int nx = 0;
    int ny = 0;
    int k = 0; /* Number of elements in the last row. */
    char *p = string;
    while ((*p) != 0)
      { if ((*p) == ';')
          { /* End of a row. */
            if (ny == 0) 
              { nx = k; } 
            else 
              { demand(nx == k, "inconsistent row length"); }
            ny++; k = 0;
          }
        else
          { k++; }
        p++;
      }
    demand(k == 0, "missing line terminator ';'");
    (*nxP) = nx;
    (*nyP) = ny;
  }

tmaze_pattern_t tmaze_make_pattern
  ( int nx, 
    int ny, 
    bool_t torus,
    bool_t sub,
    int val
  )
  { 
    demand(nx > 0, "maze width must be positive");
    demand(ny > 0, "maze height must be positive");

    /* Alocate and fill the header record: */
    tmaze_pattern_t M;
    M.nx = nx;
    M.ny = ny;
    M.torus = torus;
    M.sub = sub;
    M.val = val;
    
    /* Compute number of cells {nt}: */
    int nt = nx*ny;  /* Number of cells in maze. */
    
    /* Allocate the cell array: */
    M.tset = notnull(malloc(nt*sizeof(tmaze_tileset_t)), "no mem");
    tmaze_fill_pattern(&M, tmaze_tileset_NONE);
    return M;
  }

void tmaze_fill_pattern(tmaze_pattern_t *M, tmaze_tileset_t S)
  {
    int x; /* Cell x index (column), in {0..nx-1}. */
    int y; /* Cell y index (row), in {0..ny-1}. */
    for (y = 0; y < M->ny; y++)
      { for (x = 0; x < M->nx; x++)
          { int k = x + M->nx * y;
            M->tset[k] = S;
          }
      }
  }

void tmaze_random_pattern(tmaze_pattern_t *M, int seed, tmaze_tileset_t S)
  { srandom(seed);
    int x; /* Cell x index (column), in {0..nx-1}. */
    int y; /* Cell y index (row), in {0..ny-1}. */
    for (y = 0; y < M->ny; y++)
      { for (x = 0; x < M->nx; x++)
          { int k = x + M->nx * y;
            M->tset[k] = (1 << tmaze_pick_tile(S));
          }
      }
  }

void tmaze_make_dir(char *dir, char *sub)
  {
    mkdir(dir, 0777);
    if ((sub != NULL) && (strlen(sub) != 0) && (strcmp(sub, ".") != 0))
      { char *dirsub = NULL;
        asprintf(&dirsub, "%s/%s", dir, sub);
        mkdir(dirsub, 0777);
        free(dirsub);
      }
  }
   
void tmaze_plot_maze
  ( char *dir,
    char *name,
    tmaze_t *M,
    bool_t cells,  /* TRUE paints the background of each tile. */
    bool_t grids,  /* TRUE paints the cell grid. */
    tmaze_plot_tile_proc_t plot_tile, /* Tile plotter, or NULL, */
    int val,       /* Valency of tiles, for ports plotting (0 to omit ports). */
    double tsize,  /* Cell size (mm) */
    double margin  /* Margin width (mm). */
  )
  {
    /* Assemble the directory and name into the figure name and prefix: */
    
    /* Create the EPS stream {eps}: */
    double mm = epswr_pt_per_mm;  /* One millimeter in Postscript points. */
    double hsize = tsize * mm * M->nx; /* H size of plotting area (pt). */
    double vsize = tsize * mm * M->ny; /* V size of plotting area (pt). */
    double mrg = margin*mm; /* Margin around the plotting area (pt) */
    bool_t verbose = TRUE;
    epswr_figure_t *eps = epswr_new_named_figure
      ( dir, NULL, name, -1, NULL,
        hsize, vsize, mrg, mrg, mrg, mrg,
        verbose
      );
    
    /* Set the window mapping so that cells are 1 unit wide: */
    epswr_set_window
      ( eps,
        mrg, mrg+hsize,      mrg, mrg+vsize, 
        FALSE,
        0.0, (double)M->nx,  0.0, (double)M->ny
      );
      
    /* Determine whether tiles should be printed: */
    bool_t tiles = (plot_tile != NULL);
    
    /* Check whether there are any ports to be printed: */
    bool_t ports = (val != 0);
    
    /* Plot the tiles: */
    int layer; /* Layer of plot. */
    for (layer = tmaze_layer_MIN; layer <= tmaze_layer_MAX; layer++) 
      { /* Set the default pen and fill color for each layer: */
        switch(layer)
          {
          case tmaze_layer_HOLE:
            if (cells) 
              { epswr_set_fill_color(eps, 1.0000,1.0000,1.0000); }
            break;

          case tmaze_layer_GRID: 
            if (grids) 
              { epswr_set_pen(eps, 0.8000,0.9000,1.000, 0.10, 0,0); }
            break;

          case tmaze_layer_ROAD: 
          case tmaze_layer_CURB: 
            if (tiles) 
              { epswr_set_fill_color(eps, 1.0000,0.9500,0.6000);
                epswr_set_pen(eps, 0.0000,0.0000,0.0000, 0.15, 0,0);
              } 
            break;

          case tmaze_layer_DOTS: 
            if (ports) 
               { epswr_set_fill_color(eps, 1.0000,0.0000,0.0000);
                 epswr_set_pen(eps, 0.0000,0.0000,0.0000, 0.10, 0,0);
               }
            break;
          }
        
        /* Paint all cells, all layers: */
        int x; /* Cell x index (column), in {0..nx-1}. */
        int y; /* Cell y index (row), in {0..ny-1}. */
        for (y = 0; y < M->ny; y++)
          { for (x = 0; x < M->nx; x++)
              { 
                switch(layer)
                  {
                    case tmaze_layer_HOLE:
                      if (cells) { tmaze_plot_cell(eps, x, y, TRUE, FALSE); }
                      break;

                    case tmaze_layer_GRID: 
                      if (grids) { tmaze_plot_cell(eps, x, y, FALSE, TRUE); } 
                      break;

                    case tmaze_layer_ROAD: 
                      if (tiles) { plot_tile(eps, M, x, y, TRUE, FALSE); }
                      break;

                    case tmaze_layer_CURB: 
                      if (tiles) { plot_tile(eps, M, x, y, FALSE, TRUE); } 
                      break;

                    case tmaze_layer_DOTS: 
                      if (ports) 
                        { tmaze_plot_cell_ports(eps, x, y, TRUE, val, 0.25, TRUE, TRUE); } 
                      break;
                  }
              }
          }
      }
    epswr_end_figure(eps);
  }

void tmaze_plot_pattern
  ( char *dir,
    char *name,
    tmaze_pattern_t *M, 
    bool_t cells,  /* TRUE paints the background of each tile. */
    bool_t grids,  /* TRUE paints the cell grid. */
    tmaze_plot_tileset_proc_t plot_tileset, /* Tileset plotter, or NULL, */
    bool_t ports,  /* TRUE paints the ports as dots. */
    double tsize,  /* Cell size (mm) */
    double margin  /* Margin width (mm). */
  )
  {
    /* Create the EPS stream {eps}: */
    double mm = epswr_pt_per_mm;  /* One millimeter in Postscript points. */
    double hsize = tsize * mm * M->nx; /* H size of plotting area (pt). */
    double vsize = tsize * mm * M->ny; /* V size of plotting area (pt). */
    double mrg = margin*mm; /* Margin around the plotting area (pt) */
    bool_t verbose = TRUE;
    epswr_figure_t *eps = epswr_new_named_figure
      ( dir, NULL, name, -1, NULL,
        hsize, vsize, mrg, mrg, mrg, mrg,
        verbose
      );
    
    /* Set the window mapping so that cells are 1 unit wide: */
    epswr_set_window
      ( eps,
        mrg, mrg+hsize,      mrg, mrg+vsize, 
        FALSE,
        0.0, (double)M->nx,  0.0, (double)M->ny
      );
      
    /* Determine whether tiles should be printed: */
    bool_t tiles = (plot_tileset != NULL);
    
    /* Check whether there are any ports to be printed: */
    if (M->val == 0) { ports = FALSE; }
    
    /* Paint the cells: */
    int layer;
    for (layer = tmaze_layer_MIN; layer <= tmaze_layer_MAX; layer++) 
      { 
        /* Set the default pen and fill color for each layer: */
        switch(layer)
          {
          case tmaze_layer_HOLE:
            if (cells) 
              { epswr_set_fill_color(eps, 1.0000,1.0000,1.0000); }
            break;

          case tmaze_layer_GRID: 
            if (grids) 
              { epswr_set_pen(eps, 0.8000,0.9000,1.000, 0.10, 0,0); }
            break;

          case tmaze_layer_ROAD: 
          case tmaze_layer_CURB: 
            if (tiles) 
              { epswr_set_fill_color(eps, 1.0000,0.9500,0.6000);
                epswr_set_pen(eps, 0.0000,0.0000,0.0000, 0.15, 0,0);
              } 
            break;

          case tmaze_layer_DOTS: 
            if (ports) 
               { epswr_set_fill_color(eps, 1.0000,0.0000,0.0000);
                 epswr_set_pen(eps, 0.0000,0.0000,0.0000, 0.10, 0,0);
               }
            break;
          }
      
        /* Paint all cells, all layers: */
        int x; /* Cell x index (column), in {0..nx-1}. */
        int y; /* Cell y index (row), in {0..ny-1}. */
        for (y = 0; y < M->ny; y++)
          { for (x = 0; x < M->nx; x++)
              { switch(layer)
                  {
                    case tmaze_layer_HOLE:
                      if (cells) { tmaze_plot_cell(eps, x, y, TRUE, FALSE); }
                      break;

                    case tmaze_layer_GRID: 
                      if (grids) { tmaze_plot_cell(eps, x, y, FALSE, TRUE); } 
                      break;

                    case tmaze_layer_ROAD: 
                      if (tiles) { plot_tileset(eps, M, x, y, TRUE, FALSE); }
                      break;

                    case tmaze_layer_CURB: 
                      if (tiles) { plot_tileset(eps, M, x, y, FALSE, TRUE); } 
                      break;

                    case tmaze_layer_DOTS: 
                      if (ports) 
                        { tmaze_plot_cell_ports(eps, x, y, TRUE, M->val, 0.25, TRUE, TRUE); } 
                      break;
                  }
              }
          }
      }
    epswr_end_figure(eps);
  }

void tmaze_plot_cell
  ( epswr_figure_t *eps,
    int x, 
    int y,
    bool_t fill, 
    bool_t draw
  )
  {
    int nc = 4; /* Number of corners in the cell's outline. */
    double xp[nc], yp[nc]; /* X and Y coordinates of corners of cell. */
    tmaze_path_cell_make(xp, yp, nc);
    tmaze_path_translate(xp, yp, nc, x + 0.5, y+ 0.5);
    epswr_polygon(eps, TRUE, xp, yp, nc, fill, draw, TRUE);

  }
  
void tmaze_plot_cross
  ( epswr_figure_t *eps,
    int x, 
    int y,
    double r,
    bool_t fill, 
    bool_t draw
  )
  {
    int nc = 12; /* Number of corners in the cross's outline. */
    double xp[nc], yp[nc]; /* X and Y coordinates of corners of cell. */
    tmaze_path_cross_make(xp, yp, nc, r);
    tmaze_path_translate(xp, yp, nc, x + 0.5, y+ 0.5);
    epswr_polygon(eps, TRUE, xp, yp, nc, fill, draw, TRUE);
  }
  
void tmaze_plot_cell_ports
  ( epswr_figure_t *eps,
    int x, 
    int y,
    bool_t once,  /* If TRUE avoids paint the same port twice. */
    int val,      /* Valency (number of ports per side). */
    double r,     /* Dot radius in millimeters. */
    bool_t fill,  /* If TRUE fills the dots with the current fill color. */
    bool_t draw   /* If TRUE draws the dot outlines with the current pen. */
  )
  {
    int id;
    for (id = 0; id < 4; id++)
      { tmaze_dir_t dir = (tmaze_dir_t)id;
        /* Avoid painting twice the same ports: */
        bool_t ok_X = ((x == 0) || (dir != tmaze_dir_W));
        bool_t ok_Y = ((y == 0) || (dir != tmaze_dir_S));
        if (ok_X && ok_Y)
          { tmaze_plot_cell_side_ports(eps, x, y, dir, val, r, TRUE, TRUE); }
      }
  }

void tmaze_plot_cell_side_ports
  ( epswr_figure_t *eps,
    int x, 
    int y,
    tmaze_dir_t dir,
    int val, 
    double r,
    bool_t fill, 
    bool_t draw
  )
  {
    double xc, yc;
    int k;
    for (k = 0; k < val; k++)
      { /* Position of port {k} along the side {dir} of the cell: */
        double s = (k + 0.5) / val;
        /* Compute coordinates {xc,yc} of port: */
        switch(dir)
          {
          case tmaze_dir_S: xc = x + s; yc = y;     break;
          case tmaze_dir_N: xc = x + s; yc = y + 1; break;
          case tmaze_dir_W: xc = x;     yc = y + s; break;
          case tmaze_dir_E: xc = x + 1; yc = y + s; break;
          default: assert(FALSE);
          }
        /* Paint the dot: */
        epswr_dot(eps, xc, yc, r, fill, draw);
      }
  }

void tmaze_path_cell_make(double xp[], double yp[], int nc)
  { /* This procedure generates 8 corners, so: */
    assert(nc == 4); 
    int kc = 0; /* Counts corners appended to the T-road's outline. */
    xp[kc] = -0.5; yp[kc] = -0.5; kc++;
    xp[kc] = +0.5; yp[kc] = -0.5; kc++;
    xp[kc] = +0.5; yp[kc] = +0.5; kc++;
    xp[kc] = -0.5; yp[kc] = +0.5; kc++;
    assert(kc == nc);
  }

void tmaze_path_cross_make(double xp[], double yp[], int nc, double r)
  { /* This procedure generates 12 corners, so: */
    assert(nc == 12); 
    double a = M_SQRT2 * r;
    double b = 2*a;
    demand(b <= 0.40, "cross is too big");
    int kc = 0; /* Counts corners appended to the T-road's outline. */
    xp[kc] = 00; yp[kc] = -a; kc++;
    xp[kc] = +a; yp[kc] = -b; kc++;
    xp[kc] = +b; yp[kc] = -a; kc++;
    xp[kc] = +a; yp[kc] = 00; kc++;
    xp[kc] = +b; yp[kc] = +a; kc++;
    xp[kc] = +a; yp[kc] = +b; kc++;
    xp[kc] = 00; yp[kc] = +a; kc++;
    xp[kc] = -a; yp[kc] = +b; kc++;
    xp[kc] = -b; yp[kc] = +a; kc++;
    xp[kc] = -a; yp[kc] = 00; kc++;
    xp[kc] = -b; yp[kc] = -a; kc++;
    xp[kc] = -a; yp[kc] = -b; kc++;
    assert(kc == nc);
  }

void tmaze_path_rotate(double xp[], double yp[], int nc, int k)
  { 
    /* Reduce {k} mod 4 to {0..3}: */
    k = ((k % 4) + 4) % 4;
  
    /* Now rotate each vertex {k*90} degrees: */
    int i;
    for (i = 0; i < nc; i++) 
      { switch (k)
          {
          case 0:
            { /* Normal orientation: */ break; }
          case 1:
            { /* Rot {+90}: */ double t = xp[i]; xp[i] = -yp[i]; yp[i] = t; break; }
          case 2:
            { /* Rot {180}: */ xp[i] = -xp[i]; yp[i] = -yp[i]; break; }
          case 3:
            { /* Rot {-90}: */ double t = xp[i]; xp[i] = yp[i]; yp[i] = -t; break; }
          default:
            { /* Should never happen: */ assert(FALSE); }
          }
      }
  }
  
void tmaze_path_scale(double xp[], double yp[], int nc, double sx, double sy)
  { int i;
    for (i = 0; i < nc; i++) { xp[i] *= sx; yp[i] *= sy; }
  }

void tmaze_path_translate(double xp[], double yp[], int nc, double dx, double dy)
  { int i;
    for (i = 0; i < nc; i++) { xp[i] += dx; yp[i] += dy; }
  }
       
void tmaze_path_reverse(double xp[], double yp[], int nc)
  { int i, j;
    for (i = 0, j = nc-1; i < j; i++, j--) 
      { double t;
        t = xp[i]; xp[i] = xp[j]; xp[j] = t;
        t = yp[i]; yp[i] = yp[j]; yp[j] = t;
      }
  }

#define SIN_PI_OVER_EIGHT (0.38268343236508977172)
#define COS_PI_OVER_EIGHT (0.92387953251128675613)
#define TAN_PI_OVER_EIGHT (0.41421356237309504879)

void tmaze_path_quarter_circle_poly_make(double xp[], double yp[], int nc)
  { /* This procedure generates 4 corners, so: */
    assert(nc == 4); 
    double t = TAN_PI_OVER_EIGHT;
    int kc = 0; /* Counts corners appended to the road's outline. */
    xp[kc] = +1; yp[kc] = 00; kc++;
    xp[kc] = +1; yp[kc] = +t; kc++;
    xp[kc] = +t; yp[kc] = +1; kc++;
    xp[kc] = 00; yp[kc] = +1; kc++;
    assert(kc == nc);
  }

void tmaze_path_quarter_circle_bezier_make(double xp[], double yp[], int nc)
  { /* This procedure generates a single bezier arc, so: */
    assert(nc == 4); 
    double q = (4*M_SQRT2 - 4)/3;
    int kc = 0; /* Counts corners appended to the road's outline. */
    xp[kc] = +1; yp[kc] = 00; kc++;
    xp[kc] = +1; yp[kc] = +q; kc++;
    xp[kc] = +q; yp[kc] = +1; kc++;
    xp[kc] = 00; yp[kc] = +1; kc++;
    assert(kc == nc);
  }
