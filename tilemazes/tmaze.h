#ifndef tmaze_H
#define tmaze_H

/* Basic tools for tiled mazes and maze patterns. */
/* Last edited on 2023-02-03 23:31:59 by stolfi */

/* !!! Rename {...cell_and...} to {...cell_position_and...} to distinguish from cell_index !!! */
/* !!! Unify tmaze_t and tmaze_pattern_t !!! */

#include <stdint.h>
#include <limits.h>
#include <stdio.h>

#include <bool.h>
#include <epswr.h>

/* MAZE TILES */

#define tmaze_MAX_TILE_TYPES 8
  /* Maximum number of distinct tile types. */
  
typedef enum
  { tmaze_tile_0 = 0,
    tmaze_tile_1 = 1,
    tmaze_tile_2 = 2,
    tmaze_tile_3 = 3,
    tmaze_tile_4 = 4,
    tmaze_tile_5 = 5,
    tmaze_tile_6 = 6,
    tmaze_tile_7 = 7
  } tmaze_tile_t;
  /* A {tmaze_tile_t} is a code that identifies a tile (specifically, a
    tile type) in a tile set. Each maze family should give better
    names for its tiles. */

/* TILE SETS */

typedef uint8_t tmaze_tileset_t;
  /* A set of tiles ({tmaze_tile_t} values). */

#define tmaze_tileset_NONE (0)
  /* The empty set of tiles. */

#define tmaze_tileset_ALL (255)
  /* The set of all possible tile types. */

bool_t tmaze_tile_is_in_set(tmaze_tile_t t, tmaze_tileset_t S);
  /* TRUE iff the tile {t} is in the tileset {S}. */

bool_t tmaze_tileset_contains(tmaze_tileset_t S, tmaze_tileset_t P);
  /* TRUE iff the tileset {S} contains (or is equal to) the tileset {P}. */

tmaze_tile_t tmaze_min_tile(tmaze_tileset_t S);
  /* Returns the tile of set {S} (which must not be empty) with
    minimum numeric value. */

tmaze_tile_t tmaze_pick_tile(tmaze_tileset_t S);
  /* Returns a random tile picked from the set {S} (which must not be empty). */

/* CELL INDEXING */

typedef uint32_t tmaze_cell_count_t;
  /* A count of cells in a maze. */

typedef uint32_t tmaze_cell_index_t;
  /* The tiles in a maze with {nx} columns and {ny} rows
    of tiles are indexed with the integers {0..nt-1}
    where {nt = ntx*nty}. */

typedef struct tmaze_t
  { int nx;                /* Number of cell columns. */
    int ny;                /* Number of cell rows. */
    bool_t torus;          /* TRUE = torus, FALSE = open rectangle. */
    tmaze_tile_t *tile;    /* Tile type in each cell. */
    /* Derived from the above: */
    tmaze_cell_count_t nt; /* Number of cells in maze. */
    int nxWE;   /* Number of vertical (East:West) tile sides per row. */
    int nySN;   /* Number of horizontal (North:South) tile sides per column. */
  } tmaze_t;
  /* A tile maze with {nx} columns and {ny} rows of cells.
    
    The rows are indexed from bottom (South) to top (North), the
    columns from left (West) to right (East), both starting at 0.
    
    If {torus} is true, the opposite sides of the board are implictly
    identified, so that the maze has the topology of a torus. In that
    case every cell has four neighbors, and there are exactly {nx}
    vertical sides per row, {ny} horizontal sides per column.
    
    If {torus} is false, the opposite edges or the board are not
    identified, and certain cells have no neighbors in certain
    directions. In that case there are {nx+1} vertical sides per row,
    and {ny+1} horizontal sides per column.
    
    The {tile} fields points to an array of {nt = nx*ny} elements
    such that {tile[k]}, where {k = tmaze_cell_index(x,y,nx,ny)},
    is the tile that occupies cell {k} in column {x} and row {y}.  */

tmaze_t tmaze_make(int nx, int ny, bool_t torus);
  /* Creates a maze {M} with the specified size and
    {M.torus} attribute, and fills all its cells with {tmaze_tile_0}. */
   
void tmaze_fill(tmaze_t *M, tmaze_tile_t t);
  /* Fills all cells of {M} with the tile {t}. */

void tmaze_random(tmaze_t *M, int seed, tmaze_tileset_t S);    
  /* Intializes the random generator with {srandom(seed)} and then
    fills {M} with random tiles chosen from {S} (which must not
    be empty). */

tmaze_cell_index_t tmaze_cell_index(int x, int y, int nx, int ny);
  /* The index of the tile that is in column {x} and row {y} of 
    an array with {nx} columns and {ny} rows. */
    
void tmaze_cell_position(tmaze_cell_index_t k, int nx, int ny, int *xP, int *yP);
  /* The inverse of {tmaze_cell_index}: given the index {k} of a tile
    in an array with {nx} columns and {ny} rows, returns the column
    index in {*xP} and the row index in {*yP}. */
    
typedef enum { 
    tmaze_dir_E = 0,
    tmaze_dir_N = 1, 
    tmaze_dir_W = 2, 
    tmaze_dir_S = 3 
  } tmaze_dir_t;
  /* A {tmaze_dir_t} value specifies one of the four cardinal
    directions. They must have numeric values {0..3} in
    counterclockwise order so that rotating a direction {dir} by {+90}
    degrees produces the direction {(dir+1) % 4}. */

#define tmaze_row_col_index_NONE (-1)

int tmaze_cell_row_col_inc(int z, int d, int n, bool_t torus);
  /* Index of the tile row (or column) {z+d} in a maze with {n} rows
    (resp. columns). Requires {z} to be in {0..n-1}; the increment {d}
    may be negative. If {torus} is true the result is {(z+d)%n}. If
    {torus} is false the result is {z+d} if that integer is in the
    range {0..n-1}, otherwise {tmaze_row_col_index_NONE} when the
    requested row (resp. column) does not exist. */

tmaze_cell_index_t tmaze_cell_neighbor_index
  ( tmaze_cell_index_t k, 
    tmaze_dir_t dir,
    int nx,
    int ny,
    bool_t torus
  );
  /* Index of the neighbor in direction {dir} of the tile with index
    {k} of a maze with {nx} columns and {ny} rows. If {torus} is true
    the neighbor always exists. Otherwise the result is
    {tmaze_cell_index_NONE} when the requested neighbor does not
    exist. */

void tmaze_size_from_string(char *string, int *nxP, int *nyP);
  /* Assumes that the {string} is a representation of the maze
    where each tile is a character, row-by-row, with each row
    terminated by a ';' character.  Returns in {*nxP} and {*nyP}
    the number of columns and rows, respectively. */

/* MAZE PATTERNS */

typedef struct tmaze_pattern_t
  { int nx;        /* Number of cell columns. */
    int ny;        /* Number of cell rows. */
    bool_t torus;  /* Topology of board */
    bool_t sub;    /* Assume it is a sub-board of a larger board. */
    int val;       /* Number of ports to draw on each side. */
    tmaze_tileset_t *tset;    /* Tiles allowed in each cell. */
  } tmaze_pattern_t;
  /* A maze pattern with {nx} columns and {ny} rows of cells, with the
     specified topology and valency. If {sub} is true, then {torus}
     must be false. The cell in column {x} and row {y} can be filled
     with any of the tiles in the set {.tset[k]} where 
     {k = tmaze_cell_index(x,y,.nx,.ny)}. */

tmaze_pattern_t tmaze_make_pattern
  ( int nx, 
    int ny, 
    bool_t torus,
    bool_t sub,
    int val
  );
  /* Creates a pattern {M} with the specified size and attributes
    {M.torus,M.sub,M.val}, and fills all its cells with the empty
    tileset. */
   
void tmaze_fill_pattern(tmaze_pattern_t *M, tmaze_tileset_t S);    
  /* Fills {M} with the tileset {S}. */
   
void tmaze_random_pattern(tmaze_pattern_t *M, int seed, tmaze_tileset_t S);    
  /* Intializes the random generator with {srandom(seed)} and then
    fills {M} with random singleton susbsets of {S} (which must not
    be empty). */

/* COMPONENT STATISTICS
  
   A /component/ of a maze {M} is a connected component of its
   underlying graph {G}.
   
   The `size' of a maze and of its components is family-dependent
   (count of tiles, vertices, edges, etc.). However, the size must be
   the same for any assignment of tiles to {M}, and the sum of the
   sizes of all components must be equal to the size of the maze. */

typedef int tmaze_size_t;
  /* Size of a maze or maze component. */

typedef int tmaze_comp_count_t;
  /* Number of components in a maze. */
    
/* FILE TOOLS */

void tmaze_make_dir(char *dir, char *sub);
  /* Creates a directory "{dir}/{sub}" if it does not exist. No-op if it exists.
    Bombs out otherwise. */

/* PLOTTING MAZES */

typedef void tmaze_plot_tile_proc_t
  ( epswr_figure_t *eps,
    tmaze_t *M, 
    int x, 
    int y,
    bool_t fill, 
    bool_t draw
  ); 
  /* Type of a procedure that plots into {p} the tile in column {x} and 
    row {y} of a maze {M}.  The procedure should assume that
    the cell is a square of unit size whose lower left corner has coordinates {(x,y)}.
    
    If {fill} is true, the procedure should first fill the tile's figures
    (roads, symbols, etc.), either with specific colors or with the 
    current fill color of {eps}.
    
    Then, if {draw} is true, the procedure should draw the outlines of 
    those figures that need outlines, either with specific pens or with the 
    current fill color of {eps}.  
    
    The {fill} and {draw} parameters must be honored only for graphic
    items that touch the boundary of the cell. The reason for this
    requirement is that lines that are suposed to connect with lines
    of adjacent cells must extend beyond the boundary, by the pen's
    radius, to join seamlessly. Ditto for outlines of colored regions
    that touch the cell's boundary. Therefore, to paint a multi-cell
    maze, the colored regions of all cells must be filled first, then
    all the outlines mut be drawn; otherwise the filling of one cell
    may erase parts of the strokes of neighboring cells.
    
    However, for parts of the plot that do not get too close to the
    tile's boundary, the procedure may fill and/or stroke them in any
    order it needs to --- including stroke some of them when {fill} is
    true, or fill some of them when {draw} is true. `Not too close'
    means that the boundary of any filled area and the midpath of any
    stroked line stay at least one pen-radius away from the cell's
    boundary. */

void tmaze_plot_maze
  ( char *dir,
    char *name,
    tmaze_t *M,
    bool_t cells,  /* TRUE paints the background of each tile. */
    bool_t grids,  /* TRUE paints the cell grid. */
    tmaze_plot_tile_proc_t *plot_tile, /* Tile plotter, or NULL, */
    int val,       /* Valency of tiles, for ports plotting (0 to omit ports). */
    double tsize,  /* Cell size (mm) */
    double margin  /* Margin width (mm). */
  );
  /* Generates an EPS figure of the maze {M}. The file will be
    called "{dir}/{name}.eps". The directory "{dir}" must exist.
    
    If {grids} is true, the procedure paints all cells and draws the
    grid lines with appropriate colors.
    
    Then, if the {plot_tile} function is not NULL, the procedure calls
    {plot_tile(eps, M, x, y, fill, draw)} where {x} and {y} range
    over all columns and rows of the board. It does one pass over all
    cells with {fill=TRUE} and {draw=FALSE}, then a second pass with
    the other way around.
    
    Then, if {ports} is true and {M->val > 0}, the procedure draws
    {M->val} dots along each side of each cell of {M}. */
  
/* PLOTTING MAZE PATTERNS */

typedef void tmaze_plot_tileset_proc_t
  ( epswr_figure_t *eps,
    tmaze_pattern_t *M, 
    int x, 
    int y,
    bool_t fill, 
    bool_t draw
  ); 
  /* A procedure that plots a tileset of a maze pattern.
    The parameters are the same as those of a {tmaze_plot_tile_proc_t},
    except that {M} is a maze pattern instead of a maze instance.
    !!! Unify with {tmaze_plot_tile_proc_t} !!! */

void tmaze_plot_pattern
  ( char *dir, 
    char *name, 
    tmaze_pattern_t *M, 
    bool_t cells,  /* TRUE paints the background of each tile. */
    bool_t grids,  /* TRUE paints the cell grid. */
    tmaze_plot_tileset_proc_t *plot_tileset, /* Tleset plotter, or NULL, */
    bool_t ports,  /* TRUE paints the ports as dots. */
    double tsize,  /* Cell size (mm) */
    double margin  /* Margin width (mm). */
  );
  /* Generates an EPS figure of the maze pattern {M}. 
    The parameters are the same as those of {tmaze_plot_maze};
    except that {M} is a {tmaze_pattern_t}, and {plot_tileset} 
    is a {tmaze_plot_tileset_proc_t}.
    !!! Unify with {tmaze_plot_maze} !!! */

typedef enum
  { tmaze_layer_HOLE = 0, /* Fill the cell's background. */
    tmaze_layer_GRID = 1, /* Draw the cell outline. */
    tmaze_layer_ROAD = 2, /* Paint the colored pathways in the cell. */
    tmaze_layer_CURB = 3, /* Draw the outlines of the cell's pathway. */
    tmaze_layer_DOTS = 4  /* Mark the ports with dots on the cell's boundary. */
  } tmaze_layer_t;
  /* A {tmaze_layer_t} value specifies one layer of a drawing. 
    The numeric values must be contiguous and in increasing order
    from background to foreground. */
  
#define tmaze_layer_MIN tmaze_layer_HOLE
#define tmaze_layer_MAX tmaze_layer_DOTS
  /* The first and last layer in plotting order. */

void tmaze_plot_cell
  ( epswr_figure_t *eps,
    int x, 
    int y,
    bool_t fill, 
    bool_t draw
  );
  /* Plots into {eps} the cell in column {x} and row {y} of a maze
    board. First, if {fill} is true, fills the cell with the
    current fill color of {eps}. Then, if {draw} is true, draws the
    outline with the current pen setting. */

void tmaze_plot_cross
  ( epswr_figure_t *eps,
    int x, 
    int y,
    double r,
    bool_t fill, 
    bool_t draw
  );
  /* Plots a diagonal pentomino cross centered on the cell at column
    [x} and row {y} of the maze. First, if {fill} is true, fills the
    cross with the current fill color of {eps}. Then, if {draw} is
    true, draws the outline with the current pen setting. */

void tmaze_plot_cell_ports
  ( epswr_figure_t *eps,
    int x, 
    int y,
    bool_t once,  /* If TRUE avoids paint the same port twice. */
    int val,      /* Valency (number of ports per side). */
    double r,     /* Dot radius in millimeters. */
    bool_t fill,  /* If TRUE fills the dots with the current fill color. */
    bool_t draw   /* If TRUE draws the dot outlines with the current pen. */
  );
  /* Plots into {eps} the {val} ports on each sides of the cell in
    column {x} and row {y} of a maze board. The ports are depicted as
    equally spaced dots drawn with {epswr_dot(eps,x,y,r,fill,draw)}.
    
    If {once} is true, tries to avoid plotting twice the same port.
    Namely, plots the ports on the South side only if {y == 0}, and
    the ports on the West side only if {x == 0}. */

void tmaze_plot_cell_side_ports
  ( epswr_figure_t *eps,
    int x, 
    int y,
    tmaze_dir_t dir, /* Which side to draw the dots on. */
    int val,         /* Valency (number of ports on this side). */
    double r,        /* Dot radius in millimeters. */
    bool_t fill,     /* If TRUE fills the dots with the current fill color. */
    bool_t draw      /* If TRUE draws the dot outlines with the current pen. */
  );
  /* Plots into {eps} the {val} ports on side {dir} of the cell in
    column {x} and row {y} of a maze board. The ports are depicted
    as equally spaced dots with {epswr_dot(eps,xc,yc,r,fill,draw)}. */

/* LOW_LEVEL PLOTTING TOOLS */

void tmaze_path_cell_make(double xp[], double yp[], int nc);
  /* Stores into {(xp[0..nc-1],yp[0..nc-1]} the corners of the outline
    of a cell of the maze. Assumes that the cell is a unit square
    centered at the origin. Fails if {nc} is not the correct number of
    corners (currently 4). */

void tmaze_path_cross_make(double xp[], double yp[], int nc, double r);
  /* Stores into {(xp[0..nc-1],yp[0..nc-1]} the corners of a diagonal
    pentomino cross with bars of half-width {r}. Assumes that the cell
    is a unit square centered at the origin. Fails if {nc} is not the
    correct number of corners (currently 12). */

void tmaze_path_rotate(double xp[], double yp[], int nc, int k);
  /* Rotates the path {(xp[0..nc-1],yp[0..nc-1])} by {k*90} degrees
    counterclockwise. In particular, if {k} is a multiple of 4,
    does nothing. */
    
void tmaze_path_scale(double xp[], double yp[], int nc, double sx, double sy);
  /* Scales the path {(xp[0..nc-1],yp[0..nc-1])} by multiplying the
    X coordinates by {sx} and the Y coordinates by {sy}. */
    
void tmaze_path_translate(double xp[], double yp[], int nc, double dx, double dy);
  /* Translates the path {(xp[0..nc-1],yp[0..nc-1])} by adding
    {(dx,dy)} to every corner. */
    
void tmaze_path_reverse(double xp[], double yp[], int nc);
  /* Translates the path {(xp[0..nc-1],yp[0..nc-1])} by adding
    {(dx,dy)} to every corner. */
    
void tmaze_path_quarter_circle_poly_make(double xp[], double yp[], int nc);
  /* Stores into {x[0..nc-1],yp[0..nc-1]} a four-corner polygonal
    approximation of a quarter-circle arc. The arc starts vertically
    at point {(1,0)} and ends horizontally at point {(0,1)}. Join 4 of
    these arcs to make a regular octagon with inscribed radius 1. */

void tmaze_path_quarter_circle_bezier_make(double xp[], double yp[], int nc);
  /* Stores into {x[0..nc-1],yp[0..nc-1]} the four control points
    of a Bézier arc that approximates a quarter-circle. The arc starts vertically
    at point {(1,0)} and ends horizontally at point {(0,1)}. Join 4 of
    these arcs to make an approximate circle with radius 1. */

/* LIMITS */

#define tmaze_MAX_SIZE ((INT32_MAX-1)/4)
  /* Max size of a maze, hence of any connected component thereof.
    The definition of `size' is family-dependent. */

#define tmaze_MAX_CELL_COUNT ((UINT32_MAX-1)/4)
  /* Max number of cells in maze. This limit ensures that there aren't
    too many vertices, even if the maze is a 1-wide strip with open
    topology. */
     
#define tmaze_MAX_CELL_INDEX (tmaze_MAX_CELL_COUNT - 1)
  /* The maximum valid {tmaze_cell_index_t} value in any maze. */

#define tmaze_cell_index_NONE (tmaze_MAX_CELL_COUNT)
  /* A {tmaze_cell_index_t} value that means `no such cell'. */

#endif
