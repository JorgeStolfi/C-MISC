#ifndef tmaze_ba_H
#define tmaze_ba_H

/* Topology and geometry of the Brasilia Airport maze. */
/* Last edited on 2024-12-21 11:31:12 by stolfi */

#include <stdint.h>
#include <limits.h>
#include <stdio.h>

#include <bool.h>
#include <epswr.h>
#include <jsrandom.h>

#include <tmaze.h>

/* BRASILIA AIRPORT MAZES

  This interface describes a family of mazes inspired by a tile
  mural at the Brasilia Airport.  
  
  A maze from this family is a rectangular array of square /cells/
  which are filled with four possible /tiles/ (actually the same tile
  in four different orientations). Painted on each tile is a
  'T'-shaped road that connects the mid-sections of three of the four
  sides of the tile. The T-road in each tile can be oriented in any of
  the four possible ways. The T-roads of tiles in adjacent cells of
  the maze may or may not be connected across the shared side,
  depending on the orientation of the two tiles. */

typedef enum
  { tmaze_ba_tile_E = tmaze_dir_E, 
    tmaze_ba_tile_N = tmaze_dir_N, 
    tmaze_ba_tile_W = tmaze_dir_W, 
    tmaze_ba_tile_S = tmaze_dir_S 
  } tmaze_ba_tile_t;
  /* The four tile types in the Brasilia Airport family are copies of
    the same tile (a T-shaped road) rotated in the four possible ways.
    
    The name of each tile type indicates the /orientation/ of the tile,
    namely the cardinal direction of the tile's side that is NOT
    served by the T-road. Namely,
    
      {tile == tmaze_ba_tile_N} : not connected to cell {x,y+1}.
      {tile == tmaze_ba_tile_S} : not connected to cell {x,y-1}.
      {tile == tmaze_ba_tile_W} : not connected to cell {x-1,y}.
      {tile == tmaze_ba_tile_E} : not connected to cell {x+1,y}.
      
    The numerical value of each type must be the same as the
    {tmaze_dir_t} of its orientation. */
    
#define tmaze_ba_tileset_FOUR (15)
  /* The set of all four Brasilia Airport tile types. */

tmaze_t tmaze_ba_make(int32_t nx, int32_t ny, bool_t torus);
  /* Creates a Brasilia Airport tile maze with {nx} columns and {ny}
    rows of cells. The maze has toroidal topology iff {torus} is true.
    All cells will be filled with {tmaze_ba_tile_E} tiles */
    
void tmaze_ba_random(tmaze_t *M, int32_t seed);    
  /* Intializes the random generator with {srandom(seed)} and then
    fills {M} with randomly chosen tiles. */

tmaze_t tmaze_ba_from_string(char *string, bool_t torus);
  /* Creates a Brasilia Airport tile maze {M} with given {M.torus} attribute
    whose size and tiles are specified by {string}. The cells are
    filled row-by-row, from TOP TO BOTTOM and left to right. Codes:
    'N', 'S', 'E', 'W' (tiles), ';' (row terminator). */

tmaze_pattern_t tmaze_ba_pattern_from_string(char *string, bool_t torus, bool_t sub, int32_t val);
  /* Creates a Brasilia Airport tile maze pattern {M} with given attributes {M.torus,M.sub,M.val}
    whose size and tilesets are specified by {string}. The cells are
    filled row-by-row, from TOP TO BOTTOM and left to right. Codes:
    'N', 'S', 'E', 'W' (specific tiles), '*' (any), 'x' (none), 
    ';' (row terminator). */
    
/* PRINTING */

#define tmaze_ba_tile_chars "NWSE"
  /* Maps a {tmaze_ba_dir_t} value to a character codes. */

void tmaze_ba_print(FILE *wr, tmaze_t *M);
  /* Prints the maze {M} to {wr}, showing the orientation of
    each tile by one of the letters 'N','W','S','E'. */

/* MAZE AND COMPONENT SIZES */

tmaze_size_t tmaze_ba_tot_size(tmaze_t *M);
  /* Returs the total size of any Brasilia maze
    with the dimensions and topology of {M}.  That is simply the
    number of cells in {M}. */

tmaze_size_t tmaze_ba_max_comp_size(tmaze_t *M);
  /* Returs the maximum possible size of any component of
    any Brasilia maze with the dimensions and topology of {M}.
    That is simply the number of cells in {M}. */

/* POSTSCRIPT PLOTTING */

void tmaze_ba_plot_maze
  ( char *dir,
    char *name,
    tmaze_t *M,
    bool_t cells,  /* TRUE paints the background of each tile. */
    bool_t grids,  /* TRUE draws the cell gridlines. */
    double tsize   /* Cell size (mm) */
  );
  /* Writes a picture of maze {M} as an Encapsulated
    Postscript file, called "{dir}/{name}.eps". */
     
void tmaze_ba_plot_maze_tile
  ( epswr_figure_t *eps,
    tmaze_t *M, 
    int32_t x, 
    int32_t y,
    bool_t fill, 
    bool_t draw
  );
  /* Plots into {eps} the tile assigned to the cell on column {x} and
    row {y} of maze {M}. First, if {fill} is true, fills the tile's
    representation with the current fill color of {eps}. Then, if
    {draw} is true, draws the outline with the current pen setting.
    
    The representation of the tile is its T-road. Leaves the 
    ends of the road open or closed depending on the neighboring 
    tiles. */
   
void tmaze_ba_plot_pattern_tileset
  ( epswr_figure_t *eps,
    tmaze_pattern_t *M, 
    int32_t x, 
    int32_t y,
    bool_t fill,    /* If TRUE, paints the T-road with the current fill color. */ 
    bool_t draw     /* If TRUE, strokes the T-road's outline with the current pen. */
  );
  /* Plots into {eps} the tileset assigned to the cell on column {x}
    and row {y} of maze pattern {M}. First, if {fill} is true, fills
    the tile's representation with the current fill color of {eps}.
    Then, if {draw} is true, draws the outline with the current pen
    setting.
    
    The representation of the empty tileset is an X-shaped spot. The
    representation of a singleton tileset is the tile's T-road. The
    representation of {tmaze_tileset_FOUR} and its supersets is
    nothing. For any other tileset the procedure fails. */

void tmaze_ba_plot_tile
  ( epswr_figure_t *eps,
    tmaze_ba_tile_t tile,
    int32_t x,
    int32_t y,
    bool_t open[],  /* Which T-road ends should be left open. */
    bool_t fill,    /* If TRUE, paints the T-road with the current fill color. */
    bool_t draw     /* If TRUE, strokes the T-road's outline with the current pen. */
  );
  /* Plots onto stream {eps} the T-road of the tile on column {x} and
    row {y} of a maze. The tile will span a square of unit size with
    lowest corner {(x,y)}, in client plotting units. Assumes that the
    tile has type {tile}.

    If {fill} is true, fills the tile's T-road with the current fill
    color. Then, if {draw} is true, draws the bounding edges of the
    T-road with the current pen style. For each {d} in {0..3}, the
    bounding edge of the T-road that lies on the cell's {d}-side 
    (if any) is stroked only if {open[d]} is FALSE. */

void tmaze_ba_plot_cross
  ( epswr_figure_t *eps,
    int32_t x,
    int32_t y,
    bool_t fill,    /* If TRUE, paints the cross with the current fill color. */
    bool_t draw     /* If TRUE, strokes the cross's outline with the current pen. */
  );
  /* Plots onto stream {eps} a diagonal cross centered in cell {x,y} of
    a maze, style-wise compatible with the T-road of normal tiles. */

/* LOW_LEVEL PLOTTING TOOLS */

void tmaze_ba_path_tee_make
  ( double xp[], 
    double yp[], 
    int32_t nc, 
    double r
  );
  /* Stores into {(xp[0..nc-1],yp[0..nc-1]} the corners of the outline
    of a tile's T-road, in the upright ({tmaze_dir_N}) orientation.
    The outline starts with corner 0 on the South side of the tile, at
    the lower right corner of the stem of the `T', and proceeds
    counterclockwise until the lower left corner of the stem,
    again on the South side.
    
    Assumes that the tile is an axis-aligned unit square centered at
    the origin. The three arms of the T-road will have half-width {r}.
    Fails if {nc} is not the correct number of corners (currently 8). */

void tmaze_ba_get_open_ends_of_maze_tile
  ( tmaze_t *M, 
    tmaze_cell_index_t k, 
    bool_t open[]
  );
  /* Determines if the outline of the T-road in cell {k} of maze {M}
    should be left open (curb-less) where it touches
    the sides of the cell.  
    
    More precisely, for each direction {dir} in {0..3}, the procedure
    sets {open[dir]} to FALSE iff the leg or the T-road that points in
    that direction should be closed off with a stroke of the pen. The
    end will be closed if and only if the cell {k} has no
    {dir}-neighbor, or when that neighbor's tile is oriented so that
    it cannot connect with the tile in cell {k}. (Note that the result
    does not depend on the tile assigned to cell {k}.) */

void tmaze_ba_get_open_ends_of_pattern_tileset
  ( tmaze_pattern_t *M, 
    tmaze_cell_index_t k, 
    bool_t open[]
  );
  /* Analogous to {tmaze_ba_get_open_ends_of_maze_tile}, but for 
    cell {k} in a maze pattern {M} instead of a plain maze.  
    
    The procedure will set {open[dir]} to FALSE if an only if the cell
    {k} has no {dir}-neighbor, or when that neighbor has an empty
    tileset, or when that neighbor's tileset has a single tile that is
    oriented so that it cannot connect with any tile that one could
    place in cell {k}. (Note that the result does not depend on the
    tileset assigned to cell {k}.) */

#endif
