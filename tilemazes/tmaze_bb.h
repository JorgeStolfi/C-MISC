#ifndef tmaze_bb_H
#define tmaze_bb_H

/* Topology and geometry of the Blip-Blop maze. */
/* Last edited on 2023-02-03 23:44:20 by stolfi */

#include <stdint.h>
#include <limits.h>
#include <stdio.h>

#include <bool.h>
#include <epswr.h>
#include <frgb.h>

#include <tmaze.h>

/* BLIP-BLOP MAZES

  This interface describes a family of mazes nicknamed "Blip-Blop".
  
  A maze from this family is a rectangular array of square /cells/
  which are filled with two possible /tiles/ (actually the same tile
  in two different orientations). Painted on each tile is a diagonal
  line. In /blip/ tiles the wall runs NW to SE, in /blop/ tiles it
  runs SW to NE.
  
  The diagonal line can be interpreted as a wall that separates the
  ports on two adjacent sides from the ports on the other two sides.
  So, for example, in a blip tile the only connections are N--E and
  W--S, while in a blop tile they are N--W and E--S.
  
  The line can be seen also as a mirror, so that any light entering
  the cell perpendicularly through its North side will be refected by
  the mirror, either leftwards (blip) or rightwards (blop). Thus a
  photon launched into this maze, along some cardinal direction, will
  follow an angular path consisting of unit-length cardinal steps
  alternating with {-90} degree or {+90} degree turns. */

typedef enum
  { tmaze_bb_tile_BLIP = 0, /* Connections are N--E, S--W. */
    tmaze_bb_tile_BLOP = 1  /* Connections are N--W, S--E. */
  } tmaze_bb_tile_t;
  /* The two tile types in the Blip-Blop family are copies of
    the same tile, rotated 90 degrees. */
    
#define tmaze_bb_tileset_BOTH (3)
  /* The set of both Blip-Blop tile types. */
 
tmaze_bb_tile_t tmaze_bb_get_tile_from_dirs(tmaze_dir_t dir1, tmaze_dir_t dir2);
  /* Returns the tile that has a road connecting the side {dir1} to the side {dir2}. */

tmaze_t tmaze_bb_make(int32_t nx, int32_t ny, bool_t torus);
  /* Creates a Blip Blop tile maze with {nx} columns and {ny}
    rows of cells. The maze has toroidal topology iff {torus} is true.
    All cells will be filled with blip tiles. */
    
void tmaze_bb_random(tmaze_t *M, int32_t seed);    
  /* Intializes the random generator with {srandom(seed)} and then
    fills {M} with randomly chosen blip or blop tiles. */

tmaze_t tmaze_bb_from_string(char *string, bool_t torus);
  /* Creates a Blip Blop tile maze {M} with given {M.torus} attribute
    whose size and tiles are specified by {string}. The cells are
    filled row-by-row, from TOP TO BOTTOM and left to right. Codes:
    '>' (blip), '<' (blop), ';' (row terminator). */

tmaze_pattern_t tmaze_bb_pattern_from_string(char *string, bool_t torus, bool_t sub, int32_t val);
  /* Creates a Blip Blop tile maze pattern {M} with given attributes
    {M.torus,M.sub,M.val} whose size and tilesets are specified by
    {string}. The cells are filled row-by-row, from TOP TO BOTTOM and
    left to right. Codes: '>' (blip), '<' (blop), '*' (any), 
    'x' (none), ';' (row terminator). */
 
/* MAZE AND COMPONENT SIZES */

tmaze_size_t tmaze_bb_tot_size(tmaze_t *M);
  /* Returs the total size of any Blip-Blop maze with the dimensions
    and topology of {M}. This is simply the number of vertices of the
    graph {G} of {M}, which is the number of tile sides (counting each
    shared side as 1). */

tmaze_size_t tmaze_bb_max_comp_size(tmaze_t *M);
  /* Returs the maximum possible size of any component of the graph {G} 
    of any Blip-Blop maze with the dimensions and topology of {M}.
    This is the same as {tmaze_bb_tot_size}. */
 
/* PRINTING */

#define tmaze_bb_dir_chars "><"
  /* Maps a {tmaze_bb_tile_t} value to a character codes. */

void tmaze_bb_print(FILE *wr, tmaze_t *M);
  /* Prints the maze {M} to {wr}, showing the orientation of
    each tile by one of the letters '>' (blip) or '<' (blop). */

/* POSTSCRIPT PLOTTING */

void tmaze_bb_plot_maze
  ( char *dir,
    char *name,
    tmaze_t *M,
    bool_t cells,  /* TRUE paints the background of each tile. */
    bool_t grids,  /* TRUE paints the cell grid. */
    double tsize,  /* Cell size (mm) */
    bool_t curved  /* If TRUE uses curved roads, FALSE uses polygonal ones. */
  );
  /* Writes a picture of maze {M} as an Encapsulated
    Postscript file, called "{dir}/{name}.eps". */
     
void tmaze_bb_plot_maze_tile
  ( epswr_figure_t *eps,
    tmaze_t *M, 
    int32_t x, 
    int32_t y,
    bool_t curved,  /* If TRUE uses curved roads, FALSE uses polygonal ones. */
    frgb_t *rgb_S,  /* Fill color to use for rod that ends on South side. */
    frgb_t *rgb_N,  /* Fill color to use for rod that ends on North side. */
    bool_t fill,    /* If TRUE, paints the road pairs with the current fill color. */ 
    bool_t draw     /* If TRUE, strokes the road pairs's outline with the current pen. */
  );
  /* Plots into {eps} the tile assigned to the cell on column {x} and
    row {y} of maze {M}. First, if {fill} is true, fills the tile's
    representation with the proper fill color. Then, if
    {draw} is true, draws the outline with the current pen setting.
    
    The representation of the tile is two roads connecting the pairs
    of sides on the same side of the wall. Each road is filled with
    {rgb_S} or {rgb_N} depending on which horizontal side it is
    connected to. (If both colors are NULL, use the current fill color
    of {eps}.) Leaves the ends of the roads open unless it is an open
    (non-toroidal) board but not a sub-board. */

void tmaze_bb_plot_tile
  ( epswr_figure_t *eps,
    tmaze_bb_tile_t tile,
    int32_t x,
    int32_t y,
    bool_t curved,  /* If TRUE uses curved roads, FALSE uses polygonal ones. */
    frgb_t *rgb_S,  /* Fill color to use for rod that ends on South side. */
    frgb_t *rgb_N,  /* Fill color to use for rod that ends on North side. */
    bool_t open[],  /* Which road ends should be left open. */
    bool_t fill,    /* If TRUE, paints the road pairs with the current fill color. */
    bool_t draw     /* If TRUE, strokes the road pairs's outline with the current pen. */
  );
  /* Plots onto stream {eps} the road pairs of the tile on column {x} and
    row {y} of a maze. The tile will span a square of unit size with
    lowest corner {(x,y)}, in client plotting units. Assumes that the
    tile has type {tile}.

    If {fill} is true, fills the tile's pair of roads with the current fill
    color. Then, if {draw} is true, draws the bounding edges of the
    road pairs with the current pen style. For each {d} in {0..3}, the
    bounding edge of the road pairs that lies on the cell's {d}-side 
    (if any) is stroked only if {open[d]} is FALSE. */
   
void tmaze_bb_plot_pattern_tileset
  ( epswr_figure_t *eps,
    tmaze_pattern_t *M, 
    int32_t x, 
    int32_t y,
    bool_t curved,  /* If TRUE uses curved roads, FALSE uses polygonal ones. */
    bool_t fill,    /* If TRUE, paints the road pairs with the current fill color. */ 
    bool_t draw     /* If TRUE, strokes the road pairs's outline with the current pen. */
  );
  /* Plots into {eps} the tileset assigned to the cell on column {x}
    and row {y} of maze pattern {M}. First, if {fill} is true, fills
    the tile's representation with the current fill color of {eps}.
    Then, if {draw} is true, draws the outline with the current pen
    setting.
    
    The representation of the empty tileset is an X-shaped spot. The
    representation of a singleton tileset is the tile's pair of roads. The
    representation of {tmaze_tileset_BOTH} and its supersets is
    nothing. For any other tileset the procedure fails. */

void tmaze_bb_plot_cross
  ( epswr_figure_t *eps,
    int32_t x,
    int32_t y,
    bool_t fill,    /* If TRUE, paints the cross with the current fill color. */
    bool_t draw     /* If TRUE, strokes the cross's outline with the current pen. */
  );
  /* Plots onto stream {eps} a diagonal cross centered in cell {x,y} of
    a maze, style-wise compatible with the road pairs of normal tiles. */

/* LOW_LEVEL PLOTTING TOOLS */

void tmaze_bb_path_elbow_poly_make(double xp[], double yp[], int32_t nc, double r);
  /* Stores into {(xp[0..nc-1],yp[0..nc-1]} the corners of the outline
    of a polygonal road connecting the East and North sides of a cell.
    The roads of a blip cell are this road and a copy rotated 180 degrees.
    The roads of a blop cell are the same two roads, rotated 90 degrees. 
    
    Assumes that the tile is an axis-aligned unit square centered at
    the origin. The road will have half-width {r} and will be
    contained in the North-East half of the cell. The first corner is
    on the North side, West of center. Fails if {nc} is not the
    correct number of corners (currently 8). */

void tmaze_bb_path_elbow_bezier_make(double xp[], double yp[], int32_t nc, double r);
  /* Same as {tmaze_bb_path_elbow_poly_make}, but uses two Bezier arcs
    for the sides of the road. The vectors {x,y} are suitable for
    {epswr_bezier_polygon}. The curved sides (suitable for
    {epswr_curve}) are defined by the control points {0..3} and {4..7}.
    The road meets the East tile side at segment defined by points 3
    and 4, and the North side at points 7 and 0. Fails if {nc} is not
    the correct number of control points (currently 8). */

void tmaze_bb_get_open_ends_of_maze_tile
  ( tmaze_t *M, 
    tmaze_cell_index_t k, 
    bool_t open[]
  );
  /* Determines if the outline of any road in cell {k} of maze {M}
    should be left open (curb-less) where it touches
    the sides of the cell.  
    
    More precisely, for each direction {dir} in {0..3}, the procedure
    sets {open[dir]} to FALSE iff the road that ends on side {dir}
    should be closed off with a stroke of the pen. The end will be
    closed if and only if the cell {k} has no {dir}-neighbor. (Note
    that the result does not depend on the tile assigned to cell
    {k}.) */

void tmaze_bb_get_open_ends_of_pattern_tileset
  ( tmaze_pattern_t *M, 
    tmaze_cell_index_t k, 
    bool_t open[]
  );
  /* Analogous to {tmaze_bb_get_open_ends_of_maze_tile}, but for 
    cell {k} in a maze pattern {M} instead of a plain maze.  
    
    The procedure will set {open[dir]} to FALSE if an only if the cell
    {k} has no {dir}-neighbor, or when that neighbor has an empty
    tileset. (Note that the result does not depend on the
    tileset assigned to cell {k}.) */

#endif
