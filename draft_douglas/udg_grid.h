/* Abstract unidimensional dyadic grids. */
/* Last edited on 2014-05-15 22:38:35 by stolfilocal */

#ifndef udg_grid_H
#define udg_grid_H

#define _GNU_SOURCE
#include <stdint.h>
#include <stdio.h>

#include <interval.h>
#include <bool.h>
#include <sign.h>
#include <vec.h>
#include <psp_grid.h>

/*
  THE UNIDIMENSIONAL DYADIC MULTIGRID

  An (/infinite/) dyadic uni-dimensional multigrid/ {G*} is the union
  of an infinite collection of unidimensional circular grids -- the
  /levels/ or /layers/ of {G*} -- each identified by a natural number,
  its /rank/. 

  By definition, level 0 has a single cell, the /root cell/. Every
  cell in level {r} contains exactly two cells of rank {r+1}, its
  /children/; and, except for the root, every cell is contained in
  exactly one cell of rank {r-1}, its /parent/. The two children of a
  cell are /sisters/ to each other, and the root cell is its own
  sister by definition.
  
  The children of a cell are obtained by splitting {C} in half. The
  half-cell that lies towards {-oo} along that axis is the /low/,
  /inferior/, or /0-/child; and the other half is the /high/,
  /superior/, or /1-/child.
  
  Note that all levels of {G*} cover the same
  domain {D = \U G*}, namely the union of the root cell plus its
  inferior vertex, with circular topology. */

typedef int32_t udg_rank_t;
  /* Identifies a layer of {G*}. It may also be a relative rank (i.e.
    difference of two ranks). */

sign_t udg_rank_cmp(udg_rank_t *x, udg_rank_t *y);
  /* A comparison function for ranks. */

/*
  SIZE OF LEVELS */
 
typedef psp_grid_size_t udg_grid_size_t;
  /* Number of cells in some level of {G*}. */

udg_grid_size_t udg_grid_size(udg_rank_t r);
  /* The number of cells in level {r} of {G*}, namely {2^r}. */

/*
  CELL POSITIONS INSIDE EACH LEVEL
  
  Each cell of {G*} can be identified by its rank {r} and a
  position {p} within that level. Because of the assumed circular
  topolgy, any two positions that differ by a multiple of {2^r} are
  equivalent. */
 
typedef psp_grid_pos_t udg_grid_pos_t;
  /* Position of a cell within some level of {G*}
    Also a relative displacement between two cells
    in the same level. */
    
udg_grid_pos_t udg_max_grid_pos(udg_rank_t r);
  /* Maximum position of a cell in level {r} of {G*}, namely {2^r-1}. */
    
#define udg_rank_MAX (64 - 1)
  /* Maximum rank of any level in a multigrid. The limit 
    is such that a {udg_cell_id_t} (below) will fit in one 64-bit word.
    
    With this definition, the finest level (rank 62) has {2^62 =
    4.6×10^18} cells, which means 10 picometer (sub-atomic) resolution
    over the Earth's circumference, or 1 millisecond resolution in a
    150,000-year simulation period. That should be enough for most
    practical purposes. */

/*
  CELL IDENTIFIERS
  
  We can uniquely identify a cell of a dyadic multigrid {G*}, at any
  level, by a a single /cell identifier/ {id}. By definition, the root cell
  in layer 0 has identifier 1; and the children of cell {id} are the cells
  with identifiers {2*id} (low) and {2*id + 1} (high). The parent of cell
  {id} has therefore identifier {floor(id/2)}. Note that these formulas are
  independent of the grid's dimension {d}.
  
  These rules assign consecutive identifiers to all cells of each layer of
  {G*}; in layer {r} the identifiers range from {2^r} to {2^{r+1}-1},
  inclusive. Therefore, the position of the leading 1 bit in the identifier
  of a cell is its rank. The remaining bits are the bits of the cell's
  position within that layer. */
    
typedef uint64_t udg_cell_id_t;
  /* Identifies a cell in some multigrid. */

#define udg_NO_CELL ((dg_cell_id_t)0)
  /* A cell identifier that means /no such cell/. */

#define udg_ROOT_CELL ((dg_cell_id_t)1)
  /* The identifier of the root cell of any multigrid. */

#define udg_cell_min_id(r) (udg_ROOT_CELL << (r))
#define udg_cell_max_id(r) ((udg_ROOT_CELL << (r+1)) - ((udg_cell_id_t)1))
  /* Cells in level {r} have conscutive identifiers in the range
    {[udg_cell_min_id(r) .. udg_cell_max_id(r)]}. The minimum identifier is
    the base cell of layer {r}, the maximum identifier is the cell
    adjacent to the superior corner of the domain. */

udg_rank_t udg_cell_rank(udg_cell_id_t id);
  /* The rank of the cell whose identifier is {id}. */
  
#define udg_check_cell_rank(id,r) (((id)>>(r)) == ((dg_cell_id_t)1))
  /* TRUE iff cell {id} belongs to level {r}, */

udg_grid_pos_t udg_cell_position(udg_cell_id_t id);
  /* Returns the position {pos} of cell {id} within its layer. */
    
udg_cell_id_t udg_cell_from_position(udg_rank_t r, udg_grid_pos_t pos);
  /* Returns the identifier of the cell of level {r} whose position is {pos}. */
    
udg_cell_id_t udg_cell_shift(udg_cell_id_t id, udg_grid_pos_t dp);
  /* The identifier of the cell that is displaced from cell {id} by {dp}
    grid steps. Takes into account the circular topology of each layer. */

/* VECTORS OF CELL IDENTIFIERS */

vec_typedef(udg_cell_id_vec_t,udg_cell_id_vec,udg_cell_id_t);
  /* A vector of {udg_cell_id_t}. */

/* 
  ROOT-RELATIVE ARGUMENT
  
  A point {x} of {D} can identified by its /root-relative position/,
  which is 0 at the low end and 1 at the high end of the root cell.
  Root-relative coordinates that differ by an integer are implicitly
  equivalent, and therefore they are often reduced modulo 1 to real
  fractions in the range {[0_1)}. */

void udg_cell_interval_root_relative(udg_cell_id_t id, interval_t *B);
  /* Stores in {B[0..d-1]} the inferior and superior root-relative 
    coordinates of the cell {id}, in a {d}-dimensional grid. The
    inferior bounds will be in the range {[0_1)}; the superior bounds will
    be strictly greater than the inferior bounds, in the range
    {(0_1]}. */

/* 
  PRINTOUT */

void udg_cell_print(FILE *wr, udg_cell_id_t id);
  /* Prints the cell {id} in the format "{RANK}:{POS}" */

#endif
