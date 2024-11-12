#ifndef dg_grid_H
#define dg_grid_H

/* Abstract multi-dimensional dyadic grids */
/* Last edited on 2014-05-15 22:36:54 by stolfilocal */

#define _GNU_SOURCE
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>

#include <bool.h>
#include <sign.h>
#include <interval.h>
#include <psp_basic.h>
#include <box.h>
#include <vec.h>
#include <dg_basic.h>
#include <udg_grid.h>
  
/*
  THE DYADIC MULTIGRID

  An (/infinite dyadic/) /{d}-dimensional multigrid/ {G*} is the union
  of an infinite collection of {d}-dimensional toroidal grids -- the
  /levels/ or /layers/ of {G*} -- each identified by a natural number,
  its /rank/. 

  By definition, level 0 has a single cell, the /root cell/,
  the {d}-dimensional cube {U^d = (0_1)^d}. Every cell in level {r} contains
  exactly two cells of rank {r+1}, its /children/; and, except for the
  root, every cell is contained in exactly one cell of rank {r-1}, its
  /parent/.  The two children of a cell are /sisters/ to each other,
  and the root cell is its own sister by definition.
  
  The children of a cell are obtained by splitting {C} in half,
  perpendicularly to axis {r \bmod d}. The half-cell that lies towards
  {-oo} along that axis is the /low child/ or /0-/child, and the other
  half is the /high child/ or /1-child/.
  
  Note that all levels of a multigrid {G*} cover the same domain {D =
  \U G* = [0_1)^d}, namely the union of the root cell {U^d} plus all
  its inferior faces, with toroidal topology. */

typedef udg_rank_t dg_rank_t;
  /* Identifies a layer in a dyadic multigrid, or the depth of a node
    in a binary tree, etc. It may also be a relative rank (i.e.
    difference of two ranks). */

sign_t dg_rank_cmp(dg_rank_t *x, dg_rank_t *y);
  /* Returns {-1} if {*x < *y}, {0} if {*x == *y}, {+1} if {*x > *y}.
    Useful for sorting things by rank. */ 

dg_axis_t dg_split_axis(dg_dim_t d, dg_rank_t r);
  /* The axis perpendicular to the cut that divides a {d}-dimensional
    cell of level {r} into its children cells; namely, axis {r \bmod d}. */

  
typedef udg_grid_size_t dg_grid_size_t;
  /* Size of a toroidal grid along an axis, in multiples of 
    the corresponding grid step. */

void dg_grid_size(dg_dim_t d, dg_rank_t r, dg_grid_size_t gsz[]);
  /* Stores in {gsz[0..d-1]} the number of grid cells along axis {i} in
    level {r} of the {d}-dimensional multigrid. */

/* 
   CELL POSITIONS IN A LAYER */

typedef udg_grid_pos_t dg_grid_pos_t;
  /* Relative position of a cell within a 
    level of {G*}, along an axis, in multiples of 
    the corresponding grid step.  Also relative displacement
    between two cells. */
    
dg_grid_pos_t dg_max_grid_pos(dg_dim_t d, dg_rank_t r, dg_axis_t i);
  /* Maximum position along axis {i} of level {r} in a
    {d}-dimensional multigrid. */
    
#define dg_rank_MAX (udg_rank_MAX)
  /* Maximum rank of any level in a multigrid. The limit 
    is such that a {dg_cell_id_t} (below) will fit
    in one 64-bit word.

    In a 4-dimensional domain (for example, 3D space + time), the
    finest level has at least {2^15 = 32768} cells along each axis.
    That is enough for 1m resolution in a field 30 km wide, together
    with 5 min resolution in a 100-day simulation period.

    In a 3-dimensional domain (3D space, or 2D space + time), the
    finest level has {2^21 = 2097152} cells along each axis, which
    means 5 cm resolution in a field 100 km wide, together with 1 min
    resolution in a 1400-day simulation period.
    
    In a 2-dimensional domain (2D space, or 1D space + time), the
    finest level has at least {2^31 = 2147483648} cells along each
    axis, which means 1 mm resolution in a field 2000 km wide,
    together with 1 sec resolution in a 68-year simulation period.
    
    In a 1-dimensional domain (space or time), the finest level (rank
    62) has {2^62 = 4.6×10^18} cells, which means 10 picometer
    (sub-atomic) resolution over the Earth's circumference, or 1
    millisecond resolution in a 150,000-year simulation period. That
    should be enough for most practical purposes. */

/*
  CELL IDENTIFIERS
  
  We can uniquely identify a cell of a fixed multigrid {G*}, at any
  level, by a a single /cell identifier/ {id}. By definition, the root cell
  in layer 0 has identifier 1; and the children of cell {id} are the cells
  with identifiers {2*id} (low) and {2*id + 1} (high). The parent of cell
  {id} has therefore identifier {floor(id/2)}. Note that these formulas are
  independent of the grid's dimension {d}.
  
  These rules assign consecutive identifiers to all cells of each layer of
  {G*}; in layer {r} the identifiers range from {2^r} to {2^{r+1}-1},
  inclusive. Therefore, the position of the leading 1 bit in the identifier
  of a cell is its rank. The remaining bits are the bits of the cell's
  position vector, interleaved: that is, bit {2^j} of {pos[i]} is bit
  {2^{d*j + i}} of the cell identifier {id}, for each axis {i} and each bit
  position {j} such that {d*j + i < r}.  */
    
typedef udg_cell_id_t dg_cell_id_t;
  /* Identifies a cell in some multigrid. */

#define dg_NO_CELL (udg_NO_CELL)
  /* A cell identifier that means /no such cell/. */

#define dg_ROOT_CELL ((dg_cell_id_t)1)
  /* The identifier of the root cell of any multigrid. */

#define dg_cell_min_id(r) (dg_ROOT_CELL << (r))
#define dg_cell_max_id(r) ((dg_ROOT_CELL << (r+1)) - ((dg_cell_id_t)1))
  /* Cells in level {r} have conscutive identifiers in the range
    {[dg_cell_min_id(r) .. dg_cell_max_id(r)]}. The minimum identifier is
    the base cell of layer {r}, the maximum identifier is the cell
    adjacent to the superior corner of the domain. */

dg_rank_t dg_cell_rank(dg_cell_id_t id);
  /* The rank of the cell whose identifier is {id}. */
  
#define dg_check_cell_rank(id,r)   (((id)>>(r)) == ((dg_cell_id_t)1))
  /* TRUE iff cell {id} belongs to level {r}, */

void dg_cell_position(dg_dim_t d, dg_cell_id_t id, dg_grid_pos_t pos[]);
  /* For each axis {i} in {0..d-1}, stores in {pos[i]} the 
    position of cell {id} within its layer. */
    
dg_cell_id_t dg_cell_from_position(dg_dim_t d, dg_rank_t r, dg_grid_pos_t pos[]);
  /* Returns the identifier of the cell of level {r} whose position along axis {i}
    is {pos[i]}, for {i} in {0..d-1}. */
    
dg_cell_id_t dg_cell_axis_shift(dg_dim_t d, dg_cell_id_t id, dg_axis_t ax, dg_grid_pos_t dp);
  /* The identifier of the cell that is displaced from cell {id} by {dp}
    grid steps along axis {ax}. Takes into account the toroidal
    topology of each layer. */
    
dg_cell_id_t dg_cell_shift(dg_dim_t d, dg_cell_id_t id, dg_grid_pos_t dp[]);
  /* Returns the identifier of the cell in a {d}-dimensional multigrid
    whose position vector is {dp[0..d-1]} relative to the inferior
    corner of cell {id}, within the layer {r} that contains {id}. Takes
    into account the toroidal topology of each layer. */
    
dg_cell_id_t dg_cell_id_add(dg_dim_t d, dg_cell_id_t ida, dg_cell_id_t idb);
  /* The identifier of the cell whose position vector is the sum of the
    position vectors corresponding to identifiers {ida} and {idb}. The 
    two cells should have the same rank. */

dg_cell_id_t dg_cell_neighbor(dg_dim_t d, dg_cell_id_t id, box_face_index_t fi);
  /* The cell {id'} of rank as {id}, such that the highest-dimensional
    face shared by {id} and {id'} is face {fi} of {id}. In particular, if
    {fi == 0}, returns {id} itself. */

/* 
  ROOT-RELATIVE COORDINATES 
  
  The coordinates of a point {x} of {D} are said to be
  /root-relative/, because in them the inferior and superior corners of the
  root cell {U^d} are {(0,..0)} and {(1,..1}). Root-relative
  coordinates that differ by an integer are implicitly equivalent, and
  therefore they are often reduced modulo 1 to real fractions in the
  range {[0_1)}. */

void dg_cell_box_root_relative(dg_dim_t d, dg_cell_id_t id, interval_t B[]);
  /* Stores in {B[0..d-1]} the inferior and superior coordinates of the cell
    {id}, relative to the root cell, in a {d}-dimensional grid. The
    inferior bounds will be in the range {[0_1)}; the superior bounds will
    be strictly greater than the inferior bounds, in the range
    {(0_1]}. */

/*
  THE CANONICAL GEOMETRY

  In the /canonical geometry/ of the dyadic hierarchy, 
  the grid of level 0 includes the /canonical root cell/ {C0},
  whose projection along axis {i} is the interval {C0_i = [0 _ s^i]}
  where {s = (1/2)^(1/d)}. It follows that every cell of any rank {r}
  is a copy of {C0}, uniformly scaled by a factor {s^r}, and with its
  axes cyclically permuted {r} times.
  
  The main purpose of the canonical geometry is to justify and
  visualize the topological features of the grid, in particular the
  cyclic permutation of the splitting axis as one descends in the
  multigrid. It is not necessarily the geometry used in an actual
  application, where the grid may be arbitrary distorted by some
  `shape function'. */

const double *dg_cell_root_canonical_sizes(dg_dim_t d);
  /* Returns a pointer to a (statically allocated) vector {sz[0..d-1]}
    such that {sz[i]} is the extent along axis {i} of the canonical
    root cell of the {d}-dimensional grid. Clients should treat that
    vector as read-only. */

void dg_cell_box_canonical(dg_dim_t d, dg_cell_id_t id, interval_t B[]);
  /* Stores in {B[0..d-1]} the inferior and superior coordinates
    of the cell {id} in a {d}-dimensional grid. */
 
/* PRINTOUT */

void dg_cell_print(FILE *wr, dg_dim_t d, dg_cell_id_t id);
  /* Prints the cell {id} in the format 
    "{RANK}:({POS[0]},..{POS[d-1]})" */
   
/* VECTORS OF CELL IDENTIFIERS */

vec_typedef(dg_cell_id_vec_t,dg_cell_id_vec,dg_cell_id_t);
     /* A vector of {dg_cell_id_t}. */

#endif

