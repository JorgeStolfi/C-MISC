/* SOGrid.h -- Abstract dyadic grids in arbitrary dimension */
/* Last edited on 2007-01-04 00:20:23 by stolfi */

#ifndef SOGrid_H
#define SOGrid_H

#include <dg_tree.h>
#include <dg_grid.h>

#include <vec.h>
#include <stdlib.h>

/* THE DYADIC CELL GRID */

/* The `dyadic cell tree' is an infinte binary tree of cells in 
  some {d}-dimensional Cartesian space {R^d}, whose axes 
  are numbered from {0} to {d-1}. 
  
  By definition, the nominal geometry of the root cell is an
  axis-aligned box, whose projection along axis {i} is the interval
  {[0 _ s^i]} where {s = (1/2)^(1/d)}. The two children of a box {b}
  are obtained by splitting {b} in half, perpendicularly to the axis
  of longest extent.  
  
  As we descend in the tree, the splitting axis changes cyclically.
  All cells with the same rank {r} constitute a regular rectangular
  grid, and have their longest extent along axis {r mod d}. The
  canonical cell of rank {r} is a copy of the root cell, scaled down
  by {s^r} and with its axes permuted cyclically {r} times.
  
  This nominal geometry is only a handy way to define the grid's
  topology, in particular to give meaning to concepts such as `axis of
  longest extent' and thus determine the direction of split at each
  rank. It is not necessarily the geometry used in an actual
  application, where the grid may be arbitrary distorted by some
  `shape function'. */
 
typedef unsigned char SOGrid_Degree;
  /* Degree of a polynomial, spline, etc.. */

/* FINITE BINARY TREES */
  
/* A `finite binary tree' is an any finite initial segment of the
  infinite binary tree, with the condition that every node has either
  zero or two children. The data structure below is meant to represent
  a finite binary tree. */

typedef struct SOGrid_Tree
  { dg_dim_t d;        /* Dimension of the grid's domain (positive). */
    dg_tree_node_t *root;   /* Root node of the tree. */
  } SOGrid_Tree;
  /* Header for a finite tree of {dg_tree_node_t}s.

    If we interpret each node as a dyadic cell, then the tree
    represents a hierarchical partition of the root cell into dyadic
    cells of various sizes. */
    
#define SOGRID_MAX_INDEX (1152921504606846975L)
#define SOGRID_MAX_TREE_HEIGHT (60)
#define SOGRID_MAX_RANK (SOGRID_MAX_TREE_HEIGHT-1)
  /* This package is designed to handle finite binary trees whose node
    indices are in the range {0 .. 2^60-1}. Thus the maximum height of
    an actual tree is 60, and the maximum node rank is 59. */

#define SOGRID_MAX_AXIS_HEIGHT(d) (SOGRID_MAX_TREE_HEIGHT/(d))
#define SOGRID_MAX_AXIS_RANK(d) (SOGRID_MAX_AXIS_HEIGHT(d)-1)
  /* If the tree is used to represent a grid in 4-dimensional space
    (for instance, 3D space + time), the finest achievable grid has
    {2^15 = 32768} cells along each axis. That is enough for 1m
    resolution in an field 30 km wide, together with 10 min resolution
    in a 200-day simulation period.

    In three dimensions (3D space, or 2D space + time), the finest
    achievable grid has {2^20 = 1048576} cells along each axis, which
    means 10 cm resolution in a field 100 km wide, together with 1 min
    resolution in a 700-day simulation period. */

SOGrid_Tree *SOGrid_Tree_new(dg_dim_t d);
  /* Creates a new {d}-dimensional tree with a single node (the root). */

void SOGrid_Tree_free(SOGrid_Tree *t);
  /* Recursively reclaims the whole tree {t}, including the header. */

void SOGrid_Tree_write(FILE *wr, SOGrid_Tree *t);
  /* Writes tree {t} to {wr}. */

SOGrid_Tree *SOGrid_Tree_read(FILE *rd);
  /* Reads an {SOGrid_Tree} from a file {rd} which was created by
    {SOGrid_Tree_Write}. */

#endif

