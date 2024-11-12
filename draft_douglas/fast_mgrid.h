#ifndef fast_mgrid_H
#define fast_mgrid_H

/* Fast implementation of multigrid traversal with stencils. */
/* Last edited on 2014-05-19 16:51:19 by stolfilocal */

#include <stdint.h>

#include <bool.h>
#include <interval.h>
#include <utable.h>
#include <dg_basic.h>
#include <dg_grid.h>

typedef struct fast_mgrid_node_t { } fast_mgrid_node_t;
  /* Tree node; probably not needed. */

typedef struct fast_mgrid_level_t
  { dg_dim_t d;    /* Domain dimension. */
    interval_t *D; /* Domain box, spans the interval {D[i]} on axis {i}, for {i} in {0..d-1}. */
    dg_rank_t k;   /* Level index (0 = root). Divisible by {d} if quad-grid. */
    int nv;        /* Number of problem variables per cell. */   
    int hw;        /* Half-width of buffer zone (in cells). */
    int hw_max;    /* Max value of {hw} ({ixw} is allocated with {hw_max+2} elements). */
    int *ixw;      /* Number of cells at each grid distance from the active cells, indexed {0..hw+1}. */
    /* Active cells */
    int na_max;    /* Number of allocated entries in {det,node}. */
    fast_mgrid_node_t **node;  /* Tree nodes of active cells, indexed {0..na-1} where {na=ixw[1]}. */
    double *det;   /* Details of active cells, indexed {0..nv*na-1}. */
    /* Active and guard cells: */
    int nr_max;    /* Number of allocated entries in {val,ic_to_id}. */
    double *val;   /* Values problem variables in relevant cells, indexed {0..nv*nr-1} where {nr=ixw[hw+1]}. */
    uint64_t *ix_to_id;  /* Maps local indices of relevant cells to cell identifiers; indexed {0..nr-1}. */
    utable_t *id_to_ix;  /* Hash table that maps relative identifiers of relevant cells to local cell indices. */    
  } fast_mgrid_level_t;
  /* 
    Representation of level {k} of an adaptive hierarchical
    multi-resolution grid (dyadic or quad).
    
    The /global identifier/ of a cell in level {k} varies from {2^k} to
    {2^{k+1}-1}. The root has identifier 1; zero or negative is an
    invalid identifier.  The cells of the full grid of level {k} have
    identifiers {2^k .. 2^{k+1}-1}. The /relative identifier/ of a cell
    is its global identifier minus {2^k}.
    
    The /local index/ of a cell is a number in {0..nr-1} that designates
    the cell within the set {R} of relevant cells.
    
    The /active cells/ are an arbitrary subset of the complete level
    {k}.
    
    The /position/ of a cell {C} is a vector {pos[0..d-1]} of natural
    numbers such that {pos[i]} is the number of cells with lower
    coordinate along the unidimensional slice of the complete level that
    contains {C}.
    
    The /grid distance/ between two cells is the maximum absolute value
    of their positions along any axis.
    
    The /relevant cells/ include all the cells of the complete level {k}
    whose grid distance from the nearest active cell is at most {hw}.
    They include the active cells, and are the union of all the
    {hw}-packs of all the active cells, where the /{hw}-pack/ of a cell
    is a sub-array of the cells of the complete level {k+dk} with
    {2*hw+1} cells along each axis, with toroidal wrap-around.
    
    More specifically, the cells that are at grid distance {t} from the
    nearest active cell, for each {t} in {0..hw}, have local indices in
    the range {ixw[t]..ixw[t+1]-1}. In particular, {ixw[0] = 0}, all
    relevant cells have local indices in {0..nr-1} where {nr =
    ixw[hw+1]}, and the active cells have indices in the range {0..na-1}
    where {na = ixw[1]}.
    
    Let {C} be a relevant cell with local index {j} (in {0..nr-1}). The
    global identifier of that cell is {ix_to_id[j]}, and the value of the
    problem variable with index {r} is {val[j*nv + r]}, for {r} in
    {0..nv-1}. If {C} is active (that is, {j < na}), then its tree node
    is {node[j]}, and its detail values are {det[j*nv + r]}.

    The relevant cells presumably include all the cells necessary to
    apply the numerical methods for the active cells, including the
    interpolation filter and the problem-specific discrete differential
    operator. */

uint64_t fast_mgrid_level_get_cell_index(fast_mgrid_level_t *L, dg_cell_it_t id);
  /* Returns the index {ix} in {L} of the cell whose global identifier if {id}.
    If {ix == utable_NULL}, the cell is not in {L} (not even as a guard cell). 
    Otherwise, the values of that cell are in {L->val[ix*L->nv + k]}, for {k} in {0..nv-1}.
    If {ix < L->ixw[0]}, then the cell is active, and its details are in 
    {L->det[ix*L->nv + k]}, for {k} in {0..nv-1}. */

fast_mgrid_level_t *fast_mgrid_level_create(dg_dim_t d, interval_t *D, dg_rank_t k, int nv, int hw, int na_max, int nr_max);
  /* Allocates a {fast_mgrid_level_t} structure for some level {k} of a multiscale mesh.
    The structure will have space for maximum {na_max} active cells, and
    {nr_max} total cells (active cells plus guard cells).
    The guard cells can be at most {hw} away from the active cells, in any axis.
    Initially it has no cells. */

fast_mgrid_level_t *fast_mgrid_refine(fast_mgrid_level_t *L0, dg_rank_t dk, int hw1);
  /* Given the representation {L0} of some level {k} for an adaptive
    hierarchical multigrid structure, creates a representation {L1} for
    the level {k+dk}; where {dk=1} for the dyadic grid, and {dk=d} for a
    quad-grid, where {d} is the dimension of the domain.
    
    The active cells of {L1} are all the children of the active cells of
    {L0} that have at least one child in the tree. The relevant cells of
    {L1} will be all cells that are at grid distance {hw1} from the
    nearest active cell. 
    
    The value of {hw1} must be at most {2*hw}. */

int32_t fast_mgrid_level_add_cell(fast_mgrid_level_t *L, uint64_t id, fast_mgrid_node_t *node);
  /* Adds a cell with global identifier {id} to the level {L} of a multiscale mesh.
    Assumes that the cell is in the last layer {L->hw} of guard cells -- namely, 
    at grid distance {L->hw} from the nearest active cell.
    The problem variables are all set to {NAN}. If the cell is active
    (i.e. {L->hw == 0}), then the details are set to {NAN} and the tree node
    pointer is set to {node}.  Otherwise {node} must be {NULL}. */

void fast_mgrid_level_add_buffer_layer(fast_mgrid_level_t *L);
  /* Adds to level {L} of an incomplete multiscale grid a new layer of guard celss, namely
    all cells that are at grid distance {L->hw+1} from the active cells. 
    Fails if {L->hw} is already equal to {L->hw_max}. */

void fast_mgrid_level_add_cell_neighbors(fast_mgrid_level_t *L, dg_grid_size_t gsz[], uint64_t id);
  /* Adds to layer {L} of an adaptive multiscale mesh 
    the buffer cells that are chess-king neighbors of (that is, at grid distance 1 from)
    the cell with global identifier {id} and are not already in {L}. 
    
    Assumes that the cell {id} is in layer {L->hw-1} (the next-to-last layer) of buffer
    cells, and that all buffer cells in layers {0..L->hw-1} are already in {L}.
    Therefore, any cell that is a neighbor of cell {id} and is not yet in {L}
    must be in the last layer {L->hw}. */

#endif
