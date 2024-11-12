/* Last edited on 2014-05-19 16:46:59 by stolfilocal */

#include <stdint.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#include <bool.h>
#include <interval.h>
#include <utable.h>
#include <affirm.h>
#include <dg_grid.h>

#include <fast_mgrid.h>

fast_mgrid_level_t *fast_mgrid_refine(fast_mgrid_level_t *L0, dg_rank_t dk, int hw1)
  {
    dg_dim_t d = L0->d;
    dg_rank_t k0 = L0->k;        /* Rank (level index) of {L0} */
    int na0 = L0->ixw[1];        /* Number of active cells of {L0} */
    int nr0 = L0->ixw[L0->hw+1]; /* Number of cells in all layers. */
    
    demand(dk > 0, "rank increment must be positive");
    dg_rank_t k1 = (dg_rank_t )(k0 + dk);     /* Rank (level index) of {L1}. */
    demand((dk >= 1) && (dk <= d), "cannot refine by more than {d} levels");
    int nc = (1 << dk);  /* Number of children in {L1} of a cell of {L0}. */
    
    /* Create level {L1}, initially empty. */
    fast_mgrid_level_t *L1 = fast_mgrid_level_create(d, L0->D, k1, L0->nv, hw1, nc*na0, nc*nr0);
    
    /* Add the children of active cells, with zero details for now: */
    assert(L0->ixw[0] == 0);
    assert(L0->ixw[1] == na0);
    int na1 = 0;
    int ixa0;
    for (ixa0 = 0; ixa0 < na0; ixa0++)
      { uint64_t ida0 = L0->ix_to_id[ixa0]; /* Identifier of active cell of {L0}. */
        /* Add the children of cell {ida0}: */
        int ic;
        for (ic = 0; ic < nc; ic++) 
          { /* Get global identifier of child number {ic} of cell {ida0}: */
            uint64_t ida1 = ida0 * nc + ic;
            int32_t ixa1 = fast_mgrid_level_add_cell(L1, ida1, NULL);
            assert(ixa1 == na1); /* The children must be all distinct. */
            na1++;
          }
      }
    assert(L1->ixw[1] == na1);
    assert(utable_num_entries(L1->id_to_ix) == na1);
      
    /* Add the buffer cells, layer by layer, leaving their values undefined: */
    int r1;
    for (r1 = 1; r1 < hw1; r1++)
      { assert(L1->hw == r1-1);
        fast_mgrid_level_add_buffer_layer(L1);
        assert(L1->hw == r1);
      }
    assert(L1->hw == hw1);
    int nr1 = L1->ixw[hw1+1];
    assert(utable_num_entries(L1->id_to_ix) == nr1);  

    return L1;
  }

fast_mgrid_level_t *fast_mgrid_level_create(dg_dim_t d, interval_t *D, dg_rank_t k, int nv, int hw_max, int na_max, int nr_max)
  {
    fast_mgrid_level_t *L = notnull(malloc(sizeof(fast_mgrid_level_t)), "no mem");
    L->d = d;
    L->D = D;
    L->k = k;
    L->nv = nv;
    L->hw = 0;
    L->hw_max = hw_max;
    /* Index of first cell in each layer: */
    L->ixw = notnull(malloc((hw_max+2)*sizeof(int)), "no mem");
    int h;
    for (h = 0; h <= L->hw+1; h++) { L->ixw[h] = 0; }
    
    /* Active cells */
    L->na_max = na_max;
    L->node = notnull(malloc(na_max*sizeof(fast_mgrid_node_t *)), "no mem");
    L->det = notnull(malloc(nv*na_max*sizeof(double)), "no mem");
    
    /* Active and guard cells: */
    L->nr_max = nr_max;
    L->val = notnull(malloc(nv*nr_max*sizeof(double)), "no mem");
    L->ix_to_id = notnull(malloc(nr_max*sizeof(uint64_t)), "no mem");
    uint64_t nfull = (1 << k); /* Number of cells in complete grid of this level. */
    L->id_to_ix = utable_new(nfull); 
    
    return L;
  }

int32_t fast_mgrid_level_add_cell(fast_mgrid_level_t *L, uint64_t id, fast_mgrid_node_t *node)
  {
    int h = L->hw;  /* Index of layer where the cell will be inserted. */
    assert(h <= L->hw_max);
    int32_t ix = L->ixw[h+1]; /* Next local index after last cell of layer {h}. */
    demand(ix < L->nr_max, "no space for new cell");
    L->ixw[h+1]++;
    int iv;
    if (h == 0)
      { demand(ix < L->na_max, "no space for new active cell");
        L->node[ix] = node;
        for (iv = 0; iv < L->nv; iv++) { L->det[ix*L->nv + iv] = NAN; }
      }
    else
      { demand(node == NULL, "guard cells should not have tree nodes"); }
    for (iv = 0; iv < L->nv; iv++) { L->val[ix*L->nv + iv] = NAN; }
    L->ix_to_id[ix] = id;
    uint64_t id_min = (1 << L->k); /* Global identifier of first cell in full grid. */
    utable_set(L->id_to_ix, id - id_min, ix);
    return ix;
  }

void fast_mgrid_level_add_buffer_layer(fast_mgrid_level_t *L)
  { 
    dg_dim_t d = L->d;
    int h0 = L->hw; /* Index of last layer. */
    int h1 = h0 + 1; /* Index of new layer. */
    demand(h1 <= L->hw_max, "too many layers already");
    L->hw++;
    L->ixw[h1+1] = L->ixw[h1]; /* New layer starts empty. */
    /* Enumerate all cells of layer {h0}: */
    int ix0;
    dg_grid_size_t gsz[d]; /* Size of the full grid, along each axis. */
    dg_grid_size(d, L->k, gsz);
    for (ix0 = L->ixw[h0]; ix0 < L->ixw[h0+1]; ix0++)
      { uint64_t id0 = L->ix_to_id[ix0]; /* Global identifier of cell. */
        fast_mgrid_level_add_cell_neighbors(L, gsz, id0);
      }
  }
  
uint64_t fast_mgrid_level_get_cell_index(fast_mgrid_level_t *L, dg_cell_it_t id)
  {
    uint64_t id_min = (1 << L->k); /* Global identifier of first cell in full grid. */
    uint64_t ix = utable_get(L->id_to_ix, id - id_min);
    return ix;
  }

void fast_mgrid_level_add_cell_neighbors(fast_mgrid_level_t *L, dg_grid_size_t gsz[], uint64_t id)
  {
    demand(L->hw > 0, "active cells must be all present");
    dg_dim_t d = L->d;
    uint64_t id_min = (1 << L->k); /* Global identifier of first cell in full grid. */
    /* Enumerate the neighbors of the cell: */
    dg_grid_pos_t dp[d]; /* Position increments from cell {id} to its neighbors. */
    int ia;
    for (ia = 0; ia < d; ia++) { dp[ia] = -1; }
    while (TRUE) 
      { /* Get the global cell identifier {id1} of {dp}-neighbor: */
        uint64_t id1 = dg_cell_shift(d, id, dp);
        /* Try to find the cell with global identifier {id1} in {L}: */
        uint64_t ix1 = utable_get(L->id_to_ix, id1 - id_min);
        if (ix1 == utable_NULL)
          { /* This is a new cell: */
            ix1 = fast_mgrid_level_add_cell(L, id1, NULL);
          }
        /* Generate the next increment vector: */
        ia = d-1;
        while ((ia >= 0) && (dp[ia] == +1)) { dp[ia] = -1; ia--; }
        if (ia < 0) { break; }
        dp[ia] ++;
      }
  }
