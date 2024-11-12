/* See dg_grid.h */
/* Last edited on 2014-05-15 22:43:14 by stolfilocal */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
  
#include <affirm.h>
#include <bool.h>
#include <ix.h>

#include <udg_grid.h>
#include <dg_grid.h>

box_axis_set_t dg_axis_set_complement(dg_dim_t d, dg_axis_set_t A)
  {  return set32_difference(set32_range(0,d-1), A); }

/* THE DYADIC MULTIGRID */

dg_axis_t dg_split_axis(dg_dim_t d, dg_rank_t r)
  { return (dg_axis_t)(r % ((dg_rank_t)d)); }
  
#define dg_GRID_POS_BITS(d,r,j)  (((r)+(d)-1-(j))/(d))
  /* Number of bits needed to represent a grid position within level {r}
    along axis {j} of a {d}-dimensional multigrid. */

void dg_grid_size(dg_dim_t d, dg_rank_t r, dg_grid_size_t gsz[])
  { int i;
    for (i = 0; i < d; i++) { gsz[i] = (((dg_grid_size_t)1) << dg_GRID_POS_BITS(d,r,i)); }
  }

dg_grid_pos_t dg_max_grid_pos(dg_dim_t d, dg_rank_t r, dg_axis_t i)
  { return (dg_grid_pos_t)(((dg_grid_size_t)1) << dg_GRID_POS_BITS(d,r,i)) - ((dg_grid_size_t)1); }

sign_t dg_rank_cmp(dg_rank_t *x, dg_rank_t *y)
  { if ((*x) < (*y))
      { return -1; }
    else if ((*x) > (*y))
      { return +1; }
    else
      { return 0; }
  }
  
/* CELL IDENTIFIERS */
  
dg_rank_t dg_cell_rank(dg_cell_id_t id)
  { int r = 0;
    while (id > (((dg_cell_id_t)1) << 16)) { id >>= 16; r+= 16; }
    while (id > (((dg_cell_id_t)1) << 4)) { id >>= 4; r+= 4; }
    while (id > ((dg_cell_id_t)1)) { id >>= 1; r+= 1; }
    return r;
  }

void dg_cell_position(dg_dim_t d, dg_cell_id_t id, dg_grid_pos_t pos[])
  { dg_cell_id_t t = ((((dg_cell_id_t)1) << dg_cell_rank(id)) >> 1);
    int ax;
    for (ax = 0; ax < d; ax++) { pos[ax] = 0; }
    ax = 0;
    while (t != 0)
      { pos[ax] <<= 1; if ((id & t) != 0) { pos[ax] |= 1;  }
        ax++; if (ax >= d) { ax = 0; }
        t >>= 1;
      }
  }
  
dg_cell_id_t dg_cell_from_position(dg_dim_t d, dg_rank_t r, dg_grid_pos_t pos[])
  { /* First cell of layer {r}: */
    dg_cell_id_t id = (((dg_cell_id_t)1) << r);
    /* Limit for bit values that can be added to {id}: */
    dg_cell_id_t tlim = (((dg_cell_id_t)1) << r);
    int ax;
    for (ax = 0; ax < d; ax++) 
      { dg_grid_pos_t posax = pos[ax]; 
        /* Add {posax} to {id}, starting at bit {ax} with stride {d}: */
        dg_axis_t ay = dg_split_axis(d, r + d - ax);
        uint64_t t = (((dg_cell_id_t)1) << ay);
        while ((posax != 0) && (t < tlim))
          { if (posax & 1) 
              { if ((id & t) != 0) 
                  { id &= ~t; posax += 2; }
                else
                  { id |= t; }
              }
            posax >>= 1; t <<= d;
          }
      }
    return id;
  }

dg_cell_id_t dg_cell_axis_shift(dg_dim_t d, dg_cell_id_t id, dg_axis_t ax, dg_grid_pos_t dp)
  { dg_rank_t rank = dg_cell_rank(id);
    /* Limit for bit values that can be added to {id}: */
    dg_cell_id_t tlim = (((dg_cell_id_t)1) << rank);
    /* Add {dp} to {id}, starting at bit {ax} with stride {d}: */
    dg_axis_t ay = dg_split_axis(d, rank + d - ax);
    dg_cell_id_t t = (((dg_cell_id_t)1) << ay);
    while ((dp != 0) && (t < tlim))
      { if (dp & 1) 
          { if ((id & t) != 0) 
              { id &= ~t; dp = dp + 2; }
            else
              { id |= t; }
          }
        dp >>= 1; t <<= d;
      }
    return id;
  }
  
dg_cell_id_t dg_cell_shift(dg_dim_t d, dg_cell_id_t id, dg_grid_pos_t dp[])
  { dg_rank_t rank = dg_cell_rank(id);
    /* Limit for bit values that can be added to {id}: */
    dg_cell_id_t tlim = (((dg_cell_id_t)1) << rank);
    int ax;
    for (ax = 0; ax < d; ax++) 
      { dg_grid_pos_t dpax = dp[ax];
        /* Add {dpax} to {id}, starting at bit {ax} with stride {d}: */
        dg_axis_t ay = dg_split_axis(d, rank + d - ax);
        uint64_t t = (((dg_cell_id_t)1) << ay);
        while ((dpax != 0) && (t < tlim))
          { if (dpax & 1) 
              { if ((id & t) != 0) 
                  { id &= ~t; dpax += 2; }
                else
                  { id |= t; }
              }
            dpax >>= 1; t <<= d;
          }
      }
    return id;
  }

dg_cell_id_t dg_cell_id_add(dg_dim_t d, dg_cell_id_t ida, dg_cell_id_t idb)
  { dg_rank_t ra = dg_cell_rank(ida);
    dg_rank_t rb = dg_cell_rank(idb);
    affirm(ra == rb, "unequal ranks");
    /* Limit for bit values that can be added to {ida}: */
    dg_cell_id_t tlim = (((dg_cell_id_t)1) << ra);
    /* Adds the interleaved cell positions. */
    /* Variable {carry} contains the {d} carry bits, rotated. */
    uint64_t t = 1;
    int carry = 0, cd = (1 << d);
    /* Discard leading bit of {idb}: */
    idb -= tlim;
    while (((idb > 0) || (carry != 0)) && (t < tlim))
      { if (idb & carry & 1) 
          { carry |= cd; }
        else if ((idb | carry) & 1)
          { if ((ida & t) != 0) 
              { ida &= ~t; carry |= cd; }
            else
              { ida |= t; }
          }
        idb >>= 1; carry >>= 1; t <<= 1;
      }
    return ida;
  }

void dg_cell_box_root_relative(dg_dim_t d, dg_cell_id_t id, interval_t B[])
  { dg_grid_pos_t pos[dg_dim_MAX];
    dg_cell_position(d, id, pos);
    
    dg_grid_size_t gsz[dg_dim_MAX];
    dg_grid_size(d, dg_cell_rank(id), gsz);
    
    int i;
    for (i = 0; i < d; i++)
      { double scale = 1.0/((double)(gsz[i]));
        LO(B[i]) = scale*((double)pos[i]);
        HI(B[i]) = scale*((double)(pos[i] + 1));
      }
  }

/* Sizes of canonical cells for various dimensions: */
static const double dg_root_sizes_1[1] = 
  { 1.0 };
  
static const double dg_root_sizes_2[2] = 
  { 1.0, 0.7071067811865475 };

static const double dg_root_sizes_3[3] = 
  { 1.0, 0.7937005259840997, 0.6299605249474366 };

static const double dg_root_sizes_4[4] = 
  { 1.0, 0.8408964152537145, 0.7071067811865475, 0.5946035575013605 };

static const double *dg_root_sizes[dg_dim_MAX] = 
  { dg_root_sizes_1,
    dg_root_sizes_2,
    dg_root_sizes_3,
    dg_root_sizes_4
  };
  /* The pointer {dg_root_sizes[d]} is the address of a vector which gives the
    extent of the cannonical cell of dimension {d} along each axis. */

const double *dg_cell_root_canonical_sizes(dg_dim_t d)
  { if (d == 0) { return NULL; }
    demand(d <= dg_dim_MAX, "invalid dimension");
    return dg_root_sizes[d-1];
  }

void dg_cell_box_canonical(dg_dim_t d, dg_cell_id_t id, interval_t B[])
  { dg_cell_box_root_relative(d, id, B); 
    int i;
    const double *sz = dg_root_sizes[d-1];
    for (i = 0; i < d; i++)
      { double szi = sz[i];
        LO(B[i]) *= szi;
        HI(B[i]) *= szi;
      }
  }
  
/* PRINTOUT */

void dg_cell_print(FILE *wr, dg_dim_t d, dg_cell_id_t id)
  { 
    int i;
    dg_grid_pos_t idp[dg_dim_MAX];
    dg_rank_t idr = dg_cell_rank(id);
    dg_cell_position(d, id, idp);
    fprintf(wr, "%d:(", (int)idr);
    for (i = 0; i < d; i++) 
      { fprintf(wr, "%s%d", (i == 0 ? "" : ","), (int)idp[i]); }
    fprintf(wr, ")");
  }
  
/* VECTORS OF CELL IDENTIFIERS */

vec_typeimpl(dg_cell_id_vec_t,dg_cell_id_vec,dg_cell_id_t);

