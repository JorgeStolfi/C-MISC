/* See udg_grid.h */
/* Last edited on 2014-05-15 22:44:39 by stolfilocal */

#define _GNU_SOURCE
#include <affirm.h>
#include <bool.h>
#include <ix.h>
#include <stdlib.h>
#include <assert.h>
  
#include <udg_grid.h>

udg_grid_size_t udg_grid_size(udg_rank_t r)
  { return (((udg_cell_id_t)1) << r);  }

udg_grid_pos_t udg_max_grid_pos(udg_rank_t r)
  { return (((udg_cell_id_t)1) << r) - 1; }

sign_t udg_rank_cmp(udg_rank_t *x, udg_rank_t *y)
  { if ((*x) < (*y))
      { return -1; }
    else if ((*x) > (*y))
      { return +1; }
    else
      { return 0; }
  }
  
/* CELL IDENTIFIERS */
  
udg_rank_t udg_cell_rank(udg_cell_id_t id)
  { int r = 0;
    while (id > (((udg_cell_id_t)1) << 16)) { id >>= 16; r+= 16; }
    while (id > (((udg_cell_id_t)1) << 4)) { id >>= 4; r+= 4; }
    while (id > 1) { id >>= 1; r+= 1; }
    return r;
  }

udg_grid_pos_t udg_cell_position(udg_cell_id_t id)
  { udg_cell_id_t t = ((((udg_cell_id_t)1) << udg_cell_rank(id)) >> 1);
    return id - t;
  }
  
udg_cell_id_t udg_cell_from_position(udg_rank_t r, udg_grid_pos_t pos)
  { /* First cell of layer {r}: */
    udg_cell_id_t id0 = (((udg_cell_id_t)1) << r);
    /* Number of cells {sz} in layer {r}: */
    udg_cell_id_t sz = (((udg_cell_id_t)1) << r);
    /* Displacement {t} from {id0} is {pos % sz}: */
    uint64_t t = pos & (sz - 1);
    /* Compute index: */
    return id0 + t;
  }

udg_cell_id_t udg_cell_shift(udg_cell_id_t id, udg_grid_pos_t dp)
  { udg_rank_t r = udg_cell_rank(id);
    /* First cell of layer {r}: */
    udg_cell_id_t id0 = (((udg_cell_id_t)1) << r);
    /* Number of cells {sz} in layer {r}: */
    udg_cell_id_t sz = (((udg_cell_id_t)1) << r);
    /* Displacement {t} from {id0} is {(id - id0 + dp) % sz}: */
    uint64_t t = (id - id0 + dp) & (sz - 1);
    /* Compute index: */
    return id0 + t;
  }

void udg_cell_interval_root_relative(udg_cell_id_t id, interval_t *B)
  { udg_grid_pos_t pos = udg_cell_position(id);
    udg_grid_size_t gsz = udg_grid_size(udg_cell_rank(id));
    double scale = 1.0/((double)gsz);
    LO((*B)) = scale*((double)pos);
    HI((*B)) = scale*((double)(pos + 1));
  }

/* PRINTOUT */

void udg_cell_print(FILE *wr, udg_cell_id_t id)
  { 
    udg_grid_pos_t idp = udg_cell_position(id);
    udg_rank_t idr = udg_cell_rank(id);
    fprintf(wr, "%d:", (int)idr);
    fprintf(wr, "%d", (int)idp);
  }
  
/* VECTORS OF CELL IDENTIFIERS */

vec_typeimpl(udg_cell_id_vec_t,udg_cell_id_vec,udg_cell_id_t);
