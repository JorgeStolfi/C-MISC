/* See {csc_pack.h}. */
/* Last edited on 2020-11-06 22:27:00 by jstolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <values.h>
#include <string.h>
#include <assert.h>

#include <sheet_cut.h>

#include <vec.h>
#include <argparser.h>
#include <r2.h>
#include <bool.h>
#include <fget.h>
#include <jsfile.h>
#include <affirm.h>

#include <csc_stock_sheet.h>
#include <csc_pack.h>

vec_typeimpl(csc_pack_state_vec_t, csc_pack_state_vec, csc_pack_state_ref_t);

double csc_pack_state_used_area_m2(csc_pack_state_t *st) 
  { 
    double area1 = sheet_cut_area_m2( (r2_t){{ st->cur_px, st->cur_ht }} );
    double area2 = sheet_cut_area_m2( (r2_t){{ st->max_size.c[0] - st->cur_px, st->cur_py }} );
    return area1 + area2;
  }

csc_pack_state_vec_t csc_pack_nodes_into_sheets
  ( sheet_cut_node_vec_t *itv,
    csc_stock_sheet_vec_t *shv,
    double scrap_mrg,
    double cut_wd
  )
  {
    bool_t debug = FALSE;
    
    /* Create a vector of packing states, one for eeach sheet: */
    int32_t nsh = shv->ne; /* Number of sheets available, and of their packing states. */
    csc_pack_state_vec_t stv = csc_pack_state_vec_new(nsh);
    for (uint32_t ish = 0;  ish < nsh; ish++)
      { csc_stock_sheet_t shi = shv->e[ish];
        csc_pack_state_t *sti = csc_pack_state_new(shi, scrap_mrg, cut_wd);
        fprintf(stderr, "starting sheet %d\n", ish);
        csc_pack_state_print(stderr, 2, sti);
        stv.e[ish] = sti;
      }
    
    /* Pack each item into the first compatible stock sheet where it fits: */
    int32_t nit = itv->ne; /* Number of items to pack. */
    for (uint32_t kit = 0;  kit < nit; kit++)
      { sheet_cut_node_t *itk = itv->e[kit];
        assert(itk != NULL);
        assert((itk->pos.c[0] == 0.0) &&(itk->pos.c[1] == 0.0));
        if (debug) 
          { fprintf(stderr, "trying to place %s", itk->tag); 
            fprintf(stderr, " size =  (%.1f,%.1f)\n", itk->size.c[0], itk->size.c[1]); 
          }
        /* Try to fit in each stock sheet in turn: */
        bool_t ok = FALSE; /* Set to true when item is successfully placed. */
        int32_t ish; /* Index of sheet where placement succeded. */
        for (ish = 0; ish < nsh; ish++)
          { csc_pack_state_t *sti = stv.e[ish];
            csc_stock_sheet_t shi = shv->e[ish];
            if (debug) 
              { fprintf(stderr, "  trying to place on sheet %d\n", ish);
                csc_pack_state_print(stderr, 2, sti);
              }
            /* Check is node material matches sheet's: */
            if ((strcmp(itk->mat, shi.mat) != 0) || (itk->thk != shi.thk)) { continue; };
            /* Try stuffing the node in the sheet's packing state: */
            ok = csc_place_node_in_sheet(itk, sti);
            if (ok) { break; }
          }
        if (! ok) 
          { fprintf(stderr, "** node %s cannot be placed in any stock sheet\n", itk->tag); }
        else
          { assert(ish < nsh);
            csc_stock_sheet_t shi = shv->e[ish];
            double xk = itk->pos.c[0], yk = itk->pos.c[1];
            fprintf(stderr, "placing %s on sheet %d=%s at (%.1f,%.1f)\n", itk->tag, ish, shi.tag, xk, yk);
          }
      }
      
    return stv;
  } 

csc_pack_state_t* csc_pack_state_new(csc_stock_sheet_t sheet, double scrap_mrg, double cut_wd)
  { 
    csc_pack_state_t* st = notnull(malloc(sizeof(csc_pack_state_t)), "no mem");
    st->blk = NULL;
    st->sheet = sheet;
    st->cur_py = 0.0;
    st->cur_ht = 0.0;
    st->cur_px = 0.0; 

    r2_t scrap_size = (r2_t){{ 2*scrap_mrg, 2*scrap_mrg }};
    r2_sub(&(sheet.size), &scrap_size, &(st->max_size));
    st->cut_wd = cut_wd;
    
    return st;
  }
  
bool_t csc_place_node_in_sheet(sheet_cut_node_t *it, csc_pack_state_t *st)
  { 
    if ((st->blk != NULL) && (! sheet_cut_same_mat_thk(st->blk, it)))
      { /* Cannot pack because of material or thickness change: */
        return FALSE;
      }
    
    /* Item size: */
    double it_dx = it->size.c[0];
    double it_dy = it->size.c[1];
    assert((it->pos.c[0] == 0.0) &&(it->pos.c[1] == 0.0));

    /* Canvas size: */
    double sh_dx = st->max_size.c[0];
    double sh_dy = st->max_size.c[1];
    
   /* Boundary of the used area in the st's canvas: */
    double sh_ht = st->cur_ht;
    double sh_py = st->cur_py;
    double sh_px = st->cur_px;
    
    /* Cutting safety margin width: */
    double cut_wd = st->cut_wd;
    
    /* Note that the step boundary may overshoot the canvas boundary by {cur_wd}. */
    assert((0.0 <= sh_py) && (sh_py <= sh_ht) && (sh_ht <= sh_dy + cut_wd));
    assert((0.0 <= sh_px) && (sh_px <= sh_dx + cut_wd));
    
    assert((sh_px == 0.0) || (sh_py < sh_ht));  /* Staircase if cur row is not empty. */
    
    /* Set {(it_x,it_y)} to the chosen placement, or to {(NAN,NAN)} if failure: */
    double it_x = NAN, it_y = NAN;

    if ((sh_dx - sh_px >= it_dx) && (sh_dy - sh_py >= it_dy))
      { /* Place {it} at the end of the current row: */
        it_x = sh_px; it_y = sh_py;
        st->cur_px = sh_px + it_dx + cut_wd;
        st->cur_ht = fmax(sh_ht, sh_py + it_dy + cut_wd);
      }
    else if ((sh_dx >= it_dx) && (sh_dy - sh_ht >= it_dy))
      { /* Start a new row of items, place {it} in it: */
        it_x = 0.0; it_y = sh_ht;
        st->cur_py = sh_ht;
        st->cur_ht = sh_ht + it_dy + cut_wd;
        st->cur_px = it_dx + cut_wd;
      }
    
    if (! isnan(it_x))
      { /* Succeeded: */
        if ((it_x + it_dx > sh_dx + 0.1) || (it_y + it_dy > sh_dy + 0.1))
          { /* Desperation placement: */
            fprintf(stderr, "!! warning: item extends outside canvas boundaries.\n");
          }

        /* Place the piece: */
        it->pos = (r2_t){{ it_x, it_y }};

        /* Add piece to block: */
        st->blk = sheet_cut_add_child(st->blk, it);
        return TRUE;
      }
    else
      { /* Failed to place: */
        return FALSE;
      }
  }
    
void csc_pack_state_print(FILE *wr, int32_t ind, csc_pack_state_t *st)
  { fprintf(wr, "%*ssheet %s", ind, "", st->sheet.tag);
    fprintf(wr, " mat = %s thk = %.1f", st->sheet.mat, st->sheet.thk);
    fprintf(wr, " cut_wd = %.1f\n", st->cut_wd);
    fprintf(wr, "%*stotal size = (%.1f,%.1f)\n", ind, "", st->sheet.size.c[0], st->sheet.size.c[1]);
    fprintf(wr, "%*susable size = (%.1f,%.1f)\n", ind, "", st->max_size.c[0], st->max_size.c[1]);
    fprintf(wr, "%*sused boundary = (%.1f,%.1f,%.1f)\n", ind, "", st->cur_ht, st->cur_px, st->cur_py);
    if (st->blk != NULL) { sheet_cut_print_node(wr, ind, st->blk, FALSE); }
  }
