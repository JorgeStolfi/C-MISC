#ifndef csc_pack_H
#define csc_pack_H
/* Procedures for packing plates and blocks into stock sheets. */
/* Last edited on 2020-01-01 02:14:05 by jstolfi */

/* Copyright © 2020 by the State University of Campinas (UNICAMP). */
/* See the copyright, authorship, and warranty notice at end of file. */

#define csc_pack_INFO \
  "AUTHORS\n" \
  "  Created 2020-01-01 By Jorge Stolfi, IC-UNICAMP." \
  "\n" \
  "HISTORY\n" \
  "  2019-12-01 Created, from the packing logic from {pack_pieces.gawk} (J.Stolfi).\n" \
  "  2019-12-01 Added {scrapMargin} and {cutWidth} options. (J.Stolfi)." \
  "  2019-12-31 Added individual stock sheet description. (J.Stolfi)." \
  "  2019-12-31 Separated the grouping logic from the packing logic. (J.Stolfi)." \
  "  2020-01-01 Added stock sheet labels. (J.Stolfi)."

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>

#include <sheet_cut.h>

#include <vec.h>
#include <argparser.h>
#include <r2.h>
#include <bool.h>

#include <csc_stock_sheet.h>

typedef struct csc_pack_state_t
  { csc_stock_sheet_t sheet; /* Properties of the sheet where the pieces are to be placed. */
    sheet_cut_node_t *blk;   /* The plates already grouped and placed. */
    r2_t max_size;    /* The the bounding box of all plates in {blk} must fit here. */
    double cut_wd;    /* Assumed width of cuts to leave between plates. */
    double cur_py;    /* Baseline of current row of children. */
    double cur_px;    /* Next available X position on the current row of children. */
    double cur_ht;    /* Topline of current row of children. */
  } csc_pack_state_t;
  /* A node of this type describes the state of one packing task or sub-task.
    It escribes a block node {.blk} to which new /items/ 
    (plates or blocks) are to be added as children,
    and the area still available for them.
    
    If {.blk} is {NULL}, the state is /pristine/: no items yet have
    yet been added.  
    
    The current and future items are constrained to fit in a rectangular
    /canvas/ with low corner at {(0,0)} and size {.max_size}. The free
    area of the canvas is the top part of that rectangle, above a
    step-like line that starts horizontally at {.cur_th} on the left
    edge of the canvas, ends horizontally at {.cur_py} on the right
    edge, and drops vertically between the two at the abscissa
    {.cur_px}.
    
    The step boundary of the used area already includes thin strips of
    material that will be destroyed during the physical cutting and
    finishing procedures. Additional cutting safety areas are supposed
    to surround the canvas rectangle, just on the outside of it.
    Therefore, any new items can be placed flush against this boundary
    or against the top and side edges of the canvas.
    
    As each new item is placed, the free area is reduced to account for
    its position and size, plus a safety margin of width {.cut_wd} on any
    side that is not flush against the boundaries of the free area. */
  
typedef csc_pack_state_t* csc_pack_state_ref_t;
vec_typedef(csc_pack_state_vec_t,csc_pack_state_vec,csc_pack_state_ref_t);

csc_pack_state_t* csc_pack_state_new(csc_stock_sheet_t sheet, double scrap_mrg, double cut_wd);
  /* Creates a new packing state record.
    
    Sets the field {.blk} to {NULL} to indicate that there 
    are no pieces yet. The stock sheet properties
    (dimensions, thickness, and material)
    are taken from {sheet}.  The field {.max_size}
    is set to the sheet's dimension minus a margin of {scrap_mrg}
    all around.  The field {.cut_wd} is set as specified.
    
    Sets the used area to empty, namely
    {st.cur_ht = st.cur_py = st.cur_px = 0}.*/

double csc_pack_state_used_area_m2(csc_pack_state_t *st);
  /* Returns the area in {m^2} of the used part of the packing state {st};
    namely, the part of the rectangle from {(0,0)} to {st->max_size}
    that is below the staircase defined by {st->cur_px}, {st->cur_py},
    and  {st->cur_ht}. */

csc_pack_state_vec_t csc_pack_nodes_into_sheets
  ( sheet_cut_node_vec_t *itv,
    csc_stock_sheet_vec_t *shv,
    double scrap_mrg,
    double cut_wd
  );
  /*  
    The procedure receives a list {itv} of zero or more
    input items ({sheet_cut_node_t} records) which may be plates or blocks,
    and packs those input items into the stock sheets
    {shv.e[0..nsh-1]} where {nsh = shv.ne}. 
    
    The procedure leaves a margin of 
    scrap material of with {scrap_mrg} all around the sheet, and a cutting safety margin
    of width at least {cut_wd} between any two plates.
    
    Each item is packed into the first stock sheet of compatible material and thickness
    where it will fit.  Items in each sheet are placed by rows (bottom up, left to right). */
    
bool_t csc_place_node_in_sheet(sheet_cut_node_t *it, csc_pack_state_t *st);
  /* Tries to add the item (plate or block) {it} to the 
    block {st->blk}.  If it succeeds, updates {*st} and returns {TRUE}.
    Otherwise returns FALSE without changing {*st}.
    
    Currently, packs items greedily by rows. Namely, tries to fit the
    item {it} on the "down" part of the step boundary of {*st}, namely at
    {(cur_px,cur_py)}.  If it fits, increments {cur_px} and updates
    {cur_ht}. Otherwise finishes off the current row and places the
    item at {(0,cur_ht)}. */

void csc_pack_state_print(FILE *wr, int32_t ind, csc_pack_state_t *st);
  /* Writes the attributes of state {st} to {wr}, with lines indented by {ind}. */

#endif
