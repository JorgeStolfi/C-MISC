#ifndef csc_group_H
#define csc_group_H

/* Prodecures for grouping plates and blocks into larger blocks. */
/* Last edited on 2020-01-01 02:35:55 by jstolfi */
/* Copyright © 2020 by the State University of Campinas (UNICAMP). */
/* See the copyright, authorship, and warranty notice at end of file. */

#define csc_group_INFO \
  "AUTHORS\n" \
  "  Created 2020-01-01 By Jorge Stolfi, IC-UNICAMP." \
  "\n" \
  "HISTORY\n" \
  "  2019-12-06 Created inside {choose_sheet_cuts.c} (J.Stolfi).\n" \
  "  2020-01-01 Separated {csc_group.{h,c}} from {choose_sheet_cuts.c}. (J.Stolfi)." \
  "  2019-12-31 Separated the grouping logic from the packing logic. (J.Stolfi)."

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>

#include <sheet_cut.h>

#include <vec.h>
#include <r2.h>
#include <bool.h>

sheet_cut_node_vec_t csc_group_read_items
  ( FILE *rd,
    double cut_wd,
    sheet_cut_node_vec_t *plv,
    int32_t level
  );
  /* Reads from the file {rd} a description of
    how to group, orient, and pack plates, in the format 
    described in {csc_group_file_INFO}.
    
    More specifically, the procedure reads from the file zero or more
    input items ({sheet_cut_node_t} records) which may be plates or blocks.
    The procedure returns a vector of the (top-level) items in the file
    
    When the procedure reads a plate tag, it locates the 
    corresponding {sheet_cut_node_t} in the vector {plv},
    and uses that record instead.  That entry is then
    set to {NULL} (so that each plate can be used only once).
    
    When the procedure reads a '[', it calls itself recursively to parse
    the grouping data up to the matching ']'. That call returns a list
    of {sheet_cut_node_t} records (plates or blocks). These nodes are
    made into children of a block item, stacked horizontally, with
    {cut_wd} of space between them.
    
    The {level} is the recursion depth, used only for debugging purposes. */

#define csc_group_file_INFO \
  "The command codes are either '[', ']', or '@'.\n" \
  "\n" \
  "  The commands '[' and ']' delimit" \
  " one or more plate-like items that should be packed together" \
  " horizontally, with proper spacing to account for non-zero cut" \
  " width.  These '[...]' groups can be nested.  Each" \
  " group of items will be treated as a single plate for" \
  " packing purposes, whose dimensions are the common bounding" \
  " box of all its items.  \n" \
  "\n" \
  "  The command '@' indicates that the following plate or" \
  " group should be implicitly transposed (that is, should have" \
  " its X and Y dimensions swapped) before being added to the group.\n" \
  "\n" \
  "  Any text from a '#' character to the end of the" \
  " same line is considered a comment and ignored.  Blank" \
  " lines are ignored too.\n"

#endif
