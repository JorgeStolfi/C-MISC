#ifndef csc_stock_sheet_H
#define csc_stock_sheet_H

/* Tools for packing plates and blocks into stock sheets. */
/* Last edited on 2020-01-01 02:44:01 by jstolfi */
/* Copyright © 2020 by the State University of Campinas (UNICAMP). */
/* See the copyright, authorship, and warranty notice at end of file. */

#define choose_sheet_cuts_pack_INFO \
  "AUTHORS\n" \
  "  Created 2020-01-01 By Jorge Stolfi, IC-UNICAMP." \
  "\n" \
  "HISTORY\n" \
  "  2019-12-31 Created inside {choose_sheet_cuts.c} (J.Stolfi).\n" \
  "  2020-01-01 Added stock sheet labels. (J.Stolfi)."

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>

#include <vec.h>
#include <argparser.h>
#include <r2.h>
#include <bool.h>

typedef struct csc_stock_sheet_t 
  { char *tag;     /* A tag identifiying the sheet. */
    r2_t size;     /* Size of stock sheet. */
    char *mat;     /* Material of stock sheets. */
    double_t thk;  /* Thickness of stock sheet. */
  } csc_stock_sheet_t;
  /* Attributes of a sheet of stock material. */
  
vec_typedef(csc_stock_sheet_vec_t,csc_stock_sheet_vec,csc_stock_sheet_t);

csc_stock_sheet_vec_t csc_stock_sheet_parse_options(argparser_t *pp);
  /* Parses the "-stockSheets" command line arguments, returning a vector of the sheet specs.
    The format and semantics is defined by {csc_stock_sheet_option_HELP} and
    {csc_stock_sheet_option_INFO} below. */
  
#define csc_stock_sheet_option_HELP \
  "-stockSheets {TAG} {N} {MAT} {THK} {SZX} {SZY}"
  
#define csc_stock_sheet_option_INFO \
  "  " csc_stock_sheet_option_HELP "\n"\
  "    This command line argument specifies that there are {N} stock sheets"\
  " of material {MAT} and thickness {THK} available.  This argument"\
  " can be repeated multiple times, each instance appending"\
  " more sheets to the list.  The seets will be"\
  " identified by the tags \"{TAG}.1\", \"{TAG}.2\", etc."

#endif
