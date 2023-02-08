/* See {choose_sheet_cuts.h}. */
/* Last edited on 2020-10-27 15:32:18 by jstolfi */

#define PROG_NAME "choose_sheet_cuts"
#define PROG_DESC "fits rectangular plates on rectangular sheets and draws them"
#define PROG_VERS "1.1"

/* Copyright © 2019 by the State University of Campinas (UNICAMP). */
/* See the copyright, authorship, and warranty notice at end of file. */

#define PROG_HELP \
  "  " PROG_NAME " \\\n" \
  "  { " csc_stock_sheet_option_HELP " }.. \\\n" \
  "  -scrapMargin {SMRG} \\\n" \
  "  -cutWidth {CWD} \\\n" \
  "  -inPrefix {INPREFIX} \\\n" \
  "  -outPrefix {OUTPREFIX} "

#define PROG_INFO \
  "NAME\n" \
  "  " PROG_NAME " - " PROG_DESC "\n" \
  "\n" \
  "SYNOPSIS\n" \
  PROG_HELP "\n" \
  "\n" \
  "DESCRIPTION\n" \
  "  Reads from file \"{INPREFIX}_plates.txt\"a description of a" \
  " set of /plates/: rectangles to be cut out of a list of rectangular stock" \
  " sheets of specified material, and thicknesses.\n" \
  "\n" \
  "  Then reads" \
  " from file \"{INPREFIX}_grouping.txt\"instructions" \
  " on which plates to consider and in which order and orientation; and" \
  " how to group some of those plates into blocks, where each" \
  " block will be handled as a single plate.\n" \
  "\n" \
  "  Then takes from" \
  " the command line the specifications of one or more available" \
  " stock sheets, including their sizes, materials, and" \
  " thicknesses.\n" \
  "\n" \
  "  Packs the plates and blocks into those sheets.\n" \
  "\n" \
  "  Writes one or more" \
  " Encapsulated Postscript (EPS) files with drawings of" \
  " the sheets, with the plates placed onto the sheets.\n" \
  "\n" \
  "  Each plate or block is added to the first stock sheet" \
  " of the proper material and thickness where it fits.  In" \
  " each stock sheet, the plates and blocks are added by" \
  " rows, bottom to top, and left to right within each row.\n" \
  "\n" \
  "  The figures will be actual scale.\n" \
  "\n" \
  "INPUT FILES\n" \
  "\n"\
  "  The file with plate data is \"{INPREFIX}_plates.txt\", and contains plate" \
  " data, one per line.  " csc_PLATE_LINE_INFO "\n" \
  "\n" \
  "  The file with grouping specs is \"{INPREFIX}_grouping.txt\", and contains" \
  " a series of plate tags interspersed with groupling and placing" \
  " commands.  " csc_group_file_INFO "\n" \
  "\n" \
  "OUTPUT FILES\n" \
  "  The EPS files are called \"{OUTPREFIX}_{SHEET}.eps\" where" \
  " {SHEET} is a three-digit sheet number starting from 1.\n" \
  "\n" \
  "  The program also writes to \"{OUTPREFIX}_layout.txt\" a textual" \
  " description of the plates and their layout, one plate per line.\n" \
  "\n" \
  "  " sheet_cut_write_plate_INFO "\n" \
  "\n" \
  "OPTIONS\n" \
  "\n"\
  csc_stock_sheet_option_INFO "\n"\
  "\n"\
  "AUTHORS\n" \
  "  Created 2019-11-28 by Jorge Stolfi, IC-UNICAMP." \
  "\n" \
  "HISTORY\n" \
  "  2019-11-28 Created (J.Stolfi).\n" \
  "  2019-12-01 Included the packing logic from {pack_pieces.gawk}. (J.Stolfi).\n" \
  "  2019-12-01 Added {scrapMargin} and {cutWidth} options. (J.Stolfi)." \
  "  2019-12-06 Added recursive groups and rotation codes. (J.Stolfi)." \
  "  2019-12-07 Spinned off some code as library {libsheetcut}. (J.Stolfi)." \
  "  2019-12-08 Split input into plate info and grouping info. (J.Stolfi)." \
  "  2019-12-31 Added options to specify the available stock sheets. (J.Stolfi)." \
  "  2019-12-31 Separated the grouping logic from the packing logic. (J.Stolfi)." \
  "  2020-01-01 Separated {csc_stock_sheet.{h,c}} and {csc_pack.{h,c}}. (J.Stolfi)." \
  "  2020-01-01 Added stock sheet labels. (J.Stolfi)."

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <values.h>
#include <string.h>
#include <assert.h>

#include <sheet_cut.h>
#include <sheet_cut_plot.h>
#include <sheet_cut_write.h>

#include <vec.h>
#include <epswr.h>
#include <argparser.h>
#include <r2.h>
#include <bool.h>
#include <fget.h>
#include <jsfile.h>
#include <affirm.h>

#include <csc_stock_sheet.h>
#include <csc_group.h>
#include <csc_pack.h>

typedef struct csc_options_t
  { csc_stock_sheet_vec_t stockSheets;  /* Properties of stock sheets available. */
    double scrapMargin; /* Width of materials to be discarded around edges of sheet. */
    double cutWidth;    /* Assumed width of cuts. */
    char *inPrefix;     /* Prefix for input file names. */
    char *outPrefix;    /* Prefix for output file names. */
  } csc_options_t;
  /* Arguments from command line. */

/* INTERNAL PROTOTYPES */

csc_options_t *csc_parse_options(int argc, char **argv);
  /* Parses the command line arguments. */

int main(int argc, char **argv);

sheet_cut_node_vec_t csc_read_all_plates(FILE *rd);
  /* Reads from the file {rd} a sequence of plate data.  Builds
    a plate node for each plate, with {.pos = (0,0)}.  Returns the vector
    of such plates. */

sheet_cut_node_t* csc_read_plate(FILE *rd);
  /* Reads the next line from the plate input file {rd}, 
    in the format {csc_PLATE_LINE_INFO} below.
    Creates a new {sheet_cut_node_t} record for it, and returns its address.
    
    if it runs into end-of-file, returns {NULL}.  */
                
#define csc_PLATE_LINE_INFO \
  "The description of one plate consists of 5 fields, all in a line by themselves:\n" \
  "\n" \
  "    \"{MAT} {TAG} {DX} {DY} {THK}\"\n" \
  "\n" \
  "  where\n" \
  "\n" \
  "    {MAT} is a string with no blanks that describes the sheet's material,\n" \
  "    {TAG} is a string with no blanks that identifies the plate,\n" \
  "    {DX} and {DY} are the dimensions of a plate, and\n" \
  "    {THK} is the sheet thickness.\n" \
  "\n" \
  "  All dimensions and coordinates are floats in arbitrary units (but" \
  " this documentation assumes millimeters).\n" \

FILE* csc_open_file(char *prefix, char *compl, bool_t out);
  /* Opens the file char "{inPrefix}{compl}" for reading
    (if {out=FALSE}) or writing (if {out=TRUE}). */

/* IMPLEMENTATIONS */
 
int main (int argc, char **argv)
  {
    csc_options_t *o = csc_parse_options(argc, argv);
    
    /* Read the plate data without packing: */
    FILE *plate_rd = csc_open_file(o->inPrefix, "_plates.txt", FALSE);
    sheet_cut_node_vec_t plv = csc_read_all_plates(plate_rd);
    fclose(plate_rd);
    fprintf(stderr, "read %d plates\n", plv.ne);
    
    /* Process the  grouping file, get a list of the top-level items: */ 
    FILE *pack_rd = csc_open_file(o->inPrefix, "_grouping.txt", FALSE);
    sheet_cut_node_vec_t shv = csc_group_read_items(pack_rd, o->cutWidth, &plv, 0);
    if (! feof(pack_rd)) 
      { fprintf(stderr, "**error: spurious input, possibly mismatched ']'\n"); exit(1); }
    fclose(pack_rd); 
    
    /* Check for unused plates: */
    int32_t plate_ct = 0; /* Plates that were referenced by the grouping file. */
    for (int32_t k = 0; k < plv.ne; k++) 
      { if (plv.e[k] != NULL)
          { fprintf(stderr, "!! warning: plate %s was not used\n", plv.e[k]->tag); }
        else
          { plate_ct++; }
      }
    fprintf(stderr, "groupings list %d plates\n", plate_ct);
    
    /* Pack the items into sheets: */
    csc_pack_state_vec_t stv = 
      csc_pack_nodes_into_sheets(&(shv), &(o->stockSheets), o->scrapMargin, o->cutWidth);
    
    /* Draw the sheets and write them out: */
    FILE* wr_placed = csc_open_file(o->outPrefix, "_layout.txt", TRUE);
    for (int32_t ist = 0; ist < stv.ne; ist++)
      { csc_pack_state_t *sti = stv.e[ist];
        if (sti->blk != NULL)
          { /* Plot the layout: */
            fprintf(stderr, "starting eps figure for sheet %d...\n", ist);
            epswr_figure_t *eps = sheet_cut_plot_new_figure
              ( o->outPrefix, 
                sti->sheet.tag, sti->sheet.mat, sti->sheet.thk, sti->sheet.size, 
                o->scrapMargin,
                sti->cur_px, sti->cur_py, sti->cur_ht
              );
            sheet_cut_plot_all_nodes(eps, sti->blk, (r2_t){{0,0}});
            fprintf(stderr, "end eps figure.\n");
            sheet_cut_plot_end_figure(eps);
            
            /* Write the plates to {stdout}: */
            int32_t sheet_plate_ct = 0;
            sheet_cut_write_all_plates(wr_placed, ist, &sheet_plate_ct, sti->blk, (r2_t){{0,0}});
            fflush(wr_placed);

            /* Write sheet summary on {stderr}: */
            double area_total = sheet_cut_area_m2(sti->sheet.size);
            double area_usable = sheet_cut_area_m2(sti->max_size);
            double area_used = csc_pack_state_used_area_m2(sti);
            double pct_used = 100*area_used/area_usable/
            fprintf(stderr, "sheet %d (%s %.1f mm):\n", ist, sti->sheet.mat, sti->sheet.thk);
            fprintf(stderr, "  %d plates\n", sheet_plate_ct);
            fprintf(stderr, "  total area = %.2f m^2\n", area_total);
            fprintf(stderr, "  usable area = %.2f m^2\n", area_usable);
            fprintf(stderr, "  used area = %.2f m^2 (%.1f%%)\n", area_used, pct_used);
          }
      }
    fclose(wr_placed);

    fprintf(stderr, "groupings list %d plates on %d sheets\n", plate_ct, stv.ne);
     
    return 0;
  }
  
FILE* csc_open_file(char *prefix, char *compl, bool_t out)
  { 
    char *fname = NULL;
    asprintf(&fname, "%s%s", prefix, compl);
    FILE *ff;
    if (out)
      { ff = open_write(fname, TRUE); }
    else
      { ff = open_read(fname, TRUE); }
    free(fname);
    return ff;
  }
    
sheet_cut_node_vec_t csc_read_all_plates(FILE *rd)
  { 
    /* Vector of plates read: */
    sheet_cut_node_vec_t plv = sheet_cut_node_vec_new(500); 
    int32_t npl = 0; /* Plates read so far are {plv.e[0..npl-1]}. */

    /* Loop on input plates until EOF: */
    while (TRUE)
      { sheet_cut_node_t *pc = csc_read_plate(rd);
        if (pc == NULL) { break; }
        sheet_cut_node_vec_expand(&plv, npl);
        plv.e[npl] = pc; npl++;
      }
    sheet_cut_node_vec_trim(&plv, npl);
    return plv;
  }
  
sheet_cut_node_t* csc_read_plate(FILE *rd)
  { 
    bool_t debug = FALSE;
    fget_skip_formatting_chars(rd);
    if (feof(rd))
      { /* No more input, return: */
        return NULL;
      }
    else
      { /* {SHEET} {MAT} {THK} {PX} {PY} {DX} {DY} {TAG} */
        char *mat = fget_string(rd);
        char *tag = fget_string(rd);
        if (debug) { fprintf(stderr, "  mat = %s tag = %s", mat, tag); }
        r2_t size;
        size.c[0] = fget_double(rd);  
        size.c[1] = fget_double(rd);  
        double thk = fget_double(rd); 
        if (debug) { fprintf(stderr, "  dx = %.1f dy = %.1f thk = %.1f\n", size.c[0], size.c[1], thk); }
        sheet_cut_node_t *pc = sheet_cut_new_plate(mat, thk, size, tag);
        /* Check rest of line: */
        fget_comment_or_eol(rd, '#');
        return pc;
      }
  }
    
csc_options_t *csc_parse_options(int argc, char **argv)
  {
    /* Initialize argument parser: */
    argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, PROG_NAME " version " PROG_VERS ", usage:\n" PROG_HELP);
    argparser_set_info(pp, PROG_INFO);
    argparser_process_help_info_options(pp);

    /* Allocate the command line argument record: */
    csc_options_t *o = notnull(malloc(sizeof(csc_options_t)), "no mem");

    /* Parse keyword parameters: */

    o->stockSheets = csc_stock_sheet_parse_options(pp);

    argparser_get_keyword(pp, "-scrapMargin");
    o->scrapMargin = argparser_get_next_double(pp, 0.0, 100.0);

    argparser_get_keyword(pp, "-cutWidth");
    o->cutWidth = argparser_get_next_double(pp, 0.1, 10.0);

    argparser_get_keyword(pp, "-outPrefix");
    o->outPrefix = argparser_get_next_non_keyword(pp);

    argparser_get_keyword(pp, "-inPrefix");
    o->inPrefix = argparser_get_next_non_keyword(pp);

    /* Check for spurious arguments: */
    argparser_finish(pp);

    return o;
  }

