#define PROG_NAME "tmaze_figs"
#define PROG_DESC "create figures for the paper on tile mazes"
#define PROG_VERS "1.0"

#define tmaze_figs_C_COPYRIGHT \
  "Copyright © 2009 by the State University of Campinas (UNICAMP)"

/* Last edited on 2009-11-08 23:02:02 by stolfi */

#define PROG_HELP \
  "  " PROG_NAME " \\\n" \
  "    " argparser_help_info_HELP

#define PROG_INFO \
  "NAME\n" \
  "  " PROG_NAME " - " PROG_DESC "\n" \
  "\n" \
  "SYNOPSIS\n" \
  PROG_HELP "\n" \
  "\n" \
  "DESCRIPTION\n" \
  "  The program writes figures called \"out/{DIR}/{NAME}.eps\".\n" \
  "\n" \
  "OPTIONS\n" \
  "\n" \
  "DOCUMENTATION OPTIONS\n" \
  argparser_help_info_HELP_INFO "\n" \
  "\n" \
  "SEE ALSO\n" \
  "  mirrormaze(1).\n" \
  "\n" \
  "AUTHOR\n" \
  "  Created 2009-01-22 by Jorge Stolfi, IC-UNICAMP.\n" \
  "\n" \
  "MODIFICATION HISTORY\n" \
  "  Option bla bla added by J. Stolfi, IC-UNICAMP.\n" \
  "\n" \
  "WARRANTY\n" \
  argparser_help_info_NO_WARRANTY "\n" \
  "\n" \
  "RIGHTS\n" \
  "  " tmaze_figs_C_COPYRIGHT ".\n" \
  "\n" \
  argparser_help_info_STANDARD_RIGHTS

#define _GNU_SOURCE
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <errno.h>

#include <bool.h>
#include <jsfile.h>
#include <jsstring.h>
#include <argparser.h>

#include <tmaze.h>
#include <tmaze_ba.h>
#include <tmaze_bb.h>
#include <tmaze_figs_ba.h>
#include <tmaze_figs_bb.h>

/* TYPES AND CONSTANTS */

#define tmaze_tileset_FOUR (15)
  /* The set of the four tile types {0,1,2,3}. */

/* COMMAND-LINE OPTIONS */

typedef struct options_t
  { 
  } options_t;

/* PROTOTYPES */

void tmaze_figs_plt_empty_board(char *dir, char *name, int nx, int ny, bool_t torus, int val);
void tmaze_figs_plt_tile_graph(char *dir, char *name, int val, bool_t directed);
void tmaze_figs_plt_tile_roads(char *dir, char *name, int val);
  /* Each of these procedures writes an EPS file called "{dir}/{name}.eps"
    with a figure for the tile mazes paper (generic mazes). */

options_t *parse_options(int argc, char **argv);
  /* Parses the command line arguments and packs them as an {options_t}. */

void show_options(options_t *o);
  /* Prints the options {o} to {stderr}. */

int main(int argc,char** argv);

/* IMPLEMENTATIONS */

int main(int argc, char** argv)
  {
    options_t *o = parse_options(argc, argv);
    show_options(o);
    
    /* Make sure that directories exist: */
    mkdir("out", 0777);
    mkdir("out/basic", 0777);
    mkdir("out/brasilia", 0777);
    mkdir("out/blipblop", 0777);
    
    /* Illustration of an empty board: */
    tmaze_figs_plt_empty_board("basic","board", 7,4,TRUE, 0);
    
    /* Illustration of an empty board: */
    tmaze_figs_plt_empty_board("basic","board-ports-3", 7,4,TRUE,3);
    
    /* Illustration of a tile as directed and undirected graph: */
    tmaze_figs_plt_tile_graph("basic","tile-dgraph", 3,TRUE);
    tmaze_figs_plt_tile_graph("basic","tile-ugraph", 3,FALSE);
    
    /* Illustration of a tile with undirected graph as a road network: */
    tmaze_figs_plt_tile_roads("basic","tile-roads", 3);
    
    /* Brasilia Airport maze figures: */
    
    /* Illustrations of the maze tiles: */
    tmaze_figs_plt_ba_tileset("brasilia","tile-E", (1 << tmaze_ba_tile_E));
    tmaze_figs_plt_ba_tileset("brasilia","tile-W", (1 << tmaze_ba_tile_W));
    tmaze_figs_plt_ba_tileset("brasilia","tile-N", (1 << tmaze_ba_tile_N));
    tmaze_figs_plt_ba_tileset("brasilia","tile-S", (1 << tmaze_ba_tile_S));
    
    /* Illustrations of a sample maze: */
    tmaze_figs_plt_ba_sample("brasilia","maze-sma",  8,  4, 4561, FALSE);
    tmaze_figs_plt_ba_sample("brasilia","maze-med", 12,  7, 4651, FALSE);
    tmaze_figs_plt_ba_sample("brasilia","maze-big", 35, 41, 4643,  TRUE);
    
    /* Illustration of a conjunctive pattern: */
    tmaze_figs_plt_ba_conjuntive_pat("brasilia","cj-pattern");
    
    /* Illustration of patterns for connected components: */
    tmaze_figs_plt_ba_component_pattern("brasilia","pattern", 0, 0);
    tmaze_figs_plt_ba_component_pattern("brasilia","pattern", 0, 1);
    
    tmaze_figs_plt_ba_component_pattern("brasilia","pattern", 1, 0);
    
    tmaze_figs_plt_ba_component_pattern("brasilia","pattern", 2, 0);
    tmaze_figs_plt_ba_component_pattern("brasilia","pattern", 2, 1);
    tmaze_figs_plt_ba_component_pattern("brasilia","pattern", 2, 2);
    tmaze_figs_plt_ba_component_pattern("brasilia","pattern", 2, 3);
    
    tmaze_figs_plt_ba_component_pattern("brasilia","pattern", 3, 3);
    tmaze_figs_plt_ba_component_pattern("brasilia","pattern", 3, 7);
    
    tmaze_figs_plt_ba_component_pattern("brasilia","pattern", 4, 0);
    tmaze_figs_plt_ba_component_pattern("brasilia","pattern", 4, 1);
    
    /* Blip Blop maze figures: */
    
    /* Illustrations of the maze tiles: */
    tmaze_figs_plt_bb_tileset("blipblop","tile-blip", (1 << tmaze_bb_tile_BLIP));
    tmaze_figs_plt_bb_tileset("blipblop","tile-blop", (1 << tmaze_bb_tile_BLOP));
    
    /* Illustrations of a sample maze: */
    tmaze_figs_plt_bb_sample("blipblop","maze-sma",  6,  5, 4165, FALSE);
    tmaze_figs_plt_bb_sample("blipblop","maze-med", 12,  7, 4615, FALSE);
    tmaze_figs_plt_bb_sample("blipblop","maze-big", 35, 41, 4635,  TRUE);
    
    /* Illustration of a conjunctive pattern: */
    tmaze_figs_plt_bb_conjuntive_pat("blipblop","cj-pattern");

    /* Illustration of patterns for connected components: */
    tmaze_figs_plt_bb_component_pattern("blipblop","pattern", 4, 0);
    
    tmaze_figs_plt_bb_component_pattern("blipblop","pattern", 8, 0);
    tmaze_figs_plt_bb_component_pattern("blipblop","pattern", 8, 1);
    
    tmaze_figs_plt_bb_component_pattern("blipblop","pattern",12, 0);
    tmaze_figs_plt_bb_component_pattern("blipblop","pattern",12, 1);
    tmaze_figs_plt_bb_component_pattern("blipblop","pattern",12, 2);
    
    tmaze_figs_plt_bb_component_pattern("blipblop","pattern",16, 0);
    tmaze_figs_plt_bb_component_pattern("blipblop","pattern",16, 1);
    tmaze_figs_plt_bb_component_pattern("blipblop","pattern",16, 2);
    tmaze_figs_plt_bb_component_pattern("blipblop","pattern",16, 3);
    tmaze_figs_plt_bb_component_pattern("blipblop","pattern",16, 4);
    tmaze_figs_plt_bb_component_pattern("blipblop","pattern",16, 5);
    
    tmaze_figs_plt_bb_component_pattern("blipblop","pattern",20, 0);
    
    return 0;
  }
  
void tmaze_figs_plt_empty_board(char *dir, char *name, int nx, int ny, bool_t torus, int val)
  {
    tmaze_pattern_t M = tmaze_make_pattern(nx, ny, torus, FALSE, val);
    tmaze_fill_pattern(&M, tmaze_tileset_NONE);
    tmaze_plot_tileset_proc_t *proc = NULL;
    tmaze_plot_pattern(txtcat("out/",dir), name, &M, TRUE, TRUE, proc, (val > 0), 8.0, 2.0);
  }
  
void tmaze_figs_plt_tile_graph(char *dir, char *name, int val, bool_t directed)
  {
    fprintf(stderr, "** %s not implemented yet\n", __FUNCTION__);
  }
  
void tmaze_figs_plt_tile_roads(char *dir, char *name, int val)
  {
    fprintf(stderr, "** %s not implemented yet\n", __FUNCTION__);
  }

options_t *parse_options(int argc, char **argv)
  {
    /* Initialize argument parser: */
    argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, PROG_NAME " version " PROG_VERS ", usage:\n" PROG_HELP);
    argparser_set_info(pp, PROG_INFO);
    argparser_process_help_info_options(pp);
    
    /* Allocate the command line argument record: */
    options_t *o = (options_t *)malloc(sizeof(options_t)); 
    
    /* Parse the keyword parameters: */

    /* Parse the positional parameters: */
    
    argparser_skip_parsed(pp);
    
    /* Check for spurious arguments: */
    
    argparser_finish(pp);
    
    return o;
  }

void show_options(options_t *o)
  {
    fprintf(stderr, "-- command line arguments --------------------------------\n");
    fprintf(stderr, "----------------------------------------------------------\n");
  }

