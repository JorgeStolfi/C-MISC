#define PROG_NAME "tmaze_stats"
#define PROG_DESC "analysis of tile mazes"
#define PROG_VERS "1.0"

#define tmaze_stats_C_COPYRIGHT \
  "Copyright © 2009 by the State University of Campinas (UNICAMP)"

/* Last edited on 2023-02-03 23:46:39 by stolfi */

#define PROG_HELP \
  "  " PROG_NAME " \\\n" \
  "    -size {NX} {NY} \\\n" \
  "    -family {FAMILY} \\\n" \
  "    [ -topology { open | torus} } ] \\\n" \
  "    [ -seed {SEED} ] \\\n" \
  "    [ -trials {TRIALS} ] \\\n" \
  "    " argparser_help_info_HELP

#define PROG_INFO \
  "NAME\n" \
  "  " PROG_NAME " - " PROG_DESC "\n" \
  "\n" \
  "SYNOPSIS\n" \
  PROG_HELP "\n" \
  "\n" \
  "DESCRIPTION\n" \
  "  The program computes and writes out the count of" \
  " components of a given size in a maze with specified" \
  " dimensions, family, and and topology, averaged over" \
  " several random instances of such mazes.  Also writes" \
  " out the observed deviation of each count and its" \
  " theoretically expected value.\n" \
  "\n" \
  "  All output file names begin with the same {PREFIX}, namely" \
  " \"out/{FAMILY}/{NX}-{NY}-{SEED}-{TOPTAG}-{TRIALS}\" where {TOPTAG}" \
  " is 'o' if {TOPOLOGY} is \"open\", 't' if it is \"torus\".\n" \
  "\n" \
  "  In particular, the table of average component counts per size" \
  " are written to the files \"{PREFIX}-large.act\" and" \
  " \"{PREFIX}-small.act\"  The first set considers only" \
  " the largest components in each trial, and the second" \
  " set considers all cmponents except the largest ones.\n" \
  "\n" \
  "OPTIONS\n" \
  "  -size {NX} {NY}\n" \
  "    This mandatory parameter specifies the number of cell columns {NX} and" \
  " cell rows {NY} in the maze, respectively.\n" \
  "\n" \
  "  -family {FAMILY}\n" \
  "    This mandatory parameter specifies the maze family.  The valid" \
  " options currently are \"brasilia\" and \"blipblop\".\n" \
  "\n" \
  "  -trials {TRIALS}\n" \
  "    This optional parameter specifies the number of" \
  " instances that should be generated and analyzed.  The" \
  " default is 2 instances.\n" \
  "\n" \
  "  -topology { open | torus }\n" \
  "    This optional parameter specifies what happens at the edges of the" \
  " cell array.  The \"open\" option specifies that" \
  " paths cannot cross those edges, so that cells along them" \
  " have no neighbors.  The \"torus\"option specifies that" \
  " opposite edges of the maze are implicltly" \
  " connected with toroidal topology, so that" \
  " every cell has four neighbors.  The default is \"-topology open\".\n" \
  "\n" \
  "  -seed {SEED}\n" \
  "    This optional parameter specifies the (integer) seed for the random" \
  " number generator for the first instance.  For subsequent" \
  " instances, the seed is incremented by a fixed constant.  The" \
  " default is \"-seed 1\".\n" \
  "\n" \
  "DOCUMENTATION OPTIONS\n" \
  argparser_help_info_HELP_INFO "\n" \
  "\n" \
  "SEE ALSO\n" \
  "  mirrormaze(1).\n" \
  "\n" \
  "AUTHOR\n" \
  "  Created 2009-01-15 by Jorge Stolfi, IC-UNICAMP.\n" \
  "\n" \
  "MODIFICATION HISTORY\n" \
  "  2009-01-15 created {tmaze_stats.c} for Brasilia mazes [J.Stolfi]\n" \
  "\n" \
  "  2009-02-04 renamed {tmaze_stats.c} to {tmaze_ba_stats.c} [J.Stolfi]\n" \
  "\n" \
  "  2009-02-04 created {tmaze_bb_stats.c} for Blip-Blop mazes [J.Stolfi]\n" \
  "\n" \
  "  2009-11-08 Merged again the two into {tmaze_stats.c}. Added" \
  " the \"-trials\" option and changed printout to averages. [J.Stolfi]\n" \
  "\n" \
  "  2009-11-09 Moved the main statistics-gathering loop to {tmaze_gen.h}. Moved" \
  " the non-stats actions (plot, etc.) to {tmaze_test.c}.  Provided for output" \
  " to file (instead of stdout), standard deviations, and separate tallying" \
  " of small and large components. [J.Stolfi]\n" \
  "\n" \
  "WARRANTY\n" \
  argparser_help_info_NO_WARRANTY "\n" \
  "\n" \
  "RIGHTS\n" \
  "  " tmaze_stats_C_COPYRIGHT ".\n" \
  "\n" \
  argparser_help_info_STANDARD_RIGHTS

#define _GNU_SOURCE
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>

#include <bool.h>
#include <jsfile.h>
#include <argparser.h>
#include <dgraph.h>

#include <tmaze.h>
#include <tmaze_gen.h>

/* COMMAND-LINE OPTIONS */

typedef struct options_t
  { int nx;                 /* Number of cell columns in maze. */
    int ny;                 /* Number of cell rows in maze. */
    tmaze_family_t family;  /* Maze family to analize.*/
    bool_t torus;           /* TRUE use toroidal topology. */
    int trials;             /* Number of instances to generate and analyze. */
    int seed;               /* Seed for random generator. */
  } options_t;

options_t *parse_options(int argc, char **argv);
  /* Parses the command line arguments and packs them as an {options_t}. */

int main(int argc,char** argv);

/* IMPLEMENTATIONS */

int main(int argc, char** argv)
  {
    options_t *o = parse_options(argc, argv);
    
    /* Get the family's name as a string: */
    char *fam_name = tmaze_gen_family_name_lc(o->family);

    /* Make sure that the output directory exists: */
    tmaze_make_dir("out", fam_name);

    /* Create a maze with the requested size and topology: */
    tmaze_t M = tmaze_make(o->nx, o->ny, o->torus);
    
    /* Compose the output filename prefix: */
    char *out_name = jsprintf("out/%s/%06d-%06d-%010d-%c-%06d", fam_name, M.nx, M.ny, o->seed, "ot"[M.torus], o->trials);
     
    /* Get the total size of the maze {tot_size} and the max component size {max_size}: */
    int max_size = tmaze_gen_max_comp_size(&M, o->family);
    
    /* Get the statistics: */
    tmaze_size_t ms; /* Max component size in stats vectors. */
    double *ct_avg[2] = { NULL, NULL }; /* Average over all trials of comp count by class and size. */
    double *ct_dev[2] = { NULL, NULL }; /* Deviation over all trials of comp count by class and size. */
    tmaze_gen_comp_stats_gather
      ( &M, o->family, o->trials, o->seed, 
        &ms, &(ct_avg[0]), &(ct_dev[0]), &(ct_avg[1]), &(ct_dev[1])
      );
    assert(ms <= max_size);

    /* Write out the distribution: */
    int class; /* 0 = small, 1 = large */
    for (class = 0; class <= 1; class++)
      { char *dfile_name = NULL;
        char *dfile_name = jsprintf("%s-%s.act", out_name, (class == 0 ? "small" : "large"));
        FILE *dfile = open_write(dfile_name, TRUE);
        tmaze_gen_print_comp_size_distr(dfile, &M, o->family, ms, ct_avg[class], ct_dev[class]);
        fclose(dfile);
        free(dfile_name);
      }

    /* Free storage (good manners): */
    free(M.tile);
    return 0;
  }

#define MIN_MAZE_SIZE (2)
#define MAX_MAZE_SIZE (65535)
  /* Maze size range in both dimensions. The minimum width and minimum
    height are defined as 2, since the graph of a 1-wide or 1-tall
    maze has parallel edges, which are not well-supported by the
    representation. */

#define MAX_SEED (1073741823L)

#define MAX_TRIALS (1000000)

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
    
    argparser_get_keyword(pp, "-size");
    o->nx = (int32_t)argparser_get_next_int(pp, MIN_MAZE_SIZE, MAX_MAZE_SIZE);
    o->ny = (int32_t)argparser_get_next_int(pp, MIN_MAZE_SIZE, MAX_MAZE_SIZE);
    
    argparser_get_keyword(pp, "-family");
    if (argparser_keyword_present_next(pp, "brasilia")) 
      { o->family = tmaze_family_BRASILIA; }
    else if (argparser_keyword_present_next(pp, "blipblop"))
      { o->family = tmaze_family_BLIPBLOP; }
    else
      { argparser_error(pp, "unknown maze family"); }

    if (argparser_keyword_present(pp, "-topology"))
      { if (argparser_keyword_present_next(pp, "torus")) 
          { o->torus = TRUE; }
        else if (argparser_keyword_present_next(pp, "open"))
          { o->torus = FALSE; }
        else
          { argparser_error(pp, "unknown topology"); }
      }
    else
      { o->torus = TRUE; }

    if (argparser_keyword_present(pp, "-seed"))
      { o->seed = (int32_t)argparser_get_next_int(pp, 0, MAX_SEED); }
    else
      { o->seed = 1; }

    if (argparser_keyword_present(pp, "-trials"))
      { o->trials = (int32_t)argparser_get_next_int(pp, 1, MAX_TRIALS); }
    else
      { o->trials = 2; }

    /* Parse the positional parameters: */
    
    argparser_skip_parsed(pp);
    
    /* Check for spurious arguments: */
    
    argparser_finish(pp);
    
    return o;
  }
