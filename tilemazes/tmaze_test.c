#define PROG_NAME "tmaze_test"
#define PROG_DESC "analysis of tile mazes"
#define PROG_VERS "1.0"

#define tmaze_test_C_COPYRIGHT \
  "Copyright © 2009 by the State University of Campinas (UNICAMP)"

/* Last edited on 2023-02-03 23:33:56 by stolfi */

#define PROG_HELP \
  "  " PROG_NAME " \\\n" \
  "    -size {NX} {NY} \\\n" \
  "    -family {FAMILY} \\\n" \
  "    [ -topology { open | torus} } ] \\\n" \
  "    [ -seed {SEED} ] \\\n" \
  "    [ -print ] \\\n" \
  "    [ -graph ] \\\n" \
  "    [ -comps ] \\\n" \
  "    [ -plot {TILESIZE} ] \\\n" \
  "    [ -cells ] [ -grid ] \\\n" \
  "    [ -style { straight | curved } ] \\\n" \
  "    " argparser_help_info_HELP

#define PROG_INFO \
  "NAME\n" \
  "  " PROG_NAME " - " PROG_DESC "\n" \
  "\n" \
  "SYNOPSIS\n" \
  PROG_HELP "\n" \
  "\n" \
  "DESCRIPTION\n" \
  "  The program creates a random maze of a specified family," \
  " and optionally prints it out and/or plots it and/or writes the" \
  " maze's graph out to disk.\n" \
  "\n" \
  "  All output file names begin with the same {PREFIX}, namely" \
  " \"out/{FAMILY}/{NX}-{NY}-{SEED}-{TOPTAG}\" where {TOPTAG}" \
  " is 'o' if {TOPOLOGY} is \"open\", 't' if it is \"torus\".\n" \
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
  " number generator for the first instance.  The" \
  " default is \"-seed 1\".\n" \
  "\n" \
  "  -print\n" \
  "    This optional parameter requests a text printout of" \
  " the maze to file \"{PREFIX}.prt\".  The tiles" \
  " are printed *from top to bottom*,  each row in a separate line.  Each" \
  " cell will be printed as an ascii character, in a family-dependent" \
  " code,as defined by {tmaze_gen_print}.  If the maze has" \
  " toroidal topology, the edge wrap-around" \
  " is indicated by '|' and '-' characters.\n" \
  "\n" \
  "  -comps\n" \
  "    This optional parameter requests a table giving the" \
  " number of maze components of each size.  The table is written" \
  " to the file \"{PREFIX}.cts\", in the format" \
  " of {tmaze_print_comp_size_counts}.\n" \
  "\n" \
  "  -graph\n" \
  "    This optional parameter requests the program to write the underlying" \
  " graph of the first maze to the file \"{PREFIX}.grf\", in the format" \
  " of {dgraph_print}.\n" \
  "\n" \
  "  -plot {CELLSIZE}\n" \
  "    This optional parameter requests an Encapsulated Postscript plot of the" \
  " first maze, which is written to a file \"{PREFIX}.eps\". Each" \
  " cell will be plotted as a square with {CELLSIZE} by {CELLSIZE}" \
  " millimeters.\n" \
  "\n" \
  "  -cells\n" \
  "    This optional parameter requests that the cell backgrounds be painted.  The" \
  " default is to leave them transparent.\n" \
  "\n" \
  "  -grid\n" \
  "    This optional parameter requests that the cell outlines be stroked.  The" \
  " default is not to stroke the cell outlines.\n" \
  "\n" \
  "  -style { straight | curved }\n" \
  "    This optional parameter requests that the roads be plotted" \
  " with either straight or curved sides.  The default is curved.\n" \
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
  "  2009-11-09 Split the non-stats actions (plot, etc.) to" \
  " {tmaze_test.c},leving the stats in {tmaze_stats.c}.  General cleanup. [J.Stolfi]\n" \
  "\n" \
  "WARRANTY\n" \
  argparser_help_info_NO_WARRANTY "\n" \
  "\n" \
  "RIGHTS\n" \
  "  " tmaze_test_C_COPYRIGHT ".\n" \
  "\n" \
  argparser_help_info_STANDARD_RIGHTS

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
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
  { int32_t nx;                    /* Number of cell columns in maze. */
    int32_t ny;                    /* Number of cell rows in maze. */
    tmaze_family_t family;     /* Maze family to analize.*/
    bool_t torus;              /* TRUE use toroidal topology. */
    int32_t seed;                  /* Seed for random generator. */
    bool_t print;              /* TRUE requests an ascii printout of the maze. */
    bool_t comps;              /* TRUE requests a count of connected components. */
    bool_t graph;              /* TRUE requests a dump of the graph. */
    bool_t plot;               /* TRUE requests a plot of the maze. */
    /* Sub-options (logically) of "-plot": */
    double plot_cell_size;     /* Size of each cell in plot (mm), or 0 if no plot. */
    bool_t plot_cells;         /* TRUE to paint the cell backgrounds. */
    bool_t plot_grid;          /* TRUE to stroke the cell outlines. */
    bool_t plot_curved;        /* TRUE to draw curved roads in plot. */
  } options_t;

options_t *parse_options(int32_t argc, char **argv);
  /* Parses the command line arguments and packs them as an {options_t}. */

int32_t main(int32_t argc,char** argv);

/* IMPLEMENTATIONS */

int32_t main(int32_t argc, char** argv)
  {
    options_t *o = parse_options(argc, argv);

    /* Get the family's name as a string: */
    char *fam_name = tmaze_gen_family_name_lc(o->family);
    
    /* Make sure that the output directory exists: */
    tmaze_make_dir("out", fam_name);

    /* Create a maze with the requested size and topology: */
    tmaze_t M = tmaze_make(o->nx, o->ny, o->torus);

    /* Compose the output filename prefix: */
    char *out_name = NULL;
    asprintf(&out_name, "out/%s/%06d-%06d-%010d-%c", fam_name, M.nx, M.ny, o->seed, "ot"[M.torus]);

    /* Get the maze's size and the max component size: */
    int32_t max_size = tmaze_gen_max_comp_size(&M, o->family);

    /* Fill the maze with random tiles: */
    tmaze_gen_random(&M, o->family, o->seed);

    /* Get the underlying graph {G}: */
    dgraph_t G = tmaze_gen_graph_make(&M, o->family);
    assert(G.rows == G.cols);

    /* Print out the maze, if so requested: */
    if (o->print)
      { char *pfile_name = NULL;
        asprintf(&pfile_name, "%s.grf", out_name);
        FILE *pfile = open_write(pfile_name, TRUE);
        tmaze_gen_print(pfile, &M, o->family);
        fclose(pfile);
        free(pfile_name);
      }

    /* Plot the maze, if so requested: */
    if (o->plot_cell_size > 0)
      { tmaze_gen_plot_maze
          ( ".", out_name, &M, o->family,
            o->plot_cells, o->plot_grid, o->plot_cell_size, o->plot_curved
          );
      }

     /* Dump the graph, if so requested: */
    if (o->graph)
      { char *gfile_name = NULL;
        asprintf(&gfile_name, "%s.grf", out_name);
        FILE *gfile = open_write(gfile_name, TRUE);
        dgraph_write(gfile, &G);
        fclose(gfile);
        free(gfile_name);
      }

    /* Print the component counts, if so requested: */
    if (o->comps)
      { /* Get the connected components as a function of size: */
        tmaze_size_t ms;  /* Max component size in count vector. */
        tmaze_comp_count_t *ct; /* Count of components by size. */
        tmaze_gen_graph_count_components_by_size (&G, &M, o->family, &ms, &ct);
        assert(ms <= max_size);

        /* Write the counts: */
        char *cfile_name = NULL;
        asprintf(&cfile_name, "%s.cts", out_name);
        FILE *cfile = open_write(cfile_name, TRUE);
        tmaze_gen_print_comp_size_counts(cfile, &M, o->family, ms, ct);
        free(ct);
        fclose(cfile);
        free(cfile_name);
      }

    /* Free storage (good manners): */
    free(M.tile);
    free(out_name);
    dgraph_trim(&G, 0);
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

options_t *parse_options(int32_t argc, char **argv)
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

    o->print = argparser_keyword_present(pp, "-print");

    o->comps = argparser_keyword_present(pp, "-comps");

    o->graph = argparser_keyword_present(pp, "-graph");

    if (argparser_keyword_present(pp, "-plot"))
      { o->plot = TRUE;
        o->plot_cell_size = argparser_get_next_double(pp, 0, 100.0);
      }
    else
      { o->plot = FALSE;
        o->plot_cell_size = 0;
      }

    o->plot_cells = argparser_keyword_present(pp, "-cells");

    o->plot_grid = argparser_keyword_present(pp, "-grid");

    if (argparser_keyword_present(pp, "-style"))
      { if (argparser_keyword_present_next(pp, "straight"))
          { o->plot_curved = FALSE; }
        else if (argparser_keyword_present_next(pp, "curved"))
          { o->plot_curved = TRUE; }
        else
          { argparser_error(pp, "unrecognized \"-style\" option"); }
      }
    else
      { o->plot_curved = TRUE; }

    /* Parse the positional parameters: */

    argparser_skip_parsed(pp);

    /* Check for spurious arguments: */

    argparser_finish(pp);

    return o;
  }

