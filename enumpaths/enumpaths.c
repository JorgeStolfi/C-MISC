#define PROG_NAME "enumpaths"
#define PROG_DESC "enumerates and counts paths in a directed graph"
#define PROG_VERS "" ?? ""

#define PROG_C_COPYRIGHT \
  "Copyright © 2009 by the State University of Campinas (UNICAMP)"

/* Last edited on 2024-12-21 11:56:34 by stolfi */

#define PROG_HELP \
  "  " PROG_NAME " \\\n" \
  "    -graph {GRAPH_FILE} \\\n" \
  "    {OTHER_OPTIONS} \\\n" \
  "    " argparser_help_info_HELP ""

#define PROG_INFO \
  "NAME\n" \
  "  " PROG_NAME " - " PROG_DESC "\n" \
  "\n" \
  "SYNOPSIS\n" \
  PROG_HELP "\n" \
  "\n" \
  "DESCRIPTION\n" \
  "  Reads a directed graph. Enumerates and counts the paths of various lengths.\n" \
  "\n" \
  "OPTIONS\n" \
  "  -graph {GRAPH_FILE}\n" \
  "    Specifies the name of the file containing the graph.\n" \
  "\n" \
  "DOCUMENTATION OPTIONS\n" \
  argparser_help_info_HELP_INFO "\n" \
  "\n" \
  "SEE ALSO\n" \
  "  emacs(1), a great editor.  Don't leave /home without it.\n" \
  "\n" \
  "AUTHOR\n" \
  "  Created 2009-11-01 by Jorge Stolfi, IC-UNICAMP.\n" \
  "\n" \
  "MODIFICATION HISTORY\n" \
  "  2009-11-01 Created [JS].\n" \
  "\n" \
  "WARRANTY\n" \
  argparser_help_info_NO_WARRANTY "\n" \
  "\n" \
  "RIGHTS\n" \
  "  " PROG_C_COPYRIGHT ".\n" \
  "\n" \
  argparser_help_info_STANDARD_RIGHTS

#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <math.h>

#include <bool.h>
#include <argparser.h>
#include <dgraph.h>
#include <jsfile.h>
#include <jsrandom.h>

/* COMMAND-LINE OPTIONS */



typedef struct dgraph_t
  { int nv;
    bool_t *adj;
  } dgraph_t;
  /* A directed graph with vertices {0..nv-1}.  The {.adj}
    is an {nv} by {nv} adjacency matrix, linearized by rows;
    so {adj[nv*u + v]} is true iff there is a directed edge
    from {u} to {v}. */

typedef int vertex_t;

typedef struct options_t
  { char* graph;      /* Name of graph file. */
    vertex_t root;    /* Root vertex, or {dgraph_vertex_NONE} if any. */
  } options_t;

/* INTERNAL PROTOTYPES */

int main(int argc,char** argv);
options_t *enp_parse_options(int argc, char **argv);
  /* Parses the command line arguments and packs them as an {options_t}. */

dgraph_t *enp_graph_read(char *fname);
  /* Reads a graph from {fname}. The format is defined by {spmat_io.h} */

dgraph_t *enp_graph_throw(int nv, int ne, bool_t val);
  /* Generates a random loopfree directed graph with {nv} vertices
    and {ne} edges (if {val=TRUE}) or {ne} non-edges (if {val=FALSE}).
    uses {int32_abrandom} as a source of randomness. */

/* IMPLEMENTATIONS */

int main(int argc, char** argv)
  {
    options_t *o = enp_parse_options(argc, argv);
    graph_t *G = dgraph_read(o->graph);
    word_counter_t *K = count_words(G);
    print_word_counts(stdout, K);
    return 0;
  }
  
dgraph_t *enp_graph_read(char *fname)
  {
    FILE *rd = open_read(fname, TRUE);
    dgraph_t *G = dgraph_read(rd);
    return G;
  }

dgraph_t *enp_graph_throw(int nv, int ne, bool_t val)
  {
    max_edges = nv*(nv-1);
    demand(ne <= max_edges, "too many edges");
    if (ne > max_edges/2) { return enp_graph_throw(nv, max_edges - ne, (! val)); }
    int nv2 = nv*nv;
    bool_t *G = notnull(malloc(nv2*sizeof(bool_t)), "no mem");  /* Graph,expanded. */
    int uv;
    for (uv = 0; uv < nv2; uv++) { G[uv] = (! val); }
    int me = 0; /* Edges currently in {G}. */
    while (me < ne)
      { assert(nv >= 2);
        /* Pick a random vertex {u}: */
        int u = int32_abrandom(0,nv-1);
        /* Pick a random vertex {v} distinct from {u}: */
        int v = int32_abrandom(0,nv-2); if (v >= u) { v++; }
        int uv = nv*u + v;
        if (G[uv] != val)
          { /* Add edge and update count: */
            G[uv] = val;
            me++;
          }
      }
    /* Compactify the graph: */
    
    
        
      
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
    
    /* Parse keyword parameters: */
    
    if (argparser_keyword_present(pp, "-graph"))
      { o->graph = argparser_get_next(pp); }
    else
      { o->graph = "NONE"; }

    /* Parse positional arguments: */
    argparser_skip_parsed(pp);
    
    if (argparser_next(pp) != NULL)
      { o->infile = argparser_get_next(pp); }
    else
      { o->infile = "-"; }

    if (argparser_next(pp) != NULL)
      { o->outfile = argparser_get_next(pp); }
    else
      { o->outfile = "-"; }

    /* Check for spurious arguments: */
    argparser_finish(pp);
    
    return o;
  }
