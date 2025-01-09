#define PROG_NAME "gdr_test_draw"
#define PROG_VERS "2.0"

/* Draws the evolution of asexual traits. */

/* Last edited on 2023-06-25 19:33:05 by stolfi */

#define gdr_test_draw_C_COPYRIGHT \
  "Duh?"

#define PROG_HELP \
  " " PROG_NAME " \\\n" \
  "    -yStart {yStart} \\\n" \
  "    -yStop {yStop} \\\n" \
  "    -iniSize {iniSize} \\\n" \
  "    -finSize {finSize} \\\n" \
  "    -outPrefix {outPrefix} \\\n" \
  "    -tag {tag} \\\n" \
  "    -demoParms 0 {demoParms[0]} \n" \
  "    -demoParms 1 {demoParms[1]} \n" \
  "  where each {demoParms[s]} is \\\n" \
  "    " gdr_demo_options_HELP ""

#define PROG_INFO \
  "SYNOPSIS\n" \
  "  " PROG_HELP "\n" \
  "\n" \
  "DESCRIPTION\n" \
  "  Simulates the history of a asexually reproduced traits in" \
  " a bisexual population over several years, where each individual" \
  " has a variable number of children with a specific probability" \
  " distribution depending on age and sex.  Produces a graphical" \
  " representation of the evolution.\n" \
  "\n" \
  "  The main parameters are the initial population size, the" \
  " demographic parameters of each sex, and the overall" \
  " population growth rate.\n" \
  "\n" \
  "  The simulation spans a total of {ny} years, numbered" \
  " from {yStart} (most remote in the past) to {yStop} (most recent).  Individuals" \
  " are identified by unique /indices/ starting from 0, in chrono order of birth.\n" \
  "\n" \
  "OUTPUTS\n" \
  "  The main output file is\n" \
  "\n" \
  "    \"{outPrefix}-{tag}-{s}-{ID}-plot.eps\" for each sex {s} (0 or 1)\n" \
  "      " gdr_plot_state_file_INFO "\n" \
  "\n" \
  "OPTIONS\n" \
  "  -yStart {yStart}\n" \
  "  -yStop {yStop}\n" \
  "    These mandatory arguments specify the initial and final year of the simulation. The values may be negative.\n" \
  "\n" \
  "  -iniSize {iniSize}\n" \
  "  -finSize {finSize}\n" \
  "    These mandatory arguments specify the maximum cohort" \
  " size (number of individuals born) on years {yStart} and" \
  " {yStop}, respectively. The max cohort size at intermediate" \
  " years will be interpolated exponentially between" \
  " these two values.\n" \
  "\n" \
  "  -outPrefix {outPrefix}\n" \
  "    This mandatory argument specifies a common prefix for all output file names.\n" \
  "\n" \
  "  -tag {tag}\n" \
  "    This mandatory argument specifies a tag that can be used to" \
  " distinguish separate runs of the program.  It must be a single" \
  " upercase letter. It is also used in the names of macros" \
  " written to the TeX parameters file \"{outPrefix}-{tag}-parms.tex\".\n" \
  "\n" \
  "  -demoParms {s} {demoParms[s]} \n" \
  "    This mandatory argument specifies the key parameters" \
  " of the demographic parameters of each sex. It must appear" \
  " exactly twice, with {s} 0 (male) and 1 (female).  The" \
  " syntax of the parameters is\n" \
  "\n" \
  gdr_demo_options_INFO "\n" \
  "\n" \
  "Duh9?"

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#include <bool.h>
#include <cmp.h>
#include <jsmath.h>
#include <jsfile.h>
#include <jsrandom.h>
#include <argparser.h>
#include <in.h>
#include <rn.h>
#include <vec.h> 
#include <epswr.h> 

#include <gdr_demo.h>
#include <gdr_sim.h>
#include <gdr_lineage.h>
#include <gdr_draw.h>

typedef gdr_demo_parms_t gdr_chops_t; /* For short. */

typedef struct gdr_test_draw_options_t
  { int32_t iniSize;            /* Number of individuals in initial population. */
    int32_t finSize;            /* Desired number of individuals in final population. */
    int32_t yStart;             /* First year covered by simulation. */
    int32_t yStop;              /* Last year covered by simulation. */
    gdr_chops_t *demoParms[2];  /* Demographic parameters from command line. */
    char *outPrefix;            /* Prefix for output files. */
    char *tag;                  /* Tag identifying the run of the program. */
  } gdr_test_draw_options_t;
  /* The command line options.  */
   
/* DATA STRUCTURES AND NOTATIONS

  See also the comments in {gdr_demo.h}, {gdr_sim.h}, {gdr_lineage.h}.

  SEXES 
  
  The sex of an individual is represented by an integer either 0
  (sometimes called 'male') and 1 ('female'). 
  
  YEARS
  
  Externally the years are numbered {o.yStart..o.yStop}, but internally
  they are numbered {0..ny-1}. So internal year {y} means external year
  {yStart+y} in files. Unless said otherwise, in the following comments
  the internal numbering is used. */

#define gdr_cohort_size_MAX 20
#define gdr_num_years_MAX 500
  /* Limits for a reasonable plot file. */

void gdr_test_draw_check_limits(int32_t ny, int32_t nbt[]);
  /* Checks whether {ny,nbt[0..ny-1]} are not too big. */
  
void gdr_test_draw_plot_state
  ( char *outPrefix, 
    char *tag,
    int32_t s,
    int32_t id, 
    gdr_sim_state_t *st,
    int32_t y0,
    int32_t y1, 
    int32_t yRef
  );
  /* Calls {gdr_draw_state_plot(prefix,st,y0,y1,yRef)} where {prefix}
    is {gdr_test_draw_file_prefix(outPrefix, tag, s, id)} */

char *gdr_test_draw_file_prefix(char *outPrefix, char *tag, int32_t s, int32_t id);
  /* Returns the string "{outPrefix}-{tag}-{s}-{ID}.eps",
    where {ID} is {id} formatted to two digits, zero-filled.  The "-{s}" and/or "-{ID}" parts are
    omitted if the parameters are negative. */ 

void gdr_test_draw_dump_counts(char *name, int32_t s, int32_t ny, int32_t count[]);
  /* Writes the table {count[0..ny-1]} to {stderr} with a title that includes {name} and sex {s}. */
  
gdr_test_draw_options_t *gdr_parse_options(int32_t argc, char **argv);
  /* Parses the command line arguments. */
  
int32_t main(int32_t argc, char **argv);

/* IMPLEMENTATIONS */

int32_t main(int32_t argc, char **argv)
  {
    bool_t debug = FALSE;
    if (debug) { fprintf(stderr, "> %s\n", __FUNCTION__); }
    
    gdr_test_draw_options_t *o = gdr_parse_options(argc, argv);
    srand(417); srandom(417);
    
    int32_t ny = o->yStop - o->yStart + 1;  /* Number of years in history. */

    /* Compute the desired size {nbt[0..ny-1]} of each year: */
    int32_t *nbt = gdr_sim_compute_target_year_cohort_sizes(ny, o->iniSize, o->finSize);
    
    /* Sanity limits check: */
    gdr_test_draw_check_limits(ny, nbt);

    /* Population and lineage counts for a single run: */
    int32_t pop[ny];     /* {pop[y]} is the population for year {y}, for some sex {s} and run {r}. */ 
    int32_t nsl[ny];     /* {nsl[y]} is the surv lineage count for year {y}, for some sex {s} and run {r}. */ 
    
    for (uint32_t s = 0;  s < 2; s++) 
      { 
        fprintf(stderr, "  ### sex %d ##################################################\n", s);

        gdr_demo_parms_t *dmp = o->demoParms[s];
        double_vec_t cProb;  /* Nominal child probabilities for sex {s}. */

        cProb = gdr_demo_compute_child_count_distr(dmp);
        gdr_sim_state_t *st = gdr_sim_state_new
          (dmp->fMin, dmp->fMax, &(cProb), dmp->cPrec, ny, nbt);

        gdr_demo_show_distr("target child count", s, st->cProb, NULL);
    
        /* Run the simulation once for {st.ny} years: */
        if (debug) { fprintf(stderr, "    running the simulation\n"); }
        gdr_sim_run(st);
 
        gdr_sim_compute_populations(st, pop);
        if (debug) { gdr_test_draw_dump_counts("pop", s, ny, pop); }
        
        gdr_lineage_compute_surviving_counts(st, nsl);
        if (debug) { gdr_test_draw_dump_counts("nsl", s, ny, nsl); }
       
        /* Draw state showing lineages at several years: */
        int32_t yRef0 = ny-1;
        gdr_test_draw_plot_state(o->outPrefix, o->tag, s, 1, st, 0, ny-1, yRef0);
        /* int32_t yRef1 = yRef0 - 100;                                 */
        /* gdr_test_draw_plot_state(o->outPrefix, o->tag, s, 2, st, 0, ny-1, yRef1); */
        /* int32_t yRef2 = yRef1 - 5;                                   */
        /* gdr_test_draw_plot_state(o->outPrefix, o->tag, s, 3, st, 0, ny-1, yRef2); */
        
        /* Get actual child count count frequencies {cFreq} for sex {s}: */
        int32_t ybrMin = 2*st->fMax;            /* Try to avoid the initial transient. */  
        int32_t ybrMax = ny - 1 - st->fMax - 1; /* Avoid final period with unrestricted births. */
        int64_vec_t cCount = gdr_sim_compute_child_count_histogram(st, ybrMin, ybrMax);
        double_vec_t cFreq = gdr_demo_freqs_from_counts(cCount.ne, cCount.e);
        gdr_demo_show_distr("child count", s, st->cProb, &(cFreq));
        
        free(cCount.e);
        free(cProb.e);
        free(cFreq.e);
        gdr_sim_state_free(st);
        
        fprintf(stderr, "  ############################################################\n");
      }

    for (uint32_t s = 0;  s < 2; s++) 
      { 
        
      }
    if (debug) { fprintf(stderr, "< %s\n", __FUNCTION__); }
    return 0;
  }

char *gdr_test_draw_file_prefix(char *outPrefix, char *tag, int32_t s, int32_t id)
  { 
    char *sx = NULL; if (s >= 0) { char *sx = jsprintf("-%d", s); } else { sx = ""; }
    char *idx = NULL; if (id >= 0) { char *idx = jsprintf("-%02d", id); } else { idx = ""; }
    char *prefix = jsprintf("%s-%s%s%s", outPrefix, tag, sx, idx);
    if (s >= 0) { free(sx); }
    if (id >= 0) { free(idx); }
    return prefix;
  }

void gdr_test_draw_plot_state
  ( char *outPrefix, 
    char *tag,
    int32_t s,
    int32_t id, 
    gdr_sim_state_t *st,
    int32_t y0,
    int32_t y1,
    int32_t yRef
  )
  {
    char *prefix = gdr_test_draw_file_prefix(outPrefix, tag, s, id);
    gdr_draw_state_plot(prefix, st, y0, y1, yRef);
    free(prefix);
  }

void gdr_test_draw_dump_counts(char *name, int32_t s, int32_t ny, int32_t count[])
  {
    fprintf(stderr, "  --- %s sex %d ---\n", name, s);
    for (uint32_t y = 0;  y < ny; y++)
      { fprintf(stderr, "     %4d ", y);
        fprintf(stderr, " %8d", count[y]);
        fprintf(stderr, "\n");
      }
    fprintf(stderr, "  --------------\n");
  }

void gdr_test_draw_check_limits(int32_t ny, int32_t nbt[]) 
  {
    bool_t debug = FALSE;
    if (debug) { fprintf(stderr, "  > %s\n", __FUNCTION__); }
    
    int32_t nbtMax = INT32_MAX/ny/100;
    if (debug) { fprintf(stderr, "    ny = %d nbtMax = %d\n", ny, nbtMax); }
    demand(ny <= gdr_lineage_num_years_MAX, "too many years");
    for (uint32_t y = 0;  y < ny; y++) 
      { if (debug) { fprintf(stderr, "    nbt[%d] = %d\n", y, nbt[y]); }
        demand(nbt[y] <= nbtMax, "too many births");
      }

    if (debug) { fprintf(stderr, "  < %s\n", __FUNCTION__); }
  }

gdr_test_draw_options_t *gdr_parse_options(int32_t argc, char **argv)
  { 
    argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, PROG_HELP);
    argparser_set_info(pp, PROG_INFO);
    argparser_process_help_info_options(pp);
    
    gdr_test_draw_options_t *o = notnull(malloc(sizeof(gdr_test_draw_options_t)), "no mem");

    argparser_get_keyword(pp, "-iniSize");
    o->iniSize = (int32_t)argparser_get_next_int(pp, 1, gdr_cohort_size_MAX);

    argparser_get_keyword(pp, "-finSize");
    o->finSize = (int32_t)argparser_get_next_int(pp, 1, gdr_cohort_size_MAX);

    argparser_get_keyword(pp, "-yStart");
    o->yStart = (int32_t)argparser_get_next_int(pp, -100000, +3000);

    argparser_get_keyword(pp, "-yStop");
    int32_t yStopMax = o->yStart + gdr_num_years_MAX - 1;
    o->yStop = (int32_t)argparser_get_next_int(pp, o->yStart+2, yStopMax);

    /* The keyword "-demoParms" must appear exactly twice, once for each sex: */
    for (uint32_t s = 0;  s < 2; s++) { o->demoParms[s] = NULL; }
    for (uint32_t ks = 0;  ks < 2; ks++)
      { argparser_get_keyword(pp, "-demoParms");
        int32_t s = (int32_t)argparser_get_next_int(pp, 0, 1);
        if (o->demoParms[s] != NULL) 
          { argparser_error(pp, "multiple \"-demoParms\" for the same sex"); }
        o->demoParms[s] = gdr_demo_options_parse(pp);
      }
    if (o->demoParms[0]->cPrec != o->demoParms[1]->cPrec)
      { argparser_error(pp, "the {cPrec} parameter must be the same for both sexes"); }

    argparser_get_keyword(pp, "-outPrefix");
    o->outPrefix = argparser_get_next(pp);

    argparser_get_keyword(pp, "-tag");
    o->tag = argparser_get_next(pp);
    if (strlen(o->tag) != 1)
      { argparser_error(pp, "\"-tag\" must be just one letter"); }
    if ((strcmp(o->tag, "A") < 0) || (strcmp(o->tag, "Z") > 0)) 
      { argparser_error(pp, "\"-tag\" must be in [A-Z]"); }

    argparser_skip_parsed(pp);
    argparser_finish(pp);
    return o;
  }
