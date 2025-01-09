#define PROG_NAME "gdr_test_sim"
#define PROG_VERS "2.0"

/* Tests {gdr_sim.h} fundtions. */

/* Last edited on 2023-06-01 07:31:49 by stolfi */

#define gdr_test_sim_C_COPYRIGHT \
  "Duh?"

#define PROG_HELP \
  "  " PROG_NAME " \\\n" \
  "    -nTests {NTESTS} \\\n" \
  "    -nYears {NYEARS} \\\n" \
  "    -iniSize {INI_SIZE} \\\n" \
  "    -finSize {FIN_SIZE} \\\n" \
  "    -demoParms " gdr_demo_options_HELP ""

#define PROG_INFO \
  "SYNOPSIS\n" \
  "" PROG_HELP "\n" \
  "\n" \
  "DESCRIPTION\n" \
  "  Tests the functions {gdr_demo_parms_parse} and {gdr_demo_compute_tables}.\n" \
  "\n" \
  "OPTIONS\n" \
  "  -nTests {NTESTS} \n" \
  "    This mandatory argument specifies how many times to run the test.\n" \
  "\n" \
  "  -nYears {NYEARS} \n" \
  "    This mandatory argument specifies the number of years in the interval" \
  " of interest. Those years are numbered {0..NYEARS-1}.\n" \
  "\n" \
  "  -iniSize {INI_SIZE} \n" \
  "  -finSize {FIN_SIZE} \n" \
  "    These mandatory arguments specify the target size" \
  " of the first (year 0) and last (year {NYEARS-1}) cohorts; that" \
  " is the number of individuals that should be born in those" \
  " years.  The target size of other cohorts is interpolated exponentially" \
  " from those two limits.\n" \
  "\n" \
  "  -demoParms {PARMS} \n" \
  "    This mandatory argument specifies the key demographic parameters" \
  " for the simulation. The" \
  " syntax of {PARMS} is as follows:\n" \
  "\n" \
  "    " gdr_demo_options_INFO "\n" \
  "\n" \
  "Duh9?"

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <math.h>

#include <bool.h>
#include <jsmath.h>
#include <rn.h>
#include <jsrandom.h>
#include <argparser.h>
#include <vec.h> 

#include <gdr_demo.h>
#include <gdr_sim.h>

typedef struct gdr_tsi_options_t
  { 
    int32_t nTests;           /* Number of trials. */
    gdr_demo_parms_t *demoParms;  /* Demographic parameters from the command line. */
    int32_t nYears;           /* Number of years to simulate. */
    int32_t iniSize;          /* Target size of cohort 0. */
    int32_t finSize;          /* Target size of cohort {ny-1}. */
  } gdr_tsi_options_t;
 
/* PROTOTYPES */

gdr_tsi_options_t *gdr_tsi_options_parse(int32_t argc, char **argv);
  /* Parses the command line options */

void gdr_tsi_test_one_run
  ( gdr_demo_parms_t *dmp, 
    int32_t nt,
    int32_t ny,
    int32_t *nbt,
    double_vec_t *cFreq
  );
  /* Calls {gdr_sim_run} {nt} times after creating a state with the
    given {dmp}, {ny}, and {nbt} parameters. Records the frequencies
    {cFreq} of the child count {c}, over all runs. These vectors are
    expanded and trimmed as needed. */
  
int32_t main(int32_t argc, char **argv);

/* IMPLEMENTATIONS */

int32_t main(int32_t argc, char **argv)
  {
    bool_t debug = TRUE;
    if (debug) { fprintf(stderr, "> %s\n", __FUNCTION__); }
    srand(417); srandom(417);
    
    gdr_tsi_options_t *o = gdr_tsi_options_parse(argc, argv);
    
    /* Get parameters for sex {s}: */
    gdr_demo_parms_t *dmp = o->demoParms; /* Goal parameters. */
    gdr_demo_show_parms("specified parameters", -1, dmp);

    /* Compute the target cohort sizes table: */
    int32_t ny = o->nYears;
    int32_t *nbt = gdr_sim_compute_target_year_cohort_sizes(ny, o->iniSize, o->finSize);
    
    /* Create the simulation state: */
    double_vec_t cProb = gdr_demo_compute_child_count_distr(o->demoParms);
    gdr_sim_state_t *st = gdr_sim_state_new
      (dmp->fMin, dmp->fMax, &(cProb), dmp->cPrec, ny, nbt);

    /* Observed numbers of individuals with each child count: */
    int64_vec_t cCount = int64_vec_new(0); int32_t cMax_obs = -1;
    
    /* Test simulation, tallying child count frequencies: */
    int32_t ybrMin = 2*dmp->fMax + 1;
    int32_t ybrMax = ny-1-dmp->fMax;
    for (uint32_t t = 0;  t < o->nTests; t++)
      { gdr_sim_run(st);
        gdr_sim_accumulate_child_count_histogram(st, ybrMin, ybrMax, &(cCount), &(cMax_obs));
        if (t == 0) { gdr_sim_write_cohort_sizes("out/gdr-tsi", "A", 0, ny, st->nbr); }
      }
    double_vec_t cFreq = gdr_demo_freqs_from_counts(cMax_obs, cCount.e);
    free(cCount.e);
    
    gdr_demo_show_distr("child count", -1, st->cProb, &(cFreq));

    gdr_sim_state_free(st);
    gdr_demo_parms_free(dmp);
    free(nbt);
    free(cFreq.e);

    return 0;
  }
  
gdr_tsi_options_t *gdr_tsi_options_parse(int32_t argc, char **argv)
  { 
    argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, PROG_HELP);
    argparser_set_info(pp, PROG_INFO);
    argparser_process_help_info_options(pp);
    
    gdr_tsi_options_t *o = notnull(malloc(sizeof(gdr_tsi_options_t)), "no mem");

    argparser_get_keyword(pp, "-nTests");
    o->nTests = (int32_t)argparser_get_next_int(pp, 1, 100);
    
    argparser_get_keyword(pp, "-nYears");
    o->nYears = (int32_t)argparser_get_next_int(pp, 1, 20000);
    
    argparser_get_keyword(pp, "-demoParms");
    o->demoParms = gdr_demo_options_parse(pp);

    argparser_get_keyword(pp, "-iniSize");
    o->iniSize = (int32_t)argparser_get_next_int(pp, 1, 1000000);
    
    argparser_get_keyword(pp, "-finSize");
    o->finSize = (int32_t)argparser_get_next_int(pp, 1, 1000000);
    
    argparser_skip_parsed(pp);
    argparser_finish(pp);
    return o;
  }
