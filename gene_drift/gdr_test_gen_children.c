#define PROG_NAME "gdr_test_gen_children"
#define PROG_VERS "2.0"

/* Tests the {gdr_demo_throw_children} function of {gdr_demo.h} tables. */

/* Last edited on 2023-06-01 07:32:58 by stolfi */

#define gdr_test_gen_children_C_COPYRIGHT \
  "Duh?"

#define PROG_HELP \
  "  " PROG_NAME " \\\n" \
  "    -nTests {NTESTS} \\\n" \
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
  "  -demoParms {PARMS} \n" \
  "    This mandatory argument specifies the key demographics parameters" \
  " for the simulations. The" \
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

typedef struct gdr_tgc_options_t
  { 
    int32_t nTests;               /* Number of trials. */
    gdr_demo_parms_t *demoParms;  /* Demographic parameters from the command line. */
  } gdr_tgc_options_t;
 
/* PROTOTYPES */

gdr_tgc_options_t *gdr_tgc_options_parse(int32_t argc, char **argv);
  /* Parses the command line options */

void gdr_tgc_test_throw_children
  ( gdr_demo_parms_t *dmp, 
    double_vec_t *cProb, 
    int32_t nt,
    double_vec_t *cFreq
  );
  /* Calls {gdr_demo_throw_children} {nt} times with the given {dmp} and
    {cProb}, recording the frequencies {cFreq} of the child count {c}.
    The vector is expanded and trimmed as needed. */

int32_t main(int32_t argc, char **argv);

/* IMPLEMENTATIONS */

int32_t main(int32_t argc, char **argv)
  {
    bool_t debug = TRUE;
    if (debug) { fprintf(stderr, "> %s\n", __FUNCTION__); }
    
    gdr_tgc_options_t *o = gdr_tgc_options_parse(argc, argv);
    srand(417); srandom(417);
    
    /* Get parameters for sex {s}: */
    gdr_demo_parms_t *dmp = o->demoParms; /* Goal parameters. */
    gdr_demo_show_parms("specified parameters", -1, dmp);

    /* Compute the goal child count probabilities for sex {s}: */
    double_vec_t cProb = gdr_demo_compute_child_count_distr(dmp);
    demand(cProb.ne == dmp->cMax + 1, "inconsistent {cMax,cProb}");

    /* Test child generation, tallying child count frequencies: */
    double_vec_t cFreq = double_vec_new(10);
    gdr_tgc_test_throw_children(dmp, &cProb, o->nTests, &cFreq);
    gdr_demo_show_distr("child count", -1, &(cProb), &(cFreq));

    gdr_demo_parms_free(dmp);
    free(cProb.e);
    free(cFreq.e);

    return 0;
  }

void gdr_tgc_test_throw_children
  ( gdr_demo_parms_t *dmp, 
    double_vec_t *cProb, 
    int32_t nt,
    double_vec_t *cFreq
  )
  {
    int32_t cMax = dmp->cMax; /* Max num of children. */
    double_vec_expand(cFreq, cMax);
    
    int32_t cAge[cMax+1]; /*Children of an individual. */
    for (int32_t it = 0; it < nt; it++)
      { int32_t cNum;
        gdr_demo_throw_children(cProb, dmp->fMin, dmp->fMax, &cNum, cAge);
        
        assert((cNum >= 0) && (cNum <= cMax));
        cFreq->e[cNum] += 1.0;
      }
      
    int32_t cFreqPrec = dmp->cPrec + 2;
    double fr_unit = pow(0.1, cFreqPrec); /* Rounding unit. */

    rn_scale(cMax+1, 1.0/nt, cFreq->e, cFreq->e);
    gdr_demo_distr_round_off(fr_unit, cFreq);
    int32_t cMaxF = cFreq->ne - 1;
    while ((cMaxF > 0) && (cFreq->e[cMaxF] == 0)) { cMaxF--; }
    double_vec_trim(cFreq, cMaxF+1);
  }
  
gdr_tgc_options_t *gdr_tgc_options_parse(int32_t argc, char **argv)
  { 
    argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, PROG_HELP);
    argparser_set_info(pp, PROG_INFO);
    argparser_process_help_info_options(pp);
    
    gdr_tgc_options_t *o = notnull(malloc(sizeof(gdr_tgc_options_t)), "no mem");

    argparser_get_keyword(pp, "-nTests");
    o->nTests = (int32_t)argparser_get_next_int(pp, 1, 1000000);
    
    argparser_get_keyword(pp, "-demoParms");
    o->demoParms = gdr_demo_options_parse(pp);

    argparser_skip_parsed(pp);
    argparser_finish(pp);
    return o;
  }
