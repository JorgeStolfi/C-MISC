#define PROG_NAME "gdr_test_demo"
#define PROG_VERS "2.0"

/* Tests the functions of {gdr_demo.h}. */

/* Last edited on 2023-06-01 07:33:56 by stolfi */

#define gdr_test_demo_C_COPYRIGHT \
  "Duh?"

#define PROG_HELP \
  "  " PROG_NAME " \\\n" \
  "    -demoParms " gdr_demo_options_HELP " \\\n" \
  "    -demoParms " gdr_demo_options_HELP ""

#define PROG_INFO \
  "SYNOPSIS\n" \
  "" PROG_HELP "\n" \
  "\n" \
  "DESCRIPTION\n" \
  "  Tests the functions {gdr_demo_parms_parse} and {gdr_demo_compute_tables}.\n" \
  "\n" \
  "OPTIONS\n" \
  "  -demoParms {PARMS} \n" \
  "    This mandatory argument specifies the key demographic parameters" \
  " for the simulation. It must appear" \
  " exactly twice, referring to sexes 0 and 1, respectively.  The" \
  " syntax of each {PARMS} is as follows:\n" \
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
#include <jsrandom.h>
#include <argparser.h>
#include <vec.h> 

#include <gdr_demo.h>

typedef struct gdr_tdm_options_t
  { 
    gdr_demo_parms_t *demoParms[2]; /* Dmographic parameters from the command line */
  } gdr_tdm_options_t;
 
#define gdr_demo_PREC 3

/* PROTOTYPES */

gdr_tdm_options_t *gdr_tdm_options_parse(int32_t argc, char **argv);
  /* Parses the command line options */

gdr_demo_parms_t *gdr_tdm_create_fake_parms(gdr_demo_parms_t *dmp);
  /* Creates fake "actual frequencies" tables by perturbing {dmp} a little. */
  
int32_t main(int32_t argc, char **argv);

/* IMPLEMENTATIONS */

int32_t main(int32_t argc, char **argv)
  {
    bool_t debug = TRUE;
    if (debug) { fprintf(stderr, "> %s\n", __FUNCTION__); }
    
    gdr_tdm_options_t *o = gdr_tdm_options_parse(argc, argv);
    srand(417); srandom(417);
    
    /* Test parenting distribution functions for both sexes: */
    for (uint32_t s = 0;  s < 2; s++) 
      { 
        /* Get parameters for sex {s}: */
        gdr_demo_parms_t *dmp_goal = o->demoParms[s]; /* Goal parameters. */
        gdr_demo_show_parms("specified parameters", s, dmp_goal);

        /* Compute the goal child count probabilities for sex {s}: */
        double_vec_t cProb = gdr_demo_compute_child_count_distr(dmp_goal);
        demand(cProb.ne == dmp_goal->cMax + 1, "inconsistent {cMax,cProb}");
        gdr_demo_show_distr("goal child count", s, &(cProb), NULL);
        
        /* Create some perturbed parameters: */
        gdr_demo_parms_t *dmp_fake = gdr_tdm_create_fake_parms(dmp_goal);
        gdr_demo_show_parms("fake parameters", s, dmp_fake);

        /* Compute the actual child count frequencies: */
        double_vec_t cFreq = gdr_demo_compute_child_count_distr(dmp_fake);
        demand(cFreq.ne == dmp_fake->cMax + 1, "inconsistent {cMax,cProb}");
        gdr_demo_show_distr("goal and fake child count", s, &(cProb), &(cFreq));

        /* Write TeX table with child count probs: */
        char *outPrefix = "out/gdr-tcp";
        char *tag = "A";
        gdr_demo_write_child_distr_tex_table(outPrefix, tag, s, &(cProb), dmp_goal->cPrec, &(cFreq), dmp_fake->cPrec);
        
        gdr_demo_parms_free(dmp_goal);
        gdr_demo_parms_free(dmp_fake);
        free(cProb.e);
        free(cFreq.e);
      }

    return 0;
  }

gdr_demo_parms_t *gdr_tdm_create_fake_parms(gdr_demo_parms_t *dmp)
  {
    gdr_demo_parms_t *dmp_fake = gdr_demo_parms_new(); /* Perturbed parameters. */
    (*dmp_fake) = (*dmp);
    
    dmp_fake->cMax = dmp->cMax + 2;
    dmp_fake->cAlpha = 0.80*dmp->cAlpha;
    dmp_fake->cPrec = dmp->cPrec + 2;
    dmp_fake->fMin = dmp->fMin - 2;
    dmp_fake->fMax = dmp->fMax + 2;
    
    return dmp_fake;
  }

gdr_tdm_options_t *gdr_tdm_options_parse(int32_t argc, char **argv)
  { 
    argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, PROG_HELP);
    argparser_set_info(pp, PROG_INFO);
    argparser_process_help_info_options(pp);
    
    gdr_tdm_options_t *o = notnull(malloc(sizeof(gdr_tdm_options_t)), "no mem");

    /* The keyword "-demoParms" must appear exactly twice, once for each sex: */
    for (uint32_t s = 0;  s < 2; s++)
      { argparser_get_keyword(pp, "-demoParms");
        o->demoParms[s] = gdr_demo_options_parse(pp);
      }

    argparser_skip_parsed(pp);
    argparser_finish(pp);
    return o;
  }
