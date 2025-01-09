#define PROG_NAME "tcoptim"
#define PROG_DESC "exploring the test-and-complement model of computation"
#define PROG_VERS "1.1"

/* Last edited on 2024-12-21 11:30:53 by stolfi */

#define tcoptim_C_COPYRIGHT "Copyright © 2006 by the State University of Campinas (UNICAMP)"
#define tcoptim_C_AUTHORS "Created 2006-may-20 by Jorge Stolfi, IC-UNICAMP"

#define PROG_HELP \
  "  " PROG_NAME " \\\n" \
  "  -size {M} {N} \\\n" \
  "  -maxBifuns {MAXBIFUNS} \\\n" \
  "  -maxGates {MAXGATES} \\\n" \
  "  [ -seed {SEED} ] \\\n" \
  "  -outName {OUTNAME}"

#define PROG_INFO \
  "NAME\n" \
  "  " PROG_NAME " - " PROG_DESC "\n" \
  "\n" \
  "SYNOPSIS\n" \
  PROG_HELP "\n" \
  "\n" \
  "DESCRIPTION\n" \
  "  The program bla bla bla" \
  " bla bla bla bla {X+Y} bla bla" \
  " bla {INFILE} bla \"foobar.ppm\" bla bla bla\n" \
  "\n" \
  "  Beware that bla bla bla BLEBBLE BLOB bla" \
  " bla bla bla bla.\n" \
  "\n" \
  "OPTIONS\n" \
  "  -size {M} {N}\n" \
  "    This mandatory switch specifies the number of input bits {M} and the number of" \
  " internal state bits (TC gate width) {N}.\n" \
  "\n" \
  "  -maxBifuns {MAXBIFUNS}\n" \
  "    This optional switch specifies the maximum number of bifuns to generate.  The" \
  " default is " stringify(MAX_BIFUNS) ".\n" \
  "\n" \
  "  -maxGates {MAXGATES}\n" \
  "    This optional switch specifies the maximum number of gates to consider in any circuit.  The" \
  " default is no limit.\n" \
  "\n" \
  "  -seed {SEED}\n" \
  "    This optional switch specifies the seed for the random number generator.  The default is 4615.\n" \
  "\n" \
  "  -outname {OUTNAME}\n" \
  "    This mandatory switch specifies the common prefix for all output file names.\n" \
  "\n" \
  "DOCUMENTATION OPTIONS\n" \
  argparser_help_info_HELP_INFO "\n" \
  "\n" \
  "SEE ALSO\n" \
  "  ls(1).\n" \
  "\n" \
  "AUTHOR\n" \
  "  " tcoptim_C_AUTHORS ".\n" \
  "\n" \
  "MODIFICATION HISTORY\n" \
  "  2008-02-11 revised by J. Stolfi, IC-UNICAMP.\n" \
  "\n" \
  "WARRANTY\n" \
  argparser_help_info_NO_WARRANTY "\n" \
  "\n" \
  "RIGHTS\n" \
  "  ." tcoptim_C_COPYRIGHT "\n" \
  "\n" \
  argparser_help_info_STANDARD_RIGHTS

#define stringify(x) strngf(x)
#define strngf(x) #x

/* Must define _GNU_SOURCE in order to get {asprintf} */
#define _GNU_SOURCE
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <values.h>
#include <math.h>

#include <bool.h>
#include <affirm.h>
#include <jsfile.h>
#include <argparser.h>
#include <jsrandom.h>

#include <tcgates.h>
#include <tc_opt.h>

#define MAX_BIFUNS 10000000
  /* Max number of bifuns to generate, if user does not specify otherwise. */

typedef struct options_t 
  { int seed;        /* Seed for random number gnerators. */
    int m;           /* Number of input bits. */
    int n;           /* Number of state bits (TC circuit width). */
    bifun_ct_t maxBifuns;  /* Max number of bifuns to generate. */
    int maxGates;          /* Max number of TC gates in circuit. */
    char *outName;   /* Prefix for output file names. */
  } options_t;

int main (int argc, char **argv);
options_t *get_options(int argc, char **argv);

int main (int argc, char **argv)
  { 
    options_t *o = get_options(argc, argv);
    srandom(o->seed);
    
    tc_opt_table_t *tb = tc_opt_table_build(o->m, o->n, o->maxBifuns, o->maxGates);
    
    char *fname = jsprintf("%s-%d-%d-opt.txt", o->outName, o->m, o->n);
    FILE *wr = open_write(fname, TRUE);
    
    tc_opt_table_print_table(wr, tb);

    fclose(wr);
    free(fname);

    return 0;
  }

options_t *get_options(int argc, char **argv)
  {
    argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, PROG_NAME " version " PROG_VERS ", usage:\n" PROG_HELP);
    argparser_set_info(pp, PROG_INFO);
    argparser_process_help_info_options(pp);
    
    options_t *o = (options_t *)malloc(sizeof(options_t)); 
    
    argparser_get_keyword(pp, "-size");
    o->m = argparser_get_next_int(pp, 0, m_max);
    o->n = argparser_get_next_int(pp, o->m, n_max);
    
    if (argparser_keyword_present(pp, "-maxBifuns"))
      { o->maxBifuns = argparser_get_next_uint(pp, 1, ULLONG_MAX); }
    else
      { o->maxBifuns = MAX_BIFUNS; }
    
    if (argparser_keyword_present(pp, "-maxGates"))
      { o->maxGates = argparser_get_next_uint(pp, 1, ULLONG_MAX); }
    else
      { o->maxGates = INT_MAX; }
    
    if (argparser_keyword_present(pp, "-seed"))
      { o->seed = argparser_get_next_int(pp, 1, INT_MAX); }
    else
      { o->seed = 4615; }
    
    argparser_get_keyword(pp, "-outName");
    o->outName = argparser_get_next(pp);
    
    argparser_finish(pp);
    
    return o;
  }

