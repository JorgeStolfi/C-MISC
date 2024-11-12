/* Last edited on 2013-10-31 02:17:14 by stolfilocal */

#define PROG_NAME "bxor"
#define PROG_DESC "Reads two files, outputs their bitwise XOR."
#define PROG_VERS "2013-10-31"

#define bxor_C_copyright \
  "Copyright © 2013 Jorge Stolfi, State University of Campinas"

#define PROG_HELP \
  PROG_NAME " [ -verbose | --verbose ] {FILE_A} {FILE_B} > {OUTFILE}"

#define PROG_INFO \
  "NAME\n" \
  "  " PROG_NAME " - " PROG_DESC "\n" \
  "\n" \
  "SYNOPSIS\n" \
  PROG_HELP "\n" \
  "\n" \
  "DESCRIPTION\n" \
  "  Reads two files, {FILE_A} and {FILE_B}, and writes" \
  " to the standard output their bitwise XOR (exclusive OR).\n" \
  "\n" \
  "  If {FILE_B} is longer than {FILE_A}, the excess bytes" \
  " in {FILE_B} are ignored.  Conversely, if {FILE_A} is" \
  " longer than {FILE_B}, fails with an error message and" \
  " a non-zero exit code.\n" \
  "\n" \
  "  Either file name may be \"-\" to mean 'read from standard input'.\n" \
  "\n" \
  "OPTIONS\n" \
  "  -verbose\n" \
  "    If present, prints some information to standard error" \
  " output, such as the number of bytes written.\n" \
  "\n" \
  "DOCUMENTATION OPTIONS\n" \
  argparser_help_info_HELP_INFO "\n" \
  "\n" \
  "SEE ALSO\n" \
  "  bmix(1).\n" \
  "\n" \
  "AUTHOR\n" \
  "  Created October 2013 by Jorge Stolfi, IC-UNICAMP.\n" \
  "\n" \
  "MODIFICATION HISTORY\n" \
  "  2013-10-31:\n" \
  "    Added \"--verbose\" option. [J.Stolfi].\n" \
  "\n" \
  "WARRANTY\n" \
  argparser_help_info_NO_WARRANTY "\n" \
  "\n" \
  "RIGHTS\n" \
  "  " bxor_C_copyright ".\n" \
  "\n" \
  argparser_help_info_STANDARD_RIGHTS

#define _GNU_SOURCE
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <errno.h>

#include <affirm.h>
#include <argparser.h>

#include <cryptoy_file.h>

typedef struct options_t 
  { char *file_a;   /* Name of first input file. */
    char *file_b;   /* Name of second input file. */
    bool_t verbose; /* The "-verbose" flag. */
  } options_t;
  /* Arguments parsed from the command line. */
 
options_t *get_options(int argc, char **argv);
  /* Parses the command line arguments. */
 
FILE *get_in_file(char *name, bool_t verbose);
  /* If {name} is "-", returns {stdin}, otherwise returns 
    an open file descriptor on file {name}.  in this second case,
    prints a message to {stderr} and returns NULL if the open failed. */

int main(int argc, char **argv);

int main(int argc, char **argv)
  {
    options_t *o = get_options(argc, argv);

    /* Open the input files: */
    FILE *a = get_in_file(o->file_a, o->verbose);
    FILE *b = get_in_file(o->file_b, o->verbose);
    if ((a == NULL) || (b == NULL)) { return 1; }
    if (a == b) { fprintf(stderr, "cannot read both files from same file descriptor"); return 1;  }
    
    if (o->verbose) { fprintf(stderr, "writing to standard output\n"); }
    
    /* Process the data: */
    int64_t nw = cryptoy_file_xor(a, b, stdout);

    /* Check outcome: */
    if (nw < 0) 
      { fprintf(stderr, "aborted.\n"); return 1; }
    else
      { fflush(stdout);
        if (o->verbose) { fprintf(stderr, "wrote %ld bytes\n", nw); }
        return 0;
      }
   }

FILE *get_in_file(char *name, bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "opening \"%s\" for input...\n", name); }
    if (strlen(name) == 0)
      { fprintf(stderr, "%s: ** empty file name\n", PROG_NAME); return NULL; }
    else if (strcmp(name, "-") == 0)
      { return stdin; }
    else
      { FILE *d = fopen(name, "rb");
        if (d == NULL) { perror(PROG_NAME); }
        return d;
      }
  }

options_t *get_options(int argc, char **argv)
  { 
    /* Initialize argument parser: */
    argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, PROG_NAME " version " PROG_VERS ", usage:\n" PROG_HELP);
    argparser_set_info(pp, PROG_INFO);
    argparser_process_help_info_options(pp);

    /* Allocate the return record: */
    options_t *o = notnull(malloc(sizeof(options_t)), "no mem");

    /* Parse the keyword options: */
    o->verbose = 
      argparser_keyword_present(pp, "-verbose") ||
      argparser_keyword_present(pp, "--verbose");
    
    /* Parse the positional arguments after options: */
    argparser_skip_parsed(pp);
    o->file_a = argparser_get_next_non_keyword(pp);
    o->file_b = argparser_get_next_non_keyword(pp);
    
    /* Check for extraneous parameters: */
    argparser_finish(pp);
    
    return o;
  }
