/* Last edited on 2013-10-31 02:19:12 by stolfilocal */

#define PROG_NAME "bmix"
#define PROG_DESC "Reads and writes two files, swapping their bits according to a key file."
#define PROG_VERS "2013-10-31"

#define bmix_C_copyright \
  "Copyright © 2013 Jorge Stolfi, State University of Campinas"

#define PROG_HELP \
  PROG_NAME " [ -verbose | --verbose ] \\\n" \
  "    {INFILE_A} {INFILE_B} \\\n" \
  "    {KEYFILE} \\\n" \
  "    {OUTFILE_A} {OUTFILE_B}"

#define PROG_INFO \
  "NAME\n" \
  "  " PROG_NAME " - " PROG_DESC "\n" \
  "\n" \
  "SYNOPSIS\n" \
  PROG_HELP "\n" \
  "\n" \
  "DESCRIPTION\n" \
  "  Reads two data files, named {INFILE_A} and {INFILE_B}, and a" \
  " control file, named {KEYFILE}.  For every bit of the latter" \
  " that is '1', swaps the corresponding bits of {INFILE_A} and" \
  " {INFILE_B}.  Then writes the resulting byte streams to" \
  " files {OUTFILE_A} and {OUTFILE_B}, respectively..\n" \
  "\n" \
  "  The two input files must have exacly the same" \
  " length, and the {KEYFILE} must be at least as long" \
  " as them, otherwise the program fails with an error" \
  " message and a non-zero exit code.  If the {KEYFILE}" \
  " is longer than the input files, its excess bytes" \
  " are ignored.\n" \
  "\n" \
  "  Any input file name may be \"-\" to mean 'read from" \
  " standard input'.  Any output file name may be \"-\" to" \
  " mean 'write to standard output'.\n" \
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
  "  bxor(1).\n" \
  "\n" \
  "AUTHOR\n" \
  "  Created October 2013 by Jorge Stolfi, IC-UNICAMP.\n" \
  "\n" \
  "WARRANTY\n" \
  argparser_help_info_NO_WARRANTY "\n" \
  "\n" \
  "RIGHTS\n" \
  "  " bmix_C_copyright ".\n" \
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
  { char *infile_a;   /* Name of first input file. */
    char *infile_b;   /* Name of second input file. */
    char *keyfile;    /* Name of key file. */
    char *outfile_a;  /* Name of first output file. */
    char *outfile_b;  /* Name of second output file. */
    bool_t verbose;   /* The "-verbose" flag. */
  } options_t;
  /* Arguments parsed from the command line. */
 
options_t *get_options(int argc, char **argv);
  /* Parses the command line arguments. */
 
FILE *get_in_file(char *name, bool_t verbose);
  /* If {name} is "-", returns {stdin}, otherwise returns 
    a file descriptor on file {name}, open for reading.  In this second case,
    prints a message to {stderr} and returns NULL if the open failed. */
 
FILE *get_out_file(char *name, bool_t verbose);
  /* If {name} is "-", returns {stdout}, otherwise returns 
    a file descriptor on file {name}, open for writing.  In this second case,
    prints a message to {stderr} and returns NULL if the open failed. */

int main(int argc, char **argv);

int main(int argc, char **argv)
  {
    options_t *o = get_options(argc, argv);
    
    /* Open the input files: */
    FILE *a = get_in_file(o->infile_a, o->verbose);
    FILE *b = get_in_file(o->infile_b, o->verbose);
    FILE *x = get_in_file(o->keyfile, o->verbose);
    if ((a == NULL) || (b == NULL) || (x == NULL)) { return 1; }
    if ((a == b) || (b == x) || (b == x))
      { fprintf(stderr, "cannot read two files from the same descriptor"); return 1;  }
    
    /* Open the output files: */
    FILE *r = get_out_file(o->outfile_a, o->verbose);
    FILE *s = get_out_file(o->outfile_b, o->verbose);
    if ((r == NULL) || (s == NULL)) { return 1; }
    if (r == s) 
      { fprintf(stderr, "cannot write two files into the same descriptor"); return 1;  }

    /* Process the data: */
    int64_t nw = cryptoy_file_mix(a, b, x, r, s);
    
    /* Check outcome: */
    if (nw < 0) 
      { fprintf(stderr, "aborted.\n"); return 1; }
    else
      { if (o->verbose) { fprintf(stderr, "wrote %ld bytes\n", nw); }
        /* Close/flush output files: */
        if (r == stdout) { fflush(r); } else { fclose(r); }
        if (s == stdout) { fflush(s); } else { fclose(s); }
        return 0;
      }
   }

FILE *get_in_file(char *name, bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "opening \"%s\" for input... \n", name); }
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

FILE *get_out_file(char *name, bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "opening \"%s\" for output... \n", name); }
    if (strlen(name) == 0)
      { fprintf(stderr, "%s: ** empty file name\n", PROG_NAME); return NULL; }
    else if (strcmp(name, "-") == 0)
      { return stdout; }
    else
      { FILE *d = fopen(name, "wb");
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
    o->infile_a = argparser_get_next_non_keyword(pp);
    o->infile_b = argparser_get_next_non_keyword(pp);
    o->keyfile = argparser_get_next_non_keyword(pp);
    o->outfile_a = argparser_get_next_non_keyword(pp);
    o->outfile_b = argparser_get_next_non_keyword(pp);
    
    /* Check for extraneous parameters: */
    argparser_finish(pp);
    
    return o;
  }
