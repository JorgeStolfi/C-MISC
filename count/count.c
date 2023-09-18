#define PROG_NAME "count"
#define PROG_DESC "Outputs a count sequence. A sane replacement for GNU's {seq}"
#define PROG_VERS "2023-09"

#define ct_C_COPYRIGHT \
  "Copyright © 2023 by the State University of Campinas (UNICAMP)"
/* Last edited on 2023-09-18 00:37:30 by stolfi */

#define PROG_HELP \
  "  " PROG_NAME " \\\n" \
  "    [ -plus ] \\\n" \
  "    [ -zeros ] \\\n" \
  "    " argparser_help_info_HELP " \\\n" \
  "    FIRST LAST [ STEP [ WIDTH [ PREC ] ] ]" \

#define PROG_INFO \
  "NAME\n" \
  "  " PROG_NAME " - " PROG_DESC "\n" \
  "\n" \
  "SYNOPSIS\n" \
  PROG_HELP "\n" \
  "\n" \
  "DESCRIPTION\n" \
  "  The program writes to standard output a sequence of numbers in" \
  " arithmetic progression, starting with {FIRST} up to but not" \
  " beyond {LAST}, with the given {STEP}.\n" \
  "\n" \
  "  The {STEP} may be negative (but not zero). If {STEP} is not" \
  " given, it defaults to 1.  The program ends without printing" \
  " anything if {LAST} is less than {FIRST} and {STEP} is" \
  " positive, or {LAST} is greater than {FIRST} and \n{STEP} is negative.\n" \
  "\n" \
  "  Any or all of {FIRST}, {LAST}, and {STEP} may be" \
  " fractional, but then the end of the count may be" \
  " affected by roundoff errors.  Exponential" \
  " notation (as in \"3.14E10\" or \"5e-10\") or embedded" \
  " blanks (as in \"+ 3\") are not allowed.  Integer" \
  " arguments up to 2^50 or so in absolute value will be exact." \
  "\n" \
  "NUMBER FORMAT\n" \
  "  Each number is printed on a separate line, with at least {WIDTH} total" \
  " bytes.  If {PREC} is nonzero, it must be positive, and then the output value" \
  " includes a period \".\" and {PREC} decimal fraction digits.  A negative value" \
  " will always be printed with a '-' sign, even if the specified {PREC} is such" \
  " that the printed value is rounded to zero.  Optionally, the value" \
  " may be padded with blanks or zeros, and positive values may be" \
  " preceded by a sign '+'; see the \"-plus\" and \"-zeros\" options below.\n" \
  "\n" \
  "  The default {WIDTH} is the maximum number of characters" \
  " in any of {FIRST}, {LAST}, or {STEP}.  Leading blanks are" \
  " considered for this purpose; so, for" \
  " example, in \"count ' 1' 10 ' 001.0'\" the {WIDTH} defaults to 6. \n" \
  "\n" \
  "  If any of {FIRST}, {LAST}, or {STEP} includes a" \
  " period '.', at least one digit must be present after" \
  " the period, and then the default {PREC} is the max number of" \
  " decimal fraction digits among those values. Otherwise {PREC} defaults to zero.\n" \
  "\n" \
  "  The program always assumes the \"C\" locale conventions, ignoring" \
  " the external locale settings.  In particular, the decimal point is always '.' and no" \
  " thousands-separators are used.\n" \
  "\n" \
  "OPTIONS\n" \
  "  -plus\n" \
  "    If this option is specified, positive values will be printed" \
  " with an explicit '+' sign, and a value of exactly zero will be" \
  " printed with an extra leading zero.  The sign '+' will be printed" \
  " even if the specified {PREC} is such that the value rounds to zero.  This" \
  " option is implicitly assumed if any of {FIRST}, {LAST}, or {STEP} has a '+' sign" \
  " before the first digit, as in \"+3\", \"+02.5\" or \"+0\".\n" \
  "\n" \
  "  -zeros\n" \
  "    If this option is specified, numbers will be padded to {WIDTH} bytes" \
  " with zeros (after the sign, if any).  If omitted, blanks (before the" \
  " sign, if any) are used instead.  This option is implicitly assumed if any" \
  " of {FIRST}, {LAST}, or {STEP} begins with a redundant zero digit, as" \
  " in \"012\", \"-03.1\", \"+005\", or \"00.000\".  However, this" \
  " default is not triggered if the integer part is a single '0' or omitted, as" \
  " in \"0\", \"+0\", \"-0\", \"+0.002\", or \".05\".\n" \
  "\n" \
  "DOCUMENTATION OPTIONS\n" \
  argparser_help_info_HELP_INFO "\n" \
  "\n" \
  "SEE ALSO\n" \
  "  seq(1).\n" \
  "\n" \
  "AUTHOR\n" \
  "  Created on 2023-09-17 by Jorge Stolfi, IC-UNICAMP.\n" \
  "\n" \
  "MODIFICATION HISTORY\n" \
  "  Changes are by J.Stolfi unless noted.\n" \
  "    2023-09-17 by J.Stolfi: created from an earlier private version.\n" \
  "\n" \
  "WARRANTY\n" \
  argparser_help_info_NO_WARRANTY "\n" \
  "\n" \
  "RIGHTS\n" \
  "  " ct_C_COPYRIGHT ".\n" \
  "\n" \
  argparser_help_info_STANDARD_RIGHTS

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <locale.h>
#include <assert.h>

#include <bool.h>
#include <argparser.h>
#include <jsstring.h>

#define ct_strlen_MAX 1024
  /* A safety limit for argument length, {WIDTH}, etc. */

typedef struct ct_options_t
  { double FIRST;
    double LAST;
    double STEP;
    int32_t WIDTH;
    int32_t PREC;
    bool_t plus;
    bool_t zeros;
  } ct_options_t;

void ct_print(double x, int32_t WIDTH, int32_t PREC, bool_t plus, bool_t zeros);
  /* Prints the value of {x} as specified by the options {WIDTH,PREC,plus,zeros} */

ct_options_t *ct_parse_options(int32_t argc, char **argv);
  /* Parses the command line options. */

void ct_parse_seq_arg(argparser_t *pp, char *name, double *val_P, bool_t required, ct_options_t *o);
  /* Parses the next unparsed command line argument as one of the
    parameters {FIRST}, {LAST}, or {STEP}. Checks for valid syntax, and
    reasonable values. The value parsed is returned in {*val_P}. Also
    updates {o->WIDTH} and {o->PREC} to be at least the length and
    number of fraction digits of that argument. Sets {o->plus} if the
    number has an explicit '+' sign, and {o->zeros} if it has a
    superfluous leading '0' followed by another digit. */

void ct_parse_int_arg(argparser_t *pp, char *name, int32_t *val_P, int32_t max);
  /* Parses the next unparsed command line argument as one of the
    parameters {WIDTH} or {PREC}. Checks for valid syntax,
    and value range {0..max}. The value parsed is returned in {*val_P}. */

int32_t main (int32_t argc, char **argv)
  {
    setlocale(LC_ALL, "C");

    ct_options_t *o = ct_parse_options(argc, argv);

    double stepsgn = ( o->STEP == 0.0 ? 0 : ( o->STEP < 0 ? -1.0 : 1.0 ) );

    /* fprintf(stderr, "seq args = %g %g %g\n", o->FIRST, o->LAST, o->STEP); */
    /* fprintf(stderr, "fmt args = %d %d\n", o->WIDTH, o->PREC); */

    double x = o->FIRST;
    while ((stepsgn*(x - o->LAST) <= 0))
      { ct_print(x, o->WIDTH, o->PREC, o->plus, o->zeros);
        x += o->STEP;
      }
    return 0;
}

void ct_print(double x, int32_t WIDTH, int32_t PREC, bool_t plus, bool_t zeros)
  { if (plus)
      { if (x == 0)
          { /* Must print a superfluous '0' explicitly: */
            int32_t wd1 = (WIDTH == 0 ? 0 : WIDTH - 1);
            if (zeros)
              { printf("0%0*.*f", wd1, PREC, x); }
            else
              { printf("0%*.*f", wd1, PREC, x); }
          }
        else
          { /* Let {printf} take care of the '+' sign: */
            if (zeros)
              { printf("%+0*.*f", WIDTH, PREC, x); }
            else
              { printf("%+*.*f", WIDTH, PREC, x); }
          }
      }
    else
      { if (zeros)
          { printf("%0*.*f", WIDTH, PREC, x); }
        else
          { printf("%*.*f", WIDTH, PREC, x); }
      }
    putchar('\n');
  }

ct_options_t *ct_parse_options(int32_t argc, char **argv)
  {
    /* Initialize argument parser: */
    argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, PROG_NAME " version " PROG_VERS ", usage:\n" PROG_HELP);
    argparser_set_info(pp, PROG_INFO);
    argparser_process_help_info_options(pp);

    /* Allocate the command line argument record: */
    ct_options_t *o = notnull(malloc(sizeof(ct_options_t)), "no mem");

    /* Parse keyword options first: */
    o->plus = argparser_keyword_present(pp, "-plus");
    o->zeros = argparser_keyword_present(pp, "-zeros");

    /* Now parse positional parameters, and provide defaults */
    /* for \"-plus\" and \"-zeros\" is not specified. */
    o->WIDTH = 0;
    o->PREC = 0;
    o->STEP = 1;

    ct_parse_seq_arg(pp, "FIRST", &(o->FIRST), TRUE, o);
    ct_parse_seq_arg(pp, "LAST", &(o->LAST), TRUE, o);
    ct_parse_seq_arg(pp, "STEP", &(o->STEP), FALSE, o);

    if (o->STEP == 0)
      { argparser_error(pp, "the {STEP} cannot be zero"); }
    if ((o->FIRST + o->STEP == o->FIRST) || (o->LAST - o->STEP == o->LAST))
      { argparser_error(pp, "the {STEP} is too small"); }

    ct_parse_int_arg(pp, "WIDTH", &(o->WIDTH), ct_strlen_MAX);
    ct_parse_int_arg(pp, "PREC", &(o->PREC), ct_strlen_MAX);

    /* Check for spurious arguments: */
    argparser_finish(pp);

    return o;
  }

void ct_parse_seq_arg(argparser_t *pp, char *name, double *val_P, bool_t required, ct_options_t *o)
  { /* Check whether there is a next argument: */
    char *arg = argparser_next(pp);
    if (arg == NULL)
      { if (required) 
          { argparser_error(pp, txtcat3("argument ", name, " is required")); }
        else
          { /* Assume that {*val_P} is set to the proper default. */
            return;
          }
      }
    /* Get the argument and mark it parsed: */
    arg = argparser_get_next(pp); 

    /* Update the default width: */
    int64_t n = strlen(arg);
    if (n == 0) { argparser_error(pp, txtcat3("empty ", name, " not allowed")); }
    if (n > ct_strlen_MAX) { argparser_error(pp, txtcat3("argument ", name, " is too long")); }
    if (n > o->WIDTH) { o->WIDTH = (int32_t)n; }

    /* Does it have a decimal point? */
    char *p = strchrnul(arg, '.');
    assert((p != NULL) && (p > arg));
    if ((*p) != 0)
      { /* There must be at least one digit after the '.': */
        char *f = p+1;
        if (((*f) == 0) || ((*f) < '0') || ((*f) > '9'))
          { argparser_error(pp, txtcat("digits are required after '.' in ", arg)); }
        /* Update the default precision: */
        int32_t prc = (int32_t)(n - (p - arg) - 1);
        assert(p > 0);
        if (prc > o->PREC) { o->PREC = prc; }
      }
    /* Now {p} points at the decimal '.' or to the final {NUL}. */
    /* Find first digit: */
    while (p > arg)
      { char *q = p-1;
        if (((*q) < '0') || ((*q) > '9')) { break; }
        p = q;
      }
    /* Now {p} points at first digit or the '.' or final {NUL}: */
    if ((*p) == '0')
      { char *q = p+1;
        if (((*q) >= '0') && ((*q) <= '9'))
          { /* Starts with a redundant '0', set {o->zeros}: */
            o->zeros = TRUE;
          }
      }
    if (p > arg)
      { char *q = p-1;
        if ((*q) == '+')
          { /* Redundant '+' sign, set {o->plus}: */
            o->plus = TRUE;
          }
      }
    /* Now parse the value as a {double}: */
    char *rest = NULL;
    double val = strtod(arg, &rest);
    if ((*rest) != 0)
      { argparser_error(pp, txtcat3("invalid value \"", arg, txtcat("\" for ", name))); }
    if (! isfinite(val))
      { argparser_error(pp, txtcat3("value of ", name, " must be finite")); }
    (*val_P) = val;
  }

void ct_parse_int_arg(argparser_t *pp, char *name, int32_t *val_P, int32_t max)
  { /* Check whether there is a next argument: */
    char *arg = argparser_next(pp);
    if (arg == NULL) { return; }

    /* Get the argument and mark it parsed: */
    arg = argparser_get_next(pp); 

    /* Now parse the value as an {int32_t}: */
    char *rest = NULL;
    int64_t val = strtol(arg, &rest, 10);
    if ((*rest) != 0)
      { argparser_error(pp, txtcat3("invalid value \"", arg, txtcat("\" for ", name))); }
    if (val < 0)
      { argparser_error(pp, txtcat3("value of ", name, " cannot be negative")); }
    if (val > max)
      { argparser_error(pp, txtcat3("value of ", name, " is too large")); }
    (*val_P) = (int32_t)val;
  }



