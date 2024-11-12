/* SPOptions.h -- facilities for parsing command line arguments. */
/* Last edited on 2005-10-06 21:48:43 by stolfi */

/* Based on Params.i3 by J.Stolfi, DEC-SRC, 1988.  */

#ifndef SPOptions_H
#define SPOptions_H

/* This interface provides simple and robust tools for parsing
  the command line arguments given to a process when it is 
  started.

  NOTE: Before reading the details, check the usage example at the
  end of this interface. */

#include <r3.h>
#include <r4.h>
#include <vec.h>
#include <bool.h>
#include <nat.h>

#include <stdio.h>

typedef struct SPOptions_Parser_t /* A parser for command line arguments. */
  { string_vec_t arg;   /* Arguments given, incl. the command name {arg.e[0]}. */
    bool_vec_t parsed;  /* {parsed.e[i]} is {TRUE} if {arg.e[i]} has been parsed. */
    nat_t next;         /* The next argument to parse is {arg[next]} */
    FILE *wr;           /* File for errors */
    nat_t nusage;       /* Number of lines of error help text. */
    string_vec_t usage; /* {usage.e[0..nusage-1]}  is the help text for errors. */
  } SPOptions_Parser_t;
  
SPOptions_Parser_t *SPOptions_NewParser(FILE *wr, int argn, char **argc);
  /* Allocates the arrays {arg} and {parsed} and
   initializes them with the parameters of the current
   process.  Marks {arg[0]} as parsed, all others as unparsed, 
   and sets {next} to 1.  Any subsequent parsing errors
   will be printed out to {wr}. */

void SPOptions_SetUsage(SPOptions_Parser_t *pp, char *msg);
  /* Appends {msg} to the help text to be written to {pp->wr} 
    just before halting, in case of argument syntax error. */

/* KEYWORD PARSING */

bool_t SPOptions_IsNext(SPOptions_Parser_t *pp, char *key);
  /* Returns TRUE if and only if {arg[next]} exists, is unparsed,
    and is equal to {key}.  Does not change the state of {pp}. */

bool_t SPOptions_TestKeyword(SPOptions_Parser_t *pp, char *key);
  /*  Looks for the first unparsed argument {arg[i]}
    that is equal to {key}.  If found, marks it as parsed, 
    sets {next} to {i+1}, and returns {TRUE}.
    Otherwise returns {FALSE} and leaves {next} unchanged.
    If {local = TRUE}, only */

void SPOptions_GetKeyword(SPOptions_Parser_t *pp, char *key);
  /* Same as {SPOptions_TestKeyword}, but raises error if the 
    keyword is not found. */

bool_t SPOptions_TestKeywordNext(SPOptions_Parser_t *pp, char *key);
  /* If {IsNext(pp, key)}, marks the next argument as parsed,
    increments {next} and returns TRUE. Otherwise does none of
    these things and returns {FALSE}. */

void SPOptions_GetKeywordNext(SPOptions_Parser_t *pp, char *key);
  /* If {IsNext(pp, key)}, marks the next argument as parsed and
    increments {next}. Otherwise raises an error. */

bool_t SPOptions_TestKeywordGen(SPOptions_Parser_t *pp, char *key, bool_t local);
  /* Equal to {SPOptions_TestKeyword} if {local=FALSE}, and 
    to {SPOptions_TestKeywordNext} if {local=TRUE}. */

void SPOptions_GetKeywordGen(SPOptions_Parser_t *pp, char *key, bool_t local);
  /* Equal to {SPOptions_GetKeyword} if {local=FALSE}, and 
    to {SPOptions_GetKeywordNext} if {local=TRUE}. */

/* ARGUMENT PARSING */

char *SPOptions_GetNext(SPOptions_Parser_t *pp);
  /* Returns {arg[next]}, marks it as parsed and increments {next}.  
    Raises error if {arg[next]} does not exist or has already 
    been parsed. */

int SPOptions_GetNextInt(SPOptions_Parser_t *pp, int min, int max);
double SPOptions_GetNextDouble(SPOptions_Parser_t *pp, double min, double max);
  /* Same as {SPOptions_GetNext}, but converts the result to the approriate
    type (using {strtol} and {strtod}, respectively).  
    Raises error if the parameter is not a valid literal, or
    lies outside of the range {[min..max]}.  */

r3_t SPOptions_GetNextR3(SPOptions_Parser_t *pp, double min, double max);
r4_t SPOptions_GetNextR4(SPOptions_Parser_t *pp, double min, double max);
  /* Parses with {SPOptions_GetNextDouble} the next {N} arguments,
    which must be in the range {[min .. max]}. */

r3_t SPOptions_GetNextDir(SPOptions_Parser_t *pp);
  /* Same as {SPOptions_GetNextR3} but normalizes result to unit length. */

int_vec_t SPOptions_GetIntList(SPOptions_Parser_t *pp, char *key, int min, int max);
  /* Parses all (zero or more) unparsed occurrences of the keyword
    {key}, not necessarily in consecutive positions. Requires that each
    occurrence is immediately followed by an integer in {[min..max]}.
    Returns an array with those integers, in the order found. */

void SPOptions_Error(SPOptions_Parser_t *pp, char *msg);
  /* Prints the given message, the help text, and terminates the program. */

void SPOptions_SkipParsed(SPOptions_Parser_t *pp);
  /* Points {next} at the first unparsed argument.
    If there are parsed arguments beyond that one,
    prints a message and raises error. */

void SPOptions_Finish(SPOptions_Parser_t *pp);
  /* Checks if all parameters have been parsed; if not,
    prints a message and raises error. */

/* FIX EXAMPLE !!!! */

/* 
  In Linux and other popular operating systems, most programs expect
  their command-line arguments to consist of a string of keywords and
  keyword-labeled arguments (`options', `switches', etc.), followed by
  a list of positional arguments.

  To help the user, programs generally allow the switches and
  keyword-labeled arguments to be given in any order.  Some of those
  parameters may be optional and/or repeatable, some may be
  mandatory; some may be required or forbidden depending on the
  values of the other parameters.  Furthermore, the value of an
  argument may be just a number or a text string, or may be a cluster
  of two or more values with their own little syntax.

  This module simplifies the parsing of such command-line parameters,
  by allowing the program to scan the arguments in their {logical
  order}.  It also detects automatically many kinds of common
  mistakes---arguments that are missing, repeated, extraneous,
  malformed, or out of range---and prints the appropriate error
  messages.

  For example, here is how this module could be used by an
  hypothetical program {prt} that concatenates a bunch of files and
  prints selected pages ranges of them, possibly in reverse order, with
  several formatting options, followed by an optional output file name.

  | int MaxPages = 10000;
  |
  | // Arguments from command line:
  | int fontSize;
  | bool_t landscape;
  | char *outFile;
  |
  | typedef struct Rec { 
  |   int ini;
  |   int fin;
  |   bool_t rev;
  |   char *file;
  |   struct *Rec next;
  | } Rec;
  |
  | void ParseOptions(int argc, char **argv)
  |   { SPOptions.NewT(stderr, argc, argv);
  |     SPOptions.SetUsage(
  |       "prt \\\n"
  |       "  -fontSize <n> \\\n"
  |       "  [ -landscape | -portrait ] \\\n"
  |       "  [ -file FILENAME | -pages NUM NUM [ -reverse ] ]... \\\n"
  |       "  [ OUTFILE ]\n"
  |     );
  |     // The {-fontSize } parameter is mandatory:
  |     SPOptions_TestKeyword(pp, "-fontSize");
  |     fontSize = SPOptions_GetNextInt(pp, 1, 100);
  |
  |     // Either {-landscape} or {-portrait}, but not both:
  |     // Default is {-portrait} unless font is too big:
  |     if (SPOptions_TestKeyword(pp, "-landscape")) 
  |       { landscape = TRUE; }
  |     else if (SPOptions_TestKeyword(pp, "-portrait")) 
  |       { landscape = FALSE; }
  |     else 
  |       { landscape = (fontSize > 8); }
  |
  |     // Parse the list of files and page ranges:
  |     Rec *ranges = NULL; // List of page ranges to print.
  |     char *file = "-";   // Current file name.
  |     while (TRUE)
  |       { if (SPOptions_TestKeywordNext(pp, "-pages"))
  |           { Rec *r = malloc(sizeof(Rec));
  |             if (r == NULL) { SPOptions_Error(pp, "no mem for arg list"); }
  |             r->ini = SPOptions_GetNextInt(pp, 1, MaxPages);
  |             r->fin = SPOptions_GetNextInt(pp, r->ini, MaxPages);
  |             r->rev = SPOptions_TestKeywordNext(pp, "-reverse");
  |             r->file = file;
  |             r->next = ranges; ranges = r;
  |           }
  |         else if (SPOptions_TestKeywordNext(pp, "-file"))
  |           { file = SPOptions_GetNext(pp); }
  |         else 
  |           { break; }
  |       }
  |
  |     // Provide a default range if none was given:
  |     if (ranges == NULL)
  |       { Rec *r = malloc(sizeof(Rec));
  |         if (r == NULL) { SPOptions_Error(pp, "no mem for arg list"); }
  |         r->ini = 1; r->fin = MaxPages;
  |         r->rev = FALSE;
  |         r->file = file;
  |         r->next = ranges; ranges = r;
  |       }
  |    
  |     // Reverse the range list:
  |     Rec *r = ranges; ranges = NULL;
  |     while (r != NULL) 
  |       { Rec *t = r->next; r->next = ranges; ranges = r; r = t; }
  |
  |     // Skip to end of parsed parameters:
  |     SPOptions_SkipParsed(pp);
  |     if (pp.next < pp.arg.ne) 
  |       { outFile = SPOptions_GetNext(pp); }
  |     else
  |       { outFile = "-"; }
  |
  |     // Check for any unparsed parameters:
  |     SPOptions_Finish(pp);
  |   }

  Note that even though this code parses the parameters 
  in a fixed order, the user may give the initial keyword parameters
  in any order, and even intermix them with the range list,
  or put them after the output file name. */
#endif
