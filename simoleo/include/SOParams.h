/* SOParams.h -- facilities for parsing command line arguments. */
/* Last edited on 2007-01-04 00:20:40 by stolfi */

/* Based on ParseParams.i3 by J.Stolfi, DEC-SRC, 1988.  */

#ifndef SOParams_H
#define SOParams_H

/* This interface provides simple and robust tools for parsing
  the command line arguments given to a process when it is 
  started.

  NOTE: Before reading the details, check the usage example at the
  end of this interface. */

#include <SOBasic.h>
#include <r2.h>
#include <r3.h>
#include <r4.h>
#include <vec.h>
#include <stdio.h>

typedef struct SOParams_T /* A parser for command line arguments. */
  { string_vec_t arg;   /* Arguments given, incl. the command name {arg.e[0]}. */
    bool_vec_t parsed;  /* {parsed.e[i]} is {TRUE} if {arg.e[i]} has been parsed. */
    unsigned next;      /* The next argument to parse is {arg[next]} */
    FILE *wr;           /* File for errors */
    unsigned nusage;    /* Number of lines of error help text. */
    string_vec_t usage; /* {usage.e[0..nusage-1]}  is the help text for errors. */
  } SOParams_T;
  
SOParams_T *SOParams_NewT(FILE *wr, int argn, char **argc);
  /* Allocates the arrays {arg} and {parsed} and
   initializes them with the parameters of the current
   process.  Marks {arg[0]} as parsed, all others as unparsed, 
   and sets {next} to 1.  Any subsequent parsing errors
   will be printed out to {wr}. */

void SOParams_SetUsage(SOParams_T *pp, char *msg);
  /* Appends {msg} to the help text to be written to {pp->wr} 
    just before halting, in case of argument syntax error. */

bool_t SOParams_KeywordPresent(SOParams_T *pp, char *key);
  /*  Looks for the first unparsed argument {arg[i]}
    that is equal to {key}.  If found, marks it as parsed, 
    sets {next} to {i+1}, and returns {TRUE}.
    Otherwise returns {FALSE} and leaves {next} unchanged. */

void SOParams_GetKeyword(SOParams_T *pp, char *key);
  /* Same as {keywordPresent}, but raises error if the 
    keyword is not found. */

bool_t SOParams_IsNext(SOParams_T *pp, char *key);
  /* Returns TRUE if and only if {arg[next]} exists, is unparsed,
    and is equal to {key}. */

bool_t SOParams_TestNext(SOParams_T *pp, char *key);
  /* If {IsNext(pp, key)}, marks the next argument as parsed,
    increments {next} and returns TRUE. Otherwise does none of
    these things and returns {FALSE}. */

void SOParams_MatchNext(SOParams_T *pp, char *key);
  /* If {IsNext(pp, key)}, marks the next argument as parsed and
    increments {next}. Otherwise raises an error. */

char *SOParams_GetNext(SOParams_T *pp);
  /* Returns {arg[next]}, marks it as parsed and increments {next}.  
    Raises error if {arg[next]} does not exist or has already 
    been parsed. */

int SOParams_GetNextInt(SOParams_T *pp, int min, int max);
double SOParams_GetNextDouble(SOParams_T *pp, double min, double max);
  /* Same as {SOParams_GetNext}, but converts the result to the approriate
    type (using {strtol} and {strtod}, respectively).  
    Raises error if the parameter is not a valid literal, or
    lies outside of the range {[min..max]}.  */

r2_t SOParams_GetNextR2(SOParams_T *pp, double min, double max);
r3_t SOParams_GetNextR3(SOParams_T *pp, double min, double max);
r4_t SOParams_GetNextR4(SOParams_T *pp, double min, double max);
  /* Parses with {SOParams_GetNextDouble} the next {N} arguments
    ({N=2,3,4}), which must be numbers in the range {[min .. max]}. */

void SOParams_GetNextRN(SOParams_T *pp, int N, double min, double max, double *x);
  /* Parses with {SOParams_GetNextDouble} the next {N} arguments,
    which must be numbers in the range {[min .. max]}. */

r2_t SOParams_GetNextR2Dir(SOParams_T *pp);
r3_t SOParams_GetNextR3Dir(SOParams_T *pp);
  /* Same as {SOParams_GetNextR2} etc., but normalizes result to unit length. */

int_vec_t SOParams_GetIntList(SOParams_T *pp, char *key, int min, int max);
double_vec_t SOParams_GetDoubleList(SOParams_T *pp, char *key, double min, double max);
string_vec_t SOParams_GetList(SOParams_T *pp, char *key);
  /* Parses all (zero or more) unparsed occurrences of the keyword
    {key}, not necessarily in consecutive positions. Requires that
    each occurrence is immediately followed by a value. Returns an
    array with those values, in the order found. The first two procs
    require that each value be a number in {[min..max]}. */

void SOParams_Error(SOParams_T *pp, char *msg);
  /* Prints the given message, the help text, and terminates the program. */

void SOParams_SkipParsed(SOParams_T *pp);
  /* Points {next} at the first unparsed argument.
    If there are parsed arguments beyond that one,
    prints a message and raises error. */

void SOParams_Finish(SOParams_T *pp);
  /* Checks if all parameters have been parsed; if not,
    prints a message and raises error. */

#endif
