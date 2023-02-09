/* ERROR MESSAGES */

static char *progName;
  /* 
    Must be set by the main program before it tries to print
    any error message. 
  */

/* Last edited on 1997-11-24 00:51:47 by stolfi */
/* See the authorship and copyright notice at the end of this file. */

/* ERROR REPORTING PROCEDURES */

#ifndef jclerror_h
#define jclerror_h

void Error(char *msg);
void ParamError(char *msg, char *key, long K);
void FileError(char *msg, char *name);
void LineError(char *msg, char *name, long R);
void ElementError(char *msg, char *name, long R, long K);
  /*
    Print error message about input record "R" and element/parameter "K",
    then die with nonzero status.
  */

#endif
