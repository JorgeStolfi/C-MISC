/* Last edited on 1997-11-24 00:51:47 by stolfi */
/* See the authorship and copyright notice at the end of this file. */
/* Changes:

*/

#include <jclerror.h>
#include <jclbasic.h>
#include <stdio.h>
#include <stdlib.h>

void LineError(char *msg, char *name, long R)
  {
    fprintf(stderr, "%s: file %s, line %ld: %s\n", progName, name, R, msg);
    fflush(stderr);
    exit(1);
  }

void FileError(char *msg, char *name)
  {
    fprintf(stderr, "%s: file %s: %s\n", progName, name, msg);
    fflush(stderr);
    exit(1);
  }

void ElementError(char *msg, char *name, long R, long K)
  {
    fprintf(stderr, "%s: file %s, line %ld, field %ld: %s\n", progName, name, R, K, msg);
    fflush(stderr);
    exit(1);
  }

void Error(char *msg)
  {
    fprintf(stderr, "%s: %s\n", progName, msg);
    fflush(stderr);
    exit(1);
  }

