/* Last edited on 2007-08-15 22:41:13 by stolfi */
/* See the authorship and copyright notice at the end of this file. */

/* COMMAND LINE PARSING */

#ifndef jcloptions_h
#define jcloptions_h

#include <jclbasic.h>
#include <bool.h>

bool_t GetLongOption(
    int argc, char **argv,
    int *ip,
    char *brief,
    char *wordy,
    long min,
    long max,
    long *vp
  );
  /*
    If "argc[i]" matches either the "brief" or the "wordy" keyword, and is
    followed by one additional parameter, parses that parameter as a long
    integer, stores it in "*vp", increments "*i" by 2, and returns
    TRUE. Otherwise leaves "*vp" unchanged and returns FALSE. Bombs if
    the value is not an integer between "min" and "max".
  */
     
bool_t GetCharOption(
    int argc, char **argv,
    int *ip,
    char *brief,
    char *wordy,
    char **vp
  );
  /*
    If "argc[i]" matches either the "brief" or the "wordy" keyword, and is
    followed by one additional parameter, then makes a copy of that paramter,
    stores its addres in "*vp", increments "*i" by 2, and returns
    TRUE. Otherwise leaves "*vp" unchanged and returns FALSE.
  */
     
bool_t GetBoolOption(
    int argc, char **argv,
    int *ip,
    char *brief,
    char *wordy,
    bool_t *vp
  );
  /*
    If "argc[i]" matches either the "brief" or the "wordy" keyword,
    stores TRUE in "*vp", increments "*i", and returns TRUE. Otherwise
    leaves "*vp" unchanged and returns FALSE.
  */

#endif
