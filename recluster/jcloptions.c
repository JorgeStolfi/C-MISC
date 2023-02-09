/* Last edited on 1997-11-24 00:51:47 by stolfi */
/* See the authorship and copyright notice at the end of this file. */
/* Changes:

*/

#include <jcloptions.h>
#include <jclbasic.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <jclerror.h>

/* INTERNAL PROTOTYPES */

/* PROCEDURES */

bool_t GetLongOption(
    int argc, char **argv,
    int *ip,
    char *brief,
    char *wordy,
    long min,
    long max,
    long *vp
  )
  {
    int i = (*ip);
    int res; long val; char dum = '\000';
    if ((i+1<argc) && ((strcmp(argv[i], brief)==0) || (strcmp(argv[i], wordy)==0)))
      { res = sscanf(argv[i+1], "%ld%c", &val, &dum); 
        if ((res < 1) || (dum != '\000')) { ParamError("bad value", argv[i], i); }
        if ((val < min) || (val > max)) { ParamError("value out of range", argv[i], i); }
        (*vp) = val;
        (*ip) += 2;
        return (TRUE);
      }
    else
      { return (FALSE); }
  }
     
bool_t GetCharOption(
    int argc, char **argv,
    int *ip,
    char *brief,
    char *wordy,
    char **vp
  )
  {
    int i = (*ip);
    if ((i+1<argc) && ((strcmp(argv[i], brief)==0) || (strcmp(argv[i], wordy)==0)))
      { (*vp) = CopyString(argv[i+1]); 
        (*ip) += 2;
        return (TRUE);
      }
    else
      { return (FALSE); }
  }
  
bool_t GetBoolOption(
    int argc, char **argv,
    int *ip,
    char *brief,
    char *wordy,
    bool_t *vp
  )
  {
    int i = (*ip);
    if ((i<argc) && ((strcmp(argv[i], brief)==0) || (strcmp(argv[i], wordy)==0)))
      { (*vp) = TRUE;
        (*ip) += 1;
        return (TRUE);
      }
    else
      { return (FALSE); }
  }
     
