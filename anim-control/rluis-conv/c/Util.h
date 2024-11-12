
#ifndef Util_H
#define Util_H



/* Miscellaneous hacks */

#include <Wr.h>

void Message(char *msg);
  /*
    Writes "msg" to "stderr". */

void Error(char *msg);
  /*
    Writes "msg" to "stderr" and bombs out. */

nat Digits(nat n);
  /*
    Number of digits in "Fmt.Int(n)" */

void ResetCPUTime();
void WriteCPUTime(FILE *wr);
  /*
    WriteCPUTime writes the CPU time elapsed since the 
    last call to "ResetCPUTime" or "WriteCPUTime", or since program
    startup if there were no such calls. */
;
} /* Util */.

#endif
