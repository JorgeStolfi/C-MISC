/* Last edited on 2007-08-16 00:06:08 by stolfi */
/* See the authorship and copyright notice at the end of this file. */

/* ITEM PARSING ROUTINES */

#ifndef jclparse_h
#define jclparse_h

#include <stdio.h>
#include <jclbasic.h>

bool_t ReadItem(byte **L, long *N, long *A, long R);
  /* 
    Returns TRUE if "rd" (whose filename is "name") was at
    end-of-file.  Otherwise reads the next line of "rd" into "*L", and
    stores its length into "*N".  Assumes the buffer "*L" is currently
    "*A" bytes long; reallocates a new one if necessary.  Stops with
    error if the last line does not end with newline.
    Reports "R" as line number in case of errors.
  */
  
void WriteItem(byte *L, long N);
  /*
    Writes bytes "L[0..N-1]" to stdout, followed by newline.
  */
  
byte* CopyItem(byte *L, long N);
  /*
    Allocates an "N"-byte string and copies "L[0..N-1]" into it.
  */

/* 
  In case of error, these procedures report "*fName" as the file name,
  "R" as the line number, and "K" as the coordinate. 
*/

void ParseFreeVector(byte *L, long W, double *v, long N, char *fName, long R);
  /*
    Parses "L[0..W-1]" as "N" real numbers "v[0..N-1]" in free format.
    The values be terminated by blanks or end-of-buffer.
    Ignores any bytes of "L" beyond "v[N-1]". 
    
    Bombs out if finds less than "N" numbers, or a format error.
  */
  
void ParseFixedVector(byte *L, long W, int wd, double *v, long N, char *fName, long R);
  /*
    Parses "L[0..W-1]" as "N" real numbers "v[0..N-1]", each "wd" bytes wide.
    Ignores any bytes of "L" beyond "v[N-1]". 
    
    Bombs out if "W < N*wd", or sees a number format error.
  */
  
double ParseNum(byte *L, long W, char *fmt, char *fName, long R, long K);
  /* 
    Parses "L[0..W-1]" as a real number.  
    If the field contains only blanks, tabs, periods, or dashes,
    returns 0; else parses the field using the format "fmt".
    
    Bombs out in case of format error.
  */

#endif
