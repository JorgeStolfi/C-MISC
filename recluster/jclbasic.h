/* Last edited on 2008-07-14 22:21:15 by stolfi */
/* See the authorship and copyright notice at the end of this file. */

/* BASIC TYPES AND PROCEDURES */

#ifndef jclbasic_h
#define jclbasic_h

#include <stdio.h>
#include <bool.h>

typedef unsigned char byte;

/* LIMITS */

#define MAX_ITEMS           (2048)
/* Maximum items that can be processed */

#define MAX_VALUES          (2048)
/* Maximum coordinates per item */

#define MAX_WIDTH           (255)
/* Maximum width of a coordinate in input record */

#define MAX_LINE_LENGTH     (2048 + MAX_VALUES*MAX_WIDTH)
/* Maximum length of an input record */

#define MAX_REPEAT          (10000)
/* Maximum iterations the user can request */

#define NUM_COLORS          (216)
/* Number of colors to use in the pseudocolor map */

void* Alloc(long size);
  /*
    Like "malloc", but bombs out if runs out of memory. 
  */
  
char *CopyString(char *p);
  /* 
    Makes a newly allocated copy of a null-terminated string. 
  */

byte* CopyBytes(byte *L, long W);
  /*
    Makes a newly allocated copy of the "W"-byte array "L[0..W-1]".
  */

FILE *OpenRead(char *name);
  /* 
    Opens a file for reading; returns "stdin" if "name" is "-".
    Bombs in case of failure.
  */

FILE *OpenWrite(char *name);
  /* 
    Opens a file for writing; returns "stdout" if "name" is "-".
    Bombs in case of failure.
  */
    
bool_t ReadLine(FILE *f, byte **L, long *W, byte **B, long *A, char *name, long R);
  /* 
    If "f" was at end-of-file, returns TRUE, and sets "*L" to NULL,
    "*W" to 0.  Otherwise reads the next line of "f", up to but not
    including the newline character, into a new memory area; then
    returns its address in "*L", and its length in "*W".
    
    Uses a work buffer "*B" which is assumed to be "*A" bytes long;
    reallocates a new one, and updates "*B" and "*A", if necessary.
    Stops with error if the last line does not end with newline.
    In case of errors, reports "name" as the file name and "R"
    as the line number. 
  */
  
void WriteLine(FILE *f, byte *L, long W);
  /*
    Writes bytes "L[0..W-1]" to "f", followed by newline.
  */
  
#endif
