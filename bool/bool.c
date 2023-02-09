/* Last edited on 2013-10-02 02:46:53 by stolfilocal */ 
/****************************************************************************/
/* (C) Copyright 1992 Universidade Estadual de Campinas (UNICAMP)           */
/*                    Campinas, SP, Brazil                                  */
/*                                                                          */
/* Authors:                                                                 */
/*                                                                          */
/*   Jorge Stolfi        - CS Dept, UNICAMP <stolfi@dcc.unicamp.br>         */
/*                                                                          */
/* This file can be freely distributed, modified, and used for any          */
/*   non-commercial purpose, provided that this copyright and authorship    */
/*   notice be included in any copy or derived version of this file.        */
/*                                                                          */
/* DISCLAIMER: This software is offered ``as is'', without any guarantee    */
/*   as to fitness for any particular purpose.  Neither the copyright       */
/*   holder nor the authors or their employers can be held responsible for  */
/*   any damages that may result from its use.                              */
/****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

typedef unsigned char byte;
typedef int bool;
#define TRUE 1
#define FALSE 0

#ifdef MSDOS

#include <process.h> 
#define RMODE "rb"
#define WMODE "wb"

#else

#define RMODE "r"
#define WMODE "w"

#endif

/* INTERNAL PROTOTYPES */

int main(int argc, char *argv[]);

FILE *OpenFile(char *name);
 /*
   Opens file "name" for reading; bombs in case of failure.
   If "name" is "-", returns stdin.
 */
     
bool ReadLine(FILE *rd, char *name, byte **L, long *N, long *A, long *R);
  /* 
    Returns TRUE if "rd" (whose filename is "name") was at end-of-file.
    Otherwise reads the next line of "rd" into "*L", and stores its length 
    into "*N", and increments "*R".  Assumes the buffer "*L" is "*A" 
    bytes long; reallocates a new one if necessary. 
    Stops with error if the last line does not end with newline.
  */

void WriteLine(byte L[], long N, long *R);
  /*
    Writes bytes "L[0..N-1]" to stdout, followed by newline.
    Increments the line counter "(*R)".
  */
 
int CompareLines(byte *aL, long aN, byte *bL, long bN, bool signedBytes);
  /*
    Lexicographic comparison of lines "aL[0..aN-1]" and "bL[0..bN-1]".
    Returns -1 for "<", 0 for "=", +1 for ">".
    Uses signed byte comparisons iff "signedBytes" is TRUE.
  */
  
void GetOptions(
    int argc, char *argv[], 
    char *inName[2],    /* Input file names */
    bool *verbose,      /* TRUE to print line counts */
    bool *signedBytes,  /* TRUE to use signed byte comparisons */
    bool *only0,        /* TRUE to write lines that are only in file 0 */
    bool *both,         /* TRUE to write lines that are in both files */
    bool *only1         /* TRUE to write lines that are only in file 1 */
  );
  /* 
    Parses the command line options.
  */

void ErrorParams(char *msg);
void ErrorCannotOpen(char *name);
void ErrorOutOfOrder(char *name, long R);
void ErrorNoFinalEOL(char *name);
  /*
    Print error message and exit
  */

/* PROCEDURES */

char *progName;

int main(int argc, char *argv[])
{
  /* Command line options: */
  bool verbose;
  bool signedBytes;
  char *inName[2];     /* Input file names */

  /* Internal variables: */
  FILE *inFile[2];   /* Input files */
  byte *thisL[2];    /* Last line of each file */
  long thisN[2];     /* Used size of "thisL[i]". */
  long thisA[2];     /* Allocated size of "thisL[i]" */
  byte *prevL[2];    /* next-to-last line of each file */
  long prevN[2];     /* Used size of "prev". */
  long prevA[2];     /* Allocated size of "prevL[i]" */
  bool done[2];      /* TRUE if file[i] is exhausted. */
  long nIn[2];       /* Count input lines by file */
  bool only0;        /* TRUE --> output records in file0 - file1 */       
  bool both;         /* TRUE --> output records in file0 * file1 */       
  bool only1;        /* TRUE --> output records in file1 - file0 */       
  long nOut;         /* Counts output lines */
  int i;
  
  /* Collect options: */
  progName = argv[0];
  GetOptions(
    argc, argv,
    inName,
    &verbose,
    &signedBytes,
    &only0,
    &both,
    &only1
  );
  
  /* Open files, allocate line buffers, and read first lines */
  for (i=0; i<2; i++)
    { done[i] = FALSE;
      nIn[i] = 0;
      thisA[i] = 250;
      thisL[i] = (byte *)malloc(thisA[i]*sizeof(byte));
      prevA[i] = 250;
      prevL[i] = (byte *)malloc(prevA[i]*sizeof(byte));
      inFile[i] = OpenFile(inName[i]);
      done[i] = ReadLine(
        inFile[i], inName[i], 
        &thisL[i], &thisN[i], &thisA[i], 
        &nIn[i]
      );
    }
  nOut = 0;

  /* Merge lines and output requested ones: */
  while(! (done[0] & done[1]))
    { bool cmp;
      /* Compare next lines, one from each file: */
      if(done[0]) 
        cmp = 1;
      else if (done[1]) 
        cmp = -1;
      else 
        cmp = CompareLines(
          thisL[0], thisN[0], 
          thisL[1], thisN[1], 
          signedBytes
        );
      /* Output desired lines: */
      switch (cmp) {
        case -1: if (only0) WriteLine(thisL[0], thisN[0], &nOut); break;
        case 00: if (both)  WriteLine(thisL[0], thisN[0], &nOut); break;
        case  1: if (only1) WriteLine(thisL[1], thisN[1], &nOut); break;
      }
      /* Read next lines: */
      for (i=0; i<2; i++)
        { if (((cmp <= 0) && (i == 0)) || ((cmp >= 0) && (i == 1)))
            { byte *tl; long tn;
              tl = prevL[i]; prevL[i] = thisL[i]; thisL[i] = tl; 
              tn = prevA[i]; prevA[i] = thisA[i]; thisA[i] = tn;
              prevN[i] = thisN[i];
              assert(!done[i]); 
              done[i] = ReadLine(
                inFile[i], inName[i], 
                &thisL[i], &thisN[i], &thisA[i],
                &nIn[i]
              );
              if (!done[i]) 
                { int ord;
                  ord = CompareLines(
                    prevL[i], prevN[i], 
                    thisL[i], thisN[i], 
                    signedBytes
                  );
                  if (ord >= 0) { ErrorOutOfOrder(inName[i], nIn[i]); }
                }
            }
        }
    }

  /* Close files and print report: */
  for (i=0; i<2; i++)
    { if (inFile[i] != stdin) { fclose(inFile[i]); }
      if (verbose)
        { fprintf(stderr, "file %d: %8ld lines\n", i+1, nIn[i]); }
    }
  fclose(stdout);
  if (verbose)
    { fprintf(stderr, "output: %8ld lines\n", nOut); }
  return(0);
}
  
FILE *OpenFile(char *name)
  { FILE *f;
    if (strcmp(name, "-") == 0) return(stdin);
    f = fopen(name, RMODE);
    if (f == NULL) { ErrorCannotOpen(name); }
    return(f);
  }

bool ReadLine(FILE *rd, char *name, byte **L, long *N, long *A, long *R)
  { int c;
    byte *curL = (*L);
    long curA = (*A);
    long curN = 0;
    c = fgetc(rd);
    if (c == EOF) { (*N) = 0; return(TRUE); }
    while (c != '\n')
      { if (curN >= curA) 
          { long newA = 2*curA;
            byte *newL = (byte *)malloc(newA*sizeof(byte));
            long i;
            for (i=0; i<curN; i++) newL[i] = curL[i];
            free(curL);
            curA = newA;
            curL = newL;
          }
        curL[curN] = c;
        curN++;
        c = fgetc(rd);
        if (c == EOF) { ErrorNoFinalEOL(name); }
      }
    (*L) = curL;
    (*A) = curA;
    (*N) = curN;
    (*R)++;
    return(FALSE);
  } 

void WriteLine(byte L[], long N, long *R)
  { long i;
    for (i=0; i<N; i++) putchar(L[i]);
    putchar('\n');
    (*R)++;
  }

int CompareLines(byte *aL, long aN, byte *bL, long bN, bool signedBytes)
  { long n = ( aN < bN ? aN : bN);
    long i;
    if (signedBytes)
      { for(i=0; i<n; i++)
          { signed char ac = (signed char)(aL[i]);
            signed char bc = (signed char)(bL[i]);
            if (ac < bc) return(-1);
            if (ac > bc) return( 1);
          }
      }
    else
      { for(i=0; i<n; i++)
          { byte ac = aL[i];
            byte bc = bL[i];
            if (ac < bc) return(-1);
            if (ac > bc) return( 1);
          }
      }
    if (aN < bN) return(-1);
    if (aN > bN) return( 1);
    return(0);
  }

void GetOptions(
    int argc, char *argv[], 
    char *inName[2], 
    bool *verbose, 
    bool *signedBytes,
    bool *only0,
    bool *both,
    bool *only1
  )
  { int i = 1;
    char *op;
    /* Parse option switches: */
    (*signedBytes) = FALSE;
    (*verbose) = FALSE;
    while ((i<argc) && (argv[i][0] == '-'))
      { if ((strcmp(argv[i], "-s")==0) || (strcmp(argv[i], "--signedBytes")==0))
          { (*signedBytes) = TRUE; }
        else if ((strcmp(argv[i], "-v")==0) || (strcmp(argv[i], "--verbose")==0))
          { (*verbose) = TRUE; }
        i++;
      }
    /* Parse operation code: */
    if ((argc - i) != 3) { ErrorParams("wrong number of arguments"); }
    op = argv[i];
    if (strcmp(op, "1-2")==0)
      { (*only0) = TRUE;  (*both) = FALSE; (*only1) = FALSE; }
    else if (strcmp(op, "2-1")==0)
      { (*only0) = FALSE; (*both) = FALSE; (*only1) = TRUE;  }
    else if ((strcmp(op, "1+2")==0) || (strcmp(op, "2+1")==0))
      { (*only0) = TRUE;  (*both) = TRUE;  (*only1) = TRUE;  }
    else if ((strcmp(op, "1.2")==0) || (strcmp(op, "2.1")==0))
      { (*only0) = FALSE; (*both) = TRUE;  (*only1) = FALSE; }
    else if ((strcmp(op, "1#2")==0) || (strcmp(op, "2#1")==0))
      { (*only0) = TRUE;  (*both) = FALSE; (*only1) = TRUE;  }
    else
      { ErrorParams("bad operation"); }
    /* Parse filenames: */
    inName[0] = argv[i+1];
    inName[1] = argv[i+2];
    if((strcmp(inName[0], "-")==0) && (strcmp(inName[1], "-")==0))
      { ErrorParams("standard input can be used only once"); }
  }
  
void ErrorParams(char *msg)
  { fprintf(stderr, "%s: %s\n", progName, msg);
    fprintf(stderr, "usage: %s \\\n", progName);
    fprintf(stderr, "  [ -v | --verbose ] \\\n");
    fprintf(stderr, "  [ -s | --signedBytes] \\\n");
    fprintf(stderr, "  [ 1-2 | 2-1 | 1+2 | 1.2 | 1#2 ] \\\n");
    fprintf(stderr, "  [ file1 | - ] [ file2 | - ] \\\n");
    fprintf(stderr, "  > outfile\n");
    fflush(stderr);
    exit(1);
  }
  
void ErrorCannotOpen(char *name)
  {
    fprintf(stderr, "%s: cannot open file \"%s\"\n", progName, name);
    fflush(stderr);
    exit(1);
  }

void ErrorOutOfOrder(char *name, long R)
  {
    fprintf(stderr, 
      "%s: lines %ld and %ld of file \"%s\" are out of order\n", 
      progName, R-1, R, name
    );
    fflush(stderr);
    exit(1);
  }

void ErrorNoFinalEOL(char *name)
  {
    fprintf(stderr, 
      "%s: missing final newline in file \"%s\"\n", 
      progName, name
    );
    fflush(stderr);
    exit(1);
  }


