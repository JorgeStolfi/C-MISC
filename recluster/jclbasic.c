/* Last edited on 2007-08-15 22:35:05 by stolfi */
/* See the authorship and copyright notice at the end of this file. */
/* Changes:

*/

#include <jclbasic.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <assert.h>
#include <jclerror.h>

void* Alloc(long size)
  { void *p = malloc(size);
    if (p == NULL) { Error("out of memory"); }
    return(p);
  }
  
char *CopyString(char *p)
  { char *q = (char*)(Alloc((strlen(p)+1)*sizeof(char)));
    char *t = q;
    while(((*t) = (*p)) != 0) { p++; t++; } 
    return q; 
  }

byte* CopyBytes(byte *L, long W)
  { byte *S = (char*)(Alloc(W*sizeof(byte)));
    long i;
    for (i=0; i<W; i++) S[i] = L[i];
    return(S);
  }

FILE *OpenRead(char *name)
  { if (strcmp(name, "-") == 0) 
      { return (stdin); }
    else
      { FILE *f = fopen(name, RMODE);
        if (f == NULL) { FileError("can't open input file ", name); }
        return(f);
      }
  }

FILE *OpenWrite(char *name)
  { if (strcmp(name, "-") == 0) 
      { return (stdout); }
    else
      { FILE *f = fopen(name, WMODE);
        if (f == NULL) { FileError("can't open output file ", name); }
        return(f);
      }
  }

bool_t ReadLine(FILE *f, byte **L, long *W, byte **B, long *A, char *name, long R)
  { int c;
    byte *curB = (*B);
    long curA = (*A);
    long curW = 0;
    c = fgetc(f);
    if (c == EOF) 
      { (*L) = NULL; (*W) = 0; fclose(f); return(TRUE); }
    while (c != '\n')
      { if (curW >= curA) 
          { if (curW >= MAX_LINE_LENGTH) 
              { LineError("line too long", name, R); }
            else
              { long newA = (curA == 0 ? 256 : 2*curA);
                byte *newL = (byte *)Alloc(newA*sizeof(byte));
                long i;
                for (i=0; i<curW; i++) newL[i] = curB[i];
                if (curA > 0) free(curB);
                curA = newA;
                curB = newL;
              }
          }
        curB[curW] = c;
        curW++;
        c = fgetc(f);
        if (c == EOF) { FileError("missing NL at end of file", name); }
      }
    (*B) = curB;
    (*A) = curA;
    (*L) = CopyBytes(curB, curW);
    (*W) = curW;
    return(FALSE);
  } 

void WriteLine(FILE *f, byte *L, long W)
  { long i;
    for (i=0; i<W; i++) putc(L[i], f);
    putc('\n', f);
  }

