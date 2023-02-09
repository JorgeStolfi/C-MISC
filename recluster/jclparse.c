/* Last edited on 2007-08-16 00:05:42 by stolfi */
/* See the authorship and copyright notice at the end of this file. */
/* Changes:

*/

#include <jclbasic.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <jclerror.h>

/* INTERNAL PROTOTYPES */

/* PROCEDURES */

void ParseFreeVector(byte *L, long W, double *v, long N, char *fName, long R)
  { long i;
    byte *eb = L;
    long nb = W;
    char fmt[10];
           
    for (i=0; i<N; i++)
      { /* At this point "eb[0..nb-1]" is the unparsed part of "L" */
        byte *ef;
        long nf;
        /* Skip blanks and tabs: */
        while ((nb > 0) && (((*eb) == ' ') || ((*eb) == '\011'))) { eb++; nb--; }
        if (nb <= 0) { LineError("too few elements", fName, R); }
        /* Find next blank, tab, comma or end-of-buffer: */
        ef = eb; nf = nb;
        while ((nf > 0) && ((*ef) != ',') && ((*ef) != ' ') && ((*ef) != '\011'))
          { ef++; nf--; }
        sprintf(fmt, "%%%ldlf", nb-nf);
        /* Parse element: */
        v[i] = ParseNum(eb, nb-nf, fmt, R, i);
        /* Skip it: */
        eb = ef; nb = nf;
        /* skip comma, if any: */
        if ((nb > 0) && ((*eb) == ',')) { eb++; nb--; }
      }
  }
  
void ParseFixedVector(byte *L, long W, int wd, double *v, long N, char *fName, long R)
  { long i;
    byte *eb = L;
    long nb = W;
    char fmt[10];
    sprintf(fmt, "%%%ldlf", wd);
           
    for (i=0; i<N; i++)
      { /* At this point "eb[0..nb-1]" is the unparsed part of "L" */
        if (nb < wd) { LineError("line too short", fName, R); }
        v[i] = ParseNum(eb, wd, fmt, R, i);
        eb += wd; nb -= wd;
      }
  }
  
double ParseNum(byte *L, long W, char* fmt, long R, long K)
  { long i;
    for (i=0; i<W; i++)
      { if ((L[i] != ' ') && (L[i] != '.') && (L[i] != '-')) 
          { double val; int res;
            res = sscanf(L, fmt, &val); 
            if (res != 1) { ElementError("format error", R, K); }
            return(val);
          }
      }
    return 0;
  }
    
bool_t ReadItem(byte **L, long *N, long *A, long R)
  { int c;
    byte *curL = (*L);
    long curA = (*A);
    long curN = 0;
    c = fgetc(stdin);
    if (c == EOF) { (*N) = 0; fclose(stdin); return(TRUE); }
    while (c != '\n')
      { if (curN >= curA) 
          { if (curN >= MAX_LINE_LENGTH) 
              { LineError("line too long", "", R); }
            else
              { long newA = (curA == 0 ? 256 : 2*curA);
                byte *newL = (byte *)Alloc(newA*sizeof(byte));
                long i;
                for (i=0; i<curN; i++) newL[i] = curL[i];
                if (curA > 0) free(curL);
                curA = newA;
                curL = newL;
              }
          }
        curL[curN] = c;
        curN++;
        c = fgetc(stdin);
        if (c == EOF) { Error("missing NL at end of input file"); }
      }
    (*L) = curL;
    (*A) = curA;
    (*N) = curN;
    return(FALSE);
  } 

char* CopyItem(byte *L, long N)
  { byte *S = (char*)(Alloc(N*sizeof(byte)));
    long i;
    for (i=0; i<N; i++) S[i] = L[i];
    return(S);
  }

void WriteItem(byte *L, long N)
  { long i;
    for (i=0; i<N; i++) putchar(L[i]);
    putchar('\n');
  }


