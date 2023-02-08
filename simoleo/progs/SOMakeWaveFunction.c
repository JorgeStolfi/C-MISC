/* Creates a procedural function and writes it to a file. */
/* Last edited on 2007-01-04 00:22:33 by stolfi */

/* This program creates an {SOProcFunction} object of a 
  specified subtype, and writes it to a file. */

#include <SOBasic.h>
#include <SOFunction.h>
#include <SOWaveFunction.h>
#include <SOParams.h>

#include <stdio.h>
#include <stdlib.h>

typedef struct Options
  { int pDim;                /* Domain dimension. */
    int fDim;                /* Range dimension. */
    int fr[MAX_PDIM];        /* Frequency vector, for {SOWaveFunction}. */
    char *outName;           /* Prefix for output file name. */
  } Options; 

Options *GetOptions(int argn, char **argc);

int main(int argn, char **argc)
  { 
    Options *o = GetOptions(argn, argc);
    SOFunction *f = 
      (SOFunction *)SOWaveFunction_FromFreq(o->pDim, o->fDim, o->fr);
    FILE *wr = open_write(txtcat(o->outName, ".fun"), TRUE);
    f->m->write(f, wr);
    fclose(wr);
    return 0;
  }

#define PPUSAGE SOParams_SetUsage

Options *GetOptions(int argn, char **argc)
  {
    Options *o = (Options *)notnull(malloc(sizeof(Options)), "no mem");
    SOParams_T *pp = SOParams_NewT(stderr, argn, argc);
    int i; 
        
    PPUSAGE(pp, "SOMakeWaveFunction \\\n");
    PPUSAGE(pp, "  -pDim NUM [ -fDim NUM ] \\\n");
    PPUSAGE(pp, "  -freq NUM ... \\\n");
    PPUSAGE(pp, "  -outName NAME \n");

    SOParams_GetKeyword(pp, "-pDim");                               
    o->pDim = SOParams_GetNextInt(pp, 1, MAX_PDIM);  
       
    if (SOParams_KeywordPresent(pp, "-fDim"))
      { o->fDim = SOParams_GetNextInt(pp, 1, MAX_FDIM); }
    else
      { o->fDim = 1; }
       
    SOParams_GetKeyword(pp, "-freq");
    for (i = 0; i < MAX_DIM; i++)
      { if (i < o->pDim)
          { o->fr[i] = SOParams_GetNextInt(pp, -100000, +100000); }
        else
          { o->fr[i] = 0; }
      }
       
    SOParams_GetKeyword(pp, "-outName");                               
    o->outName = SOParams_GetNext(pp);  

    SOParams_Finish(pp);
    return o;
  }
