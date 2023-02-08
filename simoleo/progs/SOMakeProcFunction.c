/* Creates a procedural function and writes it to a file. */
/* Last edited on 2006-03-19 12:45:09 by stolfi */

/* This program creates an {SOProcFunction} object of a 
  specified subtype, and writes it to a file. */

#include <SOBasic.h>
#include <SOFunction.h>
#include <SOProcFunction.h>
#include <SOParams.h>

#include <affirm.h>

#include <stdio.h>
#include <stdlib.h>

typedef struct Options
  { char *funcName;   /* Filename of function. */
    int pDim;         /* Domain dimension. */
    int fDim;         /* Range dimension. */
    char *outName;    /* Prefix for output file name. */
  } Options;

Options *GetOptions(int argn, char **argc);

int main(int argn, char **argc)
  { 
    Options *o = GetOptions(argn, argc);
    SOFunction *f = 
      (SOFunction *)SOProcFunction_FromName(o->funcName, o->pDim, o->fDim);
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

    PPUSAGE(pp, "SOMakeProcFunction \\\n");
    PPUSAGE(pp, "  -funcName NAME \\\n");
    PPUSAGE(pp, "  -pDim NUM [ -fDim NUM ] \\\n");
    PPUSAGE(pp, "  -outName NAME \n");

    SOParams_GetKeyword(pp, "-funcName");                               
    o->funcName = SOParams_GetNext(pp);  
       
    SOParams_GetKeyword(pp, "-pDim");                               
    o->pDim = SOParams_GetNextInt(pp, 1, MAX_PDIM);  
       
    if (SOParams_KeywordPresent(pp, "-fDim"))
      { o->fDim = SOParams_GetNextInt(pp, 1, MAX_FDIM); }
    else
      { o->fDim = 1; }
       
    SOParams_GetKeyword(pp, "-outName");                               
    o->outName = SOParams_GetNext(pp);  

    SOParams_Finish(pp);
    return o;
  }
