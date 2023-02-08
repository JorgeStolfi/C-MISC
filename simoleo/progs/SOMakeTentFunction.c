/* Creates a procedural function and writes it to a file. */
/* Last edited on 2007-01-04 00:22:30 by stolfi */

/* This program creates an {SOProcFunction} object of a 
  specified subtype, and writes it to a file. */

#include <SOBasic.h>
#include <SOFunction.h>
#include <SOTentFunction.h>
#include <SOParams.h>

#include <stdio.h>
#include <stdlib.h>

typedef struct Options
  { int pDim;                /* Domain dimension. */
    dg_cell_index_t index;      /* Index of lower-numbered brick in tent's domain. */
    char *outName;           /* Prefix for output file name. */
  } Options; 

Options *GetOptions(int argn, char **argc);

int main(int argn, char **argc)
  { 
    Options *o = GetOptions(argn, argc);
    SOFunction *f = 
      (SOFunction *)SOTentFunction_FromHighBrick(o->pDim, o->index);
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
        
    PPUSAGE(pp, "SOPlotFunction \\\n");
    PPUSAGE(pp, "  -pDim NUM \\\n");
    PPUSAGE(pp, "  -index NUM \\\n");
    PPUSAGE(pp, "  -outName NAME \n");

    SOParams_GetKeyword(pp, "-pDim");                               
    o->pDim = SOParams_GetNextInt(pp, 1, MAX_PDIM);  
       
    SOParams_GetKeyword(pp, "-index");                               
    o->index = SOParams_GetNextInt(pp, 1, MAX_INT);  

    SOParams_GetKeyword(pp, "-outName");                               
    o->outName = SOParams_GetNext(pp);  

    SOParams_Finish(pp);
    return o;
  }
