/* Creates a procedural function and writes it to a file. */
/* Last edited on 2007-01-04 00:22:27 by stolfi */

/* This program creates an {SOProcFunction} object of a 
  specified subtype, and writes it to a file. */

#include <SOBasic.h>
#include <SOFunction.h>
#include <SOProcFunction.h>
#include <SOWaveFunction.h>
#include <SOTentFunction.h>
#include <SOShiftedFunction.h>
#include <SOParams.h>

#include <stdio.h>
#include <stdlib.h>

typedef struct Options
  { char *funcType;          /* Subtype of function. */
    int pDim;                /* Domain dimension. */
    int fDim;                /* Range dimension. */
    char *funcName;          /* Name of function, for {SOProcFunction}. */
    int fr[MAX_PDIM];        /* Frequency vector, for {SOWaveFunction}. */
    interval_t B[MAX_PDIM]; /* Reference box, for {SOShiftedFunction}. */
    char *outName;           /* Prefix for output file name. */
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
    int i; 
        
    PPUSAGE(pp, "SOPlotFunction \\\n");
    PPUSAGE(pp, "  -funcType TYPE \\\n");
    PPUSAGE(pp, "  -pDim NUM [ -fDim NUM ] \\\n");
    PPUSAGE(pp, "  [ -funcName NAME ] \\\n");
    PPUSAGE(pp, "  [ -freq NUM ... ] \\\n");
    PPUSAGE(pp, "  [ -box LO HI ... ] \\\n");
    PPUSAGE(pp, "  -outName NAME \n");

    SOParams_GetKeyword(pp, "-funcType");                               
    o->funcType = SOParams_GetNext(pp); 
    
    SOParams_GetKeyword(pp, "-pDim");                               
    o->pDim = SOParams_GetNextInt(pp, 1, MAX_PDIM);  
       
    if (SOParams_KeywordPresent(pp, "-fDim"))
      { o->fDim = SOParams_GetNextInt(pp, 1, MAX_FDIM); }
    else
      { o->fDim = 1; }
       
    if (strcmp(o->funcType, SOProcFunction_SubTypeId) == 0)
      { SOParams_GetKeyword(pp, "-funcName");
        o->funcName = SOParams_GetNext(pp);
      }
    else
      { o->funcName = ""; }
       
    if (strcmp(o->funcType, SOWaveFunction_SubTypeId) == 0)
      { SOParams_GetKeyword(pp, "-freq");
        for (i = 0; i < MAX_DIM; i++)
          { if (i < o->pDim)
              { o->fr[i] = SOParams_GetNextInt(pp, -100000, +100000); }
            else
              { o->fr[i] = 0; }
          }
      }
    else
      { for (i = 0; i < MAX_DIM; i++) { o->fr[i] = 0; } }
       
    if (strcmp(o->funcType, SOShiftedFunction_SubTypeId) == 0)
      { SOParams_GetKeyword(pp, "-box");
        for (i = 0; i < MAX_DIM; i++)
          { if (i < o->pDim)
              { o->LO(B[i]) = SOParams_GetNextDouble(pp, -INFTY, +INFTY);
                o->HI(B[i]) = SOParams_GetNextDouble(pp, -INFTY, +INFTY);
              }
            else
              { o->LO(B[i]) = 0.0; o->HI(B[i]) = 1.0; }
          }
      }
    else
      { for (i = 0; i < MAX_DIM; i++) { o->fr[i] = 0; } }
       
    if (strcmp(o->funcType, SOShiftedFunction_SubTypeId) == 0)
      { SOParams_GetKeyword(pp, "-box");
        for (i = 0; i < MAX_DIM; i++)
          { if (i < o->pDim)
              { o->LO(B[i]) = SOParams_GetNextDouble(pp, -INFTY, +INFTY);
                o->HI(B[i]) = SOParams_GetNextDouble(pp, -INFTY, +INFTY);
              }
            else
              { o->LO(B[i]) = 0.0; o->HI(B[i]) = 1.0; }
          }
      }
    else
      { for (i = 0; i < MAX_DIM; i++) { o->fr[i] = 0; } }
       
    SOParams_GetKeyword(pp, "-outName");                               
    o->outName = SOParams_GetNext(pp);  

    SOParams_Finish(pp);
    return o;
  }
