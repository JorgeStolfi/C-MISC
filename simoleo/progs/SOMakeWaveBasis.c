/* SOMakeWaveBasis -- Builds a test basis of SOWaveFunctios. */
/* Last edited on 2006-03-19 12:45:41 by stolfi */

#include <SOFunction.h>
#include <SOWaveFunction.h>
#include <SOParams.h>
#include <SOBasic.h>
#include <SOeps.h>
#include <SO2DPlot.h>

#include <affirm.h>

#include <stdio.h>
#include <math.h>
#include <values.h>

typedef struct Options
  { 
    char *outName;    /* Name of basis output file (minus the ".bas" extension) */
    double maxFreq;   /* Value used to control the range of frequences relative to {N} */
  } Options;
  
/* INTERNAL PROTOTYPES */

Options *GetOptions(int argn, char **argc);  
FILE *OpenPlot(char *);

/* IMPLEMENTATIONS */

int main(int argn, char **argc)
  { int i, j, nel = 0;
    Options *o = GetOptions(argn, argc);
    int maxFI = (int)ceil(o->maxFreq);
    Basis bas = SOFunctionRef_vec_new(100);
    int k[2]; /* Frequency vector. */
    
    for(i = -maxFI; i <= +maxFI; i++)
      { for(j = -maxFI; j <= +maxFI; j++)
        { double f2 = i * i + j * j;
          k[0] = i; k[1] = j;
          if(sqrt(f2) <= o->maxFreq)
            { SOWaveFunction *wf = SOWaveFunction_FromFreq(2, k);
              nel++; 
	      SOFunctionRef_vec_expand(&bas, nel);
              bas.el[nel-1] = (SOFunction *)wf;
            }
        }
      } 
    SOFunctionRef_vec_trim(&bas, nel);

    FILE *wr = open_write(txtcat(o->outName, ".bas"), TRUE);
    SOFunction_WriteBasis(wr, bas);
    printf("\n\n #### SOMakeWaveBasis-> Wave Basis Size: %d \n", bas.nel);  
    fclose(wr);
    return 0;
  }

#define PPUSAGE SOParams_SetUsage

Options *GetOptions(int argn, char **argc)
  {
    Options *o = (Options *)notnull(malloc(sizeof(Options)), "no mem");
    SOParams_T *pp = SOParams_NewT(stderr, argn, argc);
    PPUSAGE(pp, "SOMakeWaveBasis \\\n");
    PPUSAGE(pp, "  -outName BASISNAME\n");
    PPUSAGE(pp, "  -maxFreq CONTROL_VAL\n");

    SOParams_GetKeyword(pp, "-outName");                               
    o->outName = SOParams_GetNext(pp);  

    SOParams_GetKeyword(pp, "-maxFreq");                               
    o->maxFreq = SOParams_GetNextDouble(pp, 0, INFTY); 

    SOParams_Finish(pp);
    return o;
  }
