/* SPMakeTimeOnlyProcFunction -- Creates an {SPTimeOnlyProcFunction} file. */
/* Last edited on 2006-03-19 12:49:27 by stolfi */

#define PROG_NAME "SPMakeTimeOnlyProcFunction"

#include <SPOptions.h>
#include <SPTimeOnlyFunction.h> 
#include <SPTimeOnlyProcFunction.h> 
#include <SPBasic.h>
#include <js.h>
#include <stdlib.h>
#include <limits.h>
#include <string.h>
#include <stdio.h>

typedef struct Options
  { char *function;    /* Name of function (type tag, not *file*name!). */
    char *outName;     /* Prefix for output file name (minus ".tof"). */
  } Options;

Options GetOptions(int argn, char **argc);
  
void WriteFunction(SPTimeOnlyFunction *f, char *name);
  /* Writes function {f} to the file {name} plus extension ".tof". */
 
int main(int argn, char **argc)
  { Options o = GetOptions(argn, argc);
    SPTimeOnlyFunction *f = 
      (SPTimeOnlyFunction *)SPTimeOnlyProcFunction_FromName(o.function);
    WriteFunction(f, o.outName);
    return 0;
  }
  
void WriteFunction(SPTimeOnlyFunction *f, char *name)
  { FILE *wr = open_write(txtcat(name, ".tof"), TRUE);
    f->m->write(f, wr);
    fclose(wr);
  }


Options GetOptions(int argn, char **argc)
  {
    Options o;
    SPOptions_Parser_t *pp = SPOptions_NewParser(stderr, argn, argc);
          
    SPOptions_SetUsage(pp, 
      PROG_NAME " \\\n"
      "  -function NAME \\\n"
      "  -outName NAME \n"
    );

    SPOptions_GetKeyword(pp, "-function");                               
    o.function = SPOptions_GetNext(pp);  

    SPOptions_GetKeyword(pp, "-outName");                               
    o.outName = SPOptions_GetNext(pp);  
    
    SPOptions_Finish(pp);
    return o;
  }

                        
