/* Compares two matrices. */
/* Last edited on 2006-03-19 12:48:23 by stolfi */

#define PROG_NAME "SPCompareMatrices"

#include <SPOptions.h>
#include <SPMatrix.h> 
#include <SPBasic.h>
#include <bool.h>

#include <stdio.h>

typedef struct Options
  { char *inNameA;       /* Filename prefix for first input matrix {A}. */
    char *inNameB;       /* Filename prefix for second input matrix {B}. */
    char *outName;       /* Filename prefix for comparison matrix {A-B}. */
    bool_t relative;     /* TRUE writes relative differences. */
  } Options;

Options GetOptions(int argn, char **argc);

SPMatrix ReadMatrix(char *name);
  /* Reads a matrix from the file {name} plus extension ".mat". */
  
void WriteMatrix(SPMatrix M, char *name);  
  /* Writes matrix {M} to the file {name} plus extension ".mat". */

int main(int argn, char **argc)
  { Options o = GetOptions(argn, argc);
    SPMatrix A = ReadMatrix(o.inNameA);
    SPMatrix B = ReadMatrix(o.inNameB);
    SPMatrix D;
    if (o.relative)
      { D = SPMatrix_RelDiff(A, B); }
    else
      { D = SPMatrix_Mix(1.0, A, -1.0, B); }
    WriteMatrix(D, o.outName);
    return 0;
  }

SPMatrix ReadMatrix(char *name)
  { FILE *rd = open_read(txtcat(name, ".mat"), TRUE);
    SPMatrix M = SPMatrix_Read(rd);
    fclose(rd);
    return M;
  }
  
void WriteMatrix(SPMatrix M, char *name)
  { FILE *wr = open_write(txtcat(name, ".mat"), TRUE);
    SPMatrix_Write(wr, M);
    fclose(wr);
  }


Options GetOptions(int argn, char **argc)
  {
    Options o;
    SPOptions_Parser_t *pp = SPOptions_NewParser(stderr, argn, argc);

    SPOptions_SetUsage(pp, 
      PROG_NAME " \\\n"
      "  -inNames NAMEA NAMEB -outName NAME\n"
    );

    SPOptions_GetKeyword(pp, "-inNames");                               
    o.inNameA = SPOptions_GetNext(pp);  
    o.inNameB = SPOptions_GetNext(pp);  

    SPOptions_GetKeyword(pp, "-outName");                               
    o.outName = SPOptions_GetNext(pp); 
    
    o.relative = SPOptions_TestKeyword(pp, "-relative"); 

    SPOptions_Finish(pp);
    return o;
  }
