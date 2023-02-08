/* SOCompareMatrices -- Compares two matrices. */
/* Last edited on 2007-01-04 00:41:00 by stolfi */

#include <SOParams.h>
#include <SOMatrix.h> 
#include <SOBasic.h>

#include <stdio.h>

typedef struct Options
  { char *inNameA;       /* Filename prefix for first input matrix {A}. */
    char *inNameB;       /* Filename prefix for second input matrix {B}. */
    char *outName;       /* Filename prefix for comparison matrix {A-B}. */
    bool_t relative;       /* TRUE writes relative differences. */
  } Options;

Options GetOptions(int argn, char **argc);

SOMatrix ReadMatrix(char *name);
  /* Reads a matrix from the file {name} plus extension ".mat". */
  
void WriteMatrix(SOMatrix M, char *name);  
  /* Writes matrix {M} to the file {name} plus extension ".mat". */

int main(int argn, char **argc)
  { Options o = GetOptions(argn, argc);
    SOMatrix A = ReadMatrix(o.inNameA);
    SOMatrix B = ReadMatrix(o.inNameB);
    SOMatrix D;
    if (o.relative)
      { D = SOMatrix_RelDiff(A, B); }
    else
      { D = SOMatrix_Mix(1.0, A, -1.0, B); }
    WriteMatrix(D, o.outName);
    return 0;
  }

SOMatrix ReadMatrix(char *name)
  { FILE *rd = open_read(txtcat(name, ".mat"), TRUE);
    SOMatrix M = SOMatrix_Read(rd);
    fclose(rd);
    return M;
  }
  
void WriteMatrix(SOMatrix M, char *name)
  { FILE *wr = open_write(txtcat(name, ".mat"), TRUE);
    SOMatrix_Write(wr, M);
    fclose(wr);
  }

#define PPUSAGE SOParams_SetUsage

Options GetOptions(int argn, char **argc)
  {
    Options o;
    SOParams_T *pp = SOParams_NewT(stderr, argn, argc);

    PPUSAGE(pp, "SOCompareMatrices \\\n");
    PPUSAGE(pp, "  -inNames NAMEA NAMEB -outName NAME\n");

    SOParams_GetKeyword(pp, "-inNames");                               
    o.inNameA = SOParams_GetNext(pp);  
    o.inNameB = SOParams_GetNext(pp);  

    SOParams_GetKeyword(pp, "-outName");                               
    o.outName = SOParams_GetNext(pp); 
    
    o.relative = SOParams_KeywordPresent(pp, "-relative"); 

    SOParams_Finish(pp);
    return o;
  }
