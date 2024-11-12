/* Computes the Cholesky decomp of a sparse matrix. */
/* Last edited on 2006-03-19 12:48:17 by stolfi */

#define PROG_NAME "SPCholeskyTest"

#include <SPOptions.h>
#include <SPMatrix.h> 
#include <SPBasic.h>
#include <js.h>

#include <stdio.h>

typedef struct Options
  { char *inName;        /* Filename prefix for input matrix. */
    char *outName;       /* Filename prefix for output matrix. */
  } Options;

Options GetOptions(int argn, char **argc);

SPMatrix ReadMatrix(char *name);
  /* Reads a matrix from the file {name} plus extension ".mat". */
  
void WriteMatrix(SPMatrix M, char *name);  
  /* Writes matrix {M} to the file {name} plus extension ".mat". */

int main(int argn, char **argc)
  { Options o = GetOptions(argn, argc);
    SPMatrix M = ReadMatrix(o.inName);
    SPMatrix L;
    SPMatrix_Cholesky(M, 0.0, &L);
    WriteMatrix(L, o.outName);
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
      "  -inName NAME -outName NAME\n"
    );

    SPOptions_GetKeyword(pp, "-inName");                               
    o.inName = SPOptions_GetNext(pp);  

    SPOptions_GetKeyword(pp, "-outName");                               
    o.outName = SPOptions_GetNext(pp);  

    SPOptions_Finish(pp);
    return o;
  }
