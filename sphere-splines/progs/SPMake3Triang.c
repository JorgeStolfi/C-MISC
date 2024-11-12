/* Creates a triangulation where vertices have degree 5 or 7. */
/* Last edited on 2006-03-19 12:49:10 by stolfi */

#define PROG_NAME "SPMake3Triang"

#include <SPOptions.h>
#include <r3.h>
#include <SPTriang.h>
#include <SPTriangExtra.h>
#include <SPQuad.h>
#include <SPBasic.h>
#include <vec.h>
#include <js.h>
#include <stdio.h>
#include <math.h>
#include <values.h>
#include <stdlib.h>
#include <string.h>

typedef struct Options /* Parsed command line options */
  { char *triName;   /* Input triangulation (usually {icosa}). */
    char *outName;   /* Output triangulation file (minus ".tri" extension). */
  } Options;

Options GetOptions(int argn, char **argc);
  
Triangulation *ReadTriangulation(char *name);
  /* Reads a triangulation from file {name} plus extension ".tri". */

void WriteTriangulation(Triangulation *tri, char *name);
  /* Writes triangulation {tri} to the file {name} plus extension ".bas". */

int main(int argn, char **argc)
  { Options o = GetOptions(argn, argc);
    Triangulation *triIn = ReadTriangulation(o.triName);
    Triangulation *triOt = SPTriang_Refine3Triang(triIn);
    WriteTriangulation(triOt, o.outName);
    fprintf(stderr, "done.\n");
    return 0;
  }

Triangulation *ReadTriangulation(char *name)
  { FILE *rd = open_read(txtcat(name, ".tri"), TRUE);
    Triangulation *tri = SPTriang_Read(rd, 1);
    fclose(rd);
    return tri;
  }
  
void WriteTriangulation(Triangulation *tri, char *name)
  { FILE *wr = open_write(txtcat(name, ".tri"), TRUE);
    SPTriang_Write(wr, tri);
    fclose(wr);
  }


Options GetOptions(int argn, char **argc)
  {
    Options o;
    SPOptions_Parser_t *pp = SPOptions_NewParser(stderr, argn, argc);
          
    SPOptions_SetUsage(pp, 
      PROG_NAME " \\\n"
      "  -triName NAME \\\n"
      "  -outName NAME\n"
    );
          
    SPOptions_GetKeyword(pp, "-triName");                               
    o.triName = SPOptions_GetNext(pp);  

    SPOptions_GetKeyword(pp, "-outName");                               
    o.outName = SPOptions_GetNext(pp);  

    SPOptions_Finish(pp);
    return o;
  }
