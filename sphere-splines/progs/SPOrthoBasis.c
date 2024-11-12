/* SPOrthoBasis -- Semi-orthonormalize a spherical function basis. */
/* Last edited on 2009-02-10 11:20:26 by stolfi */

#define PROG_NAME "SPOrthoBasis"

#define _GNU_SOURCE
#include <stdio.h>
#include <math.h>
#include <values.h>

#include <bool.h>

#include <SPTriang.h>
#include <SPFunction.h>
#include <SPSpline.h>
#include <SPSplineSemiOrthize.h>
#include <SPIntegral.h>
#include <SPOptions.h>
#include <SPBasic.h>

#define PROG_USAGE \
  PROG_NAME " \\\n" \
  "  -basisName NAME [ -triName TRINAME ] \\\n" \
  "  [ -fillStars ] \\\n" \
  "  [ -sortBSP | -sortLex | -sortGDist ] \\\n" \
  "  [ [ -eigen NUM ] [ -gram NUM ] \\\n" \
  "    -smpOrder NUM [ -precise] \\\n" \
  "    [ -checkOrtho TOL ] \\\n" \
  "  ] \\\n" \
  "  -outName BASISNAME\n"

/* This program reads a basis {B} and applies a partial
  orthonormalization, constrained so as not to expand the support of
  any piecewise function in the basis. */

typedef struct Options
  { char *basisName;  	/* Name of basis file (minus the ".bas" extension). */
    char *triName;    	/* Optional triangulation file (minus the ".tri" extn). */
    char *outName;    	/* Name of new basis file (minus the ".bas" extension) */
    bool_t fillStars; 	/* TRUE completes elems with >= 3 pieces to vertex stars. */
    bool_t weirdToo;  	/* TRUE to orthonize even weird elements. */
    bool_t sortBSP;   	/* TRUE uses {SPSpline_BSPSortBasis}. */
    bool_t sortLex;   	/* TRUE uses {SPSpline_LexSortBasis}. */
    bool_t sortGDist; 	/* TRUE uses {SPSpline_GDistSortBasis}. */
    int eigen;        	/* Iterations of eigen-analysis when orthonizing a class. */
    int gram;           /* Iterations of Gram-Schmidt when orthonizing a class. */
    bool_t precise;     /* TRUE computes dots products even of orthogonal elems. */ 
    double checkOrtho;  /* Tolerance for orthonization errors, or {INFINITY}. */
    int smpOrder;       /* Triangle sampling order for dot product integrals. */
  } Options;
  
Options GetOptions(int argn, char **argc);  
Basis ReadBasis(char *name);
  /* Reads a function basis from the file {name} plus extension ".bas". */

Triangulation *GetTriangulation(char *name, int smpOrder, Basis F);
  /* Gets triangulation from file {triName}, or gets
    the common triangulation from {F} if {triName} is empty. */

int main(int argn, char **argc)
  { Options o = GetOptions(argn, argc);
    SPIntegral_SetDefaultSamplingOrder(o.smpOrder);
    Basis F = ReadBasis(o.basisName);
    Triangulation *tri = GetTriangulation(o.triName, o.smpOrder, F);

    if (o.fillStars)
      { /* Enlarge support of elements with >=3 faces to vertex stars: */
        int i;
        fprintf(stderr, "filling-in vertex supports...\n");
        for (i = 0; i < F.ne; i++)
          { SPSpline *fpw = SPSpline_Cast(F.e[i]);
            if (fpw != NULL) { SPSpline_FillVertexSupport(fpw); }
          }
      }
    /* Sort basis: */
    if (o.sortBSP)
      { fprintf(stderr, "sorting basis...\n");
        SPSpline_BSPSortBasis(tri, F, 0, F.ne, TRUE);
      }
    else if (o.sortLex)
      { fprintf(stderr, "sorting basis...\n");
        SPSpline_LexSortBasis(tri, F, 0, F.ne, TRUE);
      }
    else if (o.sortGDist)
      { fprintf(stderr, "sorting basis...\n");
        SPSpline_GDistSortBasis(tri, F, 0, F.ne, TRUE);
      }

    /* Semi-orthonormalize the basis: */
    { 
      /* Metrics for semi-orthogonalization: */
      auto double EvalDot(SPFunction *f, SPFunction *g);
      auto double SLapDot(SPFunction *f, SPFunction *g);
      
      double EvalDot(SPFunction *f, SPFunction *g)
        { return SPFunction_Dot(f, NoFMap, g, NoFMap, NULL, tri, FALSE); }
        
      double SLapDot(SPFunction *f, SPFunction *g)
        { return SPFunction_SLapDot(f, g, NULL, tri, FALSE); }

      bool_t orthogonalize = (o.eigen > 0) || (o.gram > 0);
      
      if (orthogonalize)
        { /* Apply semi-orthogonalization procedure: */
          fprintf(stderr, "orthogonalizing...\n");
          SPSplineSemiOrthize_Basis
            ( F, tri, o.weirdToo, EvalDot, o.precise, o.gram, o.eigen, SLapDot, TRUE );
        }
      
      /* Write basis: */
      { FILE *wr = open_write(txtcat(o.outName, ".bas"), TRUE);
        SPFunction_WriteBasis(wr, F);
        fclose(wr);
      }
      
      if (orthogonalize && (o.checkOrtho < INFINITY))
        { /* Check whether semi-orthogonalization worked: */
          double maxNormError, maxOrthoError;
          fprintf(stderr, "checking orthogonalization...\n");
          SPSplineSemiOrthize_Check
            ( F, 0, F.ne, 
              tri, FALSE, EvalDot, 
              &maxNormError, &maxOrthoError, 
              o.checkOrtho, TRUE
            );
          if ((maxNormError > o.checkOrtho) || (maxOrthoError > o.checkOrtho))
      	  { fprintf(stderr, "semi-orthogonalization failed\n");
      	    return 1;
      	  }
        }
    }
    return 0;
  }
  
Basis ReadBasis(char *name)
  { FILE *rd = open_read(txtcat(name, ".bas"), TRUE);
    Basis F = SPFunction_ReadBasis(rd);
    fclose(rd);
    return F;
  }
  
Triangulation *GetTriangulation(char *name, int smpOrder, Basis F)
  { if ((name != NULL) && (*name != '\000'))
      { FILE *rd = open_read(txtcat(name, ".tri"), TRUE);
        Triangulation *tri = SPTriang_Read(rd, smpOrder);
        fclose(rd);
        return tri;
      }
    else
      { return SPSpline_BasisTriangulation(F); }
  }

Options GetOptions(int argn, char **argc)
  {
    Options o;
    SPOptions_Parser_t *pp = SPOptions_NewParser(stderr, argn, argc);
    SPOptions_SetUsage(pp, PROG_USAGE);

    SPOptions_GetKeyword(pp, "-basisName");                               
    o.basisName = SPOptions_GetNext(pp);

    if (SPOptions_TestKeyword(pp, "-triName"))
      { o.triName = SPOptions_GetNext(pp); }
    else
      { o.triName = ""; }  

    o.fillStars = SPOptions_TestKeyword(pp, "-fillStars");

    o.sortBSP = SPOptions_TestKeyword(pp, "-sortBSP");
    o.sortLex = SPOptions_TestKeyword(pp, "-sortLex");
    o.sortGDist = SPOptions_TestKeyword(pp, "-sortGDist");
    if (o.sortBSP + o.sortLex + o.sortGDist > 1)
      { SPOptions_Error(pp, "cannot ask more than one sort order"); } 

    o.weirdToo = SPOptions_TestKeyword(pp, "-weirdToo");

    if (SPOptions_TestKeyword(pp, "-gram"))
      { o.gram = SPOptions_GetNextInt(pp, 0, INT_MAX); }
    else
      { o.gram = 0; }  

    if (SPOptions_TestKeyword(pp, "-eigen"))
      { o.eigen = SPOptions_GetNextInt(pp, 0, INT_MAX); }
    else
      { o.eigen = 0; }  
    
    if (SPOptions_TestKeyword(pp, "-smpOrder"))
      { o.smpOrder = SPOptions_GetNextInt(pp, 1, INT_MAX); }
    else
      { if ((o.eigen > 0) || (o.gram > 0)) 
          { SPOptions_Error(pp, "\"-smpOrder\" required when orthogonalizing"); }
	o.smpOrder = 10;
      }

    o.precise = SPOptions_TestKeyword(pp, "-precise");

    if (SPOptions_TestKeyword(pp, "-checkOrtho"))
      { o.checkOrtho = SPOptions_GetNextDouble(pp, 0, 1.0); }
    else
      { o.checkOrtho = INFINITY; }

    SPOptions_GetKeyword(pp, "-outName");                               
    o.outName = SPOptions_GetNext(pp);  

    SPOptions_Finish(pp);
    return o;
  }
