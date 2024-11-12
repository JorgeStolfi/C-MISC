/* Computes approximation error maps. */
/* Last edited on 2006-03-19 12:48:32 by stolfi */

#define PROG_NAME "SPComputeErrorMap"

/* Computes the worst-case approximation error map {e(p)} for the
  gauge space {\CF_{g}} of harmonic functions of a given maximum
  degree, and an arbitrary approximation space {\CA}. */

#include <SPErrorMap.h>

#include <SPSys.h>
#include <SPHarmonic.h>
#include <SPApprox.h>
#include <SPBasic.h>
#include <SPOptions.h>
#include <SPFunction.h>
#include <SPIntegral.h>
#include <vec.h>
#include <affirm.h>
#include <bool.h>
#include <stdio.h>
#include <math.h>
#include <values.h>

typedef struct Options
  { char *basisName;   /* Name of basis file for {\CA} (minus ".bas"). */
    char *matName;     /* Prefix for names of normal matrix files of {\CA} . */
    char *outName;     /* Prefix for output file names. */
    int degree;        /* Degree of gauge function space {\CF}. */
    int smpOrder;      /* Sampling order for integration. */
    /* Options to skip the subspace {\CF \inter \CA}: */
    int exclude;       /* Max degree of harmonics to exclude, or -1. */
    bool_t homogeneous;  /* TRUE excludes only even or only odd degrees. */
  } Options;

Options GetOptions(int argn, char **argc);

SPMatrix *GetMatrix(char *name);
  /* Reads a matrix from file {matName} plus extension ".mat". 
     Returns NULL if there is no such file. */

Basis ReadBasis(char *name);
  /* Reads a function basis from the file {name} plus extension ".bas". */
  
Basis MakeGaugeSpaceOrthoBasis(int degree, int exclude, bool_t homogeneous);
  /* Computes an orthonormal basis for the gauge space.
    The gauge space consists of all spherical harmonic 
    basis elements up to degree "degree", minus the subspace
    that is contained in the space of polynomial splines
    of degree "exclude" and homogeneity "homogeneous". */  

int main(int argn, char **argc)
  { 
    auto double dot(SPFunction *f, SPFunction *g);
    
    double dot(SPFunction *f, SPFunction *g)
      { return SPFunction_Dot(f, NoFMap, g, NoFMap, NULL, NULL, FALSE); }

    Options o = GetOptions(argn, argc);
    SPIntegral_SetDefaultSamplingOrder(o.smpOrder);
    Basis A = ReadBasis(o.basisName);
    /* Get normal matrix for {A}, preferably its Choleski factor: */
    SPMatrix *L = GetMatrix(txtcat(o.matName, "-ev-choL"));
    SPMatrix *B = (L == NULL ? GetMatrix(txtcat(o.matName, "-ev")) : NULL);
    Basis H = MakeGaugeSpaceOrthoBasis(o.degree, o.exclude, o.homogeneous);
    Basis G = SPFunction_GenProjectBasisOnBasis(H, A, dot, TRUE, B, L);
    SPFunction *mu = (SPFunction *)SPErrorMap_FromBases(H, G);
    FILE *wr = open_write(txtcat(o.outName, ".sfn"), TRUE);
    double muMax = 0.0, muAvg;
    mu->m->write(mu, wr);
    fclose(wr);
    SPApprox_PrintMaxAvgValues(mu, NULL, &muMax, &muAvg);
     return 0;
  }

SPMatrix *GetMatrix(char *name)
  { char *fileName = txtcat(name, ".mat");
    FILE *rd = fopen(fileName, "r");
    if (rd == NULL) 
      { fprintf(stderr, "file %s not found, skipped\n", fileName);
        return NULL;
      }
    else
      { SPMatrix *M = (SPMatrix *)notnull(malloc(sizeof(SPMatrix)), "no mem");
        fprintf(stderr, "reading %s\n", fileName);
        (*M) = SPMatrix_Read(rd);
        fclose(rd);
        return M;
      }
  }

Basis ReadBasis(char *name)
  { FILE *rd = open_read(txtcat(name, ".bas"), TRUE);
    Basis F = SPFunction_ReadBasis(rd);
    fprintf(stderr, "dimension of approx space = %d\n", F.ne);
    return F;
  }
  
Basis MakeGaugeSpaceOrthoBasis(int degree, int exclude, bool_t homogeneous)
  { int dimH = (degree + 1)*(degree+1); /* Dimension of gauge space. */
    int nH1 = 0;
    Basis H1 = Basis_new(dimH);
    int d, m;

    /* Now create the basis: */
    fprintf(stderr, "building gauge basis (harmonics, deg 0..%d)...\n", degree);
    fprintf(stderr, "dimension of gauge space = %d\n", dimH);
    if (exclude > 0)
      { char *ex = (homogeneous ? (exclude % 2 == 0 ? "even" : "odd") : "all");
        fprintf(stderr, "excluding %s harmonics up to degree %d\n", ex, exclude);
      }
    for (d = 0; d <= degree; d++)
      { if((d > exclude) || (homogeneous && ((exclude - d) % 2 == 1)))
          { for (m = -d;  m <= +d;  m++) 
              { Basis_expand(&H1,nH1);
                H1.e[nH1] = (SPFunction *)SPHarmonic_FromTerm(d, m);
                nH1++;
              }
          }
      }
    fprintf(stderr, "base size after exclusion = %d\n", nH1);
    Basis_trim(&H1, nH1);
    return H1;
  }


Options GetOptions(int argn, char **argc)
  { Options o;
    SPOptions_Parser_t *pp = SPOptions_NewParser(stderr, argn, argc);
    SPOptions_SetUsage(pp, 
      PROG_NAME " \\\n"
      "  -basisName NAME -matName NAME \\\n"
      "  -outName NAME \\\n"
      "  [ -exclude NUM [ -homogenous | -general ] ] \\\n"
      "  -degree NUM \\\n"
      "  -smpOrder NUM "
    );             

    SPOptions_GetKeyword(pp, "-basisName");                               
    o.basisName = SPOptions_GetNext(pp);  
       
    SPOptions_GetKeyword(pp, "-matName");                               
    o.matName = SPOptions_GetNext(pp);  
        
    SPOptions_GetKeyword(pp, "-outName");                               
    o.outName = SPOptions_GetNext(pp);  
                                                 
    SPOptions_GetKeyword(pp, "-degree");
    o.degree = SPOptions_GetNextInt(pp, 0, INT_MAX);
    
    if (SPOptions_TestKeyword(pp, "-exclude"))
      { o.exclude = SPOptions_GetNextInt(pp, 0, INT_MAX); }
    else
      { o.exclude = -1; }
        
    if (SPOptions_TestKeyword(pp, "-homogeneous"))
      { o.homogeneous = TRUE; }
    else if (SPOptions_TestKeywordNext(pp, "-general"))
      { o.homogeneous = FALSE; }
    else
      { if (o.exclude >= 0)
          { SPOptions_Error(pp, "must say \"-homogeneous\" or \"-general\""); }
      }

    SPOptions_GetKeyword(pp, "-smpOrder");
    o.smpOrder = SPOptions_GetNextInt(pp, 1, INT_MAX);
        
    SPOptions_Finish(pp);
    return o;
  }
