/* SPMakeHarmonicBasis -- Creates bases of spherical harmonics. */
/* Last edited on 2006-03-19 12:49:17 by stolfi */

#define PROG_NAME "SPMakeHarmonicBasis" 

/* Creates a basis consisting of spherical harmonics up to a
  given degree. */

#include <SPBasic.h>
#include <SPOptions.h>
#include <SPIntegral.h>
#include <SPFunction.h>
#include <SPHarmonic.h>

#include <vec.h>
#include <bool.h>
#include <affirm.h>

#include <stdio.h>
#include <math.h>
#include <values.h>

#define PROG_USAGE \
  PROG_NAME " \\\n" \
  "  -degree NUM [ -checkOrtho ] [ -verbose ] \\\n" \
  "  -smpOrder NUM \\\n" \
  "  -outName BASISNAME\n"

typedef struct Options
  { int degree;        /* Maximum degree of harmonics. */
    char *outName;     /* Prefix for output file names. */
    int smpOrder;      /* Sampling order for dot product integrals. */
    bool_t checkOrtho; /* TRUE checks orthonormality of elements. */
    bool_t verbose;    /* TRUE mumbles while working. */
  } Options;

Options GetOptions(int argn, char **argc);

int main(int argn, char **argc)
  { Options o = GetOptions(argn, argc);
    SPIntegral_SetDefaultSamplingOrder(o.smpOrder);
    Basis F = SPHarmonic_MakeBasis(o.degree);
    FILE *wr = open_write(txtcat(o.outName, ".bas"), TRUE);
    SPFunction_WriteBasis(wr, F);
    fclose(wr);
    if (o.checkOrtho) 
      { int nBad = SPFunction_StdCheckOrthonization
          ( F, NULL, NULL, 1000, 1.0e-11, 1000, o.verbose );
        affirm(nBad == 0, "orthogonalization test failed"); 
      }
    return 0;
  }

Options GetOptions(int argn, char **argc)
  {
    Options o;
    SPOptions_Parser_t *pp = SPOptions_NewParser(stderr, argn, argc);
    SPOptions_SetUsage(pp, PROG_USAGE);

    SPOptions_GetKeyword(pp, "-outName");                               
    o.outName = SPOptions_GetNext(pp);  
    
    o.checkOrtho = SPOptions_TestKeyword(pp, "-checkOrtho");

    o.verbose = SPOptions_TestKeyword(pp, "-verbose");

    SPOptions_GetKeyword(pp, "-degree");
    o.degree = SPOptions_GetNextInt(pp, 0, INT_MAX);

    SPOptions_GetKeyword(pp, "-smpOrder");                               
    o.smpOrder = SPOptions_GetNextInt(pp, 1, INT_MAX);  
    
    SPOptions_Finish(pp);
    return o;
  }
