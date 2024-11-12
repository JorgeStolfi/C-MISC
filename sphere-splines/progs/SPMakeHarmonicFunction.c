/* SPMakeHarmonicFunction -- Creates a spherical harmonic function file. */
/* Last edited on 2006-03-19 12:49:19 by stolfi */

#define PROG_NAME "SPMakeHarmonicFunction" 

/* Creates a function file for a real spherical harmonic term. */

/* !!! TO DO: extend to multi-term harmonic series. */

#include <SPBasic.h>
#include <SPOptions.h>
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
  "  -degree NUM -order NUM \\\n" \
  "  -outName OUTNAME\n"

typedef struct Options
  { int degree;        /* Degree of spherical harmonic. */
    int order;         /* Order of spherical harmonic. */
    char *outName;     /* Prefix for output file names. */
  } Options;

Options GetOptions(int argn, char **argc);

int main(int argn, char **argc)
  { Options o = GetOptions(argn, argc);
    SPFunction *f = (SPFunction *)SPHarmonic_FromTerm(o.degree, o.order);
    FILE *wr = open_write(txtcat(o.outName, ".sfn"), TRUE);
    f->m->write(f, wr);
    fclose(wr);
    return 0;
  }

Options GetOptions(int argn, char **argc)
  {
    Options o;
    SPOptions_Parser_t *pp = SPOptions_NewParser(stderr, argn, argc);
    SPOptions_SetUsage(pp, PROG_USAGE);

    SPOptions_GetKeyword(pp, "-degree");
    o.degree = SPOptions_GetNextInt(pp, 0, INT_MAX);

    SPOptions_GetKeyword(pp, "-order");
    o.order = SPOptions_GetNextInt(pp, -o.degree, +o.degree);

    SPOptions_GetKeyword(pp, "-outName");                               
    o.outName = SPOptions_GetNext(pp);  

    SPOptions_Finish(pp);
    return o;
  }
