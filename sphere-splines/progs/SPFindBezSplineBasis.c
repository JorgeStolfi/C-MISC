/* Constructs bases for spherical spline spaces. */
/* Last edited on 2006-03-19 12:49:06 by stolfi */

#define PROG_NAME "SPFindBezSplineBasis"

/* This program reads a triangulation {T} and computes a basis for the
  space of piecewise polynomials of a given degree over that
  triangulation -- either homogeneous or non-homegeneous. The basis
  elements are written out, in the Bezier coefficient representation
  (as an array of {SPSpline}s, whose pieces are {SPHBezFunction}s or 
  {SPNHBezFunction}s, respectively). */

#include <SPTriang.h>
#include <SPFunction.h>
#include <SPBezSplineBasis.h>
#include <SPOptions.h>
#include <SPIntegral.h>
#include <SPSpline.h>
#include <SPBasic.h>
#include <bool.h>
#include <stdio.h>
#include <values.h>

#define PROG_USAGE \
  PROG_NAME " \\\n" \
  "  -triName TRINAME \\\n" \
  "  -degree NUM -continuity NUM \\\n" \
  "  [ -homogeneous ] \\\n" \
  "  [ -oldStyle | -newStyle ] \\\n" \
  "  -outName BASISNAME\n"

typedef struct Options
  { char *triName;      /* Name of triangulation file (minus the ".tri" extension). */
    char *outName;      /* Name of basis file (minus the ".bas" extension) */
    int degree;         /* Degree of polynomials in each triangle. */
    int continuity;     /* Order of continuity across edges & vertices. */
    bool_t newStyle;    /* TRUE builds the modified ANS basis (with boat elems). */
    bool_t homogeneous; /* TRUE builds homog basis, FALSE non-homog one. */
  } Options;
  
Options GetOptions(int argn, char **argc);  

int main(int argn, char **argc)
  { Options o = GetOptions(argn, argc);
    char *triFile = txtcat(o.triName, ".tri");
    FILE *rd = open_read(triFile, TRUE);
    int smpOrder = 8;  /* Should be irrelevant. */
    SPIntegral_SetDefaultSamplingOrder(smpOrder);
    Triangulation *tri = SPTriang_Read(rd, smpOrder);
    Basis F;
    
    int deg = o.degree;
    int cont = o.continuity;
    bool_t new = o.newStyle;
    if (o.homogeneous)
      { F = SPBezSplineBasis_BuildH(deg, cont, new, triFile, tri); }
    else
      { F = SPBezSplineBasis_BuildNH(deg, cont, new, triFile, tri); }

    /* Write basis: */
    { FILE *wr = open_write(txtcat(o.outName, ".bas"), TRUE);
      SPFunction_WriteBasis(wr, F);
      fclose(wr);
    }

    return 0;
  }

Options GetOptions(int argn, char **argc)
  {
    Options o;
    SPOptions_Parser_t *pp = SPOptions_NewParser(stderr, argn, argc);
    SPOptions_SetUsage(pp, PROG_USAGE);

    SPOptions_GetKeyword(pp, "-triName");                               
    o.triName = SPOptions_GetNext(pp);  

    SPOptions_GetKeyword(pp, "-degree");
    o.degree = SPOptions_GetNextInt(pp, 0, 1000);

    SPOptions_GetKeyword(pp, "-continuity");
    o.continuity = SPOptions_GetNextInt(pp, -1, 1000);

    o.homogeneous = SPOptions_TestKeyword(pp, "-homogeneous");

    if (SPOptions_TestKeyword(pp, "-newStyle"))
      { o.newStyle = TRUE; }
    else if (SPOptions_TestKeyword(pp, "-oldStyle"))
      { o.newStyle = FALSE; }
    else
      { o.newStyle = TRUE; }

    SPOptions_GetKeyword(pp, "-outName");                               
    o.outName = SPOptions_GetNext(pp);  

    SPOptions_Finish(pp);
    return o;
  }
