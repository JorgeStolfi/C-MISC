/* SPMakePolyBasis -- Find basis for spherical polynomial (non-spline) spaces. */
/* Last edited on 2006-03-19 12:49:21 by stolfi */

/* This program outputs a basis for the space of spherical polynomials
  (global -- NOT spherical splines) of a given total degree, either
  homogeneous or non-homegeneous. The basis elements are written out
  in the Bezier coefficient representation ({SPHBezFunction}s or
  {SPNHBezFunction}s, respectively). */

#include <SPFunction.h>
#include <SPHBezFunction.h>
#include <SPNHBezFunction.h>
#include <SPDeCasteljau.h>
#include <SPOptions.h>
#include <SPBasic.h>
#include <bool.h>
#include <stdio.h>

typedef struct Options
  { char *outName;      /* Name of basis file (minus the ".bas" extension) */
    int degree;         /* Degree of polynomials in each triangle. */
    bool_t homogeneous; /* TRUE for homogeneous basis, FALSE non-homogeneous one. */
  } Options;
  
Options GetOptions(int argn, char **argc); 
void PrintElem(int basisIndex, int degree, int coeffIndex);

int main(int argn, char **argc)
  { int i, j;
    Options o = GetOptions(argn, argc);
    Basis bas;
    if (o.homogeneous)
      { int dim = SPHBezFunction_NumCoeffs(o.degree);
        SPHBezFunction *f;
        bas = SPFunctionRef_vec_new(dim);
        for(i = 0; i < dim; i++)
          { BezCoeff_vec_t c = BezCoeff_vec_new(dim);
            for(j = 0; j < dim; j++) { c.e[j] = 0.0; }
            c.e[i] = 1.0;
            f = SPHBezFunction_FromCoeffs(o.degree, c);
            PrintElem(i, o.degree, i);
            bas.e[i] = (SPFunction *)f;
          }
      }
    else
      { int dim0 = SPHBezFunction_NumCoeffs(o.degree);
        int dim1 = SPHBezFunction_NumCoeffs(o.degree - 1);
        SPNHBezFunction *f;
        bas = SPFunctionRef_vec_new(dim0 + dim1);
        for(i = 0; i < dim0; i++)
          { BezCoeff_vec_t c0 = BezCoeff_vec_new(dim0);
            BezCoeff_vec_t c1 = BezCoeff_vec_new(dim1);
            for(j = 0; j < dim0; j++) { c0.e[j] = 0.0; }
            for(j = 0; j < dim1; j++) { c1.e[j] = 0.0; }
            c0.e[i] = 1.0;
            f = SPNHBezFunction_FromCoeffs(o.degree, c0, c1);
            PrintElem(i, o.degree, i);
            bas.e[i] = (SPFunction *)f;
          }
        for(i = 0; i < dim1; i++)
          { BezCoeff_vec_t c0 = BezCoeff_vec_new(dim0);
            BezCoeff_vec_t c1 = BezCoeff_vec_new(dim1);
            for(j = 0; j < dim0; j++) { c0.e[j] = 0.0; }
            for(j = 0; j < dim1; j++) { c1.e[j] = 0.0; }
            c1.e[i] = 1.0;
            f = SPNHBezFunction_FromCoeffs(o.degree, c0, c1);
            PrintElem(dim0 + i, o.degree-1, i);
            bas.e[dim0 + i] = (SPFunction *)f;
          }
      }
    { FILE *wr = open_write(txtcat(o.outName, ".bas"), TRUE);
      SPFunction_WriteBasis(wr, bas);
      fclose(wr);
    }
    return 0;
  }

void PrintElem(int basisIndex, int degree, int coeffIndex)
  {
    BezLabel e = SPDeCasteljau_IndexToBezLabel(coeffIndex, degree);
    fprintf(stderr, "  element %04d = B(%d,%d,%d)\n", 
      basisIndex, e.e[0], e.e[1], e.e[2]
    );
  }


Options GetOptions(int argn, char **argc)
  {
    Options o;
    SPOptions_Parser_t *pp = SPOptions_NewParser(stderr, argn, argc);
    SPOptions_SetUsage(pp, 
      "SPMakePolyBasis \\\n"
      "  -degree NUM [ -homogeneous ] \\\n"
      "  -outName BASISNAME\n"
    );

    SPOptions_GetKeyword(pp, "-degree");
    o.degree = SPOptions_GetNextInt(pp, 0, 20);

    o.homogeneous = SPOptions_TestKeyword(pp, "-homogeneous");

    SPOptions_GetKeyword(pp, "-outName");                               
    o.outName = SPOptions_GetNext(pp);  

    SPOptions_Finish(pp);
    return o;
  }
