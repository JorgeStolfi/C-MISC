/* SOComputeErrorMap - Generates an Error Map of a test space over a dyadic grid. */
/* Last edited on 2006-03-19 12:42:14 by stolfi */

/*
  This program computes a least squares approximation {g}, for each 
    function {f}, of the form
  
    g = SUM { u[i]*bas[i] : i = 0..dim-1 }       (1)
    
  where {bas[0..dim-1]} is an arbitrary (but add-compatible) basis.
  The coefficients {u[i]} are found by solving the linear system
    
    G u == b    (2)
    
  where {G} is the rigidity matrix, {G[i,j] == <bas[i] | bas[j]>},
  and the vector {b} is given by {b[i] == <f | bas[i]>}.
*/
    
#include <SOGrid.h>
#include <SOMatrix.h>
#include <SOFunction.h>
#include <SOLinCombFunction.h>
#include <SOErrorFunction.h>
#include <SOApprox.h>
#include <SOParams.h>
#include <SOBasic.h>
#include <SOIntegral.h>
#include <SO2DPlot.h>

#include <affirm.h>
#include <nat.h>
#include <vec.h>

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <values.h>
#include <limits.h>

#define MAXSTR 50
#define MAX_PLOT_DEPTH 40
#define MAX_PLOT_STEPS 12

typedef struct Options
  { char *testbasisName;  /* Prefix of the test basis (wave functions) file name. */
    char *basisName;      /* Prefix of basis file names. */
    char *matName;        /* Prefix of matrix file names (default {basisName}). */
    char *outName;        /* Prefix for output file names. */
    nat gaussOrder;       /* Gaussian quadrature order for dot products. */
    /* Method used to solve the linear system: */
    bool_t cholesky;        /* TRUE uses direct method based on Cholesky factzn. */
    bool_t gaussSeidel;     /* TRUE uses Gauss-Seidel iteration. */
    bool_t conjGrad;        /* TRUE uses the conjugate gradient method. */
    /* Parameters for Gauss-Seidel iteration: */
    double omega;         /* Relaxation factor. */
    nat maxIter;          /* Maximum iterations. */
    double relTol;        /* Stopping criterion: relative change. */
    double absTol;        /* Stopping criterion: absolute change. */
    /* Plotting options: */
    //    bool_t plot;            /* TRUE to generate function and error plots. */
    //    PlotOptions plt;      /* Plot format and style options. */
  } Options;

Options *GetOptions(int argn, char **argc);

int main(int argn, char **argc)
  { int e, m, n, dim;
    SOGrid_Tree *tree;
    SOFunction *f, *w = (SOFunction *)NULL; 
    Options *o = GetOptions(argn, argc);

    char *testbasName = txtcat(o->testbasisName, ".bas");
    char *basName = txtcat(o->basisName, ".bas");

    Basis testbas = SOApprox_ReadBasis(testbasName);
    Basis bas = SOApprox_ReadBasis(basName);
    Basis appbas = SOFunctionRef_vec_new(testbas.ne);

    affirm(testbas.ne <= bas.ne, " Basis must have more elements than test basis ");

    /* Initial values for 2D plotting */
    double min[2], max[2];
    FILE *plot_file = open_write("error_map.eps", TRUE);
    char *treeFile = txtcat(o->basisName, ".tree");
    FILE *rd = open_read(treeFile, TRUE);
    tree = SOGrid_Tree_read(rd);
    min[X] = 0.0; max[X] = 1.0;
    min[Y] = 0.0; max[Y] = (double)sqrt(2.0)/2;
    /* ******************************* */

    m = testbas.e[0]->d->pDim;
    n = testbas.e[0]->d->fDim;
    dim = bas.ne;

    double_vec_t b = double_vec_new(n * dim); // to receive <basis | f>, f's results on R^fDim
    double_vec_t uj = double_vec_new(dim); // Column {j} of {u}.
    double_vec_t bj = double_vec_new(dim); // Column {j} of {b}.
    double_vec_t y = double_vec_new(dim);
    FuncMap NoMap = IdentityMap(n, m);
    SOIntegral_GaussOrder = o->gaussOrder; 

    for(e = 0; e < testbas.ne; e++)
      {

        f = testbas.e[e];    
        double_vec_t coef = double_vec_new(n * dim);

        { int i; for (i = 0; i < b.ne; i++) { b.e[i] = 0.0; } } 

        SOApprox_ComputeRightHandSide(f, &NoMap, bas, w, tree, b, FALSE); 
 
        if (o->cholesky)
          { 
            SOMatrix GL = SOApprox_GetBasisMatrix(o->matName, TRUE);
            SOApprox_CholeskySolve(GL, n, b, coef, bj, uj, y, TRUE);
          }
        else
          { affirm(FALSE , "unspecified/invalid solution method"); }
       
        appbas.e[e] = SOApprox_BuildFunction(bas, coef, o->outName, basName, f);      
      }

    double_vec_t weight = double_vec_new(testbas.ne);
    for(e = 0; e < testbas.ne; e++) weight.e[e] = 1.0;

    SOFunction *e = (SOFunction *)SOErrorFunction_Make(
      m, n, "testbas.bas", "appbas.bas", testbas, appbas, weight
    );

    //    FILE *appbas_wr = open_write("appbas.bas", TRUE);
    //    SOFunction_WriteBasis(appbas_wr, appbas);
    //    FILE *testbas_wr = open_write("testbas.bas", TRUE);
    //    SOFunction_WriteBasis(testbas_wr, testbas);

    /* 2D Plotting */
    
    !!! Deve usar SOApprox_PlotFunction, é mais fácil
    plot_file := SO2DPlot_FunctionFigure(e, plot_file, min, max, tree, -1, 6, MAX_PLOT_STEPS,0, MAX_PLOT_DEPTH);
    
    if ()

    return 0;
  }

#define PPUSAGE SOParams_SetUsage

Options *GetOptions(int argn, char **argc)
  {
    Options *o = (Options *)notnull(malloc(sizeof(Options)), "no mem");
    SOParams_T *pp = SOParams_NewT(stderr, argn, argc);

    PPUSAGE(pp, "SOComputeApprox \\\n");
    PPUSAGE(pp, "  -funcName NAME \\\n");
    PPUSAGE(pp, "  -basisName NAME [-matName NAME ] \\\n");
    PPUSAGE(pp, "  [ -cholesky | \\\n");
    PPUSAGE(pp, "    -gaussSeidel [-omega NUM] [-tol NUM] [-maxIter NUM] | \\\n");
    PPUSAGE(pp, "    -conjGrad ] \\\n");
    PPUSAGE(pp, "  [ -gaussOrder NUM ] \\\n");
    PPUSAGE(pp, "  -outName NAME \\\n");
    //    PPUSAGE(pp, "  [ -plot ] \\\n");
    //    PPUSAGE(pp, SOPlotParams_FunctionHelp " \n");

    SOParams_GetKeyword(pp, "-testbasisName");                               
    o->testbasisName = SOParams_GetNext(pp);  
       
    SOParams_GetKeyword(pp, "-basisName");                               
    o->basisName = SOParams_GetNext(pp);  
    
    if (SOParams_KeywordPresent(pp, "-matName"))
      { o->matName = SOParams_GetNext(pp); }
    else
      { o->matName = o->basisName; }
       
    SOParams_GetKeyword(pp, "-outName");                               
    o->outName = SOParams_GetNext(pp);  
                                                 
    o->cholesky = SOParams_KeywordPresent(pp, "-cholesky");
    o->gaussSeidel = SOParams_KeywordPresent(pp, "-gaussSeidel");
    o->conjGrad = SOParams_KeywordPresent(pp, "-conjGrad");
    if (! (o->cholesky || o->gaussSeidel || o->conjGrad))
      { SOParams_Error(pp, "must specify a linear solution method"); }
    else if ((o->cholesky + o->gaussSeidel + o->conjGrad) > 1)
      { SOParams_Error(pp, "please specify only one linear solution method"); }

    if (o->gaussSeidel)
      { 
        if (SOParams_KeywordPresent(pp, "-omega"))
          { o->omega = SOParams_GetNextDouble(pp, 0.0, INFTY); }
        else
          { o->omega = 1.0; }

        if (SOParams_KeywordPresent(pp, "-maxIter"))
          { o->maxIter = SOParams_GetNextInt(pp, 0, INT_MAX); }
        else
          { o->maxIter = 20; }

        if (SOParams_KeywordPresent(pp, "-absTol"))
          { o->absTol = SOParams_GetNextDouble(pp, 0.0, INFTY); }
        else
          { o->absTol = 0.0; }

        if (SOParams_KeywordPresent(pp, "-relTol"))
          { o->relTol = SOParams_GetNextDouble(pp, 0.0, INFTY); }
        else
          { o->relTol = 0.0; }
      }

    SOParams_GetKeyword(pp, "-gaussOrder");
    o->gaussOrder = SOParams_GetNextInt(pp, 1, INT_MAX);

    /* Plotting options: */

    // o->plot = SOParams_KeywordPresent(pp, "-plot");
    // o->plt = SOPlotParams_FunctionParse(pp);

    SOParams_Finish(pp);
    return o;
  }

