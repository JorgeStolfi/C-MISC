/* SOUniSolve -- Solves the screened Poisson equation by uniscale iteration. */
/* Last edited on 2023-10-15 03:38:36 by stolfi */

/*
  This program solves a differential equation on the root cell
  
    D(f)(p) == FMap(f(p), p)       (1)
    
  where {D} is a linear differential operator, {FMap} is a given
  function, and {f} is the function to be determined. For example, 
  in the screened Poisson equation on the root cell, the differential
  operator is
  
    D(f)(p) = Lap(f)(p) + c*f(p)       (2)
    
  where {Lap} is the Laplacian operator, and {c} is a given
  constant.
  
  It is assumed that an approximate solution {g} to (2) in some
  finite-dimensional function space {V} can be found by solving the
  non-linear system
  
    H a == b         (3)
    
  where 
  
    {a[i]} is coefficient {i} of the solution 
           {g} (to be determined) in the basis {bas}; 
  
    {b[i] == <F | bas[i]>}, where {F} is the   function
            {p -> FMap(g(p),p)} and {<|>} is the scalar product;
            
    {H[i,j] == <D(bas[i]) | bas[j]>}
           where {D} is the differential operator of eq.(1).
           
  For the screened Poisson equation, we have
  
     {H[i,j] == <Lap(bas[i]) | bas[j]> + c*<bas[i]|bas[j]>}
     
  If the basis functions {bas[i]} are at least C1, the scalar product
  {G[i,j] = <Lap(bas[i]) | bas[j]>} can be reduced to 
  
     {G[i,j] = <Grd(bas[i]) | Grd(bas[j])>}
     
  where {Grd} is the gradient operator. 
        
  Note that, in general, the system (3) is not linear, since the
  right-hand side {b} depends on the function {g}, and therefore on
  the coefficient vector {a}. So (3) must be solved iteratively.
*/

#include <SOGrid.h>
#include <SOMatrix.h>
#include <SOFunction.h>
#include <SOProcFunction.h>
#include <SOApprox.h>
#include <SOFuncMap.h> 
#include <SOParams.h>
#include <SOBasic.h>
#include <SO2DPlot.h>

#include <dg_grid.h>

#include <rn.h>
#include <affirm.h>
#include <nat.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <values.h>
#include <limits.h>

#define MAX_PLOT_DEPTH 40
#define MAX_PLOT_STEPS 14

typedef struct Options
  { char *rhsName;     /* Name of right-hand-side function. */
    char *basisName;   /* Prefix for basis filenames (minus extension). */
    char *matName;     /* Prefix for matrix filenames (default {basisName}). */
    char *outName;     /* Prefix for output filenames. */
    nat maxIter;       /* Maximum iterations. */
    dg_dim_t pDim;   /* Dimension of {f}'s domain. */
    dg_dim_t fDim;   /* Dimension of {f}'s values. */
    double relTol;     /* Max change in {a} coeffs between iterations. */
    double absTol;     /* Max error compared to true sol. */
    nat gaussOrder;       /* Gaussian quadrature order for dot products. */
    /* Method used to solve the linear system: */
    bool_t cholesky;     /* TRUE uses direct method based on Cholesky factzn. */
    bool_t gaussSeidel;  /* TRUE uses Gauss-Seidel iteration. */
    bool_t conjGrad;     /* TRUE uses the conjugate gradient method. */
    double omega;      /* Relaxation factor for Gauss-Seidel. */
    /* Plotting options: */
    bool_t plotFinal;    /* TRUE to plot final solution and error. */
    bool_t plotAll;      /* TRUE to plot every iteration and error. */
    //    PlotOptions plt;   /* Plot format and style options. */
  } Options;

Options GetOptions(int argn, char **argc);
SOFuncMap_Data GetRHSFunction(char *rhsName, dg_dim_t uDim, dg_dim_t pDim, dg_dim_t vDim);
SOFunction *GetTrueSolution(char *solName, dg_dim_t pDim, dg_dim_t fDim);
char *SolutionName(char *outName, int iter);

void ComputeRightHandSide
  ( SOGrid_Tree *tree,  /* SOGrid to use for integration. */
    Basis bas,          /* Basis for the approximation space. */
    SOFunction *g,      /* Current solution. */
    FuncMap *FMap,      /* Right-hand side of differential eqn. */
    nat depth,          /* Depth of sampling grid in each triangle */
    double_vec_t b        /* OUT: Right-hand-side of (3). */ 
  );
  /* Computes the right-hand side of the system (3) for the given
    solution {g} which has coordinates {a} relative to basis {bas}. */

int main(int argn, char **argc)
  { Options o = GetOptions(argn, argc);

    /* Gets basis */
    char *basName = txtcat(o.basisName, ".bas");
    Basis bas = SOApprox_ReadBasis(basName);
    nat dim = bas.nel;

    /* Gets the {SOGrid} tree */
    char *treeFile = txtcat(o.basisName, ".tree");
    FILE *rd = open_read(treeFile, TRUE);
    SOGrid_Tree *tree = SOGrid_Tree_read(rd);

    /* Initial values for 2D plotting */
    double min[2], max[2];
    FILE *sol_file = open_write("sol.eps", TRUE);
    FILE *approx_file = open_write("approx.eps", TRUE);
    min[X] = 0.0; max[X] = 1.0;
    min[Y] = 0.0; max[Y] = (double)sqrt(2.0)/2;
    /* ******************************* */

    SOFuncMap_Data FD = SOFuncMap_FromName(o.rhsName, o.fDim, o.pDim, o.fDim);
    SOFunction *sol = GetTrueSolution(FD.sol, o.pDim, o.fDim);
    
    double_vec_t a = double_vec_new(dim); /* Coeffs of current solution. */ // vDim dimensional?????
    double_vec_t b = double_vec_new(dim); /* Right-hand side vector. */     // vDim dimensional?????

    double_vec_t r = double_vec_new(dim);
    double_vec_t s = double_vec_new(dim);
    double_vec_t t = double_vec_new(dim);

    double_vec_t aNew = double_vec_new(dim);                                // vDim dimensional?????
    
    SOMatrix H = SOMatrix_Null(dim,dim);  /* The matrix of (3), if needed. */
    SOMatrix HL = SOMatrix_Null(dim,dim); /* Lower Cholesky factor, if needed. */
    
    SOFunction *g;  /* Current solution. */
    
    FILE *errWr = open_write(txtcat(o.outName, ".erp"), TRUE);
    
    double fMax[o.fDim], eMax[o.fDim];
    double d, m;
    nat iter = 0;
    
    /* Obtain the necessary matrices: */
    if (o.cholesky)
      { HL = SOApprox_GetHelmholtzMatrix(o.matName, FD.coeff, TRUE); }
    else
      { H = SOApprox_GetHelmholtzMatrix(o.matName, FD.coeff, FALSE); }

    SOApprox_GuessSol(a); 
    while (TRUE)
      { char *uName = txtcat3(o.outName, "-", fmt_int(iter, 4));

        g = SOApprox_BuildFunction(bas, a, uName, basName, sol);

	SOApprox_PrintMaxErrorValues(g, sol, fMax, eMax);

        fprintf(errWr, "%6d %16.12f\n",  iter, eMax[X]);

        if ((iter > 0) && ((d <= o.relTol*m) && (d <= o.absTol))) { break; }
        if ((eMax[X] <= o.absTol) || (iter >= o.maxIter)) { break; }
 
        SOApprox_ComputeRightHandSide(g, &FD.FMap, bas, NULL, tree, b, TRUE);
        if (o.cholesky)
          { SOApprox_CholeskySolve(HL, 1, b, aNew, r, s, t, TRUE); }
        else
          { affirm(FALSE , "no solution method?"); }

        d = SOApprox_UpdateSolution(aNew, a);
        m = rn_norm(dim, a.el);
        iter++;
      }

    fclose(errWr);

    /* /\* 2D Plotting *\/ */
    SO2DPlot_Begin_Figure(sol_file, min, max );
    SO2DPlot_Begin_Figure(approx_file, min, max );
 
    SO2DPlot_Isolines(sol, sol_file, tree, 0, 0.07, MAX_PLOT_STEPS, 0);
    SO2DPlot_Isolines(g, approx_file, tree, -0.01, 0.07, MAX_PLOT_STEPS, 0);
 
    SO2DPlot_Plot_Tree(sol_file, tree, min, max, MAX_PLOT_DEPTH );
    SO2DPlot_Plot_Tree(approx_file, tree, min, max, MAX_PLOT_DEPTH );

    SO2DPlot_End_Figure(sol_file);
    SO2DPlot_End_Figure(approx_file);

    return 0;
  }
  
SOFuncMap_Data GetRHSFunction(char *rhsName, dg_dim_t uDim, dg_dim_t pDim, dg_dim_t vDim)
  { return SOFuncMap_FromName(rhsName, uDim, pDim, vDim); }

SOFunction *GetTrueSolution(char *solName, dg_dim_t pDim, dg_dim_t fDim)
  { SOFunction *f = (SOFunction *)SOProcFunction_FromName(solName, pDim, fDim);
    if (f == NULL) 
      { fprintf(stderr, "Unknown solution = \"%s\"\n", solName);
        affirm(FALSE, "aborted");
      }
    return f;
  }

char *SolutionName(char *outName, int iter)
  {
    if ((iter > 0) && (iter < INT_MAX))
      { return txtcat3(outName, "-", fmt_int(iter,4)); }
    else
      { return outName; }
  }

#define PPUSAGE SOParams_SetUsage

Options GetOptions(int argn, char **argc)
  {
    Options o;
    SOParams_T *pp = SOParams_NewT(stderr, argn, argc);

    PPUSAGE(pp, "SOUNiSolve \\\n");
    PPUSAGE(pp, "  -rhsName NAME \\\n");
    PPUSAGE(pp, "  -basisName NAME [ -matName NAME ] \\\n");
    PPUSAGE(pp, "  [ -cholesky | \\\n");
    PPUSAGE(pp, "    -gaussSeidel [-omega NUM] | \\\n");
    PPUSAGE(pp, "    -conjGrad ] \\\n");
    PPUSAGE(pp, "  [ -maxIter NUM ] [ -relTol NUM ] [ -absTol NUM ] \\\n");
    PPUSAGE(pp, "  [ -gaussOrder NUM ] \\\n");
    PPUSAGE(pp, "  -outName NAME \\\n");
    PPUSAGE(pp, "  -pDim pdim \\\n");
    PPUSAGE(pp, "  -fDim fdim \\\n");
    PPUSAGE(pp, "  [ -plotAll ] [ -plotFinal ] \\\n");
    //    PPUSAGE(pp, SOPlotParams_FunctionHelp " \n");

    SOParams_GetKeyword(pp, "-rhsName");                               
    o.rhsName = SOParams_GetNext(pp);  

    SOParams_GetKeyword(pp, "-basisName");                               
    o.basisName = SOParams_GetNext(pp);  
       
    if (SOParams_KeywordPresent(pp, "-matName"))
      { o.matName = SOParams_GetNext(pp); }
    else
      { o.matName = o.basisName; }

    SOParams_GetKeyword(pp, "-outName");                               
    o.outName = SOParams_GetNext(pp);  
                                                 
    o.cholesky = SOParams_KeywordPresent(pp, "-cholesky");
    o.gaussSeidel = SOParams_KeywordPresent(pp, "-gaussSeidel");
    o.conjGrad = SOParams_KeywordPresent(pp, "-conjGrad");
    if (! (o.cholesky || o.gaussSeidel || o.conjGrad))
      { SOParams_Error(pp, "must specify a linear solution method"); }
    else if ((o.cholesky + o.gaussSeidel + o.conjGrad) > 1)
      { SOParams_Error(pp, "please specify only one linear solution method"); }

    if (o.gaussSeidel)
      { if (SOParams_KeywordPresent(pp, "-omega"))
          { o.omega = SOParams_GetNextDouble(pp, 0.0, INFTY); }
        else
          { o.omega = 1.0; }
      }

    if (SOParams_KeywordPresent(pp, "-maxIter"))
      { o.maxIter = SOParams_GetNextInt(pp, 0, INT_MAX); }
    else
      { o.maxIter = 20; }

    if (SOParams_KeywordPresent(pp, "-relTol"))
      { o.relTol = SOParams_GetNextDouble(pp, 0.0, INFTY); }
    else
      { o.relTol = 0.0; }

    if (SOParams_KeywordPresent(pp, "-absTol"))
      { o.absTol = SOParams_GetNextDouble(pp, 0.0, INFTY); }
    else
      { o.absTol = 0.0; }

    SOParams_GetKeyword(pp, "-gaussOrder");
    o.gaussOrder = SOParams_GetNextInt(pp, 1, INT_MAX);

    SOParams_GetKeyword(pp, "-pDim");                               
    o.pDim = SOParams_GetNextInt(pp, 1, 4);  

    SOParams_GetKeyword(pp, "-fDim");                               
    o.fDim = SOParams_GetNextInt(pp, 1, 4);

    /* Plotting options: */

    //    o.plotAll = SOParams_KeywordPresent(pp, "-plotAll");

    //    o.plotFinal = SOParams_KeywordPresent(pp, "-plotFinal");

    //    o.plt = SOPlotParams_FunctionParse(pp);
    
    SOParams_Finish(pp);
    return o;
  }

