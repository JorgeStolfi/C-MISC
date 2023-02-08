/* SOApprox.h -- tools for function approximation. */
/* Last edited on 2007-01-04 00:21:13 by stolfi */

#ifndef SOApprox_H
#define SOApprox_H

#include <SOMatrix.h>
#include <SOFunction.h>
#include <SOPlotParams.h>
#include <SOGrid.h>
#include <SOBasic.h> 

#include <dg_tree.h>

#include <vec.h> 

/* The procedures in this interface are mostly wrappers for procedures
  in other interfaces, such as {SOMatrix.h}, {SOFunction.h}, with some
  additional facilities for measuring, printing, and plotting the
  approximation errors. */ 
  
/* READING THINGS FROM FILES */

/* The procedures in this section read some object from the file
  {fileName}. The file name must include the extension; see {addext}
  in {js.h}.
  
  If {fileName == "-"}, these procedures will read from {stdin}. They
  fail with error if the file does not exist. Procedures that return
  pointers will return NULL if {fileName == ""}. */
  
SOGrid_Tree *SOApprox_ReadTree(char *fileName);
  /* Reads a dyadic grid. Returns NULL if {fileName == ""}. */

SOFunction *SOApprox_ReadFunction(char *fileName);
  /* Reads a function. Returns NULL if {fileName == ""}. */

Basis SOApprox_ReadBasis(char *fileName);
  /* Reads a function basis.*/

SOMatrix SOApprox_ReadMatrix(char *fileName);
  /* Reads an {SOMatrix}. */

/* OBTAINING RIGIDITY MATRICES */

SOMatrix SOApprox_GetBasisMatrix(char *matName, bool_t cholesky);
  /* Obtains the rigidity matrix {E[i,j] == <B[i]|B[j]>} for some
      function basis {B}. The matrix is read from the file
    {M-ev.mat}, where {M} is the given {matName}.
    
    If {cholesky} is TRUE, instead of {E} returns the Cholesky lower 
    triangular factor {EL} of {E}, which is read from file 
    {M-ev-ch.mat}. */

SOMatrix SOApprox_GetTransferMatrix(char *matName);
  /* Obtains the basis transfer matrix {F[i,j] == <B1[i]|B2[j]>} for 
    two   function bases {B1,B2}. The matrix is read from the file
    {M-ev.mat}, where {M} is the given {matName}. */

SOMatrix SOApprox_GetHelmholtzMatrix(char *matName, double c, bool_t cholesky);
  /* Obtains the differential operator matrix 
    {H[i,j] = <D(B[i]),B[j]>}, where {B} is some   
    function basis, for the Helmholtz operator 
    {D(f)(p) = SLap(f)(p) + c*f(p)} with given coefficient {c}.
    
    The matrix {H} is actually computed according to the formula
    {H[i,j] = <SGrd(B[i]) | SGrd(B[j])> + c*<B[i]|B[j]>}. where {SGrd}
    is the   gradient operator. The matrices 
    {G[i,j] = <SGrd(B[i]) | SGrd(B[j])>} and {E[i,j] = <B[i]|B[j]>}
    are read from files {M-gr.mat} and {M-ev.mat}, where {M} is the
    given {matName}.
    
    If {cholesky} is true, computes and returns the Cholesky lower 
    triangular factor {HL} of {H}, instead of {H} itself. */

/*
  SOLVING LINEAR SYSTEMS
  
  The procedures in this section are used to solve the linear system
  {A * x == y}, which occurs in the approximation of {F} by {bas}
  (where {A[i,j] = <bas[i]|bas[j]>}), or in the solution of the
  differential equation {D(f)(p) == F(p)} 
  (where {A[i,j] == <D(bas[i])|bas[j]>}.
  
  The basis elements are assumed to be scalars (range dimension = 1);
  the target function {F} may have any range dimension {n}.

  The matrix {A} (and its triangular factors, where aplicable) are
  assumed to be square of size {m × m}, where {m} is the number of
  elements in the basis. The matrices {x} and {y} have {m} rows and
  {n} columns. They are stored as vectors of size {n*m}, in row-by-row
  order.
  
*/

void SOApprox_ComputeRightHandSide
  ( SOFunction *f,     /* Function to approximate. */
    FuncMap *FMap,     /* Right-hand side functional. */
    Basis bas,         /* Basis for the approximation space. */
    SOFunction *w,     /* Weight for integrals. */
    SOGrid_Tree *tree, /* SOGrid to use for integration. */
    double_vec_t y,    /* IN/OUT: right-hand-side matrix, linearized. */ 
    bool_t verbose       /* TRUE mumbles and prints the change in {y}. */
  );
  /* (Re)computes the right-hand-side matrix {y[i,j] == <bas[i]|F[j]>},
    where {F[j](p) = FMap(f(p), p)[j]}. In particular,
    if {FMap == NULL}, reduces to {y[i,j] == <bas[i]|f[j]>}.
    
    If {verbose} is true, also computes and prints the distance
    between the newly computed value of {y} and its input value. */
    
void SOApprox_CholeskySolve
  ( SOMatrix L,       /* Cholesky factor of system's matrix. */
    int n,            /* Number of columns in soln and rhs matrices. */
    double_vec_t y,   /* Right-hand-side matrix of system. */ 
    double_vec_t x,   /* (OUT) Coordinates of new solution relative to {bas}. */
    double_vec_t r,   /* (WORK) temporary storage area. */
    double_vec_t s,   /* (WORK) temporary storage area. */
    double_vec_t t,   /* (WORK) temporary storage area. */
    bool_t verbose      /* TRUE mumbles along the way. */
  );
  /* Solves the system {A x == y} on the unknown matrix {x}, for
    a square positive-definite symmetric matrix {A}, given the
    matrix {y} and the lower triangular factor {L} of the 
    Cholesky decomposition {L * L^t} of the matrix {A}.
    The solution is returned in the argument {x}.
    
    The vectors {r}, {s}, and {t} must have size {m}, and are used as
    temporary storage areas. */
  
void SOApprox_GaussSeidelIteration(
    SOMatrix A,        /* System's matrix. */          
    int n,             /* Number of columns in soln and rhs matrices. */
    double_vec_t y,    /* Right-hand-side of system. */ 
    double omega,      /* Overrelaxation factor. */
    double_vec_t xOld, /* (IN) Current guess at the solution. */
    double_vec_t xNew, /* (OUT) New solution. */
    double_vec_t r,    /* (WORK) temporary storage area. */
    double_vec_t s,    /* (WORK) temporary storage area. */
    bool_t verbose       /* TRUE mumbles along the way. */
  );
  /* Performs one iteration of the Gauss-Seidel algorithm for solving
    the system {A x == y}.  Upon entry, {xOld} should be the coordinates
    of a guess for the solution, relative to some basis.  Upon exit,
    {xNew} will have been updated by performing one iteration of the
    Gauss-Seidel algorithm.
    
    The vectors {r} and {s} must have size {m}, and are used as
    temporary storage areas. */
  
void SOApprox_GaussSeidelSolve
  ( SOMatrix A,        /* System's matrix. */          
    int n,             /* Number of columns in soln and rhs matrices. */
    double_vec_t y,    /* Right-hand-side of system. */ 
    double omega,      /* Relaxation factor. */
    int maxIter,       /* Maximum number of iterations. */
    double absTol,     /* Stopping criterion: absolute change in solution. */
    double relTol,     /* Stopping criterion: relative change in solution. */
    double_vec_t xOld, /* (IN) Initial guess at solution. */
    double_vec_t xNew, /* (OUT) New solution. */
    double_vec_t r,    /* (WORK) temporary storage area. */
    double_vec_t s,    /* (WORK) temporary storage area. */
    bool_t verbose       /* TRUE mumbles along the way. */
  );
  /* Attempts to solve the linear system {A x == y} by multiple
    iterations of the Gauss-Seidel method.  Upon entry, {xOld} should be
    an initial guess for the solution.  Upon exit, {xNew} will be the
    tentative solution for the system.
    
    The procedure will stop when the L2 difference between
    successive trials is at most {absTol}, or the relative difference
    is less than {relTol}, or after at most {maxIter} iterations.
    
    The vectors {r} and {s} must have size {m}, and are used as
    temporary storage areas. */
  
void SOApprox_ConjugateGradientSolve
  ( SOMatrix A,         /* System's matrix. */          
    int n,              /* Number of columns in soln and rhs matrices. */
    double_vec_t y,     /* Right-hand-side of system. */ 
    double_vec_t xOld,  /* (IN) Current solution. */
    double_vec_t xNew,  /* (OUT) New solution. */
    double_vec_t r,     /* (WORK) temporary storage area. */
    double_vec_t s,     /* (WORK) temporary storage area. */
    bool_t verbose        /* TRUE mumbles along the way. */
  );
  /* Attempts to solve the linear system {A x == y} by conjugate
    gradient minimization of the squared norm of the residual matrix 
    {y - A x}. Upon entry, {xOld} should be the a guess for the solution.
    Upon  exit, {xNew} will be a tentative solution for the system.
    
    The vectors {r} and {s} must have size {m}, and are used as
    temporary storage areas. */
  
double SOApprox_UpdateSolution
  ( double_vec_t xNew,  /* (IN) New solution vector. */
    double_vec_t x      /* (I/O)Current solution vector. */
  );
  /* This procedure is handy when using iterative methods like
    GaussSeidel. It replaces the previous solution matrix {x}
    by the newly computed solution matrix {xNew}, comparing the two.
    Returns the L2 norm of the difference. */

void SOApprox_GuessSol(double_vec_t x);
  /* Fills {x} with random values. */

SOFunction *SOApprox_BuildFunction
  ( Basis bas,       /* Basis for space {V[r]}. */
    double_vec_t u,  /* Coordinates of {g} relative to {basis} */
    char *solName,   /* Name prefix for output file. */
    char *basName,   /* Name of the (Basis bas) description file */
    SOFunction *f    /* True solution (for comparisons) */
  );
  /* Returns the function {g} with coefficient matrix {u} relative 
    to the given {basis}. Also writes {g} to disk (with the 
    given {solName} extended with "-g.fun") and compares it
    to the given function {f}, with {SOApprox_PrintSomeValues}. */

/* Added on 2003-07-04 by stolfi */
void SOApprox_PlotFunction
  ( SOFunction *f, 
    SOGrid_Tree *tree,
    char *fName, 
    PlotOptions *plt,
    double *fObsMin,    /* (IN/OUT) Min value seen during plot. */ 
    double *fObsMax     /* (IN/OUT) Max value seen during plot. */ 
  );
/* Plots the function {f} alone, to a file called "{fName}.ps" 
  or "{fName}.eps", as specified by the plotting parameters in {plt}. 
  Uses {SO2DPlot_FunctionFigure}. The parameters {fObsMin} and {fObsMax}
  must be initialized by the client, and are updated with the 
  maximum and minimum values of {f} found during the plot. */ 

//SOPlot_File *SOApprox_PlotFunctionAndError
//  ( SOPlot_File *fps,        /* Postscript file, or NULL. */
//    SOPlot_PageState *pgs,   /* Page layout and state for {fps}. */
//    char *docName,           /* Document file name, minus extension. */
//    char *funcName,          /* Prefix for separate figure file names. */
//    SOFunction *g,         /* The computed approximation. */
//    SOFunction *f,         /* The target function. */
//    SOGrid *tree,      /* Reference grid. */
//    double relMeshSize,      /* Maximum step/triangle size (radians). */
//    bool_t showTriang,         /* TRUE to plot the reference grid. */ 
//    double gMax,             /* Nominal maximum of {fabs(g(p))}. */
//    double eMax,             /* Nominal maximum of {fabs(g(p)-f(p))}. */
//    SOH3_Point *obs,          /* Viewpoint. */
//    SOH3_Point *upp           /* Camera vertical reference */
//  );
  /* Plots the function {g}, and the difference {g-f}, using
    {SOPlot_BothSides} (q.v.). 
    
    The file {fps} may be NULL, or may contain other figures; the
    procedure calls {SOPlot_BeginFigure} internally. This may cause
    the file {fps} may be closed, and another one opened in its place;
    the final file, still open, is returned as result. Be sure to call
    {SOPlot_Close} on it eventually.
    
    The {docName} parameter is passed to {SOPlot_BothSides} as is. The
    {funcName} parameter gets extended with {-app} and {-err} for the
    function and the error plots, respectively (to which
    {SOPlot_BothSides} will append {-f} for front or {-v} for back).
    
    Rememeber to call {SOPlot_Close} afterwards. */

void SOApprox_PrintSomeValues(SOFunction *g, SOFunction *f);
  /* Prints values of {f(p)}, {g(p)}, and their difference 
    at a few points in the root cell of the grid. */

void SOApprox_PrintMaxErrorValues
  ( SOFunction *g,  /* Approximation. */
    SOFunction *f,  /* Target function. */
    double *gMax,     /* OUT: max abs function value, per coordinate. */
    double *eMax      /* OUT: max abs error value, per coordinate. */
  );
  /* Computes and prints the maximum absolute value {gMax[j]} of the
    function {g}, and the maximum absolute error {eMax[j]} when
    compared to {f}, for each range coordinate {j}. */

#endif
