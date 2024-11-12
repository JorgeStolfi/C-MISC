/* SPSys.h -- tools for solving linear and non-linear eq. systems */
/* Last edited on 2005-10-25 20:58:05 by stolfi */

#ifndef SPSys_H
#define SPSys_H

#include <SPMatrix.h>
#include <SPVector.h>
#include <SPBasic.h> 
#include <SPOptions.h> 
#include <vec.h> 
#include <js.h> 

/* SOLVING LINEAR SYSTEMS */

typedef enum   /* Solution methods for linear systems: */
  { SPSys_LM_NONE,
    SPSys_LM_SVD,         /* Direct method based on SVD factzn. */  
    SPSys_LM_GaussSeidel, /* Gauss-Seidel iterative method. */                  
    SPSys_LM_Cholesky,    /* Direct method based on Cholesky factzn. */  
    SPSys_LM_ConjGrad     /* Conjugate gradient method. */           
  } SPSys_LinMethod_t;

typedef struct SPSys_LinOptions_t  /* Parameters for linear system solving. */
  { SPSys_LinMethod_t mth; /* Solution method. */
    int residual;          /* Number of residual correction steps required. */
    /* Extra parameters for Gauss-Seidel method: */
    double omega;          /* Relaxation factor for Gauss-Seidel. */
    int maxIter;           /* Maximum iterations. */
    double relTol;         /* Stopping criterion: relative change. */
    double absTol;         /* Stopping criterion: absolute change. */
  } SPSys_LinOptions_t;

#define SPSys_LinOptionsHelp \
  "  [ -linearSys \\\n"\
  "      [ SVD | \\\n" \
  "        gaussSeidel MAXITER RELTOL ABSTOL [-omega NUM] | \\\n" \
  "        cholesky | \\\n" \
  "        conjGrad \\\n" \
  "      ] \\\n" \
  "      [ -residual NUM ] \\\n" \
  "  ]"

SPSys_LinOptions_t SPSys_LinOptionsParse(SPOptions_Parser_t *pp, bool_t local);
  /* Parses a {SPSys_LinOptions_t} from the command line,
    as described by {SPSys_LinOptionsHelp}. 
    
    If {local=TRUE}, the method keyword must be the next
    parameter to be parsed, otherwise it may be anywhere in the
    command line.  The "-omega" and "-residual" options, if present, 
    must appear at the specified locations. */

void SPSys_ComputeNeededFactors
  ( SPMatrix A,             /* System's matrix */
    SPSys_LinMethod_t mth,  /* Solution method that will be used. */
    double minVal,          /* Threshold for entry cleanup. */
    SPMatrix *L,            /* Left factor of {A}. */
    SPVector *D,            /* Diagonal factor of {A}. */
    SPMatrix *R             /* Right factor of {A}. */
  );
  /* Computes also the factors {L} and/or {R}, as requested by the
    method {mth}. Specifically, if {mth == SPSys_LM_SVD}, computes the
    factors {L}, {D} and {R} of the SVD factorization {L*D*R^T} of
    {A}. If {mth = SPSys_LM_Cholesky}, computes only the lower factor
    {L} of the Cholesky factorization {L*L^T} of {A}, and sets {D} to
    NULL, {R} to the null matrix. In all other cases {D} is set to
    NULL, and {L} and {R} are set to the null matrix of the correct
    size. The {minVal} is passed to {SPMatrix_Cholesky} or
    {SPMatrix_SVD}, as appropriate. */
    
void SPSys_LinSolve
  ( SPMatrix A,            /* System's matrix */
    SPMatrix L,            /* Left factor of {A}. */
    SPVector D,            /* Diagonal factor of {A}. */
    SPMatrix R,            /* Right factor of {A}. */
    SPVector y,            /* Right-hand-side vector of system. */
    SPSys_LinOptions_t *o, /* Solution mathod and specific parms. */
    SPVector t,            /* Initial guess and/or work area. */
    SPVector x,            /* (OUT) Solution vector. */
    bool_t verbose         /* TRUE mumbles while working. */           
  );
  /* Solves the linear system {A*x = y} by the method described
    by {o}. 
    
    Some methods use the vector {t} as a temporary storage area.
    For Gauss-Seidel and Conjugate Gradient, the vector {t}
    provides the starting guess, and should be initialized as 
    such by the client.

    The matrix {A} must be square. The vectors {t} and {x} must be
    allocated by the caller with {dim} elements, where {dim} is the
    number of rows and columns in {A}.

    Assumes that the matrices {A}, {L}, {D} and {R} are set as in the
    output of {SPSys_ComputeNeededFactors}.
    
    This procedure ignores the flag {o->residual}. If that flag is
    set, the client should compute a residual vector {dy} as
    appropriate, then call {SPSys_SolveResidual} to compute the
    adjustment {dx}, and finally {SPSys_Adjust} to update {x}. */

void SPSys_LinSolveResidual
  ( SPMatrix A,            /* System's matrix */                       
    SPMatrix L,            /* Left factor of {A}. */       
    SPVector D,            /* Diagonal factor of {A}. */
    SPMatrix R,            /* Right factor of {A}. */       
    SPVector dy,           /* Right-hand-side residual vector. */                
    SPSys_LinOptions_t *o, /* Solution method and specific parms. */   
    SPVector t,            /* Initial guess and/or work area. */       
    SPVector dx,           /* (OUT) Solution adjustment vector. */                
    bool_t verbose         /* TRUE mumbles while working. */           
  );
  /* Similar to {SPSys_LinSolve}, but specialized to solve the residual
    system {A*dx = dy} where {dy} is a `residual' of some sort ---
    e.g. {dy = y - A*x} where {x} is the solution returned for the
    system {A*x = Y}.
    
    The main difference is in the initial guess for {x} in the
    Gauss-Seidel or Conjugate Gradient methods (zeros, instead of a
    random vector) and the number of Gauss-Seidel iterations (2
    instead of {o.maxIter}. */

void SPSys_AdjustSolution(SPVector x, SPVector dx, bool_t verbose);
  /* Adds to {x} the increment {dx}, and reports the maximum 
    change to {stderr}.
    
    Requires {dx.ne == x.ne}. If {verbose} is true, also prints the
    old and new values of {x[i]}, and the change {dx[i]}, for all {i}. */

char *SPSys_LinMethodName(SPSys_LinMethod_t mth);
  /* Returns a string that describes the method {mth}. */

char *SPSys_LinOptionsText(SPSys_LinOptions_t *lso);
  /* Returns a newly allocated string with {lso} formatted as text. */

/* SOLVING NONLINEAR SYSTEMS */

typedef enum   /* Solution methods for linear systems: */
  { SPSys_GM_NONE,
    SPSys_GM_Iterative,   /* Simple iteration. */  
  } SPSys_GenMethod_t;


typedef struct SPSys_GenOptions_t  /* Parameters for non-linear system solving. */
  { SPSys_GenMethod_t mth; /* Solution method. */
    int maxIter;           /* Maximum iterations per time step. */
    double relTol;         /* Max rel change in vars between iterations. */
    double absTol;         /* Max abs change in vars between iterations. */
  } SPSys_GenOptions_t;

#define SPSys_GenOptionsHelp \
  "  [ -nonLinearSys iterative MAXITER RELTOL ABSTOL ]"

SPSys_GenOptions_t SPSys_GenOptionsParse(SPOptions_Parser_t *pp, bool_t local);
  /* Parses a {SPSys_GenOptions_t} from the command line,
    as described by {SPSys_GenOptionsHelp}. 
    
    If {local=TRUE}, the "-nonLinear" keyword must be the next
    parameter to be parsed, otherwise it may be anywhere in the
    command line. */

void SPSys_GuessSol(SPVector x);
  /* Fills {x} with {random} values. */
  
double SPSys_UpdateSolution(SPVector xNew, SPVector x);
  /* This procedure is handy when using iterative methods like
    Gauss-Seidel or non-linear iteration. It stores a newly computed
    solution vector {xNew} into the previous solution vector {x},
    and compares the two. It returns the L2 norm of the difference. */

char *SPSys_GenMethodName(SPSys_GenMethod_t mth);
  /* Returns a string that describes the method {mth}. */

char *SPSys_GenOptionsText(SPSys_GenOptions_t *gso);
  /* Returns a newly allocated string with {gso} formatted as text. */

/* SPECIALIZED SOLVERS */

void SPSys_SVDSolve
  ( SPMatrix L,    /* Left SVD factor of system's matrix. */          
    SPVector D,    /* Diagonal factor of system's matrix. */
    SPMatrix R,    /* Right SVD factor of system's matrix. */          
    SPVector y,    /* Right-hand-side of system. */ 
    SPVector t,    /* (WORK) temporary storage area. */
    SPVector x,    /* (OUT) Coordinates of new solution relative to {bas}. */
    bool_t verbose /* TRUE mumbles along the way. */
  );
  /* Solves the system {A x == y} on the unknown {x}, for a square
    positive-definite symmetric matrix {A}, given the vector {y} and
    the left, middle and right factors {L,D,R} of the SVD
    decomposition {L*D*R^T} of the matrix {A}. The middle factor {D},
    a diagonal matrix, should be given as a vector. The solution is
    stored in to the argument vector {x}. The argument vector {t} is
    used as a temporary storage area. */
  
void SPSys_CholeskySolve
  ( SPMatrix L,      /* Cholesky factor of system's matrix. */          
    SPVector y,      /* Right-hand-side of system. */ 
    SPVector t,      /* (WORK) temporary storage area. */
    SPVector x,      /* (OUT) Coordinates of new solution relative to {bas}. */
    bool_t verbose   /* TRUE mumbles along the way. */
  );
  /* Solves the system {A x == y} on the unknown {x}, for
    a square positive-definite symmetric matrix {A}, given the
    vector {y} and the lower triangular factor {L} of the 
    Cholesky decomposition {L * L^t} of the matrix {A}.  The solution
    is stored in to the argument vector {x}.  The argument 
    vector {t} is used as a temporary storage area. */
  
void SPSys_GaussSeidelIteration(
    SPMatrix A,     /* System's matrix. */          
    SPVector y,     /* Right-hand-side of system. */ 
    double omega,   /* Overrelaxation factor. */
    SPVector xOld,  /* (IN) Current guess at the solution. */
    SPVector xNew,  /* (OUT) New solution. */
    bool_t verbose  /* TRUE mumbles along the way. */
  );
  /* Performs one iteration of the Gauss-Seidel algorithm for solving
    the system {A x == y}.  Upon entry, {xOld} should be the coordinates
    of a guess for the solution, relative to some basis.  Upon exit,
    {xNew} will have been updated by performing one iteration of the
    Gauss-Seidel algorithm. */
  
void SPSys_GaussSeidelSolve
  ( SPMatrix A,     /* System's matrix. */          
    SPVector y,     /* Right-hand-side of system. */ 
    double omega,   /* Relaxation factor. */
    int maxIter,    /* Maximum number of iterations. */
    double absTol,  /* Stopping criterion: absolute change in solution. */
    double relTol,  /* Stopping criterion: relative change in solution. */
    SPVector xOld,  /* (IN) Initial guess at solution. */
    SPVector xNew,  /* (OUT) New solution. */
    bool_t verbose  /* TRUE mumbles along the way. */
  );
  /* Attempts to solve the linear system {A x == y} by multiple
    iterations of the Gauss-Seidel method.  Upon entry, {xOld} should be
    an initial guess for the solution.  Upon exit, {xNew} will be the
    tentative solution for the system.
    
    The preocedure will stop when the L2 difference between
    successive trials is at most {absTol}, or the relative difference
    is less than {relTol}, or after at most {maxIter} iterations. */
  
void SPSys_ConjugateGradientSolve
  ( SPMatrix A,     /* System's matrix. */          
    SPVector y,     /* Right-hand-side of system. */ 
    SPVector xOld,  /* (IN) Current solution. */
    SPVector xNew,  /* (OUT) New solution. */
    bool_t verbose  /* TRUE mumbles along the way. */
  );
  /* Attempts to solve the linear system {A x == y} by conjugate
    gradient minimization of the squared norm of the residual vector.
    Upon entry, {xOld} should be the a guess for the solution. Upon
    exit, {xNew} will be a tentative solution for the system. */

bool_t SPSys_GenStopCondition
  ( SPSys_GenOptions_t *gso, /* Non-linear system method and parameters. */
    double stopError,        /* Desired max distance from true solution. */
    double cNorm,            /* Norm (L_inf) of current coeff vector. */
    int nIter,               /* Number of global iterations already performed. */
    double error,            /* Distance from true solution. */
    double cDiff,            /* Distance (L_inf) between coeff vectors, curr and prev. */
    bool_t verbose           /* If TRUE, prints the stopping reason to {stderr}. */
  );
  /* Checks whether the outer (non-linear) iteration should be stopped.
    Returns TRUE iff any of these conditions hold:

      (0) The iteration method is null, i.e. {gso.mth = SPSys_GM_NONE}.

      (1) Too many global iterations, i.e. {nIter >= gso.maxIter}.

      (2) We are close enough to the true solution, i.e. {error <= stopError}.

      (3) There was little change in the last iteration, i.e.
           {cDiff < gso->absTol} or {cDiff < gso->relTol*cNorm}
           (provided {nIter > 0}). 
  */

#endif
