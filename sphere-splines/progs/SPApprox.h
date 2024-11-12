/* SPApprox.h -- tools for least-squares approximation */
/* Last edited on 2005-10-29 05:34:04 by stolfi */

#ifndef SPApprox_H
#define SPApprox_H

#include <SPMatrix.h>
#include <SPVector.h>
#include <SPSys.h>
#include <SPFunction.h>
#include <SPTriang.h>
#include <SPPlot.h>
#include <SPH3.h> 
#include <SPBasic.h> 
#include <vec.h> 
#include <js.h> 

/* These procedures are mostly wrappers for procedures in
  other interfaces, such as {SPMatrix.h}, {SPFunction.h}, 
  with some additional facilities for measuring, printing,
  and plotting the approximation errors. */ 

SPFunction *SPApprox_ReadFunction(char *fileName);
  /* Reads a spherical function from the given file. Bombs if no such file.*/

Basis SPApprox_ReadBasis(char *fileName);
  /* Reads a spherical function basis from the given file. Bombs if no such file.*/

SPMatrix SPApprox_ReadMatrix(char *fileName);
  /* Reads an {SPMatrix} from file {fileName}. Bombs if no such file. */

SPVector SPApprox_ReadVector(char *fileName);
  /* Reads a {SPVector} from file {fileName}. Bombs if no such file. */

void SPApprox_GetBasisMatrices
  ( char *matName,
    SPSys_LinMethod_t mth,
    SPMatrix *E,
    SPMatrix *EL,
    SPVector *ED,
    SPMatrix *ER
  );
  /* Obtains the rigidity matrix {E[i,j] == <B[j]|B[i]>} for some
    spherical function basis {B} The matrix {E} is read from the file
    "{M}-ev.mat".
    
    Also obtains the triangular factors {EL} and/or {ED} and/or {ER} of {E}, as
    needed by {SPSys_LinSolve} with method {mth}. The diagonal matrix {ED} 
    is stored as a vector.
    
    Specifically, if {mth == SPSys_LM_Cholesky}, reads only the
    Cholesky lower factor {EL} of {E} from file "{M}-ev-choL.mat" where
    {M = matName}. If {mth == SPSys_LM_SVD}, reads the SVD factors
    {EL}, {ED} and {ER} of {E} from files "{M}-ev-svdL.mat" and
    "{M}-ev-svdD.vec", and "{M}-ev-svdR.mat", respectively. Matrices that
    are not read are set to the null matrix of the correct size. */
    
void SPApprox_GetHelmholtzMatrices
  ( char *matName, 
    double c, 
    SPSys_LinMethod_t mth,
    double minVal,
    SPMatrix *H,
    SPMatrix *HL,
    SPVector *HD,
    SPMatrix *HR
  );
  /* Obtains the differential operator matrix {H[i,j] =
    <D(B[j]),B[i]>}, where {B} is some spherical function basis, for
    the Helmholtz operator {D(f)(p) = SLap(f)(p) + c*f(p)} with given
    coefficient {c}.
        
    The matrix {H} is computed according to the formula {H[i,j] =
    <SGrd(B[j]) | SGrd(B[i])> + c*<B[j]|B[i]>}, where {SGrd} is the
    spherical gradient operator. The matrices {G[i,j] = <SGrd(B[j]) |
    SGrd(B[i])>} and {E[i,j] = <B[j]|B[i]>} are read from files
    "{M}-gr.mat" and "{M}-ev.mat", where {M = matName}.
    
    The procedure computes also the factors {HL} and/or {HD} and/or
    {HR}, as needed by {SPSys_LinSolve} with method {mth}. Factors
    that are not needed are set to the null matrix of the correct
    size. See {SPSys_ComputeNeededFactors} for the meaning of {minVal}. */

SPMatrix SPApprox_GetTransferMatrix(char *matName);
  /* Obtains the basis transfer matrix {F[i,j] == <B1[j]|B2[i]>} for 
    two spherical function bases {B1,B2}. The matrix is read from the file
    {M-ev.mat}, where {M} is the given {matName}. */

void SPApprox_ComputeSystemRightHandSide
  ( SPFunction *f,      /* Function to approximate. */
    FuncMap FMap,       /* Right-hand side functional. */
    Basis bas,          /* Basis for the approximation space. */
    SPFunction *w,      /* Weight for integrals. */
    Triangulation *tri, /* Triangulation to use for integration. */
    SPVector y,         /* IN/OUT: right-hand-side vector. */ 
    bool_t verbose        /* TRUE mumbles and prints the change in {y}. */
  );
  /* (Re)computes the vector {y[i] == <F|bas[i]>}, where 
    {F(p) = FMap(f(p), p)}.  In particular, if {FMap == NULL},
    reduces to {y[i] == <f|bas[i]>}.  
    
    If {verbose} is true, also computes and prints the distance between the newly
    computed value of {y} and its input value.
    
    The vector {y} is the right-hand side of the linear system
    {A*x==y}, which occurs in the approximation of {F} by {bas} (where
    {A[i,j] = <bas[j]|bas[i]>}), or in the solution of the
    differential equation {D(f)(p) == F(p)} (where {A[i,j] ==
    <D(bas[j])|bas[i]>}. */
    
SPFunction *SPApprox_BuildFunction
  ( Basis basis,     /* Basis for space {V[r]}. */
    SPVector u,      /* Coordinates of {g} relative to {basis} */
    char *solName,   /* Name prefix for output file. */
    SPFunction *f    /* True solution (for comparisons) */
  );
  /* Returns the function {g} with coefficients {u} relative 
    to the given {basis}. Also writes {g} to disk (with the 
    given {solName} extended with "-app.sfn") and compares it
    to the given function {f}, with {SPApprox_PrintSomeValues}. */

void SPApprox_PlotFunctionAndError
  ( SPPlot_Stream *fps,      /* Postscript stream. */
    char *funcTag,           /* Prefix for page names. */
    SPFunction *g,           /* The computed approximation. */
    double f(S2Point *p),    /* The target function. */
    Triangulation *tri,      /* Reference triangulation. */
    double relMeshSize,      /* Maximum step/triangle size (radians). */
    bool_t showTriang,       /* TRUE to plot the reference triangulation. */ 
    bool_t plotTrueSol,      /* TRUE to plot the target function {f}. */ 
    double fMax,             /* Nominal maximum of {fabs(g(p))} and {fabs(f(p))}. */
    double eMax,             /* Nominal maximum of {fabs(g(p)-f(p))}. */
    SPH3_Point *obs,         /* Viewpoint. */
    SPH3_Point *upp,         /* Camera vertical reference */
    string_vec_t *caption,   /* Caption template. */
    int index,               /* Figure `index' for caption expansion. */
    double time              /* Figure `time' for caption expansion. */
  );
  /* Plots the function {g}, and the difference {g-f}, using
    {SPPlot_BothSides} (q.v.). 
    
    The file {fps} may be NULL, or may contain other figures; the
    procedure calls {SPPlot_BeginFigure} internally. This may cause
    the file {fps} to be closed, and another one opened in its place;
    the final file, still open, is returned as result. Be sure to call
    {SPPlot_Close} on it eventually.
    
    The {funcTag} parameter gets extended with {-app} and {-err} for
    the function and the error plots, respectively (to which
    {SPPlot_MultipleViews} will append {-f} for front or {-r} for back). */

void SPApprox_PrintSomeValues(SPFunction *g, SPFunction *f);
  /* Prints values of {f(p)}, {g(p)}, and their difference 
    at a few points (mostly the octahedron vertices, edge
    centers, and face centers). */

void SPApprox_PrintMaxAvgValues
  ( SPFunction *f,      /* A function. */
    Triangulation *tri, /* Triangulation to use for integration. */
    double *fMax,       /* OUT: max absolute value of {f}. */
    double *fAvg        /* OUT: mean value of {f}. */
  );
  /* Computes and prints the maximum absolute value {fMax}
    and the average value {fAvg} of a function {f}. */

void SPApprox_PrintMaxErrorValues
  ( SPFunction *g,      /* Approximation. */
    SPFunction *f,      /* Target function. */
    Triangulation *tri, /* Triangulation to use for integration. */
    double *gMax,       /* OUT: max abs value of {g}. */
    double *fMax,       /* OUT: max abs value of {f}. */
    double *eMax,       /* OUT: max abs error VALUE. */
    double *eAvg        /* OUT: root mean square error {g-f}. */
  );
  /* Computes and prints the maximum absolute values {gMax,fMax} of
    functions {g} and {f}, and the maximum and root-mean-square
    error {eMax,eAvg} when compared to {f}. */

#endif
