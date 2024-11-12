/* SPBasisMatrices.h -- Bases for spherical functions */
/* Last edited on 2005-10-29 05:37:03 by stolfi */

#ifndef SPBasisMatrices_H
#define SPBasisMatrices_H

#include <SPFunction.h>
#include <SPTriang.h>
#include <SPMatrix.h>
#include <SPVector.h>
#include <js.h>

/* BASES FOR SPHERICAL FUNCTION SPACES */

SPMatrix SPBasisMatrices_BuildMatrixGen
  ( Basis F,
    Basis G, 
    double dot(SPFunction *f, SPFunction *g), 
    Triangulation *tri,
    bool_t semiOrtho,
    bool_t orthogonal,
    bool_t symmetric,
    double minVal,
    bool_t verbose
  );
  /* Builds a general dot product matrix {M} for bases {F} and {G}, i.e. 
    {M[i,j] == dot(F[j], G[i])} (note the indices!).
    
    If {semiOrtho} is true, requires {F} and {G} to be the same basis,
    consisting of splines on the triangulation {tri}; and assumes that
    it has been semi-orthogonalized with respect to the {dot}
    product. That is, assumes that {dot(F[i],F[i]) = 1}, and
    {dot(F[i],F[j]) = 0} when {i != j} and the supports are nested.

    If {orthogonal} is true, assumes that {dot(F[j],G[i]) == 0} whenever {i != j}.
    
    If {symmetric} is true, assumes that {F == G} and {dot(f,g) == dot(g,f)}.
    
    If {minVal} is positive, off-diagonal entries {M[i,j]} which are less than
    {minVal} times {sqrt(dot(F[j],F[j]) * dot(G[i],G[i]))} are set to zero.
    If {verbose} is true, prints the non-zero elements. */

void SPBasisMatrices_RecomputeVectorGen
  ( SPFunction *f,
    Basis bas,
    Metric dot,
    SPVector y,
    bool_t verbose
  );
  /* (Re)computes the vector {y[i] == dot(f,bas[i])}. If {verbose} is
    true, also computes and prints the distance between the newly
    computed value of {y} and its input value. */


SPMatrix SPBasisMatrices_BuildMatrixEval
  ( Basis F,
    Basis G,
    SPFunction *w,
    Triangulation *tri,
    bool_t semiOrtho,
    bool_t orthogonal,
    bool_t ignoreTriang,
    double minVal,
    bool_t verbose
  );
  /* Builds the rigidity matrix {M} for bases {F} and {G}, i.e.
    {M[i,j] == SPFunction_Dot(F[j], NULL, G[i], NULL, w, tri, verbose)}.
    
    If {semiOrtho} is true, assumes {F} and {G} are the same 
    semi-orthogonalized spline basis with triangulation {tri};
    see {SPBasisMatrices_BuildMatrixGen} above for details.

    If {orthogonal} is true, assumes {M[i,j] == 0} whenever {i != j}.

    If {ignoreTriang} is true, uses {SPFunction_RawDot} instead of {SPFunction_Dot}.

    The parameters {minVal} and {verbose} are passed on to 
    {SPBasisMatrices_BuildMatrixGen}. */

SPMatrix SPBasisMatrices_BuildMatrixVelSGrd
  ( Basis F,
    Basis G,
    SPFunction *w,
    bool_t orthogonal,
    bool_t ignoreTriang,
    double minVal,
    bool_t verbose
  );
  /* Builds the Z-rotation flow matrix {M} for bases {F} and {G}, i.e.
    {M[i,j] == SPFunction_VelSGrdDot(F[j], G[i], w, NULL, verbose)}.

    If {orthogonal == TRUE}, assumes {M[i,j] == 0} whenever {i != j}.
    If {ignoreTriang == TRUE}, uses {SPFunction_RawVelSGrdDot} instead.
    Parameters {minVal} and {verbose} are passed on to {BuildMatrixGen}. */

SPMatrix SPBasisMatrices_BuildMatrixSLap
  ( Basis F,
    Basis G,
    SPFunction *w,
    bool_t orthogonal,
    bool_t ignoreTriang,
    double minVal,
    bool_t verbose
  );
  /* Builds the Laplacian matrix {M} for bases {F} and {G}, i.e.
    {M[i,j] == SPFunction_SLapDot(F[j], G[i], w, NULL, verbose)}.

    If {orthogonal == TRUE}, assumes {M[i,j] == 0} whenever {i != j}.
    If {ignoreTriang == TRUE}, uses {SPFunction_RawSLapDot} instead.
    Parameters {minVal} and {verbose} are passed on to {BuildMatrixGen}. */

void SPBasisMatrices_ChangeBasis
  ( SPVector a,
    SPMatrix FG,
    SPMatrix GGL,
    SPVector b
  );
  /* Given the coefficients {a[0..m-1]} of a function {f(p)} in terms of
    a functional basis {F[0..m-1]}, computes the coefficients
    {b[0..n-1]} of the least squares approximation {g(p)} to {f(p)} in
    some other basis {G[0..n-1]}.  Requires the matrix
    {FG[i,j] == <F[j],G[i]>} and the lower-triangular Cholesky
    factor {GGL} of the rigidity matrix {GG[i,j] == <G[j],G[i]>}. */

#endif
