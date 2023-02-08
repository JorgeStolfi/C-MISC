/* SOBasisMatrices.h -- Bases for functions */
/* Last edited on 2007-01-04 00:21:40 by stolfi */

#ifndef SOBasisMatrices_H
#define SOBasisMatrices_H

#include <SOFunction.h>
#include <SOGrid.h>
#include <SOMatrix.h>

//#include <SOFunctionOrtho.h>

#include <vec.h>

/* BASES FOR FUNCTION SPACES */

SOMatrix SOBasisMatrices_BuildMatrixGen
  ( Basis F,
    Basis G, 
    void dot(SOFunction *f, SOFunction *g, double *iv), 
    double minVal,
    bool_t symmetric,
    bool_t verbose
  );
  /* Builds a general dot product matrix {m} for bases {F} and {G}, i.e. 
    {m[i,j] == dot(F[i], G[j])}.
    
    If {minVal} is positive, entries {m[i,j]} which are less than
    {minVal} times {sqrt(dot(F[i],F[i]) * dot(G[j],G[j]))} are set to zero.
    If {symmetric} is true, assumes that {F == G} and {dot(f,g) == dot(g,f)}.
    If {verbose} is true, prints the non-zero elements. */

void SOBasisMatrices_RecomputeVectorGen
  ( SOFunction *f,
    Basis bas,
    SOFunction_fgDot dot,
    double_vec_t y,
    bool_t verbose
  );
  /* (Re)computes the vector {y[i] == dot(f,bas[i])}. If {verbose} is
    true, also computes and prints the distance between the newly
    computed value of {y} and its input value. */


SOMatrix SOBasisMatrices_BuildMatrixEval
  ( Basis F,
    Basis G,
    SOFunction *w,
    bool_t ignoreTriang,
    double minVal,
    bool_t verbose
  );
  /* Builds the rigidity matrix {m} for bases {F} and {G}, i.e.
    {m[i,j] == SOFunction_Dot(F[i], NULL, G[j], NULL, w, NULL, verbose)}.

    If {ignoreTriang == TRUE}, uses {SOFunction_RawDot} instead.
    Parameters {minVal} and {verbose} are passed on to {BuildMatrixGen}. */

SOMatrix SOBasisMatrices_BuildMatrixGrad
  ( Basis F,
    Basis G,
    SOFunction *w,
    bool_t ignoreTriang,
    double minVal,
    bool_t verbose
  );
  /* Builds the rigidity matrix {m} for bases {F} and {G}, i.e.
    {m[i,j] == SOFunction_GradDot(F[i], G[j], w, NULL, verbose)}.

    If {ignoreTriang == TRUE}, uses {SOFunction_RawGradDot} instead.
    Parameters {minVal} and {verbose} are passed on to {BuildMatrixGen}. */

void SOBasisMatrices_ChangeBasis
  ( double_vec_t a,
    SOMatrix FG,
    SOMatrix GGL,
    double_vec_t b
  );
  /* Given the coefficients {a[0..m-1]} of a function {f(p)} in terms of
    a functional basis {F[0..m-1]}, computes the coefficients
    {b[0..n-1]} of the least squares approximation {g(p)} to {f(p)} in
    some other basis {G[0..n-1]}.  Requires the matrix
    {FG[i,j] == <F[i],G[j]>} and the lower-triangular Cholesky
    factor {GGL} of the rigidity matrix {GG[i,j] == <G[i],G[j]>}. */

#endif
