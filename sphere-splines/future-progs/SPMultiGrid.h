/* SPMultiGrid.h -- Multigrid structure for point location and approximation. */
/* Last edited on 2002-12-09 13:31:53 by anamaria */

#ifndef SPMultiGrid_H
#define SPMultiGrid_H

#include <SPTriang.h>
#include <SPOverlapTable.h>
#include <SPMatrix.h>
#include <SPSpline.h>
#include <SPTriang.h>
#include <basic.h>
#include <arrays.h>
#include <stdio.h>

typedef struct GridData /*  Data for some scale. */
  { Triangulation *tri;   /* The triangular mesh for some scale. */
    char *triFile;        /* Grid filenames. */
    Basis bas;            /* A basis for the spline space {V} on {tri}. */
    ovl: OverlapTable;    /* Overlap table for coarser grid. */
  }
  
typedef struct GridDataArray { nat nel; GridData *el; } GridDataArray;

typedef GridDataArray MultiGrid; 
  /* A multigrid structure {m} consists of {NL+1} grids, {m[0..NL]},
    supposedly in order of increasing resolution. 
    
    The table {m[k].ovl} lists the triangles of grid {m.tri[k]} that
    overlap each triangle of grid {m.tri[k-1]}. The virtual grid
    {m[-1]} is the whole sphere, so {m[0].ovl} has a single row,
    {m[0].ovl[0]}, which lists all triangles of triangulation
    {m[0].tri}.

    Each level {k} of the multigrid has an associated space {V[k]}
    of spherical spline functions, defined by the basis {B[k] = m[k].bas}.

    The sparse matrix {m[k].GL} is the lower triangular Choleski
    decomposition of the rigidity matrix {G[k]} for basis {B[k]}.
    Element {[i,j]} of {G[k]}, by definition, is the dot product
    of basis functions {B[k][i]} and {B[k][j]}.  

    The sparse matrix {m.F[k]} is the connection matrix between
    spline spaces of layers {k-1} and {k}. Specifically, element
    {[i,j]} of {m.F[k]} is the dot product of {B[k-1][i]} and
    {B[k][j]}.

    The matrix {m.HL[k]} is used to solve a certain PDE 
    {D(f)(p) = FMap(f(p),p)} on the sphere. It
    is the lower triangular Choleski decomposition of the differential
    operator matrix {H[k]}, whose element {[i,j]}
    . In particular, for the Helmholtz equation
    {SLap(f)(p) + c*f(p) == F(f(p),p)}, we have {H[k] == S[k] +
    c*G[k]}, where {G[k]} is the rigidity matrix of basis {B[k]},
    defined above, and {S[k]} is the rigdity matrix of the basis
    gradients. That is, element {[i,j]} of {S[k]} is the dot product
    of {SGrd(B[k,i])} and {SGrd(B[k,j])}, where {SGrd} denotes the
    spherical gradient operator. */

typedef int Level;  /* Index {[0..m.NL]} of a layer in a multigrid {m}. */

FaceNum SPMultiGrid_ReLocate(MultiGrid *m, Level k, FaceNum fS, S2Point *p);
  /* Given a point {p} in triangle {fS} in layer {S = m.tri[k]} 
    of multigrid {m}, returns the number of the triangle in layer
    {T = m.tri[k+1]} that contains {p}. */
    
FaceNum SPMultiGrid_Locate(MultiGrid *m, S2Point *p, Level k);
  /* Returns the number of the triangle in layer {k} of {m} that contains {p}. */

void SPMultiGrid_Write(FILE *wr, MultiGrid *m);
  /* Writes a multigrid structure to file {wr}. */
    
MultiGrid Read(FILE *rd);
  /* Read a multigrid structure from file {rd}. */

#endif
