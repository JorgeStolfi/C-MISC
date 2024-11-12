/* SPBezSplineBasis.h -- The ANS basis for spherical polyn. splines */
/* Last edited on 2005-06-06 11:44:22 by stolfi */

#ifndef SPBezSplineBasis_H
#define SPBezSplineBasis_H

#include <SPTriang.h>
#include <SPFunction.h>
#include <SPSpline.h>
#include <SPDeCasteljau.h>
#include <SPBasic.h>
#include <vec.h>
#include <js.h>

Basis SPBezSplineBasis_BuildH
  ( int deg,             /* Degree of polynomials */
    int cont,            /* Continuity. */
    bool_t newStyle,       /* FALSE: old ANS, TRUE: with boat elements. */
    char *triFile,       /* Triangulation file name. */
    Triangulation *tri   /* Triangulation. */
  );
  /* Builds the ANS basis for the space of homogeneous spherical
    polynomial splines of degree {deg} and continuity {cont}, for the
    given triangulation, in the Bezier representations.
    
    If {cont == 1}, the {newStyle} parameter selects between the original
    ANS basis (FALSE, with {m+3} star elements at each vertex of order {m})
    and the improved version (TRUE, with {m} boat elements + 3 star
    elements).
    
    At present, the continuity must be either 0 (with degree in
    [2..5]??) or 1 (with degree [5..7]??). */

Basis SPBezSplineBasis_BuildNH
  ( int deg,             /* Degree of polynomials */
    int cont,            /* Continuity. */
    bool_t newStyle,       /* FALSE: old ANS, TRUE: with boat elements. */
    char *triFile,       /* Triangulation file name. */
    Triangulation *tri   /* Triangulation. */
  );
  /* Builds a basis for the space of non-homogeneous spherical splines on
    triangulation {tri}, of degree {deg} with specified continuity.
    
    The {newStyle} parameter has the same meaning as in
    {SPBezSplineBasis_BuildH}.
    
    At present, the degree must be either 0 (with degree in [3..5])
    or 1 (with degree [6..7]).  */

int HSpaceDimension(int deg, int cont, int NT);
  /* Dimension of the homogeneous spline space of degree {deg} 
    and continuity {cont}, in a triangulation with {NT}
    triangles without degeneracies (two coplanar edges
    incident to the same vertex). */

#define Spline SPSpline
#define HBezFn SPHBezFunction
#define NHBezFn SPNHBezFunction

PieceDataRef_vec_t MakeHBezPieces(int n, int deg);
  /* Creates {n} partial function records, whose {func} field is a
    homogeneous polynomial of degree {deg} in Bezier form (an
    {SPHBezFunction}), with all coefficients set to zero. */

PieceData *MakeHBezPiece(int deg);
  /* Creates one partial function record, whose {func} field is a
    homogeneous polynomial of degree {deg} in Bezier form (an
    {SPHBezFunction}), with all coefficients set to zero. */

PieceDataRef_vec_t MakeNHBezPieces(int n, int deg);
  /* Creates {n} partial function records, whose {func} field is a
    non-homogeneous polynomial of degree {deg} (an {SPNHBezFunction}),
    with all coefficients set to zero. */

void BuildFaceBasis
  ( Basis *bas,
    int *nb,
    char *triFile,
    Triangulation *tri,
    int deg,
    int cont,
    nat_vec_t fBas
  );
  /* Builds the basis elements associated with each face. Appends them
    to the basis {bas} starting at {bas[nb]}. Sets {fBas[i]} to the
    index of the first element associated with each face {i}. Also
    sets {fBas[NF}} to the index of the first element of {bas}
    following all face elements. */

void BuildEdgeBasis
  ( Basis *bas,
    int *nb,
    char *triFile,
    Triangulation *tri,
    int deg,
    int cont,
    bool_t newStyle,
    nat_vec_t eBas
  );
  /* Builds the homogeneous basis elements associated with each edge.
    Appends them to the basis {bas} starting at {bas[nb]}. Sets {eBas[i]}
    to the index of the first element associated with edge number {i}.
    Also sets {eBas[NE}} to the index of the first element of {bas}
    following all edge elements. */
  
void BuildVertBasis
  ( Basis *bas,
    int *nb,
    char *triFile,
    Triangulation *tri,
    int deg,
    int cont,
    bool_t newStyle,
    nat_vec_t vBas
  );
  /* Builds the basis elements associated with each vertex. Appends
    them to the basis {bas} starting at {bas[nb]}. Sets {vBas[i]} to
    the index of the first element associated with vertex number {i}.
    Also sets {vBas[NV}} to the index of the first element of {bas}
    following all vertex elements. */

void SetCoeffRelToEdge
  ( Triangulation *tri,  /* Triangulation. */
    int deg,             /* Degree of polynomial. */
    BezCoeff_vec_t c,    /* Coeffs vector for some face. */
    Arc e,               /* An edge of the face in question. */
    int i, int j, int k, /* Bezier coeff label relative to {e}. */
    double v             /* Coefficient value. */
  );
  /* Sets a Bezier coefficient {c_{i,j,k}} for the face {Left(e)}
    to the value {v}, where the indices {i,j,k} are given relative
    to the edge {e} --- which may not be the reference edge of 
    that triangle. */

void SetCoeffsVertex
  ( Arc e,                 /* Reference edge for basis element. */
    Triangulation *tri,    /* Triangulation. */                   
    int deg,               /* Degree of polynomial. */            
    BezCoeff_vec_t c,      /* Bezier coefficients. */             
    double *v,             /* Values to set. */                   
    int nv                 /* Number of values to set. */         
  );
  /* Stores into {c} the vertex-associated Bezier coefficients {v[k]},
    relative to the edge {e}, which is a side of the corresponding
    triangle, but not necessarily its reference edge.

    Assumes that {v[0]} is the coefficient associated to the node {Org(e)},
    {v[1]} and {v[3]} (if they exist) lie on {e}, {v[2]} and {v[5]}
    lie on {Onext(e)}, and {v[4]} between {v[3]}
    and {v[5]}. */

bool_t AllCoeffsAreZero(double *v, int nv);
  /* TRUE iff the coefficients {v[0..nv-1]} are all zero. */

void StoreBasisElem
  ( Basis *bas, 
    int *nb, 
    char *triFile, Triangulation *tri, 
    PieceDataRef_vec_t pd,
    r3_t *norm
  );
  /* Builds a piecewise function from the given pieces, stores 
     it in {bas.e[*nb]}, and increments {*nb}. The {norm} argument
     is passed to {SPSpline_ComputeSupp} (q.v.). */
     
#endif
