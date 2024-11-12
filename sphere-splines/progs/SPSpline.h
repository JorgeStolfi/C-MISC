/* SPSpline.h -- Piecewise functions on the sphere */
/* Last edited on 2006-02-28 11:50:54 by stolfi */

#ifndef SPSpline_H
#define SPSpline_H

#include <SPTriang.h>
#include <SPFunction.h>
#include <SPH3.h>
#include <SPBasic.h>

#include <r3.h>
#include <vec.h>
#include <bool.h>
#include <nat.h>

#include <stdio.h>
#include <stdlib.h>

/* PIECEWISE SPHERICAL FUNCTION OBJECT */

/* An {SPSpline} object describes a real-valued function defined 
  on {S^3} in a piecewise fashion, namely a different spherical function
  for each face of a certain triangulation of the sphere. */

typedef struct PieceData /* Data for one face of the triangulation */
  { FaceNumber face;   /* Number of triangle in triangulation. */
    SPFunction *func;  /* Partial function for this triangle. */
    bool_t bary;       /* TRUE if {func->eval} work with baryc. coords. */
  } PieceData;
    
typedef PieceData* PieceDataRef;
vec_typedef(PieceDataRef_vec_t,PieceDataRef_vec,PieceDataRef);

#define OBJ void /* Actualy {SPSpline} or a subclass */

typedef FaceNumber LocateMth(OBJ *f, S2Point *p); /* Triangle containing {p}. */

typedef struct SPSpline_Methods  /* Methods for SPSpline: */
  {
    SPFunction_Methods fn;  
      /* Methods inherited from superclass. The {write} method writes the
        partial functions but only the triangulation's file name, not
        the triangulation itself. */
      
    SPFunction_WriteMth *write;       
      /* Writes the function's def to {wr}, as an {SPSpline}. */
      
  } SPSpline_Methods;

typedef struct SPSpline_Data  /* Data fields for SPSpline: */
  { SPFunction_Data fn;       /* Data fields inherited from superclass. */
    char *triFile;            /* Name of triangulation file. */
    Triangulation *tri;       /* The reference partition of the sphere. */
    PieceDataRef_vec_t pd;    /* Data for non-zero pieces, sorted by {face} num. */ 
    SPH3_Plane supp;          /* No pieces in the negative side of this plane. */
    /* Private fields: */
    FaceNumber lastLoc;       /* Cached result of last {locate} */
  } SPSpline_Data;

typedef struct SPSpline  /* Logically a subclass of SPFunction. */
  { char *type;           /* Type identifier */
    SPSpline_Methods *m;  /* Function methods. */
    SPSpline_Data *d;     /* Function parameters. */
  } SPSpline;
  /* The data record {*d}, the piece data vector {*(d->pd.e)},
    the piece data records {*(d->pd.e[i])}, and the partial functions
    {*(d->pd.e[i]->func)} are all private to a single {SPSpline} object,
    and are reclaimed by the {free} method --- that also calls the {free}
    method of each partial function {d->pd.e[i]->func}.
    The trinagulation {*tri} and its name {*triName} may be shared by
    several {SPSpline} objects; they are NOT reclaimed by {free}. */
  
#define SPSpline_TypeId "SF.PW."

SPSpline *SPSpline_Cast(OBJ *f);
  /* If {f} (according to {f->type}) is a subclass of {SPSpline},
    returns {f} cast to that type; otherwise returns NULL. */

bool_t SPSpline_AreCompatible(SPFunction *f, SPFunction *g, Triangulation *tri);
  /* TRUE if both functions happen to be instances of {SPSpline}
    with the same underlying triangulation, and {tri} is either NULL
    or that same triangulation. */
    
bool_t SPSpline_IsCompatible(SPFunction *f, Triangulation *tri);
  /* TRUE if function {f} happens to be an instance of {SPSpline},
    and {tri} is either NULL or is the underlying triangulation of {f}. */

/* ====================================================================== */

/* The {support} of an {SPSpline} {f} is, by definition, the set
  of faces spanned by the pieces listed in {f->d->pd}. The function
  {f} is definitely zero outside its support, and arbitrary inside it.
  In particular, it *may* be identically zero in over one or more
  faces of its support.
  
  The {supporting plane} of {f} is a plane {P = f->d->supp} such that 
  no point of {f}'s support is on the negative side of {P}. The
  supporting plane is used to speed up the evaluation of {f(p)}:
  if {p} is on the negative side of {f}, then {f(p) == 0}. */

PieceData *SPSpline_PieceData_New(void);
  /* Allocates a new {PieceData} record, un-initialized. */

SPSpline *SPSpline_Read(FILE *rd);  
  /* Reads a function's description from {rd}, in the format produced by
    the {write} method. 

    The procedure obtains the triangulation from the file whose name
    is stored in the {triFile} field, with
    {SPTriang_ReadCached(triFile,smpOrder)}; so that if two
    {SPSpline}s have the same {triFile}, they will share the
    same triangulation. In any case, the triangulation recovered by
    {SPSpline_Read} must be identical to the {tri} field at the
    time the {SPSpline} was written --- including vertex, edge,
    and face numbers, and vertex coordinates; otherwise
    {SPSpline_Read} may fail, or return a malformed object. */

SPSpline *SPSpline_FromPieces
  ( char *triFile,
    Triangulation *tri,
    PieceDataRef_vec_t pd,
    r3_t *norm
  );
  /* Creates a piecewise function {f} with triangulation {tri}, which
    is assumed to come from file {triFile}. (If {tri == NULL}, the
    triangulation is automatically read from that file, with
    {SPTriang_ReadCached}.)
    
    The function {f} will be zero except in the faces of {tri} whose
    numbers are listed in {pd.e[i]->face}, for all {i}; in each of
    these triangles, it is the function {pd.e[i]->func}. If
    {pd.e[i]->bary} is true, the arguments passed to the latter will
    be corner-relative barycentric coordinates, rather than global
    Cartesian cordinates.
    
    If the {norm} vector is not null, the supporting plane of 
    the constructed function will be orthogonal to {norm}.
    Otherwise, the approximate barycenter of the faces 
    listed in {pd} will be used as {norm}. 
    
    The piece data vector {*(pd.e)}, the piece records {*(pd.e[i])},
    and the functions {*(pd.e[i]->func)} will all become private
    to the {SPSpline} object, and should not be used elsewhere.
    Moreover, all piece records and their partial functions must
    be pairwise distinct records.  The records {*tri} and {*triName}
    may be shared by several {SPSpline} objects. */

SPSpline *SPSpline_Replicate
  ( char *triFile,
    Triangulation *tri,
    bool_t bary,
    SPFunction *func
  );
  /* Creates a piecewise function {f} with triangulation {tri}, which
    is assumed to come from file {triFile}. (If {tri == NULL}, the
    triangulation is automatically read from that file, with
    {SPTriang_ReadCached}.) Each piece of this function will be a new
    {PieceData} record containing a distinct cloned copy of {func}.
    
    If {bary} is true, the arguments passed to {func} will be
    `barycentric' coordinates (relative to the corners of the face),
    rather than global Cartesian cordinates. Specifically, in a face
    with number {fn}, the barycentric coordinates of a point {p} are
    {b[0..2]} such that {p = b[0]*Dest(Onext(e)) + b[1]*Org(e) +
    b[2]*Dest(e) }, where {e = tri->side[fn]}.
    
    Since the new {SPSpline} object owns only copies of {func}, the
    latter can be reclaimed or used for other purposes after the call.
    The records {*tri} and {*triName} too may be shared by several
    {SPSpline} objects. */
    
typedef r3_t BaryCoords;       /* Triangle-relative coordinates of point */
  
R3Point SPSpline_SupportCentroid(PieceDataRef_vec_t pd, Arc_vec_t side);
  /* Computes the approximate barycenter of all vertices of the triangles
    listed in {pd}.  Note that it is a point of {R^3} not {S^2}. 
    The {side} table should map a face number to its reference arc
    (the {side} table of the triangulation). */

S2Point SPSpline_FarthestVertex(S2Point *u, PieceDataRef_vec_t pd, Arc_vec_t side);
  /* Returns the vertex {v} of the support {pd} which makes the largest
    angle with {u}.
    
    Beware that a spherical cap with apex {u} and touching {v} is
    geodesically convex only if the the angle {u,v} is at most {Pi/2}.
    Otherwise, the cap will contain all vertices of {pd}, but not
    necessarily its triangles. */
    
Sign SPSpline_PosRelPlane(PieceDataRef_vec_t pd, SPH3_Plane *Q, Arc_vec_t side);
  /* Returns 0 if the union of {pd} straddles {Q} or is empty;
    else returns the side (+1 or -1) of {Q} that contains {pd}. The
    {side} parameter should be the face-side table ({side} field) of
    the triangulation. */

r3x3_t SPSpline_SupportFrame(PieceDataRef_vec_t pd, Arc_vec_t side);
  /* Chooses an orthonormal basis {(u,v,w)} of {R^3} where {u} is
    approximately the center of the support  {pd}, {v} is parallel to its 
    longest axis, and {w} to its narrowest axis.  
    
    The three vectors are returned as the ROWS of a 3x3 matrix. The
    {side} parameter should be the face-side table ({side} field) of
    the triangulation. */
    
void FindSharedElems
  ( PieceDataRef_vec_t *pd,
    Triangulation *tri,
    nat_t *vCom, int *nvCom, 
    nat_t *eCom, int *neCom,
    nat_t *fCom, int *nfCom
  );
  /* Finds the vertices {vCom[0..nvCom-1]}, the edges {eCom[0..neCom-1]},
    and the faces {fCom[0..nfCom-1]} shared by all pieces listed in {pd}.
    The possible outcomes are:
      
      {nfCom == 1, neCom == 3, nvCom == 3} ({pd} has only one piece),
      {nfCom == 0, neCom == 1, nvCom == 2} (2 adjacent pieces),
      {nfCom == 0, neCom == 0, nvCom == 1} (>= 3 pieces around a vertex),
      {nfCom == 0, neCom == 0, nvCom == 0} (all other cases).
  */
 
SPH3_Plane SPSpline_ComputeSupp
  ( PieceDataRef_vec_t pd, 
    Arc_vec_t side, 
    r3_t *norm
  );
  /* Returns a plane {P} such that all spherical triangles {pd[i]} are
    contained in the positive half-space of {P}. The {side} table should
    map a face number to its reference arc (the {side} table of the 
    triangulation).  If the {norm} vector is not null, the plane
    will be perpendicular to it. */
    
int SPSpline_LocatePiece(SPSpline *f, S2Point *p);
  /* Finds the piece of {f->d->pd} that contains {p}. Returns a
    large number (INT_MAX) if there is no such piece. Most
    efficient if {p} is in the same triangle as the previous call. */

void SPSpline_SortPieces(PieceDataRef_vec_t pd);
  /* Sorts the pieces in the list {pd} in order of 
    increasing face number {pd[i].face}. */

Triangulation *SPSpline_CommonTriangulation(SPFunction *f, SPFunction *g);
  /* If both {f} and {g} are {SPSpline}s with the same
    triangulation {tri}, returns that triangulation.
    Else, if exactly one of {f} or {g} is an {SPSpline}
    returns its triangulation. Otherwise returns {NULL}. */
  
Triangulation *SPSpline_BasisTriangulation(Basis basis);
  /* If there is at least one element of all {basis} which is an
    {SPSpline}, and all such elements have the same
    triangulation {tri}, returns that triangulation. Otherwise returns
    {NULL}. */
    
bool_t SPSpline_WeirdSupport(SPFunction *f, Triangulation *tri);
  /* Returns TRUE iff {f} has a weird support relative to {tri};
    namely, iff {f} is NULL, is not an {SPSpline}, or its
    triangulation is not {tri}. */

bool_t SPSpline_SameSupport
  ( SPFunction *f, 
    SPFunction *g, 
    Triangulation *tri
  );
  /* Returns TRUE if {f} and {g} are both weird, relative to {tri}, or 
    are both are non-weird and their piece lists span the same
    set of faces of {tri}. */
  
bool_t SPSpline_SubSupport
  ( SPFunction *f, 
    SPFunction *g, 
    Triangulation *tri
  );
  /* Returns TRUE if {g} is weird, relative to {tri}; 
    or both functions are non-weird, and the support of {f} span 
    a subset of the faces spanned by the pieces of {g}. */
    
bool_t SPSpline_StrictSubSupport
  ( SPFunction *f, 
    SPFunction *g, 
    Triangulation *tri
  );
  /* Returns TRUE if {f} is non-weird and {g} is weird, relative to {tri}; 
    or both are non-weird, and the pieces of {f} span a strict subset
    of the faces spanned by the pieces of {g}. */
  
bool_t SPSpline_NestedSupports(SPFunction *f, SPFunction *g, Triangulation *tri);
  /* Returns TRUE if either {f} or {g} is weird, relative to {tri}; 
    or the pieces of {f} are a subset or a superset of the pieces of {g}. */
 
/* INTEGRALS AND DOT PRODUCTS OF PW FUNCTIONS */

/* The procedures in this section are analogous to their homonyms
  in {SPFunction.h}, but use the fact that the function is
  defined piecewise to speed up the computation, by skipping 
  parts where the function is zero. Also, by integrating
  with the same number of samples in each triangle, they 
  automatically adjust the accuracy to the complexity of the function.
    
  Specifically, these procedures use {SPIntegral_BySamples} on each
  triangle of the function's support, with the sample points and weights 
  stored in the triangle's {Face} record.
  
  When the weight function {w} is NULL, all these procedures assume
  that {w(p) == 1.0} for all {p}. If {verbose} is true, they print
  something every time a numerical triangle integration is needed.
  
  */

double SPSpline_IntegralPW
  ( SPSpline *f,
    FuncMap FMap, 
    SPFunction *w, 
    bool_t verbose
  );
  /* Computes the integral of {f} modified by {map} and weighted by
    {w}, that is the integral of {w(p)*FMap(f(p),p)}, over the whole
    sphere. 
    
    If {FMap == NULL}, assumes the identity map, i.e. {FMap(u,p) == u}
    for any {u,p}. Otherwise, requires {FMap.zeroPres == TRUE}. */

R3Point SPSpline_CentroidPW
  ( SPSpline *f,
    FuncMap FMap, 
    SPFunction *w, 
    bool_t verbose
  );
  /* The approximate weighted centroid of {F(p) = FMap(f(p),p)}, namely 
    {b} such that {b.c[i] == Int(p -> p.c[i]*F(p))}, for {i=0,1,2},
    where {Int(f) = SPSpline_IntegralPW(f,FMap,w)}.
    Note that the result is not an S2Point, and is not a true average
    unless {Int(f) == 1}.  
    
    If {FMap == NULL}, assumes the identity map, i.e. {FMap(u,p) == u}
    for any {u,p}. Otherwise, requires {FMap.zeroPres == TRUE}. */

double SPSpline_DotBothPW
  ( SPSpline *f,
    FuncMap FMap, 
    SPSpline *g, 
    FuncMap GMap,
    SPFunction *w, 
    bool_t verbose
  );
  /* Computes the dot product of {f} and {g} weighted by {w} and
    modified by {FMap} and {GMap} --- that is, the integral of
    {w(p)*FMap(f(p),p)*GMap(g(p),p)} over the whole sphere. 
    
    If {FMap} is NULL, then assumes the identity map, i.e. 
    {FMap(u,p) = u} for any {u,p}; otherwise, requires 
    {FMap.zeroPres == TRUE}. Requires the same condition from
    {GMap}. The procedure also requires that {f.tri == g.tri}. */

double SPSpline_DotSinglePW
  ( SPSpline *f,
    FuncMap FMap, 
    SPFunction *g, 
    FuncMap GMap,
    SPFunction *w, 
    bool_t verbose
  );
  /* Similar to {SPSpline_DotBothPW}, except that {g} may be any 
    {SPFunction}, and {GMap} needs not be zero-preserving. */
    
double SPSpline_VelSGrdDotBothPW
  ( SPSpline *f,
    SPSpline *g, 
    SPFunction *w, 
    bool_t verbose
  );
  /* Computes the weighted dot product of {v ¤ Grad(f)} and {g}, that
    is the integral of {w(p) * r3_dot(v(p),SGrd(f)(p)) * g(p))} over
    the whole sphere, where {v} is the rigid rotation velocity field
    with unit angular speed. (If {w == NULL}, assumes {w(p) == 1.0}
    for all {p}.) The procedure requires that {f.tri == g.tri}. */

double SPSpline_VelSGrdDotSinglePW
  ( SPSpline *f, 
    SPFunction *g, 
    SPFunction *w, 
    bool_t verbose
  );
  /* Similar to {SPSpline_VelSGrdDotBothPW}, except that {g} may be any 
    {SPFunction}. */


double SPSpline_SLapDotBothPW
  ( SPSpline *f,
    SPSpline *g, 
    SPFunction *w, 
    bool_t verbose
  );
  /* Computes the weighted spherical integral of {SLap(f)} and {g},
    which by Green's theorem is equal to the spherical integral of
    {w(p) * dot(SGrd(f)(p),SGrd(g)(p))}. (If {w == NULL}, assumes
    {w(p) == 1.0} for all {p}.) The procedure requires that {f.tri ==
    g.tri}. */

double SPSpline_SLapDotSinglePW
  ( SPSpline *f, 
    SPFunction *g, 
    SPFunction *w, 
    bool_t verbose
  );
  /* Similar to {SPSpline_SLapDotBothPW}, except that {g} may be any 
    {SPFunction}. */

/* SORTING AND ORTHOGONALIZING PW BASES */

/* For the procedures in this section, the `support' of a piecewise
  function {f} (an {SPSpline}) defined over the triangulation
  {tri} is the set of triangles of {tri} that appear in its piece list
  {f->d->pd} (including those pieces {p} which are listed but
  where the function happens to be zero).
  
  Both sorting procedures below produce an order that is compatible
  with the containment order of the domains.  That is, if the support
  of {bas[i]} is contained but not identical to that of {bas[j]},
  then {i < j}. Also, all elements with the same support will
  end up in consecutive positions.

  `Weird' basis elements, which are not piecewise functions
  ({SPSpline}s), or which do not use the triangulation {tri},
  are placed after all `normal' (non-weird) elements -- as if their support
  consisted of a face numbered higher than any face of {tri}. */

void SPSpline_LexSortBasis
  ( Triangulation *tri, 
    Basis bas, 
    int ini, int lim,
    bool_t verbose
  );
  /* The procedure sorts the elements {bas[ini..lim-1]} of a spline basis
    in reverse lexicographic order of the set of triangles of {tri} spanned
    by them.  
    
    Let {fn} be the largest face number occurring among the pieces of
    {f}, and let {gn} be the same for {g}. Then {f} comes before {g}
    if {fn < gn}, and after {g} if {fn > gn}. If {fn == gn}, then the
    comparison is repeaetd with that piece excluded from both.
    
    If {verbose} is true, prints the index into the sorted basis
    of the first element with each triangle set. */

void SPSpline_BSPSortBasis
  ( Triangulation *tri, 
    Basis bas, 
    int ini, int lim,
    bool_t verbose
  );
  /* Permutes the elements {bas[ini..lim-1]} in an attempt to 
    minimize fill-in of the Cholesky factors of the dot product
    matrices.
    
    Recursively, the sorted basis consists of three consecutive
    segments {N}, {P}, and {Z} such that {Z} is asymptotically 
    small --- {O(sqrt(n))} for a triangulation with {n} faces ---
    and no element from {N} overlaps any element from {P}.
    This partition is found by selecting a plane {Q} that roughly
    bisects the union of the supports, then placing into {Z}
    those elements whose supports straddle {Q}, and the other
    elements into {N} or {P} depending on which half-space 
    of {Q} they lie.
    
    If {verbose} is true, prints the index into the sorted basis
    of the first element with each triangle set. */
    
void SPSpline_GDistSortBasis
  ( Triangulation *tri, 
    Basis bas, 
    int ini, int lim,
    bool_t verbose
  );
  /* Permutes the elements {bas[ini..lim-1]} in an attempt to 
    minimize fill-in of the Cholesky factors of the dot product
    matrices.
    
    Recursively, the sorted basis consists of three consecutive
    segments {N}, {P}, and {Z} such that {Z} is asymptotically 
    small --- {O(sqrt(n))} for a triangulation with {n} faces ---
    and no element from {N} overlaps any element from {P}.
    This partition is found by collecting the set {F} of all faces that 
    occur in the support of any basis element,
    dividing {F} in two hopefully compact halves, `low' and `high',
    then classifying eac basis element as {N}, {P}, or {Z}
    according to whether its support is all-low, all-high,
    or mixed, respectively. 
    
    If {verbose} is true, prints the index into the sorted basis
    of the first element with each triangle set. */

void SPSpline_CheckBasisOrder
  ( Triangulation *tri, 
    Basis bas, 
    int ini, int lim
  );
  /* Checks whether the order of the basis elements is consistent
    with containment order, i.e. if the support of {f} is a proper
    subset of that of {g} then {f} must come before {g} in the basis.
    Caution: time {O(NB^2)} time for {NB}elements. */ 

void SPSpline_FillVertexSupport(SPSpline *f);
  /* If the support of {f} consists of a single triangle, or two
    triangles sharing an edge, does nothing. Otherwise, if the support
    triangles are all incident to a common vertex, then appends to
    {f->d->pd} any additional pieces whose triangles are also incident
    to that vertex --- so that the support of {f} becomes the entire
    star of that vertex.
    
    The procedure fails if the triangles that comprise the support of
    {f} have no common vertex.  
    
    Partial functions of new pieces are obtained by cloning the
    {func} of some existing piece and scaling the clones by 0. The
    {bary} attribute is copied too. */

void SPSpline_DescribeSupport(FILE *wr, SPFunction *f, Triangulation *tri);
  /* Prints to file {wr} a description of the support of {f}. If {f} is
    piecewise, and uses the triangulation {tri}, prints the list of
    faces in its piece list. If, moreover, that list consists of two
    faces sharing an edge, or three or more faces sharing a common
    vertex, also prints that shared element. */

void SPSpline_DescribeSupports
  ( FILE *wr, 
    Triangulation *tri, 
    Basis bas, 
    int ini, int lim
  );
  /* Prints to file {wr} a description of the supports of elements
    {bas[ini..lim-1], condensing consecutive elements which have the
    same support description. See }SPSpline_DescribeSupport{. */

#endif
