/* SPTriang.h -- Delaunay triangulation on a sphere. */
/* Last edited on 2011-12-12 11:07:05 by stolfi */ 

#ifndef SPTriang_H
#define SPTriang_H

/* Procedures that build the Delaunay triangulation of a 
  given set of points on the sphere {S^2}.  May also be 
  interpreted as the convex hull of those points. See 

   /Primitives for the Manipulation of General Subdivisions 
   and the Computation of Voronoi Diagrams/, 
   L. Guibas, J. Stolfi, ACM TOG, April 1985

*/

#include <SPQuad.h>
#include <SPBasic.h>

#include <r3x3.h>
#include <vec.h>
#include <bool.h>
#include <nat.h>

#include <stdio.h>

typedef SPQuad_Arc Arc;
typedef SPQuad_Arc_vec_t Arc_vec_t;

/* SITES (TRIANGULATION VERTICES) */

typedef nat_t SiteNumber;
typedef struct Site /* Data for a triangulation vertex */
  { 
    SiteNumber num; /* Site number (see NumberVertices). */

    S2Point pos;    /* A point on the unit sphere. */

  } Site;
  
vec_typedef(SiteRef_vec_t,SiteRef_vec,Site *);

/* TRIANGLES (FACES OF THE TRIANGULATION) */

typedef nat_t FaceNumber;
typedef struct Face
  { 
    FaceNumber num; /* Face number. */
    
    /* The following fields are set by {SPTriang_FromTopology}: */
    
    r3x3_t b2c;
    r3x3_t c2b;
      /* {r3x3_map_col(b2c,a)} computes the Cartesian coordinates
        of a point given its `barycentric' coordinates {a},
        the coordinates relative to the corners of the triangle.
        The corners are, in order, {Dest(Onext(e))},
        {Org(e)}, {Dest(e)}, where {e == side[i]}.
        {r3x3_map_col(c2b,c)} performs the inverse transformation. */
        
    S2Point_vec_t sp;
    double_vec_t wp;
      /* The spherical integral of any function {f} defined on
        this triangle can be computed by the formula 
        {SUM { wp[i]*f(sp[i]) : i = 0..sp.ne }}. */

  } Face;
  
vec_typedef(FaceRef_vec_t,FaceRef_vec,Face *);

/* EDGES (LINES OF THE TRIANGULATION) */

typedef nat_t EdgeNumber;

/* TRIANGULATION HEADER */
  
typedef struct Triangulation 
  { Arc_vec_t out;
      /* Element {out[i]} is one arc out of site number {i}. */

    Arc_vec_t arc;
      /* Elements {arc[2*i]} and {arc[2*i+1]} are
        the same edge in both orientations. */

    Arc_vec_t side;
      /* Element {side[i]} is one arc whose left face is triangle {i}. */
      
    int smpOrder;
      /* Sampling order that was used in each triangle. */

  } Triangulation;

/* QUAD-EDGE DATA POINTERS: */

Site *Org (Arc e);
Site *Dest (Arc e);

Face *Left(Arc e);
Face *Right(Arc e);

nat_t OrgDegree(Arc e);
  /* Number of arcs in the {Onext} ring of {e}. */

nat_t LeftDegree(Arc e);
  /* Number of arcs in the {Lnext} ring of {e}. */

bool_t Adjacent(Arc e, Arc f);
  /* TRUE if sites {Org(e)} and {Org(f)} are adjacent
    (the endpoints of some edge). */

/* GEOMETRIC OPS ON SITES */

bool_t SamePoint(S2Point *p, S2Point *q);
  /* True iff {p == q}. */

bool_t PositiveSide(R3Point *a, R3Point *b, R3Point *c, R3Point *d);
  /* TRUE if {d} is on the positive side of plane {Join(a,b,c)}. */

bool_t TrianglesOverlap
  ( S2Point *u, S2Point *v, S2Point *w, 
    S2Point *p, S2Point *q, S2Point *r
  );
  /* Returns TRUE if the spherical triangles {u,v,w} and {p,q,r} 
    have a non-trivial intersection. Assumes they are counterclockwise 
    oriented (i.e. {r3_det(u,v,w) > 0}, {r3_det(p,q,r) > 0} */ 

void SPTriang_ConvexCone(int *N, R3Point *v);
  /* Computes the strictly convex enclosing cone of the vectors {v[0..N-1]}.
    Upon exit, {v[0..N-1]} will be a subset of the input vectors,
    which are the extremal vectors of the cone, in counterclockwise
    order around any internal vector. If there is no strictly 
    convex cone enclosing the input vectors (that is, if they
    are not strictly contained in any half-space), returns {N == 0}. */

double TriDist(R3Point *p, r3x3_t *c2b);
  /* Computes a distance metric from {p} to the trihedron
    of the closed spherical triangle whose Cartesian-to-relative matrix 
    is {c2b}. The distance is 0 inside or on the boundary of the 
    triangle, and positive outside. */

double TriangleArea(R3Point *u, R3Point *v, R3Point *w);
  /* Area of the plane triangle {u,v,w}. */
  
double SphericalTriangleArea(S2Point *u, S2Point *v, S2Point *w);
  /* Area of the spherical triangle {u,v,w}. */
  
/* TRIANGULATION OPS */

/* The procedures in this section use only the quad-edge structure and
  the {Site} coordinates. The structure is interpreted either as a
  spherical triangulation or as a convex polyhedron. 
  
  These procedures do not use or create the {Face} records, and,
  except where stated, do not expect edges and sites to have
  meaningful {num} fields. When the structure is finished, use
  {FromTopology} below to creathe these missing fields.
  
  Beware that some of these procedures may produce a triangulation
  that is topologically or geometrically inconsistent. */

Arc SPTriang_Locate(R3Point *p, Arc e);
  /* Returns an edge {f} of the polyhedron hanging from {e} such that
    the face {Left(f)} is visible from the point {p} (which is assumed
    to be outside the polyhedron). If the point {p} coincides with an
    existing vertex, returns an edge out of that vertex. */

Arc SPTriang_Connect(Arc a, Arc b);
  /* Adds a new edge connecting {Dest(a)} to {Org(b)} across the face
    {Left(a)}, which must be the same as {Left(b)}. */

void SPTriang_Disconnect(Arc a);
  /* Disconnects the edge {a} from the triangulation,
    leaving it as a sphere with a single non-loop edge. */

void SPTriang_SwapEdge(Arc e);
  /* Let {a == Oprev(e)}, {b == Oprev(Sym(e))}. Disconnects the edge
    {e} from the triangulation, creating a four-sided face, and
    reconnects {e} as the other diagonal of that face, directed from
    {Dest(a)} to {Dest(b)}. The face records {Left(e)} and {Right(e)},
    if any, are preserved. The links {Org(e)}, {Dest(e)}, {Left(a)},
    {Left(b)} are properly updated. */

void SPTriang_InsertSiteInFace(Arc e, Site *x);
  /* Replaces the face {Left(e)}, which may have any number of sides,
     by a wheel of triangles surrounding site {x}. */

void SPTriang_RemoveSite (Arc b);
  /* Removes site {Dest(b)} and all its incident edges. */

bool_t SPTriang_EdgeIsDegnerate(Arc e);
  /* TRUE if {Left(e) == Right(e)}, or any of those two faces is a
    degenerate triangle (with a loop edge), or if their vertices
    opposite to {e} coincide. */

Sign SPTriang_FaceOrientation(Arc e);
  /* Returns +1 if the face {Left(e)}, circulated in the direction of {e}, 
    is facing outwards (away from the origin); -1 if facing inwards;
    0 if the face is coplanar with the origin. */

bool_t SPTriang_EdgeIsOutwards(Arc e);
  /* TRUE if the faces {Left(e)} and {Left(Sym(e))} are both facing 
    outwards (away from the origin). */
    
bool_t SPTriang_EdgeIsSilhouette(Arc e);
  /* TRUE if any of the the faces {Left(e)} and {Left(Sym(e))} is
    coplanar with the origin, or one of them is facing outwards 
    and the other is facing inwads. */
    
/* GENERATING SITES */

SiteRef_vec_t SPTriang_MakeSites
  ( int nSites, 
    S2Point vertex_pos(int i)
  );
  /* Creates an array of site records with numbers {0..nSites-1}
    where site {i} has coordinates {vertex_pos(i)}. */

/* CREATING TABLES */

Triangulation *SPTriang_RegularOctahedron(int smpOrder);
  /* The regular octahedron (whose faces are the cardinal octants). */

Triangulation *SPTriang_FromTopology(Arc e, int smpOrder);
  /* Enumerates all edges, vertices, and faces of the polyhedron
    hanging from {e}, assigns consecutive numbers to them, creates the
    {Face} records, computes the geometrical face data
    ({b2c,c2b,sp,wp}), then packs the tables of those elements into a
    {Triangulation} record. */

void SPTriang_CheckTopology(Triangulation *tri, bool_t verbose);
  /* Performs various consistency checks on the triangulation,
    including the vertex/edge/face tables.  If {verbose} is true,
    prints the vertex degrees. */

/* LOW-LEVEL TABLE CREATION */

/* These tools are meant to be used by other triangulation-building
  procedures (e.g. {SPTriangExtra.h}). */

Arc_vec_t SPTriang_CreateFaceRecords(Arc_vec_t arc);
  /* Creates one Face record for each face of the
    triangulation reachable from {arc[...]}, and sets the 
    Ldata fields of all edges.  The triangles
    will be numbered consecutively starting from 0.
    Returns a vector with one Arc out of each triangle.  */

void SPTriang_ComputeFaceMatrices(Arc e);
  /* (Re)computes the barycentric-to-Cartesian matrix {f->b2c} and its
    inverse {f->c2b} for the face {f = Left(e)}, assuming {e} is its
    reference side. Assumes that all {Org} fields are set. */
    
void SPTriang_ComputeFaceSamples(Arc e, int smpOrder);
  /* (Re)computes the integration sampling points {f->sp} and their
    corresponding weights {f->wp} for the face {f = Left(e)}. The {smpOrder}
    parameter is given to {SPIntegral_SampleTriangle}. Assumes that
    all {Org} fields are set. */

void SPTriang_ComputeGeometryDataOfFaces(Triangulation *tri, int smpOrder);
  /* (Re)computes all derived geometry information, from the 
    site coordinates. 
    
    Specifically, calls {SPTriang_ComputeFaceMatrices(e)} and
    {SPTriang_ComputeFaceSamples(e, smpOrder)} for the reference edges
    {e = tri->side[i]} of all faces in {tri}. Assumes that the {Org}
    and {Left} fields are properly set, and the table {tri->side} is
    up-to-date.
    
    This procedure is meant to be used after changing the site positions
    in a triangulation, provided that the topology hasn't changed. */

/* TRIANGULATION I/O */

void SPTriang_Write(FILE *wr, Triangulation *t);
  /* Writes a description of {t} into {wr}, in a format that 
    can be read back with {Read} below. Saves the site
    coordinates {pos}, but omits the face 
    matrices {c2b,b2c} and the integration sampling data
    {sp,wp}. */

Triangulation *SPTriang_Read(FILE *rd, int smpOrder);
  /* Reads a triangulation {t} from the reader {rd}. Assumes the
    format used by {Write} above. The face matrices {c2b,b2c} and the
    integration sampling data {sp,wp} are recomputed from the saved
    site corrdinates and the given {smpOrder}. */

Triangulation *SPTriang_ReadCached(char *fileName, int smpOrder);
  /* Reads a triangulation {t} from the file {fileName}, with caching. 
    It is equivalent to {SPTriang_Read(open_read(fileName), smpOrder)},
    except that the {fileName} and the resulting triangulation are
    cached internally; so that subsequent calls to same {fileName}
    will just return the cached result. */

void SPTriang_Print(Triangulation *t, char *name);
  /* Writes a file {name}".txt" with a readable description of the
    triangulation. */
  
#endif
