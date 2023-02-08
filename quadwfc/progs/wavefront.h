/* Represenatation of the wavefront. */
/* Last edited on 2013-10-02 03:08:50 by stolfilocal */

#ifndef wavefront_H
#define wavefront_H

#include <basic.h>

#include <affirm.h>
#include <bool.h>
#include <vec.h>
#include <r3.h>
#include <r4.h>
#include <interval.h>
#include <quad.h>
#include <frgb.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <values.h>

/* WAVEFRONT SAMPLE POINTS */

typedef struct sample_t /* Data for a sample point along a ray. */
  { r3_t pos;      /* Position of sample. */
    char *sgn;     /* Signature of ray remainder, starting at current mode. */
    int imd;       /* Index of current layer in geophysical model. */
    r3_t vel;      /* Wavefront velocity vector. */
  } sample_t;
  /* A {sample_t} describes a point {pos} along a traced ray,
   the signature {sgn} of the forward (future) part of the ray,
   the index {imd} of the geophysical model's layer that contains
   this point {pos}, and the local velocity vector {vel}.  
   
   The signature {sgn} is an alternating sequence of wave mode letters
   ('p' or 's'), reflector indices ('0' to '9') and interaction codes
   ('r' or 't'), starting and ending with a mode letter. A signature {sgn}
   that is NULL or points to a null character ('\000') indicates that
   the ray is `dead' (no longer being traced), e.g. because it left the 
   system or hit the wrong reflector. In that case, the position {pos}
   and other attributes  of the sample are irrelevant. */

void kill_sample(sample_t *s);
  /* Marks {s} as belonging to a `dead' ray, which is no longer being
    simulated (because e.g. it has left the system, it was totally
    absorbed, etc..) */

bool_t sample_is_dead(sample_t *s);
  /* TRUE iff {s} has been marked `dead' is zero. */

/* MESH VERTICES */

typedef struct segment_t /* Data for a mesh vertex. */
  { sample_t curr;  /* Sample of current wavefront. */
    sample_t prev;  /* Sample of previous wavefront (if any) or last interaction. */
    unsigned num;   /* Vertex number; used inside {write_mesh} procedure. */
  } segment_t;

typedef segment_t *sref_t;

/* MESH FACEDS */

typedef struct face_t 
  { bool_t omit;  /* TRUE iff the face is a hole in the mesh. */ 
    unsigned num; /* Face serial number. */
  } face_t;

typedef face_t *fref_t;

/* MESH EDGES */

typedef quad_arc_t qarc_t; /* A reference to a directed edge, primal or dual. */

#define qarc_NULL quad_arc_NULL

/* VECTORS OF MESH ELEMENTS */

typedef quad_arc_vec_t qarc_vec_t;
vec_typedef(sref_vec_t,sref_vec,sref_t);
vec_typedef(fref_vec_t,fref_vec,fref_t);
  /* An object {v} of type {XXX_vec_t} is a vector of things of type {XXX_t}. The
    number of elements is {v.ne}, the {e[0..ne-1]}, with its size.
    
    The procedure {XXX_vec_new} allocates a new {XXX_vec_t} with space
    for {ne} elements.
    
    The procedure {XXX_vec_expand} reallocates {v->e} if necessary so
    that {v->e[index]} exists. Preserves old elements and updates
    {v->ne} to the new size, strictly greater than {index}. Size is
    doubled at each reallocation to ensure total {O(N)} time for {N}
    calls.
    
    The procedure {XXX_vec_trim} reallocates {v->e} if necessary so
    that it has exactly {num} elements. Preserves the values of
    {v->e[0..ne-1]} and sets {v->ne = ne}. */

/* WAVEFRONT REPRESENTATION */

typedef struct wavefront_t /* Representation of a wavefront. */
  { qarc_t a;       /* Some edge reference of the mesh. */
    int ns;         /* Number of active samples. */
    sref_vec_t st;  /* Current nodes are {st.e[0..ns-1]}, in arbitrary order. */
  } wavefront_t;
  /* The wavefront is a 2D topological mesh, reachable from {a}. Each
    vertex of the mesh is associated to a {segment_t} record, a sample
    point on the moving wavefront. ??? ADD FACES AND HOLES. */

/* WAVEFRONT-SPECIFIC TOPOLOGICAL OPERATORS */

/* Mesh vertex pointers: */

#define quad_org(e) ((segment_t *) quad_odata(e))
#define quad_dst(e) ((segment_t *) quad_ddata(e))

#define SET_quad_org(e,v) quad_set_odata((e),(void *)(v))
#define SET_quad_dst(e,v) quad_set_ddata((e),(void *)(v))

/* Coordinates of vertex on current mesh: */

#define ORGP(e) ((quad_org(e))->curr.pos)
#define DESTP(e) ((quad_dst(e))->curr.pos)

/* Mesh face pointers: */

#define LEFT(e) ((face_t *) quad_ddata(quad_rot(e)))
#define RIGHT(e) ((face_t *) quad_odata(quad_rot(e)))

#define SET_LEFT(e,f) quad_set_ddata(quad_rot(e),(void *)(f))
#define SET_RIGHT(e,f) quad_set_odata(quad_rot(e),(void *)(f))

void create_face_records(qarc_t e);
  /* Creates one {face_t} record for each face of the triangulation
    reachable from {e} (including holes), and sets the {LEFT} fields
    of the appropriate edges. The faces will be numbered consecutively
    starting from 0. Any face records already linked by the {LEFT}
    will be ignored. */

/* TRAVERSAL */

qarc_vec_t renumber_edges(qarc_vec_t a);
  /* Enumerates undirected edges reachable from the arcs {a.e[..]}, and
    stores in their {num} fields consecutive integers.
    Returns an array {res} where {res.e[i]} is one primally reachable
    arc from the edge with number {start+i}.

    An arc {b} is /primally reachable/ from an arc {a} iff it can be obtained 
    from {a} by a finite number of {Sym} and {Onext} operations (no {Rot}s
    or {Tor}s).  An edge {e} is /reachable/ iff some of its arcs is reachable. */

qarc_vec_t arcs_from_edges(qarc_vec_t edge);
  /* Given an array {edge.e[..]} with one arc from each undirected 
    primal edge, returns an array {res} with two symmetric arcs for
    each edge.  The arcs for {edge.e[i]} are {res.e[2*i]} and 
    {res.e[2*i+1]}, in no particular order.  ??? unnecessary? */

qarc_vec_t renumber_vertices(qarc_vec_t arc);
  /* Renumbers all the origin vertices of the arcs in {arc[..]}
    consecutively from 0. Returns a vector of arcs {out} such that
    {quad_org(out[i])} is vertex number {i}, for all vertices, without
    repetitions. */

qarc_vec_t renumber_faces(qarc_vec_t arc);
  /* Renumbers all the left faces of the arcs in {arc[..]}
    consecutively from 0. Returns a vector {side} such that
    {LEFT(side[i])} is face number {i}, for all faces, without
    repetition. */

#endif
