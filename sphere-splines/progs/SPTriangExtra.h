/* SPTriangExtra.h -- creates a triangulation where vertices have degree 5 or 7. */
/* Last edited on 2003-01-17 02:01:38 by anamaria */

#ifndef SPTriangExtra_H
#define SPTriangExtra_H

#include <r3.h>
#include <SPTriang.h>
#include <SPQuad.h>
#include <SPBasic.h>
#include <vec.h>
#include <js.h>

Triangulation *SPTriang_Refine57Triang
  ( Triangulation *old, 
    r3_t *b5, 
    r3_t *b7
  );
  /* Applies the `5-7' refinement schema to the triangulation {old}.
    The new triangulation will have {NV+6*NF} vertices, {13*NE} edges,
    and {13*NF} faces, where {NV,NE,NF} are the counts for {old}. The
    old vertices are preserved (shared) with their original positions
    and degrees; all new vertices have degree 5 or 7. 
    All other elements are newly allocated.
    
    The `shape' parameters {b5} and {b7} are the barycentric
    coordinates of the new vertices relative to the corners of the
    containing face. Good (but surely not optimal) values are 
    {b5 = (0.366, 0.462, 0.173)} and {b7 = (0.664, 0.274, 0.062)}. 
    Major changes from these values may result in flipped-over triangles. */

void SPTriang_Recompute57Triang
  ( Triangulation *tri, 
    Triangulation *old, 
    r3_t *b5,
    r3_t *b7
  );
  /* Assumes that {tri} is the result of applying {Refine57Triang} to
    {old}.  Recomputes all vertex coordinates of {tri} to accunt for 
    different shape parameters {b5,b7} and different coordinates
    in {old}.  Assumes that the topology of {old} has not changed,
    and that the vertex
    numbers in {tri} are those produced by {SPTriang_Refine57Triang}. 
    Uses {tri->out} but otherwise ignores the topology of {tri}.
    Recomputes the face matrices. */

Triangulation *SPTriang_Refine3Triang(Triangulation *old);
  /* Applies the `3-fold' refinement schema to the triangulation
    {old}. The new triangulation will have {3*NV-4} vertices, {3*NE}
    edges, and {3*NF} faces, where {NV,NE,NF} are the counts for
    {old}. The old vertices are preserved (shared) with their original
    positions and degrees; all new vertices have degree 6. All other
    elements are newly allocated. */

double_vec_t SPTriang_FaceAreas(Triangulation *tri);
  /* Areas of all faces, as flat triangles. */

double_vec_t SPTriang_EdgeLengths(Triangulation *tri);
  /* Lengths of all edges, as straight segments. */
  
double_vec_t SPTriang_RelAngles(Triangulation *tri);
  /* For each directed edge {e = tri->arc[i]}, returns the spherical 
    angle {x[i]} between {e} and the next edge with same origin, namely
    between the tangent directions {dir(e)} and {dir(Onext(e))};
    divided by the `ideal' angle {2*Pi/Degree(e)}.

    This measure has the ideal value (1.0) if the directions of
    the edges out of {Org(e)} are uniformly spaced. It tends to the 
    worst value (0.0) when {Onext(e)} tends to be collinear
    with {e}. */
  
double_vec_t SPTriang_CrossRatios(Triangulation *tri);
  /* For each arc {e = tri->arc[i]}, let {p = Org(e)},
    {q = Dest(e)}, {u = Dest(Onext(e))}, {v = Dest(Oprev(e))}.
    Let {H} be the plane defined by the vectors {u,v}. 
    The procedure returns the ratio {2*dist(p,H)/(dist(p,H) + dist(q,H))}.

    This measure has the ideal value (1.0), for example,
    when the triangles adjacent to {e} are isosceles. 
    Value 0 (when {Onext(e)} and {Oprev(e)} are collinear)
    is terrible for spline basis; negative values are bad too. */

double_vec_t SPTriang_RelOppAngles(Triangulation *tri);
  /* For each arc {e = tri->arc[i]}, let {a} be the arc out of {Org(e)}
    whose starting direction {dir(a)} makes maximum angle with {dir(e)}
    but less than or equal to {Pi}. The value {x[i]} returned is the the
    angle between {-dir(e)} and {dir(a)}, in the range {[0_Pi]}, divided
    by the `ideal' angle {Pi/Degree(e)}.

    This measure has the ideal value (1.0), for example, if {Org(e)} is
    an odd vertex with uniform angles between its arcs. It has the worst
    value (0.0) if there is an arc {a} that leaves {Org(e)} in the
    direction exactly opposite to {e}. */

void SPTriang_CookValues(double_vec_t x, double tiny, double huge);
  /* Maps any values that are less than {tiny} to
    {tiny}, and any values that are larger than {huge} to 
    {huge}, in a smooth fashion. */

double SPTriang_GeomAvg(double_vec_t x);
  /* Geometric mean of the {x[i]}.  Better none of them
    be zero or negative. */

double SPTriang_Max(double_vec_t x);
double SPTriang_Min(double_vec_t x);
  /* Maximum and minimum of the {x[i]}. */

double SPTriang_RatioVar(double_vec_t x);
  /* Mean value of {(x[i]/xavg + xavg/x[i])/2 - 1}, where {xavg} is
    the geometric mean of the {x[i]}. Better none of the {x[i]} be
    zero or negative. For small deviations, it is close to the
    variance of {x[i]/xavg}. */

double SPTriang_RatioSpr(double_vec_t x);
  /* Ratio {(minx/maxx + maxx/minx)/2 - 1}, where {minx} and {maxx}
    are the extremal values of the {x[i]}. Better none of the {x[i]}
    be zero or negative. For small deviations, it is close to {((maxx
    - minx)/maxx)^2/2}. */

#endif
