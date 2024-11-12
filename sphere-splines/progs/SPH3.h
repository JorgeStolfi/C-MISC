/* SPH3.h -- Oriented projective geometry in three dimensions. */
/* Last edited on 2003-01-17 00:57:51 by anamaria */

#ifndef SPH3_H
#define SPH3_H

/* This module defines some geometric types and operations for the 
  oriented projective three-space {H^3}, using homogeneous coordinates. 
  Created 93-04-18 by Marcos C. Carrard. Based on H3.pas by J. Stolfi. */
   
#include <SPBasic.h>
#include <r3.h>
#include <r4.h>
#include <r4x4.h>
#include <r6.h>

typedef struct SPH3_Point
  { r4_t h;  /* {c[0..3]} are the homogeneous coords {[w,x,y,z]}. */
  } SPH3_Point;

#define Origin  (SPH3_Point){{{1.0, 0.0, 0.0, 0.0}}}
#define NoPoint (SPH3_Point){{{0.0, 0.0, 0.0, 0.0}}}

typedef struct SPH3_Line
  { r6_t k;  /* {k[0..5]} are the Plücker coords {[wx,wy,xy,wz,xz,yz]}. */
  } SPH3_Line;

typedef struct SPH3_Plane
  { r4_t f;  /* {f[0..3]} are the homogeneous coeffs {<W,X,Y,Z>}. */
  } SPH3_Plane;

#define Omega   (SPH3_Plane){{{1.0, 0.0, 0.0, 0.0}}}
#define NoPlane (SPH3_Plane){{{0.0, 0.0, 0.0, 0.0}}}
  
SPH3_Point SPH3_FromCartesian(r3_t *c);
  /* Given the Cartesian coordinates {c} of a point of {R^3},
    returns the corresponding point on the `hither' range 
    of {H^3} (the half-space of points with positive weight {w}). */
    
r3_t SPH3_ToCartesian(SPH3_Point *p);
  /* Cartesian coordinates of point {p} (which must be finite). */

Sign SPH3_Side(SPH3_Point *p, SPH3_Plane *Q); 
  /* Returns the position of point {p} relative to plane {Q}: 0 on the
    plane, +1 in the positive halfspace, -1 in the negative halfspace.
    May give inconsistent results for points very close to the plane. */

Sign SPH3_Orient(SPH3_Point *p, SPH3_Point* q, SPH3_Point *r, SPH3_Point *s);
  /* Returns the orientation (handedness) of the tetrahedron {p q r s}.
    Specifically, the tetrahedron {[1 0 0 0] [0 1 0 0] [0 0 1 0] [0 0 0 1]}
    has positive orientation. Returns 0 iff the points are coplanar.
    May give inconsistent results for almost-coplanar points. */
    
SPH3_Line SPH3_LineFromTwoPoints(SPH3_Point *p, SPH3_Point *q);
  /* The line from {p} to {q}. */

SPH3_Plane SPH3_PlaneFromThreePoints(SPH3_Point *p, SPH3_Point *q, SPH3_Point *r);
  /* The plane through {p}, {q}, and {r}. */
    
SPH3_Point SPH3_PointFromThreePlanes(SPH3_Plane *P, SPH3_Plane *Q, SPH3_Plane *R);
  /* The point where {P}, {Q}, and {R} meet. */
    
SPH3_Plane SPH3_PlaneFromLineAndPoint(SPH3_Line *n, SPH3_Point *r);
  /* The plane through {n} and {r}. */
    
SPH3_Line SPH3_LineFromTwoPlanes(SPH3_Plane *P, SPH3_Plane *Q);
  /* The line where {P} meets {Q}. */

SPH3_Point SPH3_PointFromLineAndPlane(SPH3_Line *n, SPH3_Plane *R);
  /* The point where {n} meets {R}. */

r3_t SPH3_Dir(SPH3_Point *frm, SPH3_Point *tto);
  /* Direction (a unit-length vector) of point {tto} seen from point {frm}.
    Works even if one of them is at infinity.  Does not work if both
    are at infinity, or coincident, or antipodal. */
    
double SPH3_Dist(SPH3_Point *a, SPH3_Point *b);
  /* Distance between {a} and {b}, which must lie in the front range. */
    
double SPH3_DistSqr(SPH3_Point *a, SPH3_Point *b);
  /* Distance squared between {a} and {b}, which must lie in the front range. */
    
r3_t SPH3_Normal(SPH3_Plane *P);
  /* The normal direction of plane {P}, which, on the hither range, points into 
    {P}'s positive halfspace.  Assumes {P} is not at infinity. */
    
/* PROJECTIVE MAPS */

typedef struct SPH3_PMap
  { r4x4_t dir; 
    r4x4_t inv;
  } SPH3_PMap;  
  /* {dir} is the map's matrix, {inv} is its inverse. */

SPH3_Point SPH3_MapPoint(SPH3_Point *p, SPH3_PMap *m);
  /* Applies projective map {m} to point {p}. */

SPH3_Point SPH3_InvMapPoint(SPH3_Point *p, SPH3_PMap *m);
  /* Applies the inverse of projective map {m} to point {p}. */

SPH3_Plane SPH3_MapPlane(SPH3_Plane *P, SPH3_PMap *m);
  /* Applies projective map {m} to plane {P} */

SPH3_Plane SPH3_InvMapPlane(SPH3_Plane *P, SPH3_PMap *m);
  /* Applies the inverse of projective map {m} to plane {P} */

SPH3_PMap SPH3_CompMap(SPH3_PMap *m, SPH3_PMap *n);
  /* Returns the composition of {m} and {n}, applied in that order */
    
SPH3_PMap SPH3_InvMap(SPH3_PMap *m);
  /* Returns the inverse of map {m}. */

SPH3_PMap SPH3_PerspMap(SPH3_Point *obs, SPH3_Point *foc, SPH3_Point *upp);
  /* Computes a perspective transformation with given viewing parameters:
      {obs} == the scenesys coords of imagesys {(0,0,d)} (the observer);
      {foc} == the scenesys coords of imagesys {(0,0,0)} (the image focus);
      {up}  == the scenesys coords of some point with 
                 imagesys {X == 0, Y > 0} (zenith reference). */

#endif
