/* SPBasic.h -- basic data types for the spherical splines project. */
/* Last edited on 2024-11-11 03:31:19 by stolfi */

#ifndef SPBasic_H
#define SPBasic_H

#include <r3.h>
#include <r6.h>
#include <js.h>
#include <nat.h>

typedef enum { negative = -1, nullary = 0, positive = 1 } Sign;

typedef char * string;

typedef r3_t R3Point;    /* Point of {R^3}, possibly not on {S^2}. */
typedef r3_t R3Gradient; /* Gradient relative to Cartesian coordinates. */
typedef r6_t R3Hessian;  /* Hessian lower triang matrix rel to Cartesian coords. */
typedef r3_t R3Vector;   /* Vector of {R^3}, possibly not tangent {S^2}. */
                         
typedef r3_t S2Point;    /* Point on the sphere {S^2}. */
typedef r3_t S2Gradient; /* Spherical gradient (tangent to {S^2}). */

typedef double ScalarField (S2Point *p);

typedef r3_t VectorField (S2Point *p);

typedef r3_t Color;

typedef signed char int8;
typedef unsigned char nat8;
    
#define PI       (3.1415926535897932384626433833)
#define TWOPI    (2*PI)
#define FOURPI   (4*PI)

#define PHI      (1.6180339887498948482)
#define SQRT3    (1.73205080756887729352)
#define SQRTHALF (0.70710678118654752440084)

/* Some useful macros: */

#define mumble(...) \
  do { if (verbose) { fprintf(stderr, __VA_ARGS__); } } while(0)

/* Some useful geometric tools: */

vec_typedef(R3Point_vec_t,R3Point_vec,R3Point);
  /* A vector of points of {R^3} ({R3Point}s). */

#define S2Point_vec_t R3Point_vec_t
  /* Arrays of points on the sphere ({S2Point}s). */
  
#define S2Point_vec_new R3Point_vec_new
#define S2Point_vec_expand R3Point_vec_expand
#define S2Point_vec_trim R3Point_vec_trim

#endif
