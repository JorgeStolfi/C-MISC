/* {stmaps.h} -- test maps for {dynstolfi.c} */
/* Luiz Henrique de Figueiredo & Jorge Stolfi -- 24 Jan 2008 15:00:00 */
/* Last edited on 2008-04-24 15:31:25 by stolfi */

#ifndef stmaps_H
#define stmaps_H

#include <gp.h>

/* The maps: */
typedef enum 
{ LOCORNER, /* Set (x,y) to {xmin,ymin} of the domain. */
  HICORNER, /* Set (x,y) to {xmax,ymax} of the domain. */
  DIRMAP,   /* Set (x,y) to the image of (x,y). */
  INVMAP    /* Set (x,y) to the inverse image of (x,y). */
} Opcode; 

typedef void Mapfunc(real* x, real* y, Opcode op);
/* Type of a procedure that defines a map
  of R^2 to R^2 (op=DIRMAP), its inverse (op=INVMAP),
  and a good domain for plotting (op=LOCORNER,op=HICORNER). */

Mapfunc *func_from_name(char *name);
/* Returns a map given its name ("henon", "sqturn" , etc.). */

void f_henon(real* x, real* y, Opcode op);
/* Hénon's map. */

void f_haoui(real* x, real* y, Opcode op);
/* A quadratic map similar to Hénon's. */

void f_sqout(real* x, real* y, Opcode op);
/* A smooth map that pushes up and rotates by 90 degrees.
  The square S = [-1_+1]×[-1_+1] is mapped to itself.
  The boundary of the square is an attractor for that domain.
  Points outside {S} are mapped to {S} by mirroring. */

void f_cirout(real* x, real* y, Opcode op);
/* A smooth map that pushes up and rotates by {1/phi} of a turn.
  The unit disk {C} is mapped to itself.
  The unit circle is an attractor for that domain.
  Points outside {C} are mapped into {C} by inversion. */

void f_cirin(real* x, real* y, Opcode op);
/* A smooth map that contracts towards the X axis and rotates {1/phi} of a turn.
  The unit disk {C} is mapped to itself.
  The origin is an attractor for that domain.
  Points outside {C} are mapped into {C} by inversion. */

void f_contract(real* x, real* y, Opcode op);
/* A trivial map that contracts {R^2} towards the origin.
  The origin is the only attractor. */

void f_sqwhorl(real* x, real* y, Opcode op);
/* A piecewise linear area-preserving `square rocambole' map.
  The square U = [0_1]×[0_1] is mapped to itself.
  The corner (1,0) is a fixed point.
  Points outside {U} are mapped to {U} by mirroring. */

void f_triwhorl(real* x, real* y, Opcode op);
/* A piecewise linear area-preserving `triangular rocambole' map.
  The upright equilateral triangle {T} with radius 1 is mapped to itself.
  The lower right corner of {T} is a fixed point.
  Points outside {T} are mapped to {T} by mirroring. */

void f_sqturn(real* x, real* y, Opcode op);
/* A simple map that rotates the plane by 90 degrees around (0.5,0.5).
  The any 4-symmetric region centered at (0.5,0.5) is mapped to itself.
  The point (0.5,0.5) is a fixed point. */

void f_rot144(real* x, real* y, Opcode op);
/* A simple map that rotates the plane by 144 degrees (2/5 of a turn).
  The origin is a fixed point. */

void f_rotphi(real* x, real* y, Opcode op);
/* A simple map that rotates the plane by {1/phi} of a turn.
  The origin is a fixed point. */

void f_trinop(real* x, real* y, Opcode op);
/* The identity map in the unit triangle.
  The upright equilateral triangle {T} with radius 1 is mapped to itself.
  The lower right corner of {T} is a fixed point.
  Points outside {T} are mapped to {T} by mirroring. */

void f_sqnop(real* x, real* y, Opcode op);
/* The identity map in the unit square.
  The square U = [0_1]×[0_1] is mapped to itself.
  Points outside {U} are mapped to {U} by mirroring. */

/* UTILITIES */

void solve_quadratic(double A, double B, double C, real *z1p, real *z2p);
/* Solves {A*z^2 + B*z + C = 0} for {z}, returns roots in {*z1p}, {*z2p},
  with {*z1p <= *z2p}. Fails if the roots are not real. */

void cart_to_bary(real x, real y, real *a, real *b, real *c);
 /* Maps Cartesian coordinates (x,y) to barycentric coordinates 
   {a,b,c} relative to the upright equilateral triangle {T}
   with height 1.5 and barycenter at the origin. */
   
void bary_to_cart(real a, real b, real c, real *x, real *y);
 /* The inverse of {cart_to_bary}, relative to the same triangle {T}. */

void reduce_to_unit_disk(real *x, real *y);
 /* Maps any point {(*x,*y)} to a point inside the unit disk
   C = {(x,y) : x^2+y^2 <= 1}, by inversion with respect to
   the boundary. */

void reduce_to_unit_square(real *x, real *y);
 /* Maps any point {(*x,*y)} to a point inside the unsigned unit square
   U = [0_1]×[0_1], by mirroring across the sides of that square. */

void reduce_to_signed_unit_square(real *x, real *y);
 /* Maps any point {(*x,*y)} to a point inside the signed unit square
   S = [-1_+1]×[-1_+1], by mirroring across the sides of that square. */

void reduce_to_unit_triangle(real *a, real *b, real *c);
 /* Maps any point with barycentric coords {(*a,*b,*c)} to a point
   inside the canonical triangle with corners (1,0,0), (0,1,0),
   (0,0,1), by mirroring across the sides of that triangle. */

#endif
