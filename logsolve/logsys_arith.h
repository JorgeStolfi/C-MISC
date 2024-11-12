#ifndef logsys_arith_H
#define logsys_arith_H

/* Last edited on 2012-12-19 23:16:43 by stolfilocal */
/* Equation systems for arithmetic. */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>

#include <bool.h>

#include <logsys.h>

logsys_t *logsys_build_mul(int nx, int ny, logsys_va_t *x[], logsys_va_t *y[], logsys_va_t *z[]);
  /* Returns a {logsys_t} for a (nx,ny,nx+ny)-bit multiplier. Returns in
    {x[0..nx-1]} and {y[0..ny-1]} the variables corresponding to the
    bits of th two multiplicands {X} and {Y}, and in
    {z[0..nx+ny-1]} the variables that correspond to the bits of the
    product {Z=X*Y}.  
    
    The arrays {x,y,z} must be allocated by the client, with the
    proper size. Bit {[0]} is the least significant in all three
    arrays. The system will contain about {6*nx*ny} internal variables
    and equations. */
     
void logsys_add_adder_slice
  ( logsys_t *S, 
    logsys_va_t *u,
    logsys_va_t *v, 
    logsys_va_t *w, 
    logsys_va_t **s, 
    logsys_va_t **c
  );
  /* Adds to {s} a set of variables and equations that implement a
    three-bit adder circuit whose inputs are {u,v,w}. The output of the
    circuit consists of two variables, a sum bit (whose handle is
    returned in {*s}) and a carry bit (whose handle is returned in {*c}.

    In general, the procedure will add five new variables and five
    equations, corresponding to the boolean formulas {s1 := (u ++ v); c1
    := (u /\ v); s := s1 ++ w; c2 := (s1 /\ w); c := c1 ++ c2}. However,
    if any of {u,v,w} is NULL, that variale is assumed to be identically
    zero, and the circuit is simplified accordingly. In particular, if
    only one argument is NULL, the procedure will add only two new
    variables and two equations. If two arguments are NULL, the
    procedure will not change {S}; the sum variable {*s} will be the
    other argument (which may be NULL), and the carry variable {*c} will
    be NULL. */

/* DEBUGGING */

void logsys_check_mul(logsys_t *S, int nx, int ny, logsys_va_t *x[], logsys_va_t *y[], logsys_va_t *z[]);
  /* Checks whether the {logsys_t} {S} indeed behaves as a 
    (nx,ny,nx+ny)-bit multiplier. */

void logsys_check_solve_mul
  ( logsys_t *S, 
    int ni, 
    logsys_va_t **iv, 
    int nx, 
    int ny, 
    logsys_va_t *x[], 
    logsys_va_t *y[], 
    logsys_va_t *z[]
  );
  /* Checks whether the {logsys_t} {S} indeed behaves as a 
    (nx,ny,nx+ny)-bit multiplier. */

#endif
