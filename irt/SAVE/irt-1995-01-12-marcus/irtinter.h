#ifndef IRTINTER_H
#define IRTINTER_H

#include <h3.h>
#include <r3.h>
#include <ia.h>
#include <aa.h>
#include <pcode.h>
#include <iaeval.h>
#include <flteval.h>

extern int irt_num_rays;   /* Total calls to irt_compute_intersection */
extern int irt_num_evals;  /* Total calls to irt_compute_intersection */

typedef enum {
    arith_ia = 0,       /* Standard interval arithmetic without derivatives */
    arith_ia_diff = 1,  /* Standard interval arithmetic with derivatives */
    arith_aa = 2        /* Affine arithmetic */
  } arith_t;
  
typedef struct { /* Description of the scene's geometry: */
    char *proc_name;           /* Name of function that defines the scene's shape */
    pcode_proc_t proc;         /* The pseudo-code for the characteristic function */
    arith_t arithmetic;        /* Type of arithmetic to use */

    /* Register and stack arrays for evaluating $proc$: */
    Interval *ia_regs, *ia_stack;      /* For interval arithmetic */
    IntervalDiff *id_regs, *id_stack;  /* Ditto, including derivative. */
    AAP *aa_regs, *aa_stack;           /* For affine arithmetic. */
    Float *flt_regs, *flt_stack;       /* For floating-point. */
    FloatDiff4 *nrm_regs, *nrm_stack;  /* Ditto, including gradient. */

  } shape_t;

void irt_compute_intersection (
    shape_t *sh,         /* The object's shape */
    h3_point_t *org,     /* Ray's origin */
    h3_point_t *dst,     /* Ray's destination */
    Interval *hit,       /* (Out) Parameter interval of surface hit */
    int *slo, int *shi,  /* (Out) Function's signs before and after hit */
    int print_ray        /* Debugging switch */
  ); 
  /*
    Computes the first proper intersection $hit$ of the given
    shape, which will be evaluated with $sh->arithmetic$.
    
    Assumes the ray is parametrized by a real $t$ in $[0 .. 1]$,
    by linear interpolation of the homogeneous coordinates`of $org$ 
    and $dst$.  
    
    Ignores improper intersections, that is, intersection intervals
    that appear to begin at $t=0$ or end at $t=1$.
    
    If the routine finds a proper intersection, it returns in $hit$
    the interval of $t$ values that appears to contain the first such
    thing.  
    
    The routine also returns in $slo$ and $shi$ the sign of the
    function (+1 or -1) before and after the intersection.
    
    If the routine finds no proper intersection, it makes $hit.lo >
    hit.hi$.  In that case, the variables $slo$ and $shi$ are set to
    the sign of the function on $[0 .. 1]$ which should be either +1
    or -1 (excluding improper roots).
    
    (The signs $slo$ and/or $shi$ may be set to zero, for instance when
    the whole of $[0 .. 1]$ is one big improper root. But in that case
    the routine will probably take forever and/or die before returning...)
    
    Setting $print_ray$ to 1 generates a trace of the root-finding algorithm.
    */

void irt_compute_surface_normal (
    shape_t *sh,        /* The object's shape */
    h3_point_t *hit,    /* The hit point. */
    r3_t *nrm           /* Out: The normal vector at $hit$. */
  );
  /*
    Computes the outward unit normal vector to the surface $sh->proc$
    at the point $hit$. */
  
void irt_debug_ray_graphically (
    shape_t *sh,         /* The shape function */
    h3_point_t *org,	 /* Ray's origin */                                  
    h3_point_t *dst,	 /* Ray's destination */                             
    char *scene_name,
    char *ray_name
  );
  /* 
    Generates a Postscript plot of the function's value along the ray from
    $org$ to $dst$.  The plot is written to a file called
    "<scene_name>_<ray_name>.ps" */

#endif

