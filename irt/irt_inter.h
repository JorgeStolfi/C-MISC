#ifndef irt_inter_H
#define irt_inter_H

/* Last edited on 2023-02-22 12:17:12 by stolfi */
/* Tools for ray-scene intersection using IA, AA, etc.. */

#define _GNU_SOURCE
#include <stdint.h>

#include <hr3.h>
#include <r3.h>
#include <ia.h>
#include <aa.h>
#include <pcode.h>
#include <epswr.h>
#include <iaeval.h>
#include <flteval.h>

#include <irt_evray.h>

extern int32_t irt_num_rays;   /* Total calls to {irt_compute_intersection} */

void irt_compute_intersection
  ( shape_t *sh,         /* The object's shape */
    arith_t arith,       /* Type of arithmetic to use */
    hr3_point_t *org,    /* Ray's origin */
    hr3_point_t *dst,    /* Ray's destination */
    Interval *hit,       /* (Out) Parameter interval of surface hit */
    int32_t *slo, int32_t *shi,  /* (Out) Function's signs before and after hit */
    bool_t print_ray,    /* TRUE to print diagnostics. */
    epswr_figure_t *epsf         /* Postscript stream for debugging plot, or NULL. */
  ); 
  /* Computes the first proper intersection {hit} of the given
    shape, which will be evaluated with {arith}.
    
    Assumes the ray is parametrized by a real {t} in {[0 .. 1]},
    by linear interpolation of the homogeneous coordinates`of {org} 
    and {dst}.  
    
    Ignores improper intersections, that is, intersection intervals
    that appear to begin at {t=0} or end at {t=1}.
    
    If the routine finds a proper intersection, it returns in {hit}
    the interval of {t} values that appears to contain the first such
    thing.  
    
    The routine also returns in {slo} and {shi} the sign of the
    function (+1 or -1) before and after the intersection.
    
    If the routine finds no proper intersection, it makes {hit.lo >
    hit.hi}.  In that case, the variables {slo} and {shi} are set to
    the sign of the function on {[0 .. 1]} which should be either +1
    or -1 (excluding improper roots).
    
    (The signs {slo} and/or {shi} may be set to zero, for instance when
    the whole of {[0 .. 1]} is one big improper root. But in that case
    the routine will probably take forever and/or die before returning...)
    
    If {print_ray} is TRUE, generates some diagnostic
    printout for this ray. If {epsf} is not NULL, writes
    to it a plot of the function along the ray.  */

#endif

