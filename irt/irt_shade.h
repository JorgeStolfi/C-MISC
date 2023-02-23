#ifndef irt_shade_H
#define irt_shade_H

/* Last edited on 2023-02-22 19:54:38 by stolfi */
/* Color and shading computations for the Interval Ray Tracer. */

#define _GNU_SOURCE
#include <stdint.h>

#include <hr3.h>
#include <r3.h>
#include <frgb.h>
#include <epswr.h>

#include <irt_scene.h>

frgb_t irt_compute_scene_color
  ( scene_t *sc,
    hr3_point_t *eye,  
    r3_t *dir,
    int32_t max_bounces,
    arith_t arith,
    bool_t print_ray,
    epswr_figure_t *epsf
  );
  /* Computes the color of the scene when seen from point {eye}
    and looking towards the direction {dir}.  The {max_bounces} parameter
    is the maximum number of times the ray is allowed to bounce off 
    mirror surfaces. Uses the arithmetic {arith} to find the first 
    ray-scene hit.
    
    It is best if {*eye} is normalized so that the largest absolute 
    coordinate is about 1. 
    
    If {print_ray} is TRUE, generates some diagnostic
    printout for this ray. If {epsf} is not NULL, writes
    to it a plot of the function along the ray.  */

#endif

