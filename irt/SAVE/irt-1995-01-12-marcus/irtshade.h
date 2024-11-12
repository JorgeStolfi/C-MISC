#ifndef IRTSHADE_H
#define IRTSHADE_H

#include <h3.h>
#include <r3.h>
#include <rgb.h>
#include <irtscene.h>

rgb_t irt_compute_scene_color (
    scene_t *sc,
    h3_point_t *eye, 
    r3_t *dir,
    int max_bounces,
    int debug_flags
  );
  /*
    Computes the color of the scene when seen from point $eye$
    and looking towards the direction $dir$.  The $max_bounces$ parameter
    is the maximum number of times the ray is allowed to bounce off 
    mirror surfaces.  
    
    It is best if *eye is normalized so that the largest absolute 
    coordinate is about 1. 
    
    Setting bit 0 (value 1) of $debug_flags$ generates some diagnostic
    printout for this ray. Setting bit 1 (value 2) generates a Postscript
    file with a plot of the function along the ray.  */

#endif

