/* Creating the initial wavefront. */
/* Last edited on 2005-09-15 09:23:10 by stolfi */

#ifndef init_wf_H
#define init_wf_H

#include <basic.h>
#include <wavefront.h>

wavefront_t init_wf 
  ( r3_t *orig, double rad, double vel, double tol,
    int imd, char *sgn
  );
  /* Creates an initial wavefront. The samples are distributed on the
    lower half of the sphere with radius {rad} centered at {orig}.
    The samples are spaced at most {tol} apart.
    
    Each sample has current velocity vector {vel} directed radially outwards.
    The current medium index is set to {imd}, and the future 
    ray's signature is set to {sgn}.

    Also creates face records. The top face (that stands for the
    missing upper hemisphere) has {omit=TRUE}, the others have
    {omit=FALSE}. */

#endif

