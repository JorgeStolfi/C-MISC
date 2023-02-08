/* Propagating the wavefront. */
/* Last edited on 2005-09-15 14:00:12 by stolfi */

#ifndef wave_prop_H
#define wave_prop_H

#include <basic.h>
#include <wavefront.h>
#include <geomodel.h>

void wave_prop
  ( wavefront_t *wf, 
    double time, 
    geomodel_t *geo, 
    double tol
  );
  /* Propagates the wavefront represented by {wf} during {time} seconds */

void ray_prop
  ( sample_t *u, 
    double time, 
    geomodel_t *geo,
    sample_t *v
  );
  /* Traces a ray from sample {u} by the specified {time}, taking care
    of interactions with the given reflectors. The ray's behavior at
    each interaction is governed by the sample's future signature
    {u.sgn}. The final sample is returned in {v}.
    
    If the sample is `dead' does nothing. If the ray dies for some
    reason (e.g. when trying to follow a 't' ray but total internal
    reflection occurs), the sample becomes dead. */

#endif
