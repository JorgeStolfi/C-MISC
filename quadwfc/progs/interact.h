/* Interaction between rays and reflectors. */
/* Last edited on 2005-09-15 18:48:58 by stolfi */

#ifndef interact_H
#define interact_H

#include <basic.h>
#include <wavefront.h>
#include <geomodel.h>

bool_t inside_bbox(r3_t *p, interval_t bb[]);
  /* TRUE iff {p} is inside the rectangular box {bb[0] × bb[1] × bb[2]}. */

sign_t reflector_side ( r3_t *p, reflector_t *rf );
  /* Returns +1, 0 or -1 depending on whether {p} lies
    above, on, or below the reflector {rf}. */

int containing_medium ( r3_t *p, geomodel_t *geo );
  /* Returns the index in {geo.md} of the layer (medium) that contains
   the point {p}. Returns -1 if {p} is outside the bounding box of {geo}.
   
   At present, assumes that the reflectors in {geo.rf} span the whole
   model's domain in X and Y, do not intersect, and are sorted by
   *decreasing* Z-coordinate. */

void first_intersection_1 
  ( r3_t *p, 
    r3_t *c, 
    double time,
    reflector_t *rf, 
    /*Out*/ 
    r3_t *q, 
    r3_t *nq, 
    double *alpha
  );
  /* Computes the first intersection of the reflector {rf} with the
    ray segment from {p} to {c}. returns the intersection point in
    {q}, the corresponding reflector normal vector in {*nq}, and the
    fraction of the segment until {*q} (a number between 0 and 1) in
    {alpha}. If there is no intersection, sets {alpha = INF}. */

void first_intersection 
  ( r3_t *p, 
    r3_t *c, 
    double time,
    geomodel_t *geo, 
    /*Out*/ 
    r3_t *q, 
    r3_t *nq, 
    double *alpha, 
    int *irf 
  );
  
void interact
  ( sample_t *u,
    r3_t *ndir,
    int irf,
    reflector_t *rf, 
    double zeps
  );
  /* Computes the continuation of a ray that hits the reflector {rf}
    at the sample point {u}, where the reflector's normal is {ndir}.
    When the ray crosses a reflector, its position may be perturbed up
    or down by {zeps}, depending on the sense of crossing, in order to
    avoid bogus double-interactions.
    
    After the interaction, the ray may follow any of four paths. The
    choice is determined by a four-byte `event signature' {ty} which
    is taken from the first four bytes of the ray's future signature,
    {u->sgn[0..3]}. The event signature {ty} must have the form
    "{M}{K}{S}{N}" where {M} is the ray's current mode ('P' or 'S'),
    {K} is the desired reflector's index (ascii decimal, 1..NRF), {S}
    is either 'r' (reflected) or 't' (transmitted), and {N} is the
    mode of the ray to be followed after the interaction ('P' or 'S').
    
    The ray may die at the interaction, e.g. if the ray hits a
    reflector with index different from {K}, or if {S} is 't' and
    there is total reflection at the interface. In that case the
    sample is marked dead, and its position and other fields are
    undefined. The procedure is a no-op if the sample is already dead
    on input. */

int snell_law
  ( r3_t *v, 
    char imode, 
    r3_t *ndir, 
    reflector_t *rf, 
    char side, 
    char omode
  );
  /* Given the velocity vector {*v} and mode {imode} of the 
    incident ray, computes the velocity vector of the reflected ray (if
    {side=='r'}) or transmitted ray (if {side=='t'}) of mode {omode}
    after interaction with the reflector {*rf}. Updates the input
    velocity {*v}. Uses the normal direction {*ndir} to the reflector
    at the hit point.  
    
    The call result is {0} if the input and output media are the same,
    otherwise it is {+1} or {-1} depending on whether the output
    medium lies below or above the input medium, respectively. */

#endif
