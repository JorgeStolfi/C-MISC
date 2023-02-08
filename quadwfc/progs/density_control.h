/* Adjusting the density of samples on the wavefront. */
/* Last edited on 2005-09-15 14:17:48 by stolfi */

#ifndef density_control_H
#define density_control_H

#include <basic.h>
#include <wavefront.h>
#include <geomodel.h>

void density_control
  ( wavefront_t *wf, 
    double tol, 
    double time, 
    geomodel_t *geo
  );
  /* Subdivides long edges and contracts short edges as needed to maintain
    the spacing between neighbors approximately equal to {tol}. 
    For each new node, the {prev} positions are interpolated from 
    those of the neighbors, and then traced forward by the specified
    {time} amount in order to get the current positions. This 
    strategy assumes that the interpolated {prev} position is in 
    the correct region of the geophysical model. */

#endif
