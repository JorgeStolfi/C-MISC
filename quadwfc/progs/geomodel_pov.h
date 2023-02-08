/* Generating POV-ray description of a geophysical model. */
/* Last edited on 2005-08-26 15:49:37 by stolfi */

#ifndef geomodel_pov_H
#define geomodel_pov_H

#include <basic.h>
#include <geomodel.h>

#include <bool.h>
#include <sign.h>
#include <r3.h>

#include <stdio.h>

void pov_geomodel(FILE *fpov, geomodel_t *geo, double radius);  
  /* Writes a POV-ray description of the geophysical model {geo}. Each
    grid cell is assumed to be two flat triangles. Edges and vertices
    are drawn as cylinders and balls with the given radius. */

#endif
