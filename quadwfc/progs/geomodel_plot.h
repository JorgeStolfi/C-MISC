/* Plotting a geophysical model. */
/* Last edited on 2023-02-12 11:45:42 by stolfi */

#ifndef geomodel_plot_H
#define geomodel_plot_H

#include <basic.h>
#include <geomodel.h>

#include <bool.h>
#include <sign.h>
#include <epswr.h>
#include <frgb.h>
#include <r3.h>
#include <r4.h>
#include <r4x4.h>
#include <hr3.h>
#include <argparser.h>

void plot_geomodel(
    epswr_figure_t *eps,
    hr3_pmap_t *map,         /* Perspective projection matrix. */
    geomodel_t *geo,         /* Mesh to plot. */
    frgb_t *color,           /* Base color of mesh. */
    r3_t *dLight,            /* Direction towards main light source. */
    double shadow            /* Amount of darkening by shadow. */
  );  
  /* Paints the geophysical model {geo}, with reflectors of color {color}.
    
    Parts of the reflectors facing away from the direction {dLight}
    are darkened slightly, and parts facing towards it are lightened.
    The {shadow} parameter specifies the maximum relative change in
    the color's intensity, in either sense. */

#endif
