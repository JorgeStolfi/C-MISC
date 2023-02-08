/* Writing and plotting wavefronts and models. */
/* Last edited on 2005-08-21 10:17:09 by stolfi */

#ifndef output_H
#define output_H

#include <basic.h>
#include <wavefront.h>
#include <geomodel.h>

void plot_wave(wavefront_t *wf, char *outName, int iter);
/* Plots the wavefront as an Encapsulated Postscript file,
   called "{outName}-{iter}.eps". */

void output_wave(wavefront_t *wf, char *outName, int iter);
/* Writes the geometry of the wavefront as a text file,
   called "{outName}-{iter}.tri". */

void output_model(geomodel_t *geo, char *outName);
/* Writes the geophysical model {geo} to file "{outName}.geo". */

#endif
