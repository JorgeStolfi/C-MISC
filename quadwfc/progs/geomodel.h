/* Geophysical model (rock layers and interfaces). */
/* Last edited on 2005-09-15 17:38:08 by stolfi */

#ifndef geomodel_H
#define geomodel_H

#include <basic.h>

typedef struct medium_t
  { int num;      /* Identifying number. */
    double v_P;   /* Primae (pressure) wave velocity. */
    double v_S;   /* Secundae (shear) wave velocity. */
  } medium_t; 

typedef struct reflector_t /* Properties of a reflecting/refracting interface */
  { medium_t *m_above;
    medium_t *m_below;
    interval_t bb[3]; /* Bounding box of all vertices. */
    int nv[2];        /* Number of vertices along each horizontal axis. */
    double *z;        /* Z values at grid points. */
  } reflector_t;
  /* The XY projection of the mesh is a rectangular grid,
    axis-aligned, with {nx=nv[0]} by {ny=nv[1]} vertices. The Z
    coordinate of the mesh vertex with indices {ix,iy} is stored in
    {p[ix + nx*iy]}, for {ix} in {0..nx-1} and {iy} in {0..ny-1}. */

typedef struct geomodel_t
  { interval_t bb[3]; /* Bounding box of domain. */
    int nmd;          /* Number of reflectors in model. */
    medium_t *md;     /* {md[0..nmd-1]} are the media. */
    int nrf;          /* Number of reflectors in model. */
    reflector_t *rf;  /* {rf[0..nrf-1]} are the reflectors. */
  } geomodel_t; 
    
void write_geomodel(FILE *wr, geomodel_t *geo);
  /* Writes a description of {geo} into {wr}, in a format that can be
    read back with {read_geomodel} below. */

geomodel_t read_geomodel(FILE *rd);
  /* Reads a geophysical model {t} from file {rd}. Assumes the
    format used by {Write} above. */

medium_t make_medium ( int num, double v_P, double v_S );
  /* Build a medium record with serial number {imd} and
    velocities {v_P,v_S}. */

reflector_t make_sinusoidal_reflector
  (
    medium_t *m_above,
    medium_t *m_below,
    int nx, 
    int ny, 
    double xinf,
    double xsup,
    double yinf,
    double ysup,
    double zinf,
    double zsup,
    double xfreq,
    double yfreq
  );
  /* Creates a reflector between media {m_above} and {m_below}.
    The reflector will span the box
    
      {[xinf _ xsup] × [yinf _ ysup] × [zinf _ zsup]}.
    
    The rectangular mesh will have {nx} by {ny} nodes, and will be the
    graph of a sinusoidal wave with frequency vector
    {(xfreq,yfreq)}. */

#endif
