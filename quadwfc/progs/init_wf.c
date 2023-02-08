/* See {init_wf.h}. */
/* Last edited on 2009-02-10 08:59:48 by stolfi */

#include <init_wf.h>

#include <basic.h>
#include <wavefront.h>
#include <triangulate.h>

wavefront_t wavefront_from_samples(sref_vec_t *s);
  /* Creates a {wavefront_t} record from the given sample points. The
    topology is that of the Delaunay triangulation of the
    stereographic projection of the points from the north pole of the
    unit sphere onto the plane {Z=0}. Assumes that the points are on
    the plane {Z = 0} or on the lower part of the unit sphere.
    The mesh has no face records. */

segment_t *make_init_segment(r3_t *pos, r3_t *dir, unsigned int num);
  /* Creates a sample node for the initial hemispherical mesh. The
     sample's current position will be {pos}, the velocity will be {dir},
     and the  previous position will be {pos - dir}. The {num} is saved in
     the sample's {num} field. The ray signature {sgn} is set NULL, and
     the current medium is set to zero. */

void project_sites(int ns, sref_t st[]);
  /* Applies a stereographic projection from the point {(0,0,1)} to
     all sites {st[0..ns-1]}, assumed to be on the lower half of the
     unit sphere. Only the X and Y coordinates are changed; the Z
     coordinate is left intact. */

void unproject_sites(int ns, sref_t st[]);
  /* Applies the inverse stereographic projection from the point
     {(0,0,1)} to all sites {st[0..ns-1]}. For each site {(X,Y,Z)},
     the procedure applies the inverse projection to the point
     {(X,Y,0)}, assuming that {Z} is already the coordinate of the
     result (on the sphere). Only the X and Y coordinates are changed;
     the Z coordinate is left intact. */

sref_vec_t make_unit_hemisphere_samples(double tol);
  /* Creates a dense set of samples on the lower half of the
    unit sphere. All samples have unit-length radial velocities
    and {(0,0,0)} as their previous positions. */

sref_vec_t make_unit_disk_samples(int nr, int nt);
  /* Creates a set of samples inside the unit circle of the plane
    {Z=0}. The samples are the vertices of a regular polar grid with
    {nr} steps along the radius and {nt} steps around. The samples
    have unit-length velocities, directed downwards. */

sref_vec_t make_unit_circle_samples(int n);
  /* Creates a set of {n} equally spaced samples around the nit circle
    uof the plane {Z=0}. All samples have unit-length velocities,
    udirected downwards. */

void scale_and_translate_segment(segment_t *s, r3_t *orig, double vel, double rad);
  /* Magnifies all the mesh samples in {s} by the scale factor {rad} and 
    translates them by the vector {orig}. The sample velocity vectors
    are scaled by {vel}. */

wavefront_t init_wf 
  ( r3_t *orig, double rad, double vel, double tol,
    int imd, char *sgn
  )
  {
    /* Create a lower hemisphere of radius 1, velocity 1, centered at (0,0,0): */
    sref_vec_t st = make_unit_hemisphere_samples(tol/rad);
    // sref_vec_t st = make_unit_disk_samples(5,17);
    // sref_vec_t st = make_unit_circle_samples(23);
    
    wavefront_t wf = wavefront_from_samples(&st);
    /* Translate and scale site positions and velocities: */
    int i;
    for (i = 0; i < wf.ns; i++)
      { segment_t *s = st.e[i];
        /* Translate and scale site positions and velocities: */
        scale_and_translate_segment(s, orig, vel, rad);
        /* Set initial signature and layer index: */
        s->curr.sgn = sgn;
        s->curr.imd = imd;
      }
    return wf;
  }

wavefront_t wavefront_from_samples(sref_vec_t *st)
  {
    project_sites(st->ne, st->e);
    qarc_t a = triangulate(st->ne, st->e);
    unproject_sites(st->ne, st->e);
    return (wavefront_t){ a, st->ne, *st };
  }
  
sref_vec_t make_unit_hemisphere_samples(double tol)
  {
    sref_vec_t st = sref_vec_new(1);  /* Pointers to all nodes. */
    int ns = 0;
    fprintf(stderr, "creating hemispherical samples (tol = %.3f)...\n", tol);

    /*  First sample (South pole): */
    { r3_t dir = (r3_t){{ 0.0, 0.0, -1.0 }};
      sref_vec_expand(&st, ns); 
      st.e[ns] = make_init_segment(&dir, &dir, ns);
      ns++;
    }

    /* Number of laTitude steps: */
    int nlat = (int)((M_PI/2)/tol + 0.5); 

    int ilat;
    for (ilat = 1; ilat <= nlat; ilat++)
      { /* Latitude in radians, from almost {-M_PI/2} to 0: */
        double lat = (ilat - nlat) * ((M_PI/2)/nlat);
        /* Number of loNgitude steps: */
        int nlon = (int)((2*M_PI*fabs(sin(lat)))/tol + 0.5); 
        fprintf(stderr, "  ilat = %2d  lat = %6.3f  nlon = %3d\n", ilat, lat, nlon);
        int ilon;
        for (ilon = 1; ilon <= nlon; ilon++)
          { /* Longitude in radians, from 0 to {2*M_PI}: */
            /* Add a little twist to avoid collinear points. */
            double lon = (ilon + 0.01*ilat)*((2*M_PI)/nlon); 
            /* Compute the unit direction at {-lat,lon}: */
            r3_t dir = (r3_t){{ cos(lon)*cos(lat), sin(lon)*cos(lat), sin(lat) }};
            /* Store sample in set {st}: */
            sref_vec_expand(&st, ns); 
            st.e[ns] = make_init_segment(&dir, &dir, ns);
            ns++;
          }
      }
    sref_vec_trim(&st, ns);
    fprintf(stderr, "created %d samples.\n", ns);
    return st;
  }

sref_vec_t make_unit_disk_samples(int nr, int nt)
  {
    sref_vec_t st = sref_vec_new(1);  /* Pointers to all nodes. */
    int ns = 0;

    /*  First sample (origin, moving down): */
    { r3_t pos = (r3_t){{ 0.0, 0.0, 0.0 }};
      r3_t dir = (r3_t){{ 0.0, 0.0, -1.0 }};
      sref_vec_expand(&st, ns); 
      st.e[ns] = make_init_segment(&pos, &dir, ns);
      ns++;
    }

    /* Other points */
    int ir, it;
    double twist = 0.5/((double)nr);
    for (ir = 1; ir <= nr; ir++)
      { double r = (double)ir/(double)nr;
        for (it = 1; it <= nt; it++)
          { double t = (it + twist*ir) * ((2*M_PI)/nt);
            double u = 0.0; // 0.002*drandom() - 0.001;
            double v = 0.0; // 0.002*drandom() - 0.001;
            r3_t pos = (r3_t){{ u + cos(t)*r, v + sin(t)*r, 0 }};
            r3_t dir = (r3_t){{ 0, 0, -1 }};
            sref_vec_expand(&st, ns); 
            st.e[ns] = make_init_segment(&pos, &dir, ns);
            ns++;
          }
      }
    sref_vec_trim(&st, ns);
    return st;
  }

sref_vec_t make_unit_circle_samples(int n)
  {
    sref_vec_t st = sref_vec_new(1);  /* Pointers to all nodes. */
    int ns = 0;
    int it;
    for (it = 0; it < n; it++)
      { /* Point on unit circle t {Z = 0}, disturbed, moving down: */
        double t = it * (2*M_PI/n);
        double r = 1.00 + 0.01*sin(2*t); /* Distance from origin to point. */
        r3_t pos = (r3_t){{ r*cos(t), r*sin(t), 0 }};
        r3_t dir = (r3_t){{ 0, 0, -1 }};
        sref_vec_expand(&st, ns); 
        st.e[ns] = make_init_segment(&pos, &dir, ns);
        ns++;
      }
    sref_vec_trim(&st, ns);
    return st;
  }

segment_t *make_init_segment(r3_t *pos, r3_t *dir, unsigned int num)
  {
    segment_t *s = (segment_t *)notnull(malloc(sizeof(segment_t)), "no mem");

    /* Velocity vector: */
    s->curr.vel = (*dir);
    /* Current and previous positions: */
    s->curr.pos = (*pos); 
    r3_sub(pos, dir, &(s->prev.pos));
    s->curr.sgn = NULL; /* Will be set later. */
    s->curr.imd = 0;
    /* Serial number: */
    s->num = num;
    return s;
  }

void scale_and_translate_segment(segment_t *s, r3_t *orig, double vel, double rad)
  { 
    r3_scale(vel, &(s->curr.vel), &(s->curr.vel)); 
    r3_mix(1.0, orig, rad, &(s->curr.pos), &(s->curr.pos)); 
    r3_mix(1.0, orig, rad, &(s->prev.pos), &(s->prev.pos)); 
  }
  
void project_sites(int ns, sref_t st[])
  {
    int i;
    for (i = 0; i < ns; i++)
      { 
        segment_t *s = st[i];
        double scale = 1.0 - s->curr.pos.c[2];
        s->curr.pos.c[0] /= scale;
        s->curr.pos.c[1] /= scale;
      }
  }
  
void unproject_sites(int ns, sref_t st[])
  {
    int i;
    for (i = 0; i < ns; i++)
      { 
        segment_t *s = st[i];
        double scale = 1.0 - s->curr.pos.c[2];
        s->curr.pos.c[0] *= scale;
        s->curr.pos.c[1] *= scale;
      }
  }
