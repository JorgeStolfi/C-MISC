/* See {wave_prop.h}. */
/* Last edited on 2009-02-10 10:06:14 by stolfi */

#include <wave_prop.h>

#include <basic.h>
#include <interact.h>
#include <wavefront.h>
#include <density_control.h>

/* Thickness of reflector, for perturbations: */
#define ZEPS 0.01

void wave_prop
  ( wavefront_t *wf, 
    double time, 
    geomodel_t *geo, 
    double tol
  )
  {
    sref_vec_t *st = &(wf->st);

    /* Apply a ray tracing step on all nodes: */
    int i;
    for (i = 0; i < wf->ns; i++)
      { 
        st->e[i]->prev = st->e[i]->curr;
        ray_prop(&(st->e[i]->prev), time, geo, &(st->e[i]->curr));
        /* !!! Must delete dead vertices !!! */
      }

    /* Adjust density as needed: */
    density_control( wf, tol, time, geo );
  }

void ray_prop
  ( sample_t *u, 
    double time, 
    geomodel_t *geo, 
    sample_t *v
  )
  {
    /* Moving sample along ray segment: */
    sample_t w = (*u);

    while (! sample_is_dead(&w))
      {
        /* Advance {w} to the next intersection or to the end of the segment: */
        
        /* Compute final position {p} along ray from {w}, ignoring reflectors: */
        r3_t p;
        r3_mix(1.0, &(w.pos), time, &(w.vel), &(p));

        /* Compute first intersection of that ray with reflectors: */
        int irf;   /* Index of reflector that was hit. */
        r3_t q;    /* Point of intersection. */
        r3_t nq;   /* Surface normal at intersection. */
        double alpha; /* Fraction until intersection. */
        first_intersection(&(w.pos), &p, time, geo, /*Out*/ &q, &nq, &alpha, &irf);
        
        /* If there are no intersections, we are done: */
        if (alpha == INF) { w.pos = p; break; }

        /* Move to the intersection point: */
        w.pos = q; 

        /* Adjust time remaining for this step: */
        time = (1-alpha)*time;

        /* Compute new velocity after interaction: */
        interact(&(w), &nq, irf, &(geo->rf[irf]), ZEPS);
      }
      
    /* Sample points that leave the system must die: */
    if (! inside_bbox(&(w.pos), geo->bb)) { kill_sample(&w); }
    
    /* Return final sample point: */
    (*v) = w;
  }

