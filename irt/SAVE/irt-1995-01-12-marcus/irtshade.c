/* Last edited on 2007-01-04 00:19:25 by stolfi */ 

#include "irtshade.h"
#include <rgb.h>
#include <ia.h>
#include <r3.h>
#include <h3.h>
#include <ioprotos.h>
#include <stdio.h>
#include <math.h>
#include <irtinter.h>

/*** INTERNAL PROTOTYPES ***/

rgb_t irt_compute_hit_color (
    scene_t *sc,
    r3_t *dir, 
    h3_point_t *npt, 
    h3_point_t *fpt, 
    int slo, int shi,
    int max_bounces,
    int debug_flags
  );
  /*
    Computes the apparent color of a point on the surface of the scene,
    when seen from very close, looking in the direction $dir$.
    
    The caller must give $npt$ and $fpt$, two points bracketing the 
    ray-surface intersection (with $npt$ nearest to the observer).
    It should also give the sign $slo$ of the function right in front of 
    $npt$, and the sign $shi$ right beyond $fpt$. */
  
rgb_t irt_compute_scene_transparency (
    scene_t *sc,
    h3_point_t *org,
    h3_point_t *dst,
    double dist,
    int print_ray
  );
  /*
    Computes the straight-line transparency between the two points $org$
    and $dst$ on the scene.  
    
    The parameter $dist$ should be the Euclidean distance between the
    two points (used to compute the attenuation inside transparent
    materials.)
    
    It is best if $org$ and $dst$ are normalized so that their coordinates
    have about the same magnitude. 
    
    Setting $print_ray$ generates a debugging printout of the ray. */
  
void irt_compute_hit_endpoints (
    h3_point_t *org,   /* Start of ray */
    h3_point_t *dst,   /* End of ray */
    Interval *hit,     /* Parameter interval containing hit point */
    h3_point_t *npt,   /* Out: "near" hit endpoint */
    h3_point_t *fpt    /* Out: "far" hit endpoint */
  );
  /*
    Computes the "near" and "far" points on the segment $org$--$dst$,
    corresponding to the low and high points of the parameter 
    interval $hit$.  */

void irt_compute_mean_hit_point (
    h3_point_t *npt,
    h3_point_t *fpt,
    h3_point_t *mpt
  );
  /*
    Computes the "mean" (nominal) hit point $mpt$, somewhere
    between $npt$ and $fpt$. */

/*** IMPLEMENTATIONS ***/

#define FUDGE 1.0e-6

rgb_t irt_compute_scene_color(
    scene_t *sc,
    h3_point_t *eye, 
    r3_t *dir,
    int max_bounces,
    int debug_flags
  )
  {
    h3_point_t inf;  /* Ray's endpoint (at infinity) */

    h3_point_t npt;  /* "Near" (towards observer) hit point */
    h3_point_t fpt;  /* "Mean" hit point */
    
    Interval hit;    /* Parameter interval containing hit point */
    int slo, shi;    /* Sign of function before and after $hit$ (see zerofind.h) */
    rgb_t color;
    
    int print_ray = ((debug_flags & 1) != 0);
    int plot_ray = ((debug_flags & 2) != 0);
    
    h3_infty(dir, &inf);

    /* Trace the ray: */
    irt_compute_intersection (&(sc->shape), eye, &inf, &hit, &slo, &shi, print_ray);

    if (plot_ray)
      { irt_debug_ray_graphically(&(sc->shape), eye, &inf, sc->name, "sample-ray"); }

    if (hit.lo <= hit.hi)  /* ray intersected object */
      { 
        irt_compute_hit_endpoints (eye, &inf, &hit, &npt, &fpt);
        if (slo == -1)
          { if (sc->looks.has_transp)
              { error("irt_compute_scene_color: transparent solids not implemented yet"); }
            else
              { fprintf(stderr, "  *** irt_compute_scene_color: observer inside opaque object\n");
                fprintf(stderr, "      eye = "); h3_print_point(stderr,  eye); fprintf(stderr, "\n");
                fprintf(stderr, "      inf = "); h3_print_point(stderr, &inf); fprintf(stderr, "\n");
                fprintf(stderr, "      hit = "); ia_print(stderr, hit); fprintf(stderr, "\n");
                fprintf(stderr, "      npt = "); h3_print_point(stderr, &npt); fprintf(stderr, "\n");
                fprintf(stderr, "      fpt = "); h3_print_point(stderr, &fpt); fprintf(stderr, "\n");
                irt_debug_ray_graphically(&sc->shape, eye, &inf, sc->name, "bogus-ray");
                error("aborted");
              }
          }
        color = irt_compute_hit_color (sc, dir, &npt, &fpt, slo, shi, max_bounces, print_ray);
      }
    else
      { color = sc->looks.background_color; }
    return(color);
  }

rgb_t irt_compute_hit_color (
    scene_t *sc,
    r3_t *dir, 
    h3_point_t *npt, 
    h3_point_t *fpt, 
    int slo, int shi,
    int max_bounces,
    int print_ray
  )
  { h3_point_t mpt;  /* "Mean" (nominal) hit point */
    
    r3_t nrm;  /* Outwards-pointing surface normal at $mpt$ */ 
    r3_t rfl;  /* Reflected viewing ray's direction (for mirror and shiny looks) */

    r3_t npc;  /* Cartesian coordinates of "near" hit point, perturbed */
    r3_t fpc;  /* Cartesian coordinates of "far" hit point, perturbed */
    
    h3_point_t npp; /* Homogeneous coordinates of "near" hit point, perturbed */
    h3_point_t fpp; /* Homogeneous coordinates of "far" hit point, perturbed */

    double cos_dir_nrm; /* Cosine of ray direction and normal */
    
    rgb_t color; /* Total apparent color of surface at hit point */
    rgb_t matte_light; /* Matte component of color */
    rgb_t shiny_light; /* Shiny or mirror component of color */
    
    int lnum;

    irt_compute_mean_hit_point(npt, fpt, &mpt);
    
    irt_compute_surface_normal (&(sc->shape), &mpt, &nrm);
    cos_dir_nrm = r3_dot(dir, &nrm);
        
    /* Perturb the "near" hit point along normal, towards the observer's side: */
    h3_to_cart (npt, &npc);
    r3_mix_in (copysign(FUDGE, -cos_dir_nrm), &nrm, &npc);
    h3_from_cart (&npc, &npp);
    h3_inf_reduce(&npp);
    
    if (sc->looks.has_shiny || sc->looks.has_mirror || (slo == -1)) 
      { /* calculate reflected ray direction */
        /* Used for highlights, mirroring, or internal reflection. */
	double k = -2.0 * r3_dot(dir, &nrm);
	rfl.c[0] = k * nrm.c[0] + dir->c[0];
	rfl.c[1] = k * nrm.c[1] + dir->c[1];
	rfl.c[2] = k * nrm.c[2] + dir->c[2];
        (void) r3_normalize(&rfl);
      }
    
    if (sc->looks.has_transp)
      {
	/* Perturb the "far" hit point along normal, away from observer's: */
        /* Used for refracted rays */
	h3_to_cart (fpt, &fpc);
	r3_mix_in (copysign(FUDGE, cos_dir_nrm), &nrm, &fpc);
	h3_from_cart (&fpc, &fpp);
	h3_inf_reduce(&fpp);
      }

    /* Initialize components: */
    color = sc->looks.ambient_color;
    rgb_black(&matte_light);
    rgb_black(&shiny_light);

    for(lnum=0; lnum < sc->num_lights; lnum++)
      { light_t *lt = &(sc->light[lnum]); 
        
        r3_t ldr;             /* Direction of light source */
        double ldist;         /* Distance to light source */
        double rdist;         /* Light's reference distance */
        
        double cos_nrm_ldr;   /* Cosine of angle between normal and light direction */
        
	/* get light's direction and distance */
        r3_sub(&(lt->position), &npc, &ldr);
        ldist = r3_normalize(&ldr);
        rdist = lt->ref_distance;

	cos_nrm_ldr = r3_dot (&nrm, &ldr);
	if(cos_nrm_ldr * cos_dir_nrm < 0.0)
          { /* Object faces light; compute matte- and shiny-scattered light: */
	    
            h3_point_t lpt;       /* Position of light source */
            rgb_t shadow_factor;  /* Transparency coeffs from light to hit point */
            double inv_sq_factor; /* Inverse square law attenuation factor */
	    rgb_t incident_light; /* Color of incident light */
            
	    /* Attenuate light's color by shadow and inverse-square factors: */
            h3_from_cart(&(lt->position), &lpt); h3_inf_reduce(&lpt);
            shadow_factor = irt_compute_scene_transparency (sc, &npp, &lpt, ldist, print_ray);
	    rgb_weigh (&(lt->color), &shadow_factor, &incident_light);
            inv_sq_factor = 2.0/((ldist*ldist)/(rdist*rdist) + 1.0);
            rgb_scale (inv_sq_factor, &incident_light, &incident_light);
            
            if (! rgb_is_almost_black(&incident_light))
              { /* Accumulate matte-scattered incident light: */
		rgb_mix_in (fabs(cos_nrm_ldr), &incident_light, &matte_light);

		if (sc->looks.has_shiny)
		  { double cos_rfl_ldr = r3_dot (&rfl, &ldr);
		    if (cos_rfl_ldr > sc->looks.shiny_spread)
		      { /* Acumulate shiny-scattered incident light */
                        double ang_coeff; /* Attenuation for angle with reflected ray */
			double grz_coeff; /* Enhancement for grazing reflection */

			ang_coeff = (cos_rfl_ldr - sc->looks.shiny_spread)/(1.0 - sc->looks.shiny_spread);
			ang_coeff = ang_coeff*ang_coeff;

			grz_coeff = (2.0 - cos_nrm_ldr*cos_nrm_ldr)/2.0;

			rgb_mix_in(ang_coeff*grz_coeff, &incident_light, &shiny_light);
		      }
		  }
              }                    
	  }
      }


    /* Mirror-like reflection */
    if ((sc->looks.has_mirror) && (max_bounces > 0))
      { shiny_light = irt_compute_scene_color(sc, &npp, &rfl, max_bounces - 1, print_ray); }

    /* Refraction */
    if (sc->looks.has_transp)
      { error ("irt_compute_hit_color: transparent solids not implemented yet"); }

    /* Put it all together: */
    rgb_weigh_in(&matte_light, &(sc->looks.matte_color), &color);
    rgb_weigh_in(&shiny_light, &(sc->looks.shiny_color), &color);

    return(color);
  }

rgb_t irt_compute_scene_transparency (
    scene_t *sc,
    h3_point_t *org,
    h3_point_t *dst,
    double dist,
    int print_ray
  )
  { rgb_t color;
    Interval hit;
    int slo, shi;
    
    irt_compute_intersection (&(sc->shape), org, dst, &hit, &slo, &shi, print_ray);

    if (hit.lo <= hit.hi) 
      { rgb_black(&color); }
    else if (slo == +1)
      { rgb_white(&color); }
    else if (sc->looks.has_transp)
      { error("irt_compute_scene_transparency: transparent solids not implemented"); }
    else
      { rgb_black(&color); }
    return(color);
  }

void irt_compute_hit_endpoints (
    h3_point_t *org,
    h3_point_t *dst,
    Interval *hit,
    h3_point_t *npt,
    h3_point_t *fpt
  )
  { h3_mix (1.0-hit->lo, org, hit->lo, dst, npt); h3_inf_reduce(npt);
    h3_mix (1.0-hit->hi, org, hit->hi, dst, fpt); h3_inf_reduce(fpt);
  }

void irt_compute_mean_hit_point (
    h3_point_t *npt,
    h3_point_t *fpt,
    h3_point_t *mpt
  )
  { double nw = 0.5*npt->c[0];
    double fw = 0.5*fpt->c[0];
    h3_mix (fw, npt, nw, fpt, mpt);
    h3_inf_reduce(mpt);
  }
