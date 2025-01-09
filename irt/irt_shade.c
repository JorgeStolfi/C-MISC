/* See {irt_shade.h}. */
/* Last edited on 2024-11-15 18:30:06 by stolfi */

#define _GNU_SOURCE
#include <stdint.h>
#include <stdio.h>
#include <math.h>

#include <ia.h>
#include <r3.h>
#include <hr3.h>
#include <bool.h>
#include <affirm.h>
#include <frgb.h>
#include <epswr.h>
#include <frgb_ops.h>

#include <irt_inter.h>

#include <irt_shade.h>

/*** INTERNAL PROTOTYPES ***/

bool_t color_is_almost_black(frgb_t *r, double tol);
  /* TRUE iff all components of {*r} are less than {tol}. */

frgb_t irt_compute_hit_color 
  ( scene_t *sc,
    r3_t *dir, 
    hr3_point_t *npt, 
    hr3_point_t *fpt, 
    int32_t slo, int32_t shi,
    int32_t max_bounces,
    arith_t arith,
    bool_t print_ray,
    epswr_figure_t *epsf
  );
  /* Computes the apparent color of a point on the surface of the scene,
    when seen from very close, looking in the direction {dir}.
    
    The caller must give {npt} and {fpt}, two points bracketing the 
    ray-surface intersection (with {npt} nearest to the observer).
    It should also give the sign {slo} of the function right in front of 
    {npt}, and the sign {shi} right beyond {fpt}. */
  
frgb_t irt_compute_scene_transparency
  ( scene_t *sc,
    hr3_point_t *org,
    hr3_point_t *dst,
    double dist,
    arith_t arith,
    int32_t print_ray,
    epswr_figure_t *epsf
  );
  /* Computes the straight-line transparency between the two points {org}
    and {dst} on the scene.  
    
    The parameter {dist} should be the Euclidean distance between the
    two points (used to compute the attenuation inside transparent
    materials.)
    
    It is best if {org} and {dst} are normalized so that their coordinates
    have about the same magnitude. 
    
    Setting {print_ray} generates a debugging printout of the ray. If
    {epsf} is not NULL, writes to it a plot of the function along the
    ray. */
  
void irt_compute_hit_endpoints
  ( hr3_point_t *org,   /* Start of ray */
    hr3_point_t *dst,   /* End of ray */
    Interval *hit,     /* Parameter interval containing hit point */
    hr3_point_t *npt,   /* Out: "near" hit endpoint */
    hr3_point_t *fpt    /* Out: "far" hit endpoint */
  );
  /* Computes the "near" and "far" points on the segment {org}--{dst},
    corresponding to the low and high points of the parameter 
    interval {hit}.  */

void irt_compute_mean_hit_point
  ( hr3_point_t *npt,
    hr3_point_t *fpt,
    hr3_point_t *mpt
  );
  /* Computes the "mean" (nominal) hit point {mpt}, somewhere
    between {npt} and {fpt}. */

/*** IMPLEMENTATIONS ***/

#define FUDGE 1.0e-6

frgb_t irt_compute_scene_color
  ( scene_t *sc,
    hr3_point_t *eye, 
    r3_t *dir,
    int32_t max_bounces,
    arith_t arith,
    bool_t print_ray,
    epswr_figure_t *epsf
  )
  {
    hr3_point_t inf;  /* Ray's endpoint (at infinity) */

    hr3_point_t npt;  /* "Near" (towards observer) hit point */
    hr3_point_t fpt;  /* "Mean" hit point */
    
    Interval hit;    /* Parameter interval containing hit point */
    int32_t slo, shi;    /* Sign of function before and after {hit} (see zerofind.h) */
    frgb_t color;
    
    inf = hr3_point_at_infinity(dir);

    if (epsf != NULL)
      { /* Print a graph of the fucntion with equal-interval enclosures: */
        irt_debug_ray_graphically(epsf, &(sc->shape), arith, eye, &inf);
      }

    /* Trace the ray: */
    if (print_ray) { fprintf(stderr, "=== begin irt_compute_scene_color ===\n"); }
    if (epsf != NULL) { epswr_set_client_window(epsf, Zero, One, Zero, One); }
    irt_compute_intersection 
      ( &(sc->shape), arith, eye, &inf, 
        &hit, &slo, &shi, print_ray, epsf
      );
    if (epsf != NULL) { epswr_frame(epsf); }
    
    if (hit.lo <= hit.hi)  /* ray intersected object */
      { irt_compute_hit_endpoints (eye, &inf, &hit, &npt, &fpt);
        if (slo == -1)
          { if (sc->looks.has_transp)
              { fatalerror("irt_compute_scene_color: transparent solids not implemented yet"); }
            else
              { fprintf(stderr, "  *** irt_compute_scene_color: observer inside opaque object\n");
                hr3_point_print(stderr, "      eye = ", eye, "%8.4f", "\n");
                hr3_point_print(stderr, "      inf = ", &inf, "%8.4f", "\n");
                fprintf(stderr, "      hit = "); ia_print(stderr, hit); fprintf(stderr, "\n");
                hr3_point_print(stderr, "      npt = ", &npt, "%8.4f", "\n");
                hr3_point_print(stderr, "      fpt = ", &fpt, "%8.4f", "\n");
                fatalerror("aborted");
              }
          }
        color = irt_compute_hit_color (sc, dir, &npt, &fpt, slo, shi, max_bounces, arith, print_ray, epsf);
      }
    else
      { color = sc->looks.background_color; }
    
    if (print_ray)
      { frgb_print(stderr, "pixel color = ", &color, 3, "%7.5f", "\n");
        fprintf(stderr, "=== end irt_compute_scene_color ===\n");
      }
    return color;
  }

frgb_t irt_compute_hit_color
  ( scene_t *sc,
    r3_t *dir, 
    hr3_point_t *npt, 
    hr3_point_t *fpt, 
    int32_t slo, int32_t shi,
    int32_t max_bounces,
    arith_t arith,
    bool_t print_ray,
    epswr_figure_t *epsf
  )
  { hr3_point_t mpt;  /* "Mean" (nominal) hit point */
    
    r3_t nrm;  /* Outwards-pointing surface normal at {mpt} */ 
    r3_t rfl;  /* Reflected viewing ray's direction (for mirror and shiny looks) */

    r3_t npc;  /* Cartesian coordinates of "near" hit point, perturbed */
    r3_t fpc;  /* Cartesian coordinates of "far" hit point, perturbed */
    
    hr3_point_t npp; /* Homogeneous coordinates of "near" hit point, perturbed */
    hr3_point_t fpp; /* Homogeneous coordinates of "far" hit point, perturbed */

    double cos_dir_nrm; /* Cosine of ray direction and normal */
    
    frgb_t color; /* Total apparent color of surface at hit point */
    frgb_t matte_light; /* Matte-scatter component of light */
    frgb_t shiny_light; /* Shiny-scatter component of light */
    frgb_t mirrr_light; /* Mirror component of light */
    
    uint32_t lnum;

    irt_compute_mean_hit_point(npt, fpt, &mpt);
    
    irt_compute_surface_normal(&(sc->shape), arith, &mpt, &nrm);
    cos_dir_nrm = r3_dot(dir, &nrm);
        
    /* Perturb the "near" hit point along normal, towards the observer's side: */
    npc = r3_from_hr3(npt);
    r3_mix_in(copysign(FUDGE, -cos_dir_nrm), &nrm, &npc);
    npp = hr3_from_r3(&npc);
    hr3_L_inf_normalize_point(&npp);
    
    if (sc->looks.has_shiny || sc->looks.has_mirror || (slo == -1)) 
      { /* calculate reflected ray direction */
        /* Used for highlights, mirroring, or internal reflection. */
	double k = -2.0 * r3_dot(dir, &nrm);
	rfl.c[0] = k * nrm.c[0] + dir->c[0];
	rfl.c[1] = k * nrm.c[1] + dir->c[1];
	rfl.c[2] = k * nrm.c[2] + dir->c[2];
        (void) r3_dir(&rfl, &rfl);
      }
    
    if (sc->looks.has_transp)
      {	/* Perturb the "far" hit point along normal, away from observer's: */
        /* Used for refracted rays */
	fpc = r3_from_hr3(fpt);
	r3_mix_in (copysign(FUDGE, cos_dir_nrm), &nrm, &fpc);
	fpp = hr3_from_r3(&fpc);
	hr3_L_inf_normalize_point(&fpp);
      }

    /* Initialize components: */
    matte_light = sc->looks.ambient_color;
    shiny_light = frgb_Zeros;
    mirrr_light = frgb_Zeros;

    for(lnum=0; lnum < sc->num_lights; lnum++)
      { light_t *lt = &(sc->light[lnum]); 
        
        r3_t ldr;             /* Direction of light source */
        double ldist;         /* Distance to light source */
        double rdist;         /* Light's reference distance */
        
        double cos_nrm_ldr;   /* Cosine of angle between normal and light direction */
        
	/* get light's direction and distance */
        r3_sub(&(lt->position), &npc, &ldr);
        ldist = r3_dir(&ldr, &ldr);
        rdist = lt->ref_distance;

	cos_nrm_ldr = r3_dot (&nrm, &ldr);
	if(cos_nrm_ldr * cos_dir_nrm < 0.0)
          { /* Object faces light; compute matte- and shiny-scattered light: */
	    
            frgb_t shadow_factor;  /* Transparency coeffs from light to hit point */
            double inv_sq_factor; /* Inverse square law attenuation factor */
	    frgb_t incident_light; /* Color of incident light */
            
	    /* Compute fudged hit point {nfp} and light position {lfp} to avoid self-shadow: */
            demand(ldist >= 6*FUDGE, "light is too close to the hit point");
            
            r3_t nfc = npc;
            r3_mix_in(+2*FUDGE, &ldr, &nfc);
            hr3_point_t nfp = hr3_from_r3(&nfc);
            hr3_L_inf_normalize_point(&nfp);
            
            r3_t lfc = lt->position;
            r3_mix_in(-2*FUDGE, &ldr, &lfc);
            hr3_point_t lfp = hr3_from_r3(&lfc);
            hr3_L_inf_normalize_point(&lfp);
            
            /* Attenuate light's color by shadow and inverse-square factors: */
            shadow_factor = irt_compute_scene_transparency (sc, &nfp, &lfp, ldist, arith, print_ray, epsf);
	    incident_light = frgb_mul(&(lt->color), &shadow_factor);
            inv_sq_factor = 2.0/((ldist*ldist)/(rdist*rdist) + 1.0);
            incident_light = frgb_scale(inv_sq_factor, &incident_light);
            
            if (print_ray) 
              { frgb_print(stderr, "  light =        ", &incident_light, 3, "%7.5f", "\n"); }
              
            if (! color_is_almost_black(&incident_light, 0.001))
              { /* Accumulate matte-scattered incident light: */
		frgb_t matte_term = frgb_scale(fabs(cos_nrm_ldr), &incident_light);
                matte_light = frgb_add(&matte_light, &matte_term);
                if (sc->looks.has_shiny)
		  { double cos_rfl_ldr = r3_dot (&rfl, &ldr);
		    if (cos_rfl_ldr > sc->looks.shiny_spread)
		      { /* Acumulate shiny-scattered incident light */
                        double ang_coeff; /* Attenuation for angle with reflected ray */
			double grz_coeff; /* Enhancement for grazing reflection */

			ang_coeff = (cos_rfl_ldr - sc->looks.shiny_spread)/(1.0 - sc->looks.shiny_spread);
			ang_coeff = ang_coeff*ang_coeff;

			grz_coeff = (2.0 - cos_nrm_ldr*cos_nrm_ldr)/2.0;

			frgb_t shiny_term = frgb_scale(ang_coeff*grz_coeff, &incident_light);
                        shiny_light = frgb_add(&shiny_light, &shiny_term);
                      }
		  }
              }                    
	  }
      }
    
    color = frgb_Zeros;
    
    if (print_ray) { frgb_print(stderr, "  matte light =  ", &matte_light, 3, "%7.5f", "\n"); }
    frgb_t matte_color = frgb_mul(&matte_light, &(sc->looks.matte_color));
    if (print_ray) { frgb_print(stderr, "  matte color =  ", &matte_color, 3, "%7.5f", "\n"); }
    color = frgb_add(&color, &matte_color);
    
    if (sc->looks.has_shiny)
      { if (print_ray) { frgb_print(stderr, "  shiny light =  ", &shiny_light, 3, "%7.5f", "\n"); }	
        frgb_t shiny_color = frgb_mul(&shiny_light, &(sc->looks.shiny_color));
        if (print_ray) { frgb_print(stderr, "  shiny color =  ", &shiny_color, 3, "%7.5f", "\n"); }	
        color = frgb_add(&color, &shiny_color);
      }

    /* Mirror-like reflection */
    if ((sc->looks.has_mirror) && (max_bounces > 0))
      { mirrr_light = irt_compute_scene_color(sc, &npp, &rfl, max_bounces - 1, arith, print_ray, epsf);
        if (print_ray) { frgb_print(stderr, "mirror light = ", &mirrr_light, 3, "%7.5f", "\n"); }
        frgb_t mirrr_color = mirrr_light; /* frgb_mul(&mirrr_light, &(sc->looks.mirrr_color)); */
        if (print_ray) { frgb_print(stderr, "mirror color = ", &mirrr_color, 3, "%7.5f", "\n"); }
        color = frgb_add(&color, &mirrr_color);
      }

    /* Refraction */
    if (sc->looks.has_transp)
      { fatalerror ("irt_compute_hit_color: transparent solids not implemented yet"); }
    
    if (print_ray) { frgb_print(stderr, "color =        ", &color, 3, "%7.5f", "\n"); }
    return(color);
  }

frgb_t irt_compute_scene_transparency
  ( scene_t *sc,
    hr3_point_t *org,
    hr3_point_t *dst,
    double dist,
    arith_t arith,
    int32_t print_ray,
    epswr_figure_t *epsf
  )
  { frgb_t color;
    Interval hit;
    int32_t slo, shi;
    
    if (print_ray) { fprintf(stderr, "=== begin irt_compute_scene_transparency ===\n"); }
    if (epsf != NULL) { epswr_set_client_window(epsf, Zero, One, Zero, One); }
    irt_compute_intersection
      ( &(sc->shape), arith, org, dst, 
        &hit, &slo, &shi, print_ray, epsf 
      );
    if (epsf != NULL) { epswr_frame(epsf); }

    if ((hit.lo > hit.hi) || (hit.lo >= 1.0))
      { color = frgb_White; }
    else if (sc->looks.has_transp)
      { fatalerror("irt_compute_scene_transparency: transparent solids not implemented");
        color = frgb_Black;
      }
    else
      { color = frgb_Black; }
    if (print_ray) { fprintf(stderr, "=== end irt_compute_scene_transparency ===\n"); }
    return color;
  }

void irt_compute_hit_endpoints
  ( hr3_point_t *org,
    hr3_point_t *dst,
    Interval *hit,
    hr3_point_t *npt,
    hr3_point_t *fpt
  )
  { *npt = hr3_point_mix(1.0-hit->lo, org, hit->lo, dst); hr3_L_inf_normalize_point(npt);
    *fpt = hr3_point_mix(1.0-hit->hi, org, hit->hi, dst); hr3_L_inf_normalize_point(fpt);
  }

void irt_compute_mean_hit_point
  ( hr3_point_t *npt,
    hr3_point_t *fpt,
    hr3_point_t *mpt
  )
  { double nw = 0.5*npt->c.c[0];
    double fw = 0.5*fpt->c.c[0];
    (*mpt) = hr3_point_mix (fw, npt, nw, fpt);
    hr3_L_inf_normalize_point(mpt);
  }

bool_t color_is_almost_black(frgb_t *r, double tol)
  { return((r->c[0] < tol) && (r->c[1] < tol) && (r->c[2] < tol)); }

