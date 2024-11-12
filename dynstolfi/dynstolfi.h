#ifndef dynstolfi_H
#define dynstolfi_H

/* {dynstolfi.h} -- definitions and prototypes for {dynstolfi.c} */
/* Last edited on 2008-04-24 15:42:10 by stolfi */

#include <dv.h>
#include <gp.h>

/* SHOWING THE MAP'S EFFECT ON A SINGLE POINT */

#define seedcolor 1
/* Color to use for the seed point(s), chosen region, etc.. */

#define locolor 2
#define hicolor 7
#define ncolors (hicolor + 1 - locolor)

int exp_to_color(int e);
/* Maps an iteration exponent of {e} to colors {locolor .. hicolor},
 cyclically. Namely maps {e = 0,1,2,...} to
 {locolor+0,locolor+1,locolor+2,..} and {e = -1,-2,...} to
 {hicolor-1,hicolor-2,...}. */
 
void show_map_effect_on_point(real x, real y, int e1, int e2);
/* Paints the point {p = (x,y)} with {seedcolor} and 
 {f^e2(f^e1(p))} with color that depends on {e1+e2}.
 
 Note that if {e1} and {e2} have opposite signs, the result of
 {f^e2(f^e1(p))} need not be identical to {f^{e1+e2}}, because
 there may be convergence (or bugs!) in {f}. */

/* SHOWING THE MAP'S EFFECT ON A POINT CLOUD */

#define winXsize (512)
#define winYsize (512)
/* Assumed window dimensions in pixels.  */
/* !!! Should get them from the window itself! */

void show_map_with_spot(real xc, real yc, int e1, int e2, real rp, int nsmp);
/* Applies {show_map_effect_on_point(x,y,e1,e2)} on a small
 spot-like cloud of {nsmp^2} sample points with center at {(xc,yc)}
 and radius {rp} (in pixels). */

void show_map_with_circle(real xc, real yc, int e1, int e2, real rr);
/* Applies {show_map_effect_on_point(x,y,e1,e2)} on a set of sample
 points chosen densely along a circle with center at {(xc,yc)} and
 radius {rr} (relative to the plot window radius). */

void show_map_with_square(real xc, real yc, int e1, int e2, real rr);
/* Applies {show_map_effect_on_point(x,y,e1,e2)} on a set of sample
 points chosen densely along a square with center at {(xc,yc)} and
 radius {r} (relative to the plot window radius). */

/* GATHERING ORBITS */

/* These procedures manipulate a dynamic set {A} of pixels.
   A pixel can be in one of the following states w.r.t this set: */
typedef enum 
{ VIABLE=0,  /* Pixel is not in {A,f(A),ff(A), ..}. */
  CHOSEN,    /* Pixel is in {A}.*/
  EXCLUDED   /* Pixel is in one of {f(A),ff(A),...}. */
} Amark;

void clearAM(void);
/* Makes the set {A} empty. */

Amark getAM(int ix, int iy);
/* Gets the status of pixel {ix,iy} w.r.t the set {A}. */

void setAM(int ix, int iy, Amark what);
/* Marks the pixel {ix,iy} with the status {what}. */

int orbit_is_new(real x0, real y0, int ngen, int nsmp);
/* Returns TRUE iff the pixel {P} that contains (x0,y0)
  is not in the chosen set {A},  does not intersect any
  forward image of {A} (as previously marked), and 
  none of its {ngen} forward images seem to intersect {A}.
  Uses {nsmp^2} sample points inside the pixel {P}. */

void mark_and_paint_orbit(real x0, real y0, int ngen, int nsmp);
/* Paints the pixel {P} that contains the point (x0,y0) with color 1,
  then paints colors 2..7 (cyclically) all pixels that seem to meet
  the {ngen} forward images of {P}. Uses {nsmp^2} sample points in the
  pixel {P}.
  
  Also marks all pixels that meet those forward images as EXCLUDED,
  and marks {P} as CHOSEN. (Presumably, {P} is VIABLE and all pixels
  that meet the {ngen} forward images of {P} are VIABLE or EXCLUDED). */

int try_to_choose_pixel_and_paint_orbit(real x0, real y0, int ngen);
/* Tries adding to the set {A} the pixel that contains the point
  (x0,y0). If it succeeds, paints that pixel and its {ngen} forward
  images, sets the {AM} table of {P} to CHOSEN, sets the {AM} table of
  all pixels that meet those images to EXCLUDED, and returns TRUE. If
  it fails, paints nothing, sets nothing, and returns FALSE. */ 

/* FUNDAMENTAL REGION */

#define CONNECT8 (0)
/* Define as (1) to use 8-neighbor connectivity. */

/* These procedures use a global queue {Q} of pixels. A pixel can be
  in the following states relative to this queue: */

/* The queue for the flood fill algorithm: */
typedef enum 
{ NOTVISITED=0, /* Pixel was never in the queue. */
  VISITED       /* Pixel is (or has been) in the queue. */
} Qmark;


void clearQM(void);
/* Makes the queue {Q} empty. */

void mark_and_stack(int ix, int iy);
/* If the pixel {ix,iy} has never been in the queue, adds it. */

void fundamental_region(real x0, real y0, int ngen);
/* Computes a connected fundamental region that
  contains point {x0,y0}.  The pixels in the region are
  added to the set {A}, and their orbits are painted as above. */

#endif
