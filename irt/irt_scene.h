/* Last edited on 2023-02-22 20:22:05 by stolfi */
/* Scene description tools for the Interval Ray Tracer. */

#ifndef irt_scene_H
#define irt_scene_H

#define _GNU_SOURCE
#include <stdint.h>
#include <stdio.h>

#include <r3.h>
#include <pcode.h>
#include <frgb.h>

#include <irt_evray.h>

#define MAXLIGHTS 4

typedef struct light_t /* Description of a light source: */
  { r3_t position;        /* Position of light source */
    frgb_t color;         /* Color of light */
    r3_t ref_point;       /* Point of scene where light has given color */
    double ref_distance;  /* Computed: dist(ref_point, position). */
  } light_t;
    
typedef struct looks_t /* Description of scene's material, etc: */
  { /* Surface looks: */
    frgb_t ambient_color;      /* Color of ambient light */
    frgb_t matte_color;        /* Matte (Lambertian) scattering coeffs */
    frgb_t shiny_color;        /* Shiny (semi-specular) scattering coeffs */
    double shiny_spread;       /* Cosine of max shiny scattering angle (1=mirror). */
    frgb_t background_color;   /* Color of outer space */
    int32_t max_bounces;           /* Maximum number of mirror reflections */

    /* Flags computed from the looks: */
    int32_t has_shiny;             /* (shiny_color != black) & (shiny_spread != 1)  */
    int32_t has_mirror;            /* (shiny_color != black) & (shiny_spread == 1) */
    int32_t has_transp;            /* FALSE for now */
    
  } looks_t;
  
typedef struct pixel_num_t { int32_t row, col; } pixel_num_t;  /* Pixel indices */

typedef struct view_t /* Viewing and imaging parameters: */
  { r3_t observer;             /* Scene coordinates of observer */
    r3_t focus;                /* Scene coordinates of image center */
    r3_t up;                   /* Appproximate scene direction of screen "up" axis */
    /* Defined based on image size: */
    int32_t image_width;       /* Image width in pixels */
    int32_t image_height;      /* Image height in pixels */
    double pixel_size;         /* Apparent size (radians) of a pixel on the screen */
    r3_t screen_x;             /* Image system X axis vector: disp from a pixel to its right neighbor */
    r3_t screen_y;             /* Image system Y axis vector: disp from a pixel to its top neighbor */
  } view_t;

typedef struct scene_t /* Description of a scene: */
  { shape_t shape;          /* The thing's geometry (a p-code function). */
    looks_t looks;          /* The thing's colors. */
    view_t view;            /* The camera parameters. */
    
    /* Light sources */
    int32_t num_lights;        /* Number of light sources. */
    light_t light[MAXLIGHTS];  /* The lights. */
  } scene_t;

void irt_read_scene(char *pcodeFile, char *parmsFile, scene_t *sc);
  /* Reads the mathematical description of a scene from the p-code file {pcodeFile} 
    and the parameter file {parmsFile},  and stores it into {*sc}.
  
    The fields {image_width}, {image_height}, {pixel_size}, {screen_x}
    and {screen_y} of {sc->view} are left undefined; their values in the
    parameter file are ignored. The caller must use
    {irt_set_image_parameters} below to define those parameters. */
    
void irt_set_image_parameters(scene_t *sc, int32_t NX, int32_t NY, double pix_size, double D);
  /*  The parameters {NX} and {NY} should be the desired image dimensions in pixels.
    The parameter {pix_size} should be the assumed screen pixel size (mm),
    and {D} should be the assumed viewer to screen distance (mm).
    
    The procedure sets the image parameters {image_width} and {image_height} of {sc->view}
    to {NX} and {NY}, and {pixel_size} to {pix_size/D}.  It then computes {screen_x} and
    {screen_y} from those values. */

#endif
