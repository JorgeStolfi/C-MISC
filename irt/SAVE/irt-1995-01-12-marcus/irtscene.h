#ifndef IRTSCENE_H
#define IRTSCENE_H

#include <stdio.h>
#include <r3.h>
#include <rgb.h>
#include <irtinter.h>
#include <pcode.h>

#define MAXLIGHTS 4

typedef struct { /* Description of a light source: */
    r3_t position;        /* Position of light source */
    rgb_t color;          /* Color of light */
    r3_t ref_point;       /* Point of scene where light has given color */
    double ref_distance;  /* Computed: dist(ref_point, position). */
  } light_t;
    
typedef struct { /* Description of scene's material, etc: */
    /* Surface looks: */
    rgb_t ambient_color;       /* Color of ambient light */
    rgb_t matte_color;         /* Matte (Lambertian) scattering coeffs */
    rgb_t shiny_color;         /* Shiny (semi-specular) scattering coeffs */
    double shiny_spread;       /* Cosine of max shiny scattering angle (1=mirror). */
    rgb_t background_color;    /* Color of outer space */
    int max_bounces;           /* Maximum number of mirror reflections */

    /* Flags computed from the looks: */
    int has_shiny;             /* (shiny_color != black) & (shiny_spread != 1)  */
    int has_mirror;            /* (shiny_color != black) & (shiny_spread == 1) */
    int has_transp;            /* FALSE for now */
    
  } looks_t;
  
typedef struct { int row, col; } pixel_num_t;  /* Pixel indices */

typedef struct { /* Viewing and imaging parameters: */
    int image_width;           /* Image width in pixels */
    int image_height;          /* Image height in pixels */
    double pixel_size;         /* Apparent size (radians) of a pixel on the screen */
    r3_t observer;             /* Scene coordinates of observer */
    r3_t focus;                /* Scene coordinates of image center */
    r3_t up;                   /* Appproximate scene direction of screen "up" axis */
    r3_t screen_x;             /* Displacement from a pixel to its right neighbor */
    r3_t screen_y;             /* Displacement from a pixel to its top neighbor */
  } view_t;

typedef struct { /* Description of a scene: */
    char *name;             /* The name. */
    shape_t shape;          /* The thing's geometry. */
    looks_t looks;          /* The thing's colors. */
    view_t view;            /* The camera. */
    
    /* Light sources */
    int num_lights;         /* Number of light sources. */
    light_t light[MAXLIGHTS];  /* The lights. */

    pixel_num_t print_ray;  /* Print a trace of this ray's evaluation. */ 
    pixel_num_t plot_ray;   /* Plot the shape function along this ray. */

  } scene_t;

void irt_read_scene (char *scene_name, scene_t *sc);
 /*
   Reads a scene description from the file <scene_name>.parm. */

#endif
