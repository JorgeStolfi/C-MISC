#ifndef MAINDECL_H
#define MAINDECL_H

#include "raymath.h"
#include "constants.h"
#include "box.h"
#include "sphere.h"
#include "superq.h"
#include "triangle.h"
#include "object.h"
#include "light.h"
#include "shade.h"

/* intersection routines */
double (*objint[])(t_3d *pos, t_3d *ray, t_object obj) 
  = { intsph, intbox, inttri, intsup };

/* normal routines */
void (*objnrm[])(t_3d *pos, t_object obj, t_3d nrm) 
  = { nrmsph, nrmbox, nrmtri, nrmsup };

/* global variables */
int 		nlight;                 /* presently in use 	*/
int 		lightlim = LIGHTS;      /* maximun declared 	*/
t_light 	light[LIGHTS];		/* array of lights  	*/
int 		nobject;                /* presently in use 	*/
int 		objectlim = OBJECTS;    /* maximun declared 	*/
t_object 	object[OBJECTS];      	/* array of objects 	*/
int 		nsurface;               /* presently in use 	*/
int 		surfacelim = SURFACES;	/* maximun declared 	*/
t_surface 	surface[SURFACES];   	/* array of surfaces 	*/
int 		sizex, sizey;           /* image sizes 		*/
t_3d 		eyep, lookp, up;        /* view definition 	*/
double 		hfov, vfov;             /* field of view 	*/
int 		level, maxlevel;        /* reflection levels 	*/
char 		*outfilename[60];       /* pixel file name 	*/
t_color 	background; 		/* background color 	*/

#endif







