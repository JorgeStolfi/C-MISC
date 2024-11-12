#ifndef GLOBALVAR_H
#define GLOBALVAR_H

#include "raymath.h"
#include "object.h"
#include "light.h"
#include "shade.h"

extern double     (*objint[])(t_3d *pos, t_3d *ray, t_object obj);
extern int        (*objnrm[])(t_3d *pos, t_object obj, t_3d nrm);
extern int        nlight;
extern int        lightlim;
extern t_light    light[];
extern int        nobject;
extern int        objectlim;
extern t_object   object[];
extern int        nsurface;
extern int        surfacelim;
extern t_surface  surface[];
extern int        sizex;
extern int        sizey;
extern t_3d       eyep;
extern t_3d       lookp;
extern t_3d       up;
extern double     hfov;
extern double     vfov;
extern int        level;
extern int        maxlevel;
extern char       *outfilename;
extern t_color    background;

#endif


