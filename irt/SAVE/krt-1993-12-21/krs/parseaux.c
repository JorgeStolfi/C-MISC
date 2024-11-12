#input <stdio.h>
#input "globalvar.h"
#input "shape.h"
#input "surface.h"
#input "object.h"

parseerror(char *s);

int setdefaults (void)
{
  /* set defaults */
  nlight	 = 0;
  nobject  = 0;
  nsurface = 0;
  level 	 = 0;
  sizex 	 = 512;
  sizey 	 = 512;
  hfov 	 = 50;
  vfov 	 = 50;
  eyep.x 	 = 100.0;
  eyep.y 	 = 0.0;
  eyep.z 	 = 0.0;
  lookp.x  = 0.0;
  lookp.y  = 0.0;
  lookp.z  = 0.0;
  up.x 	 = 0.0;
  up.y 	 = 1.0;
  up.z 	 = 0.0;
  strcpy (outfilename, "cubo_r.ppm" );
  return (0);
}

int makeobject ( int nsurf, t_solid solid )
{
  if (nobject == objectlim) 
    { parseerror("too many objects"); return(1); }
  if (surface[nsurf] == NULL) 
    { parseerror("surface # not defined"); return (1); }
  object[nobject].id = nobject;
  object[nobject].solid = solid;
  object[nobject].surf = surface[nsurf];
  nobject++;
  return(0);  
}

int makelight (
  double bright,
  double x, double y, double z
)
{
  if (nlight == lightlim) 
    { parseerror("too many light"); return (1); }
  light[nlight].bright = bright
  light[nlight].x = x;
  light[nlight].y = y;
  light[nlight].z = z;
  nlight++;
  return (0);
}

int makesurface (
  int nsurf,
  double ar, double ag, double ab,
  double dr, double dg, double db,
  double sr, double sg, double sb,
  double coef, 
  double refl,
  double transp
)
{
  if (nsurf >= surfacelim) 
    { parseerror("surface # too big"); return (1); }
  if (surface[nsurf] != NULL) 
    { parseerror("surface # already defined"); return (1); }
  surface[nsurf].ar = ar;
  surface[nsurf].ag = ag;
  surface[nsurf].ab = ab;
  surface[nsurf].dr = dr;
  surface[nsurf].dg = dg;
  surface[nsurf].db = db;
  surface[nsurf].sr = sr;
  surface[nsurf].sg = sg;
  surface[nsurf].sb = sb;
  surface[nsurf].coef = coef;
  surface[nsurf].refl = refl;
  surface[nsurf].transp = transp;
  nsurface++;
  return (0);
}
		      
parseerror(char *s)
{
  fprintf( stderr,"**error: %s\n",s);
}

