/* {dynstolfi.c} -- exploration of dynamic system, including fundamental region */
/* Last edited on 2008-04-24 15:42:34 by stolfi */
/* Luiz Henrique de Figueiredo & Jorge Stolfi */
/* 20 May 2004 -- lhf -- Created program. */
/* 24 Jan 2008 -- stolfi & lhf -- added fund. region and orbit tracing. */

/* Things that could be improved:
 * * in {fundamental_region}:
 *   * Use a priority queue sorted by Euclidean distance from seed.
 * * in {try_to_choose_pixel_and_paint_orbit}:
 *   * Consider the inverse images as well as the forward ones.
 *   * Use adaptive sampling if {f} stretches the pixel too much.
 *   * Assume that {p} is outside the domain if {f^n(p)} is {+oo,+oo} or {NaN,NaN}.
 * * in {stmaps.c}:
 *   * Return {+oo,+oo} or {NaN,NaN} if {(x,y)} is outside the domain.
 */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#include <gp.h>
#include <dv.h>
#include <stmaps.h>
#include <dynstolfi.h>

/* The function to plot: */
static Mapfunc *f = NULL;  /* To be set by main. */

/* IMPLEMENTATIONS */

int exp_to_color(int e)
{
 return locolor + (((e % ncolors) + ncolors) % ncolors);
}

void show_map_effect_on_point(real x, real y, int e1, int e2)
{gpcolor(seedcolor); gpplot(x,y);
 if ((e1 != 0) || (e2 != 0))
 {int e;
  e = abs(e1);
  while(e--) f(&x,&y,(e1 > 0 ? DIRMAP : INVMAP));
  e = abs(e2);
  while(e--) f(&x,&y,(e2 > 0 ? DIRMAP : INVMAP));
  gpcolor(exp_to_color(e1+e2)); gpplot(x,y);
 }
}

void show_map_with_spot(real xc, real yc, int e1, int e2, real rp, int nsmp)
{if (nsmp <= 0) return;
  real xmax,xmin,ymax,ymin,x,y;
 f(&xmin,&ymin,LOCORNER);
 f(&xmax,&ymax,HICORNER);
 real sx=(rp/(double)winXsize)*(xmax-xmin)/nsmp;
 real sy=(rp/(double)winYsize)*(ymax-ymin)/nsmp;
 int n2o = nsmp*nsmp;
 int n2i = nsmp*nsmp - 3*nsmp;
 int dx,dy;
 for (dx=-(nsmp-1); dx<=+(nsmp-1); dx += 2)
 {for (dy=-(nsmp-1); dy<=+(nsmp-1); dy += 2)
  {/* Use a circular butterfly cloud: */
   int d2 = dx*dx + dy*dy;
   int ok = (d2 <= n2o) && ((d2 >= n2i) || (dx*dy >= 0));
   if(ok)
   {x = xc + dx*sx;
    y = yc + dy*sy;
    show_map_effect_on_point(x,y,e1,e2);
   }
  }
 }
}

void show_map_with_circle(real xc, real yc, int e1, int e2, real rr)
{real xmax,xmin,ymax,ymin,x,y;
 f(&xmin,&ymin,LOCORNER);
 f(&xmax,&ymax,HICORNER);
 real rx=rr*(xmax-xmin)/2;
 real ry=rr*(ymax-ymin)/2;
 real rc = hypot(rx,ry);
 /* Compute pixel size: */
 real dx=(xmax-xmin)/((double)winXsize);
 real dy=(ymax-ymin)/((double)winYsize);
 real dmin = (dx < dy ? dx : dy);
 int nsmp;
 /* Compute sample count for dense coverage of cirles: */
 nsmp = (int)ceil((2*M_PI*rc)/(0.5*dmin));
 double dt = (2*M_PI)/nsmp; /* Angular step. */
 int nsub = 3;
 int i,k;
 for (k=1; k <= nsub; k++)
 {real u = ((double)k)/((double)nsub);
  for (i=0; i < nsmp; i++)
  {double t = i*dt;
   x = xc + u*rc*cos(t);
   y = yc + u*rc*sin(t);
   show_map_effect_on_point(x,y,e1,e2);
  }
 }
 /* Compute sample count for dense coverage of rays: */
 nsmp = (int)ceil(rc/(0.5*dmin));
 nsmp = 16*((nsmp+15)/16) + 1;
 int nray = 8;
 for (k=0; k < nray; k++)
 {real t = 2*M_PI*((double)k)/((double)nray);
  real C = cos(t), S = sin(t);
  for (i=0; i <= nsmp; i++)
  {if((k%2 == 0) || ((i/8)%2 == 0)) 
   {double u = ((double)i)/((double)nsmp);
    x = xc + u*rc*C;
    y = yc + u*rc*S;
    show_map_effect_on_point(x,y,e1,e2);
   }
  }
 }
}

void show_map_with_square(real xc, real yc, int e1, int e2, real rr)
{real xmax,xmin,ymax,ymin;
 f(&xmin,&ymin,LOCORNER);
 f(&xmax,&ymax,HICORNER);
 real rx=rr*(xmax-xmin)/2;
 real ry=rr*(ymax-ymin)/2;
 real rs = hypot(rx,ry)/M_SQRT2;
 /* Compute pixel size: */
 real dx=(xmax-xmin)/((double)winXsize);
 real dy=(ymax-ymin)/((double)winYsize);
 real dmin = (dx < dy ? dx : dy);
 /* Compute sample count per {rs} for dense coverage: */
 int nsmp = (int)ceil(rs/(0.5*dmin));
 int nsub = 4;
 int i,k;
 for (k=-nsub; k <= +nsub; k += 2)
 {real u = ((double)k)/((double)nsub)*rs;
  for (i=-nsmp; i <= +nsmp; i++)
  {double v = ((double)i)/((double)nsmp)*rs;
   show_map_effect_on_point(xc + u,yc + v,e1,e2);
   show_map_effect_on_point(xc + v,yc + u,e1,e2);
  }
 }
}

/* GATHERING ORBITS */

/* The chosen vertex set {A}: */
#define N 1000
static Amark AM[N][N]; /* Set of chosen pixels. */

void clearAM(void)
{
 memset(AM,0,sizeof(AM));
}   

Amark getAM(int ix, int iy)
{
 if (ix>=0 && iy>=0 && ix<winXsize && iy<winYsize) 
  return AM[ix][iy];
 else
  return EXCLUDED;
}

void setAM(int ix, int iy, Amark what)
{
 if (ix>=0 && iy>=0 && ix<winXsize && iy<winYsize) AM[ix][iy]=what;
}

int orbit_is_new(real x0, real y0, int ngen, int nsmp)
{assert(nsmp > 0);
 real x,y,px,py;
 int dx,dy;
 real sx = 0.5/((double)nsmp); /* Half-step of sub-pixel sampling. */
 real sy = 0.5/((double)nsmp); /* Half-step of sub-pixel sampling. */
 /* Get the indices [px][py] of the pixel {P} containing (x0,y0): */
 px=x0; py=y0; gpview(&px,&py);
 /* Check whether (x0,y0) is already in {A}: */
 if (getAM(px,py) != VIABLE) return 0;
 /* Check whether (x0,y0) is already in {A}: */
 int e;
 for(dx = -(nsmp-1); dx <= +(nsmp-1); dx += 2)
  for(dy = -(nsmp-1); dy <= +(nsmp-1); dy += 2)
  {/* Offset the pixel center point by a fraction of pixel: */
   px = x0; py = y0; gpview(&px,&py);
   x = px + dx*sx; y = py + dy*sy; gpunview(&x,&y);
   /* Check the orbit of (x,y). */
   e=0;
   while (e < ngen)
   {f(&x,&y,DIRMAP); e++;
    px=x; py=y; gpview(&px,&py);
    if (getAM(px,py) == CHOSEN) return 0;
   }
  }
 return 1;
}

void mark_and_paint_orbit(real x0, real y0, int ngen, int nsmp)
{assert(nsmp > 0);
 real x,y,px,py;
 int dx,dy;
 real sx = 0.5/((double)nsmp);
 real sy = 0.5/((double)nsmp);
 int e;
 Amark am;
 for(dx = -(nsmp-1); dx <= +(nsmp-1); dx += 2)
  for(dy = -(nsmp-1); dy <= +(nsmp-1); dy += 2)
  {/* Offset the pixel center point by a fraction of pixel: */
   px = x0; py = y0; gpview(&px,&py);
   x = px + dx*sx; y = py + dy*sy; gpunview(&x,&y);
   /* Paint and mark the forward orbit of (x,y). */
   e = 0;
   while (e < ngen)
   {f(&x,&y,DIRMAP); e++;
    px=x; py=y; gpview(&px,&py);
    am = getAM(px,py);
    assert(am != CHOSEN); /* Paranoia... */
    if (am == VIABLE)
    {setAM(px,py,EXCLUDED);
     gpcolor(exp_to_color(e));
     gpplot(x,y);
    }
   
   }
  }
 /* Paint and mark the pixel {P}: */
 gpcolor(seedcolor);
 gpplot(x0,y0);
 px=x0; py=y0; gpview(&px,&py);
 setAM(px,py,CHOSEN);
}

int try_to_choose_pixel_and_paint_orbit(real x0, real y0, int ngen)
{
 int nsmp = 5; /* Will use {nsmp^2} samples in pixel. */
 if(! orbit_is_new(x0,y0,ngen,nsmp)) return 0;
 mark_and_paint_orbit(x0,y0,ngen,nsmp);
 return 1;
}

/* FUNDAMENTAL REGION */
static Qmark QM[N][N]; /* Set of pixels that were enqueued. */
static int Q[N*N];
static int QB,QE;

void clearQM(void)
{
 memset(QM,0,sizeof(QM));
}   

void mark_and_stack(int ix, int iy)
{
 if (ix>=0 && iy>=0 && ix<winXsize && iy<winYsize && (QM[ix][iy]==NOTVISITED))
 {
  int k=N*ix+iy;
  /* printf("  stacked %d %d %d\n" , k, ix, iy); */
  Q[QE]=k; QE++;
  QM[ix][iy]=VISITED;
 }
}

void fundamental_region(real x0, real y0, int ngen)
{
 real px,py;
 int k;
 printf("\ndoing fundamental_region  seed = (%g %g)",x0,y0);
 /* Get indices [px][py] of pixel containing (x0,y0): */
 px=x0; py=y0; gpview(&px,&py);
 printf(" = [%g][%g]\n",px,py);
 /* Push that pixel into the queue: */
 QB=QE=0;
 mark_and_stack(px,py);
 while (QE>QB)
 {
  /* Get the next candidate from the queue: */
  /* printf("Q: %d .. %d %d\r",QB,QE,QE-QB); */
  k=Q[QB]; QB++;
  px=k/N; py=k%N;
  /* printf("  popped %d %g %g ...",k,px,py); */
  gpunview(&px,&py);
  if (try_to_choose_pixel_and_paint_orbit(px,py,ngen)) 
  {
   px=k/N; py=k%N;
   /* printf(" (chosen)"); */
   mark_and_stack(px-1,py+0);
   mark_and_stack(px+1,py+0);
   mark_and_stack(px+0,py-1);
   mark_and_stack(px+0,py+1);
#if CONNECT8
   mark_and_stack(px-1,py-1);
   mark_and_stack(px-1,py+1);
   mark_and_stack(px+1,py-1);
   mark_and_stack(px+1,py+1);
#endif
  }
  /* printf("\n"); */
 }
 printf("\nregion done Q: %d .. %d %d\n",QB,QE,QE-QB);
}

int main(int argc, char* argv[])
{
 char *fname;
 int ngen; 
 real xmin,ymin,xmax,ymax;
 
 /* Get the function: */
 fname = "henon"; 
 if (argc>1) fname = argv[1];
 fprintf(stderr, "using function \"%s\"\n", fname);
 f = func_from_name(fname);
 
 /* Get the generation count for {orbit},{fundamental_region}. */
 ngen=100; 
 if (argc>2) sscanf(argv[2],"%d",&ngen);
 
 /* Get the window of interest: */
 f(&xmin,&ymin,LOCORNER);
 f(&xmax,&ymax,HICORNER);

 /* Graphics window setup: */
 char *title = NULL;
 asprintf(&title, "%s (^%d) [%g %g]×[%g %g]", fname,ngen,xmin,xmax,ymin,ymax); 
 gpopen(title);
 gpwindow(xmin,xmax,ymin,ymax);
 gppalette(0,"black");
 gppalette(1,"white");
 gppalette(2,"red");
 gppalette(3,"orange");
 gppalette(4,"yellow");
 gppalette(5,"green");
 gppalette(6,"cyan");
 gppalette(7,"magenta");
 while (1)
 {
  real x,y;
  real spotrad =  4.0; /* Radius of showation spot (pixels). */
  real circlerad = 0.25; /* Relative radius of showation circle. */
  real squarerad = 0.25; /* Relative radius of showation square. */
  char* s=gpevent(1,&x,&y);
  
  /* Showing the map: */
  
  if (s[0]=='k' && s[1]=='d')
  {/* Show the effect of the direct map, small spot, cumulative: */
   show_map_with_spot(x,y, +1, 0, spotrad, 10);
   continue;
  }
  if (s[0]=='k' && s[1]=='i')
  {/* Show the effect of the inverse map, small spot, cumulative: */
   show_map_with_spot(x,y, -1,  0, spotrad, 10);
   continue;
  }
  if (s[0]=='k' && s[1]=='b')
  {/* Show the effect of {f(f^{-1}(x))}, small spot, cumulative: */
   show_map_with_spot(x,y, -1, +1, spotrad, 10);
   continue;
  }

  if (s[0]=='k' && s[1]=='c')
  {/* Show the effect of the direct map, big circle, cumulative, : */
   show_map_with_circle(x,y, +1, 0, circlerad);
   continue;
  }
  
  if (s[0]=='k' && s[1]=='s')
  {/* Show the effect of the direct map, big circle, cumulative, : */
   show_map_with_square(x,y, +1, 0, squarerad);
   continue;
  }
  
  if ((s[0]=='b' || s[0]=='m') && s[1]=='1')
  {/* Show the effect of the direct map, small spot, labile: */
   gpclear(0);
   if (s[0]=='b') clearAM(); clearQM();
   if (s[2]=='+') show_map_with_spot(x,y, +1, 0, spotrad, 10);
   continue;
  }
  if ((s[0]=='b' || s[0]=='m') && s[1]=='3')
  {/* Show the effect of the inverse map, small spot, labile: */
   gpclear(0);
   if (s[0]=='b') clearAM(); clearQM();
   if (s[2]=='+') show_map_with_spot(x,y, -1, 0, spotrad, 10);
   continue;
  }
  if ((s[0]=='b' || s[0]=='m') && s[1]=='2')
  {/* Show the effect of {f(f^{-1}(x))}, small spot, labile: */
   gpclear(0);
   if (s[0]=='b') clearAM(); clearQM();
   if (s[2]=='+') show_map_with_spot(x,y, -1, +1, spotrad, 10);
   continue;
  }

  /* Orbits and such: */
  
  if (s[0]=='k' && s[1]=='o')
  {/* Show the forward orbit: */
   /* gpmark(3,"B"); */
   try_to_choose_pixel_and_paint_orbit(x,y,ngen);
   /* gpmark(0,"."); */
   continue;
  }
  if (s[0]=='k' && s[1]=='r')
  {/* Find a fundamental region: */
   fundamental_region(x,y,ngen);
   continue;
  }

  /* Miscellaneous: */
  
  if (s[0]=='k' && s[1]=='z') 
  {/* Clear the screen and work tables: */
   gpclear(0); clearAM(); clearQM();
   continue;
  }
  if (s[0]=='k' && s[1]=='q') 
  {/* Quit: */
   gpclose(0);
   break;
  }
  if (s[0]=='k' && s[1]=='w')
  {/* Print the mouse coords: */
   printf("%g\t%g\n",x,y);
   continue;
  }
  
 }
 return 0;
}
