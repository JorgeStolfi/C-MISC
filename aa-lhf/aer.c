/*
* aer.c
* adaptive enumeration for rendering
* Luiz Henrique de Figueiredo
* 07 Nov 94
*/
/* Last edited on 2023-01-14 12:20:00 by stolfi */

#if 1
#define FRONTFIRST
#define SHADE
#define RENDER
#define ZBUFFER
#endif

#undef DEBUG

#define XY

#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef RENDER
#include <gp.h>
#else
#undef SHADE
#undef ZBUFFER
#include <nogp.h>
#endif

#include <arith.h>
#include <fn.h>

typedef struct
{
 real xmin,xmax,ymin,ymax,zmin,zmax;
 int x,y,z,n;
} Cell;

int putenv(char *string);

void explore(Cell* c);
int haszero(Cell* c);
int toosmall(Cell* c);
void image(Cell* c, real* fmin, real* fmax);
void plot(Cell* c);

#undef N
#define	N 1024 

static short Z[N][N];
static int level=0;
static unsigned long int Ktoosmall, Kimage, Kplot;
static unsigned long int Kimage, Ktoosmall, Kplot;

#define Explore(N,X,Y,Z,XMIN,XMAX,YMIN,YMAX,ZMIN,ZMAX)\
do					\
{					\
 Cell c; 				\
 c.n=N;					\
 c.x=X;	c.xmin=XMIN;	c.xmax=XMAX;	\
 c.y=Y;	c.ymin=YMIN;	c.ymax=YMAX;	\
 c.z=Z;	c.zmin=ZMIN;	c.zmax=ZMAX;	\
 explore(&c);				\
} while(0)

int main(int argc, char* argv[])
{
 int n;
 if (argc>1 && sscanf(argv[1],"%d",&level)==1) /* ok */; else level=6;
 n=1<<level;
 {
  int x,y;
  for (x=0; x<N; x++)
   for (y=0; y<N; y++)
    Z[x][y]=-1;
 }
 {
  char g[]="GP=1000x1000";
  sprintf(g,"GP=%dx%d",n,n);
  putenv(g);
 }
 dvopen("adaptive enumeration for rendering");
 {
  int c;
  for (c=0; c<=100; c++)
  {
   char s[]="grey100";
   sprintf(s,"grey%d",c);
   dvpalette(c,s);
  }
  for (c=0; c<=100; c++)
  {
   dvrgb(c+101,c/100.0,0,0);
  }
  dvpalette(251,"grey");
  dvpalette(252,"yellow");
  dvpalette(254,"blue");
  dvpalette(255,"red");
#ifdef PALETTE
  for (c=0; c<=255; c++)
  {
   dvcolor(c);
   dvline(c,0,c,n);
  }
  dvwait(-1);
#endif
  dvcolor(100);
  dvbox(0,n,0,n);
  dvcolor(254);
 }
 Explore(n,0,0,0,-3,3,-3,3,-3,3);
 dvclose(1);
 printf(
  "%s %s level %d visited %lu grey %lu leaves %lu\n",
  fnname(),aa_id(),level,Kimage,Ktoosmall,Kplot
 );
 return 0;
}

void explore(Cell* c)
{
 if (haszero(c))
 {
  if (toosmall(c))
   plot(c);
  else
  {
   int n=c->n/2;
   int x=c->x; real xmin=c->xmin; real xmax=c->xmax; real xmid=(xmin+xmax)/2;
   int y=c->y; real ymin=c->ymin; real ymax=c->ymax; real ymid=(ymin+ymax)/2;
   int z=c->z; real zmin=c->zmin; real zmax=c->zmax; real zmid=(zmin+zmax)/2;
 --level;
#ifdef FRONTFIRST
   Explore(n,x+n,y  ,z+n,xmid,xmax,ymid,ymax,zmid,zmax);
   Explore(n,x  ,y  ,z+n,xmin,xmid,ymid,ymax,zmid,zmax);
   Explore(n,x  ,y+n,z+n,xmin,xmid,ymin,ymid,zmid,zmax);
   Explore(n,x+n,y+n,z+n,xmid,xmax,ymin,ymid,zmid,zmax);
   Explore(n,x  ,y  ,z  ,xmin,xmid,ymid,ymax,zmin,zmid);
   Explore(n,x+n,y  ,z  ,xmid,xmax,ymid,ymax,zmin,zmid);
   Explore(n,x  ,y+n,z  ,xmin,xmid,ymin,ymid,zmin,zmid);
   Explore(n,x+n,y+n,z  ,xmid,xmax,ymin,ymid,zmin,zmid);
#else
   Explore(n,x  ,y  ,z  ,xmin,xmid,ymid,ymax,zmin,zmid);
   Explore(n,x  ,y+n,z  ,xmin,xmid,ymin,ymid,zmin,zmid);
   Explore(n,x+n,y  ,z  ,xmid,xmax,ymid,ymax,zmin,zmid);
   Explore(n,x+n,y+n,z  ,xmid,xmax,ymin,ymid,zmin,zmid);
   Explore(n,x  ,y  ,z+n,xmin,xmid,ymid,ymax,zmid,zmax);
   Explore(n,x  ,y+n,z+n,xmin,xmid,ymin,ymid,zmid,zmax);
   Explore(n,x+n,y  ,z+n,xmid,xmax,ymid,ymax,zmid,zmax);
   Explore(n,x+n,y+n,z+n,xmid,xmax,ymin,ymid,zmid,zmax);
#endif
 ++level;
  }
 }
}

int haszero(Cell* c)
{
 real fmin,fmax;
 image(c,&fmin,&fmax);
 return (fmin<=0) && (fmax>=0);
}

int toosmall(Cell* c)
{
++Ktoosmall;
 return (level<=0);
}

void image(Cell* c, real* fmin, real* fmax)
{
 AAform x,y,z,f;
++Kimage;
 aa_open();
#ifdef XY
 aa_interval(x,c->xmin,c->xmax);	/* x = [ xmin .. xmax ] */
 aa_interval(y,c->ymin,c->ymax);	/* y = [ ymin .. ymax ] */
 aa_interval(z,c->zmin,c->zmax);	/* z = [ zmin .. zmax ] */
#else
 aa_interval(x,c->xmin,c->xmax);	/* y = [ zmin .. zmax ] */
 aa_interval(y,c->zmin,c->zmax);	/* x = [ ymin .. ymax ] */
 aa_interval(z,c->ymin,c->ymax);	/* x = [ xmin .. xmax ] */
#endif
 aa_debug(x);
 aa_debug(y);
 aa_debug(z);
 fneval(x, y, z, f); aa_debug(f);
 aa_range(f,fmin,fmax);                 /* f= [ fmin .. fmax ] */
 aa_close();

#ifdef DEBUG
 fprintf(stderr,"image: (%d) %d,%d,%d %g %g\n",level,c->x,c->y,c->z,*fmin,*fmax);
#endif

#ifdef SHADE
 {
  real u=f[1]/x[1]; 
  real v=f[2]/y[2]; 
  real w=f[3]/z[3]; 
  real n=sqrt(u*u+v*v+w*w);
  int c;
#ifndef XY
  w=v;
#endif
  if (n==0.0 || w<0.0)
   c=101+0.5+(100.0*(-w)/n);
  else
   c=0.5+(100.0*w/n);
#if 0
fprintf(stderr,"%lg %lg %lg %lg %d\n",u,v,w,n,c);
#endif
  dvcolor(c);
#if 0
  if (c>30) dvcolor(251); else dvcolor(252);
#endif
 }
#endif
}

void plot(Cell* c)
{
#ifdef RENDER
 int x=c->x;
 int y=c->y;
 int z=c->z;
#ifdef ZBUFFER
 if (Z[x][y]<z)
 {
  Z[x][y]=z;
  dvplot(x,y);
 }
#else
  dvplot(x,y);
#endif
#else
#endif
#if 0
  dvflush();
  printf("%d %d %d %d\n",c->x,c->y,c->z,dvcolor(0));
#endif
++Kplot;
}
