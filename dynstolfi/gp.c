/* See {gp.h}. */
/* Luiz Henrique de Figueiredo (lhf@visgraf.impa.br) - 28 Apr 93 */
/* Last edited on 2008-04-24 15:08:12 by stolfi */

#include <gp.h>

static struct
{
 Box w,v,d;
 real ax,bx;
 real ay,by;
} gp=
{
 {0.0, 1.0, 0.0, 1.0, 1.0, 1.0},
 {0.0, 1.0, 0.0, 1.0, 1.0, 1.0},
 {0.0, 1.0, 0.0, 1.0, 1.0, 1.0},
  1.0, 0.0,
  1.0, 0.0,
};

real gpopen(char* name)
{
 gp.d=*dvopen(name);
 if (gp.d.xu>gp.d.yu)
 {
  gp.d.xu/=gp.d.yu;
  gp.d.yu=1.0;
 }
 else
 {
  gp.d.yu/=gp.d.xu;
  gp.d.xu=1.0;
 }
 gpwindow(0.0,1.0,0.0,1.0);
 gpviewport(0.0,1.0,0.0,1.0);
 gppalette(0,"white");				/* black on white */
 gppalette(1,"black");
 gpcolor(1);
 return gp.d.xu/gp.d.yu;
}

int gppoint(real x, real y)
{
 gpview(&x,&y);
 return dvpoint(x,y);
}

void gpplot(real x, real y)
{
 gpview(&x,&y);
 dvplot(x,y);
}

void gpline(real x1, real y1, real x2, real y2)
{
 gpview(&x1,&y1);
 gpview(&x2,&y2);
 dvline(x1,y1,x2,y2);
}

void gpbox(real xmin, real xmax, real ymin, real ymax)
{
 gpview(&xmin,&ymin);
 gpview(&xmax,&ymax);
 if (xmin>xmax) { real t=xmin; xmin=xmax; xmax=t; }
 if (ymin>ymax) { real t=ymin; ymin=ymax; ymax=t; }
 dvbox(xmin,xmax,ymin,ymax);
}

void gptri(real x1, real y1, real x2, real y2, real x3, real y3)
{
 gpbegin('f');
  gppoint(x1,y1);
  gppoint(x2,y2);
  gppoint(x3,y3);
 gpend();
}

void gpcircle(real x, real y, real r)
{
 real x1=x+r;
 real y1=y;
 gpview(&x,&y);
 gpview(&x1,&y1);
 if (x1>x) r=x1-x; else r=x-x1;
 dvcircle(x,y,r);
}

void gptext(real x, real y, char* s, char* mode)
{
 gpview(&x,&y);
 dvtext(x,y,s,mode);
}

char* gpevent(int wait, real* x, real* y)
{
 int ix,iy;
 char* r=dvevent(wait,&ix,&iy);
 *x=ix; *y=iy;
 gpunview(x,y);
 return r;
}

real gpwindow(real xmin, real xmax, real ymin, real ymax)
{
 gp.w.xmin=xmin;
 gp.w.xmax=xmax;
 gp.w.ymin=ymin;
 gp.w.ymax=ymax;
 gpmake();
 return (xmax-xmin)/(ymax-ymin);
}

real gpviewport(real xmin, real xmax, real ymin, real ymax)
{
 real a=(xmax-xmin)/(ymax-ymin);
 gp.v.xmin=xmin;
 gp.v.xmax=xmax;
 gp.v.ymin=ymin;
 gp.v.ymax=ymax;
 gpmake();
 xmin=gp.w.xmin; ymin=gp.w.ymin; gpview(&xmin,&ymin);
 xmax=gp.w.xmax; ymax=gp.w.ymax; gpview(&xmax,&ymax);
 if (xmin>xmax) { real t=xmin; xmin=xmax; xmax=t; }
 if (ymin>ymax) { real t=ymin; ymin=ymax; ymax=t; }
 dvclip(xmin,xmax,ymin,ymax);
 return a;
}

void gpview(real* x, real* y)
{
 *x=(int)(gp.ax*(*x)+gp.bx);
 *y=(int)(gp.ay*(*y)+gp.by);
}

void gpunview(real* x, real* y)
{
 *x=(*x-gp.bx+0.5)/gp.ax;
 *y=(*y-gp.by+0.5)/gp.ay;
}

void gpmake(void)
{
 real Ax=(gp.d.xmax-gp.d.xmin)/gp.d.xu;			/* aspect correction */
 real Ay=(gp.d.ymax-gp.d.ymin)/gp.d.yu;
 gp.ax = (gp.v.xmax-gp.v.xmin)/(gp.w.xmax-gp.w.xmin);	/* map wc to ndc */
 gp.ay = (gp.v.ymax-gp.v.ymin)/(gp.w.ymax-gp.w.ymin);
 gp.bx =  gp.v.xmin-gp.ax*gp.w.xmin;
 gp.by =  gp.v.ymin-gp.ay*gp.w.ymin;
 gp.ax = Ax*gp.ax;					/* map ndc to dc */
 gp.ay = Ay*gp.ay;
 gp.bx = Ax*gp.bx+gp.d.xmin;
 gp.by = Ay*gp.by+gp.d.ymin;
}
