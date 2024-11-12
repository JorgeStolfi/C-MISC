/* See {gp.h}. */
/* Luiz Henrique de Figueiredo (lhf@visgraf.impa.br) - 28 Apr 93 */
/* Last edited on 2016-09-13 15:40:01 by stolfilocal */

#include <gp.h>

static struct
{
 Box w,v,d;
 double ax,bx;
 double ay,by;
} gp=
{
 {0.0, 1.0, 0.0, 1.0, 1.0, 1.0},
 {0.0, 1.0, 0.0, 1.0, 1.0, 1.0},
 {0.0, 1.0, 0.0, 1.0, 1.0, 1.0},
  1.0, 0.0,
  1.0, 0.0,
};

double gpopen(char* name)
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

int gppoint(double x, double y)
{
 gpview(&x,&y);
 return dvpoint(x,y);
}

void gpplot(double x, double y)
{
 gpview(&x,&y);
 dvplot(x,y);
}

void gpline(double x1, double y1, double x2, double y2)
{
 gpview(&x1,&y1);
 gpview(&x2,&y2);
 dvline(x1,y1,x2,y2);
}

void gpbox(double xmin, double xmax, double ymin, double ymax)
{
 gpview(&xmin,&ymin);
 gpview(&xmax,&ymax);
 if (xmin>xmax) { double t=xmin; xmin=xmax; xmax=t; }
 if (ymin>ymax) { double t=ymin; ymin=ymax; ymax=t; }
 dvbox(xmin,xmax,ymin,ymax);
}

void gptri(double x1, double y1, double x2, double y2, double x3, double y3)
{
 gpbegin('f');
  gppoint(x1,y1);
  gppoint(x2,y2);
  gppoint(x3,y3);
 gpend();
}

void gpcircle(double x, double y, double r)
{
 double x1=x+r;
 double y1=y;
 gpview(&x,&y);
 gpview(&x1,&y1);
 if (x1>x) r=x1-x; else r=x-x1;
 dvcircle(x,y,r);
}

void gptext(double x, double y, char* s, char* mode)
{
 gpview(&x,&y);
 dvtext(x,y,s,mode);
}

char* gpevent(int wait, double* x, double* y)
{
 int ix,iy;
 char* r=dvevent(wait,&ix,&iy);
 *x=ix; *y=iy;
 gpunview(x,y);
 return r;
}

double gpwindow(double xmin, double xmax, double ymin, double ymax)
{
 gp.w.xmin=xmin;
 gp.w.xmax=xmax;
 gp.w.ymin=ymin;
 gp.w.ymax=ymax;
 gpmake();
 return (xmax-xmin)/(ymax-ymin);
}

double gpviewport(double xmin, double xmax, double ymin, double ymax)
{
 double a=(xmax-xmin)/(ymax-ymin);
 gp.v.xmin=xmin;
 gp.v.xmax=xmax;
 gp.v.ymin=ymin;
 gp.v.ymax=ymax;
 gpmake();
 xmin=gp.w.xmin; ymin=gp.w.ymin; gpview(&xmin,&ymin);
 xmax=gp.w.xmax; ymax=gp.w.ymax; gpview(&xmax,&ymax);
 if (xmin>xmax) { double t=xmin; xmin=xmax; xmax=t; }
 if (ymin>ymax) { double t=ymin; ymin=ymax; ymax=t; }
 dvclip(xmin,xmax,ymin,ymax);
 return a;
}

void gpview(double* x, double* y)
{
 *x=(int)(gp.ax*(*x)+gp.bx);
 *y=(int)(gp.ay*(*y)+gp.by);
}

void gpunview(double* x, double* y)
{
 *x=(*x-gp.bx+0.5)/gp.ax;
 *y=(*y-gp.by+0.5)/gp.ay;
}

void gpmake(void)
{
 double Ax=(gp.d.xmax-gp.d.xmin)/gp.d.xu;			/* aspect correction */
 double Ay=(gp.d.ymax-gp.d.ymin)/gp.d.yu;
 gp.ax = (gp.v.xmax-gp.v.xmin)/(gp.w.xmax-gp.w.xmin);	/* map wc to ndc */
 gp.ay = (gp.v.ymax-gp.v.ymin)/(gp.w.ymax-gp.w.ymin);
 gp.bx =  gp.v.xmin-gp.ax*gp.w.xmin;
 gp.by =  gp.v.ymin-gp.ay*gp.w.ymin;
 gp.ax = Ax*gp.ax;					/* map ndc to dc */
 gp.ay = Ay*gp.ay;
 gp.bx = Ax*gp.bx+gp.d.xmin;
 gp.by = Ay*gp.by+gp.d.ymin;
}
