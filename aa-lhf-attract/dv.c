/*
* dv/X11.c
* device primitives for X11
* Luiz Henrique de Figueiredo (lhf@visgraf.impa.br)
* 22 Nov 95
*/

#include <malloc.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/cursorfont.h>

#include <gp.h>

#define	vector(n,t)	( (t*) malloc((n)*sizeof(t)) )
#define	revector(p,n)	( (void*) realloc(p,(n)*sizeof(*p)) )
#define	colormapsize	256	/* was DisplayCells(display,screen) */

/* add PointerMotionMask for idle motion report */
#define	EventMask	\
  (ButtonPressMask|ButtonReleaseMask|KeyPressMask|ButtonMotionMask|ExposureMask)
#define	Idle		(-1)

typedef unsigned long	Pixel;

static Box	dv;
static Display*	display;
static int	screen;
static Window	window;
static Drawable	canvas;
static Pixmap	back=0;
static GC	gc;
static Colormap	cm;
static Cursor	waitcursor;
static Cursor	closecursor;
static XFontStruct* font=NULL;
static int	color=1;
static int	marksize=3;
static char*	marktype=".";
static Pixel*	colormap;
static int	polymode=0;
static XPoint*	point=NULL;
static int	points=0;
static int	maxpoints=0;

static void	waitinput	(Cursor cursor);
static void	waitevent	(int mask);

Box* dvopen(char* name)
{
 display=XOpenDisplay(NULL);
 if (display==NULL)
 {
  fprintf(stderr,"%s: cannot open X display [%s]\n",name,XDisplayName(NULL));
  exit(1);
 }
 screen=DefaultScreen(display);
 gc=DefaultGC(display,screen);
 cm=DefaultColormap(display,screen);
 font=XQueryFont(display,XGContextFromGC(gc));
 waitcursor=XCreateFontCursor(display,XC_icon);
 closecursor=XCreateFontCursor(display,XC_pirate);
 colormap=vector(colormapsize,Pixel);
 colormap[0]=WhitePixel(display,screen);
 colormap[1]=BlackPixel(display,screen);
 {					/* in case use color before binding */
  int c;
  for (c=2; c<colormapsize; c++)
   colormap[c]=c;			/* hope for the best */
 }
 XSetBackground(display,gc,colormap[0]);
 XSetForeground(display,gc,colormap[1]);
 XSetFillStyle(display,gc,FillSolid);
 {					/* window geometry */
  int x,y,m;
  unsigned int width,height;
  m=XParseGeometry("512x512+0+0",&x,&y,&width,&height);
  m=XParseGeometry(getenv("GP"),&x,&y,&width,&height);
  if (m&WidthValue && !((m&HeightValue))) height=width;
  if (m&XNegative) x+=DisplayWidth(display,screen)-width-2;
  if (m&YNegative) y+=DisplayHeight(display,screen)-height-2;
  window=XCreateSimpleWindow(display,RootWindow(display,screen),
	x,y,				/* origin */
	width,height,			/* size */
	1,				/* border width */
	colormap[1],			/* border color */
	colormap[0]			/* background color */
	);
  dv.xmin=0;
  dv.xmax=width-1;
  dv.ymin=height-1;
  dv.ymax=0;
  dv.xu=width;				/* should query server for aspect */
  dv.yu=height;
 }
 XStoreName(display,window,name);
 XSelectInput(display,window,EventMask);
 XMapWindow(display,window);
 waitevent(ExposureMask);		/* wait permission for drawing */
 if (DoesBackingStore(ScreenOfDisplay(display,screen)))
 {					/* lazy redraw */
  XSetWindowAttributes a;
  a.backing_store=Always;		/* even when iconified */
  a.bit_gravity=StaticGravity;		/* even when resized */
  XChangeWindowAttributes(display,window,CWBackingStore|CWBitGravity,&a);
 }
 canvas=window;
 return &dv;
}

void dvclose(int wait)
{
 if (wait) waitinput(closecursor);
 XCloseDisplay(display);
}

void dvclear(int wait)
{
 if (wait) waitinput(waitcursor);
 if (canvas==window)
  XClearWindow(display,window);
 else
 {
  int old=dvcolor(0);
  dvbox(0,dv.xu,0,dv.yu);
  dvcolor(old);
 }
}

void dvflush(void)
{
 XFlush(display);
}

void dvwait(int t)
{
 void usleep(unsigned);			/* sleep microseconds */
 dvflush();
 if (t>=0)
  usleep(1000*t);			/* sleep t miliseconds */
 else
  waitinput(waitcursor);
}

int dvpalette(int c, char* name)
{
 XColor C;
 if (c<0 || c>=colormapsize)
  return 0;
 if (XParseColor(display,cm,name,&C) && XAllocColor(display,cm,&C))
 {
  colormap[c]=C.pixel;
  if (c==0)				/* color 0 is background */
  {
   XSetWindowBackground(display,window,colormap[c]);
   XSetBackground(display,gc,colormap[c]);
   XClearWindow(display,window);
  }
  return 1;
 }
 else
  return 0;
}

#define	C(x)	round((x)*(0x0ffff))
int dvrgb(int c, real r, real g, real b)
{
 char name[]="#RRRRGGGGBBBB";
 sprintf(name,"#%04X%04X%04X",C(r),C(g),C(b));
 return dvpalette(c,name);
}
#undef C

int dvcolor(int c)
{
 if (c<0 || c>=colormapsize)
  return -colormapsize;
 else
 {
  int old=color;
  color=c;
  XSetForeground(display,gc,colormap[c]);
  return old;
 }
}

int dvfont(char* name)
{
 XFontStruct* f=XLoadQueryFont(display,name);
 if (f==NULL)
  return 0;
 else
 {
  XSetFont(display,gc,f->fid);
  font=f;
  return 1;
 }
}

void dvmark(int size, char* mark)
{
 if (size>0) marksize=size;
 if (mark!=NULL && *mark!=0) marktype=mark;
}

void dvclip(int xmin, int xmax, int ymin, int ymax)
{
 XRectangle r;
 r.x=xmin;
 r.y=ymin;
 r.width=xmax-xmin+1;
 r.height=ymax-ymin+1;
 XSetClipRectangles(display,gc,0,0,&r,1,Unsorted);
}

void dvline(int x1, int y1, int x2, int y2)
{
 XDrawLine(display,canvas,gc,x1,y1,x2,y2);
}

void dvbox(int xmin, int xmax, int ymin, int ymax)
{
 XFillRectangle(display,canvas,gc,xmin,ymin,xmax-xmin+1,ymax-ymin+1);
}

void dvtri(int x1, int y1, int x2, int y2, int x3, int y3)
{
 dvbegin('f');
  dvpoint(x1,y1);
  dvpoint(x2,y2);
  dvpoint(x3,y3);
 dvend();
}

#define	at(a,b)	(m[0]==a && m[1]==b)

void dvtext(int x, int y, char* s, char* m)
{
 int n=strlen(s);
 int w=XTextWidth(font,s,n);
 int h=font->ascent+font->descent;
 int opaque=0;
 y-=font->descent-1;
 if (m==NULL) m="sw";
 if (m[0]=='!') { opaque=1; ++m; } else opaque=0;
 if (0);
 else if (at('n', 0 )) { x-=w/2;	y+=h;	}
 else if (at('n','e')) { x-=w;		y+=h;	}
 else if (at('n','w')) { 		y+=h;	}
 else if (at('s', 0 )) { x-=w/2;		}
 else if (at('s','e')) { x-=w;			}
 else if (at('s','w')) { 			}
 else if (at('e', 0 )) { x-=w; 		y+=h/2;	}
 else if (at('w', 0 )) { 		y+=h/2;	}
 else if (at('c', 0 )) { x-=w/2;	y+=h/2;	}
 if (opaque)
  XDrawImageString(display,canvas,gc,x,y,s,n);
 else
  XDrawString(display,canvas,gc,x,y,s,n);
}

void dvcircle(int x, int y, int r)
{
 XDrawArc(display,canvas,gc,x-r,y-r,2*r,2*r,0,360*64);
}

void dvplot(int x, int y)
{
 int size=marksize;
 char* m;
 for (m=marktype; *m!=0; m++)
  switch (*m)
  {
   XSegment s[5]; XPoint p[5];
   int dx,dy;
   default:
   case '.':
    XDrawPoint(display,canvas,gc,x,y);
    break;
   case 'o':
    XDrawArc(display,canvas,gc,x-size,y-size,2*size,2*size,0,360*64);
    break;
   case 'O':
    XFillArc(display,canvas,gc,x-size,y-size,2*size,2*size,0,360*64);
    break;
   case '+':
    s[0].x1=x-size;	s[0].y1=y;	s[0].x2=x+size;	s[0].y2=y;
    s[1].x1=x;		s[1].y1=y-size;	s[1].x2=x;	s[1].y2=y+size;
    XDrawSegments(display,canvas,gc,s,2);
    break;
   case '*':
    dx=size*0.5; dy=size*0.866;
    s[0].x1=x-size;	s[0].y1=y;	s[0].x2=x+size;	s[0].y2=y;
    s[1].x1=x+dx;	s[1].y1=y+dy;	s[1].x2=x-dx;	s[1].y2=y-dy;
    s[2].x1=x+dx;	s[2].y1=y-dy;	s[2].x2=x-dx;	s[2].y2=y+dy;
    XDrawSegments(display,canvas,gc,s,3);
    break;
   case 'x':
    dx=size; dy=size;
    s[0].x1=x+dx;	s[0].y1=y+dy;	s[0].x2=x-dx;	s[0].y2=y-dy;
    s[1].x1=x+dx;	s[1].y1=y-dy;	s[1].x2=x-dx;	s[1].y2=y+dy;
    XDrawSegments(display,canvas,gc,s,2);
    break;
   case 'b':
    XDrawRectangle(display,canvas,gc,x-size,y-size,2*size,2*size);
    break;
   case 'B':
    XFillRectangle(display,canvas,gc,x-size,y-size,2*size+1,2*size+1);
    break;
   case 'd':
    p[0].x=x-size;	p[0].y=y;
    p[1].x=x;		p[1].y=y+size;
    p[2].x=x+size;	p[2].y=y;
    p[3].x=x;		p[3].y=y-size;
    p[4].x=x-size;	p[4].y=y;
    XDrawLines(display,canvas,gc,p,5,CoordModeOrigin);
    break;
   case 'D':
    p[0].x=x-size;	p[0].y=y;
    p[1].x=x;		p[1].y=y+size;
    p[2].x=x+size;	p[2].y=y;
    p[3].x=x;		p[3].y=y-size;
    p[4].x=x-size;	p[4].y=y;
    XFillPolygon(display,canvas,gc,p,5,Convex,CoordModeOrigin);
    break;
  }
}

void dvbegin(int mode)
{
 if (point==NULL) point=vector(maxpoints=16,XPoint);
 polymode=mode;
 points=0;
}

int dvpoint(int x, int y)		/* should check XMaxRequestSize */
{
 int i=points++;
 if (i>=maxpoints) point=revector(point,maxpoints*=2);
 point[i].x=x;
 point[i].y=y;
 return points;
}

void dvend(void)
{
 switch (polymode)
 {
  case 'p':				/* closed polygonal line */
   dvpoint(point[0].x,point[0].y);
  case 'l':				/* open polygonal line */
   XDrawLines(display,canvas,gc,point,points,CoordModeOrigin);
   break;
  case 'f':				/* filled polygon */
   dvpoint(point[0].x,point[0].y);
   XFillPolygon(display,canvas,gc,point,points,Complex,CoordModeOrigin);
   break;
  case 'm':				/* polymarker (only dots) */
   XDrawPoints(display,canvas,gc,point,points,CoordModeOrigin);
   break;
 }
 polymode=0;
 points=0;
}

char* dvevent(int wait, int* x, int* y)
{
 XEvent event;
 unsigned int state=0;
 if (wait)
  XWindowEvent(display,window,EventMask,&event);
 else if (!XCheckWindowEvent(display,window,EventMask,&event))
#if 0
  return NULL;
#else					/* report idle status */
 {
  Window rw,w ; int rx,ry;
  XQueryPointer(display,window,&rw,&w,&rx,&ry,x,y,&state);
  event.type=Idle;			/* fake event */
 }
#endif
 {
  static char report[]="m123+SCM";
  char* r=report;
  switch (event.type)
  {
   case Expose:
    *x=100000000+event.xexpose.x*10001+event.xexpose.width;
    *y=100000000+event.xexpose.y*10001+event.xexpose.height;
    *r++='r';
    *r++=(event.xexpose.count!=0) ? '+' : '-';	/* signal last expose */
    break;
   case ButtonPress:
   case ButtonRelease:
    state=event.xbutton.state;
    *x=event.xbutton.x;
    *y=event.xbutton.y;
    *r++='b';
    *r++='0'+event.xbutton.button;
    *r++=(event.type==ButtonPress) ? '+' : '-';
    break;
   case KeyPress:
   case KeyRelease:
    report[1]=0;
    XLookupString(&event.xkey,report+1,sizeof(report)-1,NULL,NULL);
    state=event.xkey.state;
    *x=event.xkey.x;
    *y=event.xkey.y;
    *r++='k';
    r++;				/* lazy */
    *r++=(event.type==KeyPress) ? '+' : '-';
    break;
   case MotionNotify:
#if 1
    while (XEventsQueued(display,QueuedAfterReading)>0)
    {
     XEvent ahead;
     XPeekEvent(display,&ahead);
     if (ahead.type!=MotionNotify) break;
     if (ahead.xmotion.window!=window) break;
     XWindowEvent(display,window,EventMask,&event);
    }
#endif
    state=event.xmotion.state;
#if 1
    *x=event.xmotion.x;
    *y=event.xmotion.y;
#else
{
Window rw,w ; int rx,ry;
XQueryPointer(display,window,&rw,&w,&rx,&ry,x,y,&state);
}
#endif
   case Idle:
    *r++=(event.type==MotionNotify) ? 'm' : 'i';
    if (state&Button1Mask) *r++='1';
    if (state&Button2Mask) *r++='2';
    if (state&Button3Mask) *r++='3';
    *r++='+';
    break;
  }
  if (state&ShiftMask)	 *r++='S';
  if (state&ControlMask) *r++='C';
  if (state&Mod1Mask)	 *r++='M';	/* meta is Mod1*/
  *r++=0;
  return report;
 }
}

static void waitinput(Cursor cursor)
{
 XDefineCursor(display,window,cursor);
 waitevent(ButtonReleaseMask|KeyPressMask);
 XUndefineCursor(display,window);
}

static void waitevent(int mask)
{
 XEvent event;
 XWindowEvent(display,window,mask,&event);
}

void dvdoublebuffer(int on)
{
 if (on)
 {
  back=XCreatePixmap(display,window,dv.xu,dv.yu,8);
  dvbackbuffer();
  dvclear(0);
 }
 else if (back!=0)
 {
  XFreePixmap(display,back);
  back=0;
 }
}

void dvswapbuffers(void)
{
 XCopyArea(display,back,window,gc,0,0,dv.xu,dv.yu,0,0);
#if 0
 XFlush(display);
#endif
}

void dvfrontbuffer(void)
{
 canvas=(Drawable) window;
}

void dvbackbuffer(void)
{
 canvas=(Drawable) back;
}
