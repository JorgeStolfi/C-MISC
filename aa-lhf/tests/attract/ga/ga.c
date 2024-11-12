/*
* ga.c
* approximation of global attractors
* Luiz Henrique de Figueiredo (lhf@lncc.br)
* 08 May 97
*/

#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <gp.h>
#include <ga.h>

typedef struct node Node;
struct node 
{
 real xmin,xmax,ymin,ymax;
 Node *l,*r;
};
static Node dummyNIL,*NIL=&dummyNIL;
static Node dummyOK,*OK=&dummyOK;

#define isleaf(_) (L(_)==NIL || L(_)==OK)
#define next(_) ((_)->r)
#define xmin(_) ((_)->xmin)
#define xmax(_) ((_)->xmax)
#define ymin(_) ((_)->ymin)
#define ymax(_) ((_)->ymax)
#define L(_) ((_)->l)
#define R(_) ((_)->r)
#define new(t)	malloc(sizeof(t))

int subdivide(void);
static int select(void);
void display(void);
void plot(int c, Node* t);
void show(int c, int t, real xmin, real xmax, real ymin, real ymax);
int toosmall(Node* t);
void image(Node* t, real* fxmin, real* fxmax, real* fymin, real* fymax);
Node* mknode(real xmin, real xmax, real ymin, real ymax);
void meet(Node* t, Node* s);

static real tol;
static unsigned long int Kimage, Ktoosmall, Kplot, Kmeet;
static Node* head;
static Node* top;

int main(int argc, char* argv[])
{
 if (argc>1)
  sscanf(argv[1],"%lg",&tol);
 else
  tol=-8;
 f_omega();
fprintf(stderr,"%s %s tol %s=%g ",f_id(),aa_id(),argv[1],tol);
fprintf(stderr,"visited %lu grey %lu leaves %lu meet %lu\n",
	Kimage,Ktoosmall,Kplot,Kmeet);
 return 0;
}

void enumerate(real xmin, real xmax, real ymin, real ymax)
{
 int n;
 char s[200];
 sprintf(s,"%s global attractor for %s",aa_id(),f_id());
 gpopen(s);
#if 1
 gppalette(0,"black");
 gppalette(1,"white");
 gppalette(2,"red");
 gppalette(3,"blue");
 gppalette(4,"grey");
 gppalette(5,"green");
#else
 gppalette(2,"black");
 gppalette(3,"black");
 gppalette(4,"black");
#endif
 gpwindow(xmin,xmax,ymin,ymax);
 if (tol<0)				/* negative tol means in pixels */
 {
  real x0,x1,y0,y1;
  tol=-tol;
  x0=0;   y0=0;   gpunview(&x0,&y0);
  x1=tol; y1=tol; gpunview(&x1,&y1);
  tol=min(fabs(x1-x0),fabs(y1-y0));
 }
 head=top=mknode(xmin,xmax,ymin,ymax);
 for (n=1; !toosmall(head); n++)
 {
  int k;
  fprintf(stderr,"%d",n);
  k=subdivide();
  fprintf(stderr,"\t%d",k);
  k=select();
  fprintf(stderr,"\t%d\n",k);
if(0)gpwait(-1);
 }
if(1)gpclear(0);
 display();
 gpclose(1);
}

Node* subdivide1(Node* t)
{
 real dx=xmax(t)-xmin(t);
 real dy=ymax(t)-ymin(t);
 if (dx<dy)
 {
  real ymid=(ymin(t)+ymax(t))/2.0;
  L(t)=mknode(xmin(t),xmax(t),ymin(t),ymid);
  R(t)=mknode(xmin(t),xmax(t),ymid,ymax(t));
 }
 else
 {
  real xmid=(xmin(t)+xmax(t))/2.0;
  L(t)=mknode(xmin(t),xmid,ymin(t),ymax(t));
  R(t)=mknode(xmid,xmax(t),ymin(t),ymax(t));
 }
 next(L(t))=R(t);
 return L(t);
}

int subdivide(void)
{
 int k=0;
 Node* t;
 Node pp,*p=&pp;
 for (t=head; t!=NULL; )
 {
  Node* n=next(t);
  Node* q=subdivide1(t);
  next(p)=q; p=next(q);
  t=n; ++k; ++k;
 }
 p=&pp; head=next(p);
 return k;
}

static int select(void)
{
 int k;
 Node* t;
 Node pp,*p=&pp;
 for (t=head; t!=NULL; t=next(t))
 {
  Node i;
  image(t,&i.xmin,&i.xmax,&i.ymin,&i.ymax);
if(0)show(2,'f',xmin(&i),xmax(&i),ymin(&i),ymax(&i));
  meet(top,&i);
 }
 k=0;
 next(p)=head;
 for (t=head; t!=NULL; )
 {
  Node* n=next(t);
  if (L(t)==OK)
  {
   p=t;
   ++k;
   plot(3,t);
  }
  else
  {
   next(p)=n;
   plot(0,t);
#if 0
   R(t)=NIL;
#endif
  }
  t=n;
 }
 p=&pp; head=next(p);
 return k;
}

void meet(Node* t, Node* s)
{
 if (t==NULL) return;		/* never happens? */
#if 0
 if (R(t)==NIL) return;		/* avoid deleted */
#endif
++Kmeet;
 if ((xmin(t)>xmax(s)) || (xmin(s)>xmax(t))
  || (ymin(t)>ymax(s)) || (ymin(s)>ymax(t))) return;

 if (isleaf(t)){
  L(t)=OK;
if(0)plot(5,t);
 }
 else
 {
  meet(L(t),s); meet(R(t),s);
 }
}

void display(void)
{
 Node* t;
 gpcolor(2);
 for (t=head; t!=NULL; t=next(t))
 {
#if 0
  real x=(xmin(t)+xmax(t))/2.0;
  real y=(ymin(t)+ymax(t))/2.0;
  gpplot(x,y);
#else
#if 1
  show(2,'f',xmin(t),xmax(t),ymin(t),ymax(t));
#else
  printf("%g %g %g %g B\n",xmin(t),xmax(t),ymin(t),ymax(t));
  show(1,'p',xmin(t),xmax(t),ymin(t),ymax(t));
#endif
#endif
  ++Kplot;
 }
}

int toosmall(Node* t)
{
 ++Ktoosmall;
 return (t && (xmax(t)-xmin(t))<=tol) && ((ymax(t)-ymin(t))<=tol);
 printf("%p %g %g %g %g\n",t,xmin(t),xmax(t),ymin(t),ymax(t));
}

void image(Node* t, real* fxmin, real* fxmax, real* fymin, real* fymax)
{
 AAform x,y,fx,fy;
 aa_open();
 aa_interval(x,xmin(t),xmax(t));
 aa_interval(y,ymin(t),ymax(t));
 aa_f(x,y,fx,fy);
 aa_range(fx,fxmin,fxmax);
 aa_range(fy,fymin,fymax);
#if 0
aa_trace("x",x);
aa_trace("y",y);
aa_trace("fx",fx);
aa_trace("fy",fy);
if(0)gpwait(-1);
#endif
 aa_close();
 ++Kimage;
}

void plot(int c, Node *t)
{
 show(c,'f',xmin(t),xmax(t),ymin(t),ymax(t));
 show(1,'p',xmin(t),xmax(t),ymin(t),ymax(t));
 return;
}

void show(int c, int t, real xmin, real xmax, real ymin, real ymax)
{
 gpcolor(c);
 gpbegin(t);
  gppoint(xmin,ymin);
  gppoint(xmax,ymin);
  gppoint(xmax,ymax);
  gppoint(xmin,ymax);
 gpend();
}

Node* mknode(real xmin, real xmax, real ymin, real ymax)
{
 Node* t=new(Node);
 xmin(t)=xmin; xmax(t)=xmax; ymin(t)=ymin; ymax(t)=ymax;
 L(t)=NIL;
 R(t)=NULL;
 return t;
 printf("%g %g %g %g\n",xmin(t),xmax(t),ymin(t),ymax(t));
 plot(3,t);
} 
