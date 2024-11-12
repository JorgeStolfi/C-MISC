/*
* ia.c
* interval arithmetic manager
* Luiz Henrique de Figueiredo & Jorge Stolfi
* 19 Oct 95
*/

#include <stdio.h>
#include <math.h>
#include <aa.h>

#define	MIN	0
#define	MAX	1

real aa_norm(AAform a)
{
 return fabs(a[MAX]-a[MIN]);
}

void aa_set(AAform a, AAform b)			/* a=b */
{
 a[MIN]=b[MIN];
 a[MAX]=b[MAX];
}

void aa_neg(AAform a, AAform b)			/* a=-b */
{
 a[MIN]=-b[MAX];
 a[MAX]=-b[MIN];
}

void aa_scale(AAform a, AAform b, real t)	/* a=tb */
{
 real bmin=b[MIN];
 real bmax=b[MAX];
 if (t>0)
 {
  a[MIN]=t*bmin;
  a[MAX]=t*bmax;
 }
 else
 {
  a[MIN]=t*bmax;
  a[MAX]=t*bmin;
 }
}

void aa_trans(AAform a, AAform b, real t)	/* a=b+t */
{
 a[MIN]=t+b[MIN];
 a[MAX]=t+b[MAX];
}

void aa_add(AAform a, AAform b, AAform c)	/* a=b+c */
{
 a[MIN]=b[MIN]+c[MIN];
 a[MAX]=b[MAX]+c[MAX];
}

void aa_sub(AAform a, AAform b, AAform c)	/* a=b-c */
{
 real amin=b[MIN]-c[MAX];
 real amax=b[MAX]-c[MIN];
 a[MIN]=amin;
 a[MAX]=amax;
}

#if 0
void aa_mul(AAform a, AAform b, AAform c)	/* a=b*c */
{
 real amin,amax,t;
 t=b[MIN]*c[MIN];             amin=t;                  amax=t;
 t=b[MIN]*c[MAX]; if (t<amin) amin=t; else if (t>amax) amax=t;
 t=b[MAX]*c[MIN]; if (t<amin) amin=t; else if (t>amax) amax=t;
 t=b[MAX]*c[MAX]; if (t<amin) amin=t; else if (t>amax) amax=t;
 a[MIN]=amin;
 a[MAX]=amax;
}
#endif

void aa_sqr(AAform a, AAform b)			/* a=b^2 */
{
 real bmin=b[MIN];
 real bmax=b[MAX];
 if (bmin>=0)
 {
  a[MIN]=bmin*bmin;
  a[MAX]=bmax*bmax;
 }
 else if (bmax<=0)
 {
  a[MIN]=bmax*bmax;
  a[MAX]=bmin*bmin;
 }
 else
 {
  a[MIN]=0.0;
  if (bmax<-bmin)
   a[MAX]=bmin*bmin;
  else
   a[MAX]=bmax*bmax;
 }
}

void aa_exp(AAform a, AAform b)			/* a=exp(b) */
{
 a[MIN]=exp(b[MIN]);
 a[MAX]=exp(b[MAX]);
}

void aa_log(AAform a, AAform b)			/* a=log(b) */
{
 a[MIN]=log(b[MIN]);
 a[MAX]=log(b[MAX]);
}

void aa_clear(AAform a)
{
 a[MIN]=0.0;
 a[MAX]=0.0;
}

void aa_interval(AAform a, real amin, real amax)
{
 a[MIN]=amin;
 a[MAX]=amax;
}

void aa_range(AAform a, real* amin, real* amax)
{
 *amin=a[MIN];
 *amax=a[MAX];
}

void aa_open(void)
{
}

void aa_close(void)
{
}

void aa_trace(char* name, AAform a)
{
 real amin,amax;
 aa_range(a,&amin,&amax);
 fprintf(stderr,"aa: %s= %g .. %g\n",name,amin,amax);
}

char* aa_id(void)
{
 return "ia";
}

void aa_mul(AAform a, AAform b, AAform c)	/* a=b*c */
{
 real amin,amax;
 {
  if (b[MIN]>=0.0)
  {
   if (c[MIN]>=0.0)
   {
    amin=b[MIN]*c[MIN];
    amax=b[MAX]*c[MAX];
   }
   else if (c[MAX]<=0.0)
   {
    amin=b[MAX]*c[MIN];
    amax=b[MIN]*c[MAX];
   }
   else
   {
    amin=b[MAX]*c[MIN];
    amax=b[MAX]*c[MAX];
   }
  }
  else if (b[MAX]<=0.0)
  {
   if (c[MIN]>=0.0)
   {
    amin=b[MIN]*c[MAX];
    amax=b[MAX]*c[MIN];
   }
   else if (c[MAX]<=0.0)
   {
    amin=b[MAX]*c[MAX];
    amax=b[MIN]*c[MIN];
   }
   else
   {
    amin=b[MIN]*c[MAX];
    amax=b[MIN]*c[MIN];
   }
  }
  else
  {
   if (c[MIN]>=0.0)
   {
    amin=b[MIN]*c[MAX];
    amax=b[MAX]*c[MAX];
   }
   else if (c[MAX]<=0.0)
   {
    amin=b[MAX]*c[MIN];
    amax=b[MIN]*c[MIN];
   }
   else
   {
    real t;
    amin=b[MIN]*c[MAX];
    t=b[MAX]*c[MIN];
    if (t<amin) amin=t;
    amax= b[MIN]*c[MIN];
    t=b[MAX]*c[MAX];
    if (t>amax) amax=t;
   }
  }
 }
 a[MIN]=amin;
 a[MAX]=amax;
}

void aa_inv(AAform a, AAform b)			/* a=1/b */
{
 real bmin=b[MIN];
 real bmax=b[MAX];
 if ((bmin <= 1.0e-20) && (bmax >= -1.0e-20))
 { 
  a[MIN]=-BIGR; 
  a[MAX]=BIGR; 
 }
 else if (bmax<0)
 { 
  a[MIN]=1.0/bmin;
  a[MAX]=1.0/bmax;
 }
 else
 { 
  a[MIN]=1.0/bmax;
  a[MAX]=1.0/bmin;
 }
}

void aa_sin(AAform a, AAform b)			/* a = sin(b) */
{
 real amin,amax,bmin,bmax,t;

 if ((b[MAX]-b[MIN])>=M_2PI) { a[MIN]=-1.0; a[MAX]=1.0; return; }

 bmax=b[MAX];
 bmin=b[MIN];	t=(int)(bmin/M_2PI);	if (bmin<0) --t;	t*=M_2PI;
 bmin-=t;	bmin/=M_PI_2;
 bmax-=t;	bmax/=M_PI_2;

#if 0
printf("original: %g .. %g\tconverted: %g .. %g\n",b[MIN],b[MAX],bmin,bmax);
#endif

 if ((bmin<=1 && 1<=bmax) || (bmin<=5 && 5<=bmax)) amax=+1; else amax=2;
 if (bmin<=3 && 3<=bmax) amin=-1; else amin=2;

 if (amin>1)
 {
  bmin=sin(b[MIN]); bmax=sin(b[MAX]);
  if (bmin<bmax) amin=bmin; else amin=bmax;
 }
 else bmin=2;
 if (amax>1)
 {
  if (bmin>1) { bmin=sin(b[MIN]); bmax=sin(b[MAX]); }
  if (bmin>bmax) amax=bmin; else amax=bmax;
 }
 a[MIN]=amin;
 a[MAX]=amax;

#if 0
printf("\t%g .. %g\n",amin,amax);
#endif
}

void aa_cos(AAform a, AAform b)			/* a = cos(b) */
{
 AAform t;
 aa_trans(t,b,M_PI_2);
 aa_sin(a,t);
}

void aa_sincos(AAform a, AAform b, AAform c)	/* a = sin(c), b=cos(c) */
{
 aa_sin(a,c);
 aa_cos(b,c);
}
