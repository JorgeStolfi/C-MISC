/*
* ia.c
* interval arithmetic manager
* Luiz Henrique de Figueiredo
* 24 Nov 94
*/
/* Last edited on 2023-02-22 12:34:09 by stolfi */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <arith.h>

/* PROTOTYPES */

void error (char *msg);
void assert(int test, char *msg);

/* IMPLEMENTATIONS */

#define	MIN	0
#define	MAX	1

#define BIGR (1.0e80)

void error (char *msg)
  { fprintf (stderr, "*** %s\n", msg);
    exit(1);
  }

void assert(int test, char *msg)
  {
    if (! test)
      { fprintf(stderr, "\n*** assertion is false: %s ***\n", msg);
        exit(1);
      }
  }

real aa_norm(AAform a)
{
 return fabs(a[MAX]-a[MIN]);
}

void aa_set(AAform a, AAform b)			/* a = b */
{
 a[MIN]=b[MIN];
 a[MAX]=b[MAX];
}

void aa_neg(AAform a, AAform b)			/* a = -b */
{
 a[MIN]=-b[MAX];
 a[MAX]=-b[MIN];
}

void aa_scale(AAform a, AAform b, real t)       /* a = tb */
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

void aa_trans(AAform a, AAform b, real t)	/* a = b+t */
{
 a[MIN]=t+b[MIN];
 a[MAX]=t+b[MAX];
}

void aa_add(AAform a, AAform b, AAform c)	/* a = b+c */
{
 a[MIN]=b[MIN]+c[MIN];
 a[MAX]=b[MAX]+c[MAX];
}

void aa_sub(AAform a, AAform b, AAform c)	/* a = b-c */
{
 a[MIN]=b[MIN]-c[MAX];
 a[MAX]=b[MAX]-c[MIN];
}

#if 0
void aa_mul(AAform a, AAform b, AAform c)	/* a = b*c */
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

void aa_sqr(AAform a, AAform b)			/* a = b^2 */
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

void aa_inv(AAform a, AAform b)	/* a=1/b */
{
  real bmin = b[MIN];
  real bmax = b[MAX];
  if ((bmin <= 1.0e-20) && (bmax >= -1.0e-20))
    { a[MIN] = -BIGR; a[MAX] = BIGR; }
  else if (bmax < 0) 
    { a[MAX] = 1.0/bmax;
      a[MIN] = 1.0/bmin;
    }
  else
    { a[MIN] = 1.0/bmax;
      a[MAX] = 1.0/bmin;
    }
}

void aa_sqrt(AAform a, AAform b)	/* a=sqrt(b) */
{
  assert(b[MAX] >= 0, "aa_sqrt: negative argumnt range");
  a[MIN] = (b[MIN] <= 0.0 ? 0.0 : sqrt(b[MIN]));
  a[MAX] = (b[MAX] <= 0.0 ? 0.0 : sqrt(b[MAX]));
}

void aa_tanh(AAform a, AAform b)	/* a=tanh(b) */
{
  a[MIN] = tanh(b[MIN]);
  a[MAX] = tanh(b[MAX]);
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
 fprintf(stderr, "aa: %s= %g .. %g\n",name,amin,amax);
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

