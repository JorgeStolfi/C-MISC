/*
* aa.c
* affine arithmetic manager
* Luiz Henrique de Figueiredo & Jorge Stolfi
* 04 Jun 97
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <aa.h>
#define N	AA_N

#if 0
#define	forallsymbols(i)	for (i=N-1; i>=0; i--)
#define	forallnoises(i)		for (i=N-1; i; i--)
#elif 1
#define	forallsymbols(i)	for (i=0; i<N; i++)
#define	forallnoises(i)		for (i=1; i<N; i++)
#else
#define	forallsymbols(i)	for (i=last; i>=0; i--)
#define	forallnoises(i)		for (i=last; i; i--)
#endif

static int last=0; 

static int aa_newsymbol(void)
{
 if (++last>=N)
 {
  fprintf(stderr,"aa: cannot handle more than %d symbols\n",N);
  exit(1);
 }
 return last;
}

real aa_norm(AAform a)				/* \sum_1^n a_i*/
{
 int i;
 real na=0;
 forallnoises(i)
  na+=fabs(a[i]);
 return na;
}

void aa_set(AAform a, AAform b)			/* a = b */
{
 memcpy(a,b,sizeof(AAform));
}

void aa_neg(AAform a, AAform b)			/* a = -b */
{
 int i;
 forallsymbols(i)
  a[i]=-b[i];
}

void aa_scale(AAform a, AAform b, real t)	/* a = tb */
{
 int i;
 forallsymbols(i)
  a[i]=t*b[i];
}

void aa_trans(AAform a, AAform b, real t)	/* a = b+t */
{
 aa_set(a,b);
 a[0]+=t;
}

void aa_add(AAform a, AAform b, AAform c)	/* a = b+c */
{
 int i;
 forallsymbols(i)
  a[i]=b[i]+c[i];
}

void aa_sub(AAform a, AAform b, AAform c)	/* a = b-c */
{
 int i;
 forallsymbols(i)
  a[i]=b[i]-c[i];
}

void aa_mul(AAform a, AAform b, AAform c)	/* a = b*c */
{
 int i;
 real nb=aa_norm(b);
 real nc=aa_norm(c);
 forallnoises(i)
  a[i]=b[0]*c[i]+b[i]*c[0];
 i=aa_newsymbol();
 a[i]=nb*nc;
 a[0]=b[0]*c[0];
}

void aa_sqr(AAform a, AAform b)	/* a = b^2 */
{
 int i;
 real nb=aa_norm(b); nb*=nb; nb/=2;
 forallnoises(i)
  a[i]=2*b[0]*b[i];
 i=aa_newsymbol();
 a[i]=nb;
 a[0]=b[0]*b[0]+nb;
}

void aa_exp(AAform a, AAform b)	/* a = exp(b) */
{
 int i;
 real nb=aa_norm(b);
 real eb=exp(b[0]+nb);
 real ec=(nb==0) ? eb : (eb-exp(b[0]-nb))/(2.0*nb);
 real c=log(ec);
 real E=(eb-ec*(b[0]+nb-c+1))/2.0;
 forallnoises(i)
  a[i]=ec*b[i];
 i=aa_newsymbol();
 a[i]=E;
 a[0]=ec*(b[0]-c+1)+E;
}

void aa_log(AAform a, AAform b)	/* a = log(b) */
{
 int i;
 real nb=aa_norm(b);
 real bb=b[0]+nb;
 real fb=log(bb);
 real ci=(nb==0) ? (1.0/bb) : (fb-log(b[0]-nb))/(2.0*nb);
 real bc=bb*ci;
 real E=(bc-log(bc)-1)/2.0;
 forallnoises(i)
  a[i]=ci*b[i];
 i=aa_newsymbol();
 a[i]=E;
 a[0]=ci*b[0]-log(ci)-1-E;
}

static void cheb_inv(real xlo, real xhi, real *alpha, real *gamma, real *delta)
{
 if ((xlo <= 0.0) && (xhi >= 0.00))	/* Interval straddles zero */
 {
  *alpha = 0; *gamma = 0; *delta = BIGR; 
 }
 else if (xhi == xlo) 
 { *alpha = 0; *gamma = 1.0/xlo; *delta = 0; } 
 else if (xhi < 0.0) 
 {
  cheb_inv(-xhi, -xlo, alpha, gamma, delta); *gamma = - *gamma;
 }
 else					/* 0.0 < xlo < xhi */
 {
   real dhi, dlo;
   real slo = 1.0/xlo;
   real shi = 1.0/xhi;
   if (((xhi - xlo) < 1.0) && ((shi - slo) >= (xhi - xlo)*BIGR))
   { *alpha = BIGR; }
   else
   { *alpha = (shi - slo)/(xhi - xlo); }
   
   /* Compute max and min of $1/x - \alpha x$ in $xr$: */
   if ( -(*alpha)*xhi*xhi <= 1.0 )
     { /* Difference is monotonically decreasing in $[xlo..xhi]$. */
       dlo = shi - (*alpha)*xhi;
       dhi = slo - (*alpha)*xlo;
     }
   else if ( -(*alpha)*xlo*xlo >= 1.0 )
     { /* Difference is monotonically increasing in $[xlo..xhi]$. */
       dhi = shi - (*alpha)*xhi;
       dlo = slo - (*alpha)*xlo;
     }
   else
     { /* Difference is concave in $[xlo..xhi]$. */
       dlo = 2.0 * sqrt(-(*alpha));
       { real d1 = slo - (*alpha)*xlo;
	 real d2 = shi - (*alpha)*xhi;
	 dhi = (d1 > d2 ? d1 : d2);
       }
     }
   *gamma = (dhi + dlo)/2.0;
   { real r1 = dhi - (*gamma);
     real r2 = (*gamma) - dlo;
     *delta = (r1 > r2 ? r1 : r2);
   }
 }
}

void aa_inv(AAform a, AAform b)		/* a = 1/b */
{
 int i;
 real nb = aa_norm(b);
 real alpha, gamma, delta;
 cheb_inv(b[0]-nb, b[0]+nb, &alpha, &gamma, &delta);
 a[0]=alpha*b[0] + gamma;
 forallnoises(i)
  a[i]=alpha*b[i];
 i=aa_newsymbol();
 a[i]=delta;
}

void aa_sincos(AAform s, AAform c, AAform t) /* c = cos(t), s = sin(t) */
{
 int i,j,k;
 real ct0 = cos(t[0]);
 real st0 = sin(t[0]);
 real r = aa_norm(t);

 /* Compute a good joint approximation to cos(t-t0), sin(t-t0) */
 if (r == 0)
 {
  c[0] = 1; s[0] = 0;
  forallnoises(i) { c[i] = 0; s[i] = 0; }
 }
 else
 {
  /* sine: */
  s[0] = 0;
  k=aa_newsymbol();
  if (r >= 3*M_PI_2)
  {
   forallnoises(i) { s[i] = 0; }
   s[k] = 1;
  }
  else
  {
   real sr = sin(r);
   real slap = (1+sr)/(1+r); /* an upper bound to the best answer */
   forallnoises(i) { s[i] = slap*t[i]; }
   s[k] = r*slap - sr;
  }

  /* cosine: */
  forallnoises(i) { c[i] = 0; }
  j=aa_newsymbol();
  if (r >= M_PI)
  {
   c[0] = 0; c[j] = 1;
  }
  else
  {
   real cr = cos(r);
   c[0] = (1+cr)/2;
   c[j] = (1-cr)/2;
  }
 }

 forallsymbols(i) /* including i=0 */
 {
  real ci = ct0*c[i] - st0*s[i];
  real si = st0*c[i] + ct0*s[i];
  c[i] = ci; s[i] = si;
 }
}

void aa_sin(AAform a, AAform b)			/* a = sin(b) */
{
 AAform c;
 aa_sincos(a,c,b);
}

void aa_cos(AAform a, AAform b)			/* a = cos(b) */
{
 AAform s;
 aa_sincos(s,a,b);
}

void aa_clear(AAform a)
{
 memset(a,0,sizeof(AAform));
}

void aa_interval(AAform a, real amin, real amax)
{
 int i;
 memset(a,0,sizeof(AAform));
 i=aa_newsymbol();
 a[0]=(amax+amin)/2;
 a[i]=(amax-amin)/2;
}

void aa_range(AAform a, real* amin, real* amax)
{
 real na=aa_norm(a);
 *amin=a[0]-na;
 *amax=a[0]+na;
}

void aa_open(void)
{
 last=0; 
}

void aa_close(void)
{
 last=0; 
}

void aa_stats(void)
{
 fprintf(stderr,"aa: %d symbols used\n",last);
}

char* aa_id(void)
{
 return "aa";
}

#undef forallnoises
#define	forallnoises(i)		for (i=1; i<=last; i++)

void aa_trace(char* name, AAform a)
{
 int i;
 real amin,amax;
 aa_range(a,&amin,&amax);
 fprintf(stderr,"aa: %s= %g .. %g\t",name,amin,amax);
 fprintf(stderr,"%g",a[0]);
 forallnoises(i)
  if (a[i]!=0) fprintf(stderr," %+g@%d",a[i],i);
 fprintf(stderr,"\n");
}
