/*
* aa.c
* affine arithmetic manager
* Luiz Henrique de Figueiredo
* 24 Nov 94
*/
/* Last edited on 2023-02-22 12:34:43 by stolfi */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <arith.h>
#include <string.h>

/* PROTOTYPES */

/* void bcopy(char *b1, char *b2, int length); */
/* void bzero(char *b, int length); */

void error (char *msg);
void assert(int test, char *msg);
static int aa_newsymbol(void);
void cheb_inv(real xlo, real xhi, real *alpha, real *gamma, real *delta);
void cheb_sqrt(real xlo, real xhi, real *alpha, real *gamma, real *delta);
void cheb_tanh(real xlo, real xhi, real *alpha, real *gamma, real *delta);

/* IMPLEMENTATIONS */

#define	forallsymbols(i)	for (i=0; i<N; i++)
#define	forallnoises(i)		for (i=1; i<N; i++)

#define BIGR (1.0e80)

static int last=0; 

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
 bcopy((char *)b, (char*)a, sizeof(AAform));
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

void aa_inv(AAform a, AAform b)		/* a = 1/b */
{
 int i;
 real rb = aa_norm(b)/2;
 real alpha, gamma, delta;
 cheb_inv(b[0]-rb, b[0]+rb, &alpha, &gamma, &delta);
 a[0] = alpha*b[0] + gamma;
 forallnoises(i) a[i] = alpha*b[i];
 i = aa_newsymbol();
 a[i] = delta;
}

void cheb_inv(real xlo, real xhi, real *alpha, real *gamma, real *delta)
{ 
  if ((xlo <= 0.0) && (xhi >= 0.0))
    { /* Interval straddles 0 */
      *alpha = 0; *gamma = 0; *delta = BIGR; 
    }
  else if (xhi == xlo) 
    { *alpha = 0; *gamma = 1.0/xlo; *delta = 0; } 
  else if (xhi < 0.0) 
    { cheb_inv(-xhi, -xlo, alpha, gamma, delta); *gamma = - *gamma; }
  else /* 0 < xlo < xhi */
    { real dhi, dlo;
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

      /* Sanity check: */
      { real ehi = shi - (*alpha) * xhi - *gamma;
        real elo = slo - (*alpha) * xlo - *gamma;
        assert (fabs(ehi) <= 1.0001 * (*delta), "cheb_inv bug");
        assert (fabs(elo) <= 1.0001 * (*delta), "cheb_inv bug");
      }
    }
}

void aa_sqrt(AAform a, AAform b)	/* a=sqrt(b) */
{
 int i;
 real rb = aa_norm(b)/2;
 real alpha, gamma, delta;
 cheb_sqrt(b[0]-rb, b[0]+rb, &alpha, &gamma, &delta);
 a[0] = alpha*b[0] + gamma;
 forallnoises(i) a[i] = alpha*b[i];
 i = aa_newsymbol();
 a[i] = delta;
}

void cheb_sqrt(real xlo, real xhi, real *alpha, real *gamma, real *delta)
{
  if (xhi == xlo) 
    { *alpha = 0; *gamma = sqrt(xlo); *delta = 0; } 
  else if (xhi <= 0)
    { error("cheb-sqrt: negative argument range");
      *alpha = 0; *gamma = 0; *delta = 0;
    }
  else /* xhi >= 0 */
    { real slo, shi, smd, xmd;
      if (xlo < 0) xlo = 0.0;
      slo = sqrt(xlo);
      shi = sqrt(xhi);
      *alpha = (shi - slo)/(xhi - xlo); 
      smd = 0.5*(shi + slo);
      xmd = smd*smd;
      /* Error extrema are xlo(-), xmd(+), xhi(-) */
      /* Approximation is tight (Chebyshev) */
      *gamma = (smd + slo - (*alpha) * (xmd + xlo))/2;
      *delta = (smd - slo - (*alpha) * (xmd - xlo))/2;

      /* Sanity check: */
      { real ehi = shi - (*alpha) * xhi - *gamma;
        real elo = slo - (*alpha) * xlo - *gamma;
        assert (fabs(ehi) <= 1.0001 * (*delta), "cheb_sqrt bug");
        assert (fabs(elo) <= 1.0001 * (*delta), "cheb_sqrt bug");
      }
    }
}

void aa_tanh(AAform a, AAform b)	/* a=tanh(b) */
{ 
 int i;
 real rb = aa_norm(b)/2;
 real alpha, gamma, delta;
 cheb_tanh(b[0]-rb, b[0]+rb, &alpha, &gamma, &delta);
 a[0] = alpha*b[0] + gamma;
 forallnoises(i) a[i] = alpha*b[i];
 i = aa_newsymbol();
 a[i] = delta;
}

void cheb_tanh(real xlo, real xhi, real *alpha, real *gamma, real *delta)
{
  if (xhi == xlo) 
    { *alpha = 0; *gamma = tanh(xlo); *delta = 0; } 
  else if ((xhi <= 0) || (xlo < -xhi))
    { cheb_tanh(-xhi, -xlo, alpha, gamma, delta); *gamma = -(*gamma); }
  else /* xhi > 0 & xlo >= -xhi */
    { real tlo = tanh(xlo);
      real thi = tanh(xhi);
      real tmd, xmd;
      *alpha = (thi - tlo)/(xhi - xlo); 
      tmd = sqrt(1.0 - *alpha);
      xmd = atanh(tmd);
      if ((1 - tlo*tlo) > (*alpha)) 
        { /* Error extrema are xlo(-), xmd(+), xhi(-) */
          /* Approximation is tight (Chebyshev) */
          *gamma = (tmd + tlo - (*alpha) * (xmd + xlo))/2;
          *delta = (tmd - tlo - (*alpha) * (xmd - xlo))/2;
        }
      else
        { /* Error extrema are -xmd(-), xmd(+) */
          /* This is not the Chebyshev app, but should be close enough: */
          *gamma = 0.0;
          *delta = tmd - (*alpha) * xmd;
          /* Sanity check: */
          { real ehi = thi - (*alpha) * xhi;
            real elo = tlo - (*alpha) * xlo;
            assert (fabs(ehi) <= 1.0001 * (*delta), "cheb_tanh bug");
            assert (fabs(elo) <= 1.0001 * (*delta), "cheb_tanh bug");
          }
        }
    }
}

void aa_interval(AAform a, real amin, real amax)
{
 int i;
 bzero((char*)a, sizeof(AAform));
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

void aa_trace(char* name, AAform a)
{
 int i;
 real amin,amax;
 aa_range(a,&amin,&amax);
 fprintf(stderr,"aa: %s= %g .. %g\t",name,amin,amax);
 fprintf(stderr,"%g",a[0]);
 forallnoises(i)
  if (a[i]!=0) fprintf(stderr,"%+g@%d",a[i],i);
 fprintf(stderr,"\n");
}

char* aa_id(void)
{
 return "aa";
}
