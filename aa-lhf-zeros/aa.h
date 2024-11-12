/* Last edited on 2009-02-10 12:15:36 by stolfi */

/*
* aa.h
* affine arithmetic manager
* Luiz Henrique de Figueiredo
* 14 Sep 98 10:19:53
*/

#ifndef aa_H
#define aa_H

#ifndef AA_N
#define	AA_N	32 
#endif

#ifndef NORM
#define NORM	AA_N
#define NONORM	-1.0
#endif

#ifndef real
#define	real	double
#endif

typedef real	AAform[AA_N+1];			/* full representation */

void aa_clear(AAform a);
void aa_set(AAform a, AAform b);		/* a = b */
void aa_neg(AAform a, AAform b);		/* a = -b */
void aa_scale(AAform a, AAform b, real t);	/* a = tb */
void aa_trans(AAform a, AAform b, real t);	/* a = b+t */
void aa_add(AAform a, AAform b, AAform c);	/* a = b+c */
void aa_sub(AAform a, AAform b, AAform c);	/* a = b-c */
void aa_mul(AAform a, AAform b, AAform c);	/* a = b*c */
void aa_sqr(AAform a, AAform b);		/* a = b^2 */
void aa_exp(AAform a, AAform b);		/* a = exp(b) */
void aa_log(AAform a, AAform b);		/* a = log(b) */
void aa_inv(AAform a, AAform b);		/* a = 1/b */
void aa_sin(AAform a, AAform b);		/* a = sin(b) */
void aa_cos(AAform a, AAform b);		/* a = cos(b) */
void aa_atan(AAform a, AAform b);		/* a = atan(b) */
void aa_sincos(AAform a, AAform b, AAform c);	/* a = sin(c), b=cos(c) */
void aa_condense(AAform a, int from, int at);

real aa_norm(AAform a);
void aa_interval(AAform a, real amin, real amax);
void aa_range(AAform a, real* amin, real* amax);

void aa_open(void);
void aa_close(void);
void aa_stats(void);
void aa_trace(char* name, AAform a);
char* aa_id(void);
int aa_last(void); /* Index of last noise used so far */
void aa_reset(int k);

#ifdef DEBUG
#define	aa_debug(x)	aa_trace(#x,x)
#else
#define	aa_debug(x)
#endif

/* some useful macros */

#ifndef abs
#define	abs(_)	( ((_)<0) ? -(_) : (_) )
#endif

#ifndef BIGR
#define BIGR (1.0e80)
#endif

#ifndef M_PI
#define M_2PI	6.28318530717958647692528676655900576
#define M_PI	3.14159265358979323846264338327950288
#define M_PI_2	1.57079632679489661923132169163975144
#define M_PI_4	0.78539816339744830961566084581987572
#endif
#ifndef M_2PI
#define M_2PI	6.28318530717958647692528676655900576
#endif

#endif /* aa_H */
