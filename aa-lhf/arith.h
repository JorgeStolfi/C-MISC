/*
* arith.h
* interval arihmetic
* Luiz Henrique de Figueiredo
* 24 Nov 94
*/
/* Last edited on 2023-01-14 12:21:01 by stolfi */

#ifndef real
#define	real	double
#define MAXREAL MAXDOUBLE
#endif

#ifndef abs
#define	abs(_)	( ((_)<0) ? -(_) : (_) )
#endif

#define	N	45
typedef real	AAform[N];			/* full representation */

void aa_set(AAform a, AAform b);		/* a = b */
void aa_neg(AAform a, AAform b);		/* a = -b */
void aa_scale(AAform a, AAform b, real t);	/* a = tb */
void aa_trans(AAform a, AAform b, real t);	/* a = b+t */
void aa_add(AAform a, AAform b, AAform c);	/* a = b+c */
void aa_sub(AAform a, AAform b, AAform c);	/* a = b-c */
void aa_mul(AAform a, AAform b, AAform c);	/* a = b*c */
void aa_sqr(AAform a, AAform b);		/* a = b^2 */
void aa_inv(AAform a, AAform b);		/* a = 1/b */
void aa_sqrt(AAform a, AAform b);               /* a = sqrt(b) */
void aa_tanh(AAform a, AAform b);               /* a = tanh(b) */

real aa_norm(AAform a);
void aa_interval(AAform a, real amin, real amax);
void aa_range(AAform a, real* amin, real* amax);

void aa_open(void);
void aa_close(void);
void aa_trace(char* name, AAform a);
char* aa_id(void);

#ifdef DEBUG
#define	aa_debug(x)	aa_trace(#x,x)
#else
#define	aa_debug(x)
#endif
