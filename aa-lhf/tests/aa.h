/*
* aa.h
* affine arithmetic manager
* Luiz Henrique de Figueiredo
* 06 Sep 95
*/

#ifndef real
#define	real	double
#endif

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

#define	N	31
typedef real	AAform[N];			/* full representation */

void aa_clear(AAform a);
void aa_set(AAform a, AAform b);		/* a = b */
void aa_neg(AAform a, AAform b);		/* a = -b */
void aa_scale(AAform a, AAform b, real t);	/* a = tb */
void aa_trans(AAform a, AAform b, real t);	/* a = b+t */
void aa_add(AAform a, AAform b, AAform c);	/* a = b+c */
void aa_sub(AAform a, AAform b, AAform c);	/* a = b-c */
void aa_mul(AAform a, AAform b, AAform c);	/* a = b*c */
void aa_sqr(AAform a, AAform b);		/* a = b^2 */
void aa_inv(AAform a, AAform b);		/* a = 1/b */
void aa_sin(AAform a, AAform b);		/* a = sin(b) */
void aa_cos(AAform a, AAform b);		/* a = cos(b) */
void aa_atan(AAform a, AAform b);		/* a = atan(b) */
void aa_sincos(AAform a, AAform b, AAform c);	/* a = sin(c), b=cos(c) */

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



int main()
{
  AAform S, C, x;

  aa_clear(S);aa_clear(C);aa_clear(x);

  aa_interval(x, -21.96006130973899 -4.3466462300870035,
	         -21.96006130973899 +4.3466462300870035);

  aa_sincos(S, C, x);
  aa_trace("S", S);
  aa_trace("C", C);

  printf("x=[%f, %f]   Cos(x)=[%f, %f]\n",
	 x[0]-aa_norm(x), x[0]+aa_norm(x),
	 C[0]-aa_norm(C), C[0]+aa_norm(C));


  return 0;
}
 
