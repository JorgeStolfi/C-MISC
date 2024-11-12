/*
* henon.c
* scaled henon map (1-ax^2+y/5, 5*bx)
* Luiz Henrique de Figueiredo
* 01 Dec 94
*/

#include <ga.h>

#define a	1.2
#define b	0.2

void aa_f(AAform x, AAform y, AAform fx, AAform fy)
{
 AAform s,t;
 aa_sqr(s,x);
 aa_scale(t,s,-a);
 aa_scale(s,y,0.2);
 aa_add(t,t,s);
 aa_trans(fx,t,1);
 aa_scale(fy,x,5.0*b);
}

void f_omega(void)
{
 enumerate(-2,2,-2,2);
}

char* f_id(void)
{
 return "henons";
}
