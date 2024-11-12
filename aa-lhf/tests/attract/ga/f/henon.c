/*
* henon.c
* henon map (1-ax^2+y, bx)
* Luiz Henrique de Figueiredo
* 01 Dec 94
*/

#include <ga.h>

#define a	1.4
#define b	0.3

void aa_f(AAform x, AAform y, AAform fx, AAform fy)
{
 AAform s,t;
 aa_sqr(s,x);
 aa_scale(t,s,-a);
 aa_add(t,t,y);
 aa_trans(fx,t,1);
 aa_scale(fy,x,b);
}

void f_omega(void)
{
 enumerate(-2,2,-0.5,0.5); return;
 enumerate(-2,2,-2,2);
}

char* f_id(void)
{
 return "henon";
}
