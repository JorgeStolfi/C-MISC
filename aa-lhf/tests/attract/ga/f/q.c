/*
* henon.c
* henon map (1-ax^2+y, bx)
* Luiz Henrique de Figueiredo
* 01 Dec 94
*/

#include <ga.h>

#define a	1.3
#define b	1

void aa_f(AAform x, AAform y, AAform fx, AAform fy)
{
 aa_scale(fx,y,a);
 aa_add(fx,fx,x);
 aa_scale(fy,x,b);
 return;
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
