/*
* circle.c
* circle map (cos(x+y), sin(x+y))
* Luiz Henrique de Figueiredo
* 01 Dec 94
*/

#include <ga.h>

void aa_f(AAform x, AAform y, AAform fx, AAform fy)
{
 AAform t;
 aa_add(t,x,y);
 aa_sincos(fx,fy,t);
}

void f_omega(void)
{
 enumerate(-5,5,-5,5);
}

char* f_id(void)
{
 return "circle";
}
