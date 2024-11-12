/*
* qei.h
* enumeration of implicit curves using quadtrees and affine arithmetic
* Luiz Henrique de Figueiredo
* 01 Dec 94
*/

#include <aa.h>

void enumerate(real xmin, real xmax, real ymin, real ymax);

void aa_f(AAform x, AAform y, AAform fx, AAform fy);
void f_omega(void);
char* f_id(void);
