/* See {stmaps.h}. */
/* Luiz Henrique de Figueiredo & Jorge Stolfi -- 24 Jan 2008 15:00:00 */
/* Last edited on 2008-04-24 15:31:49 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#include <stmaps.h>

#define	a_henon	(1.4)
#define	b_henon	(0.3)

void f_henon(real* x, real* y, Opcode op)
{
  real rX,rY,hX,hY,tX,tY;
 switch(op)
 {
  case LOCORNER: *x = -1.5; *y = -1.5; return;
  case HICORNER: *x = +1.5; *y = +1.5; return;
  case DIRMAP:
   rX=*x; rY=*y;
   /* Do Y shear and X squeeze: */
   hY = rY + 1 - a_henon*rX*rX;
   hX = b_henon*rX;
   /* Swap: */
   tX = hY; tY = hX;
   *x = tX; *y = tY;
   return;
  case INVMAP:
   tX = *x; tY = *y;
   /* Unswap: */
   hX = tY; hY = tX;
   /* Undo X squeeze and Y shear: */
   rX = hX/b_henon;
   rY = hY - 1.0 + a_henon*rX*rX;
   *x = rX; *y = rY;
   return;
 }
}

#define	a_haoui	(0.75)
#define	b_haoui	(0.50)
#define	c_haoui	(0.50)

void f_haoui(real* x, real* y, Opcode op)
{
 real rX, rY, sX, sY, tX, tY;
 switch(op)
 {
  case LOCORNER: *x = -2.5; *y = -2.5; return;
  case HICORNER: *x = +2.5; *y = +2.5; return;
  case DIRMAP:
   rX = *x; rY = *y;
   /* Do Y shear and X squeeze: */
   sX = b_haoui*rX;
   sY = rY + a_haoui*(1 - rX*rX) + c_haoui;
   /* Swap: */
   tX = sY; tY = sX;
   *x = tX; *y = tY;
   return;
  case INVMAP:
   tX = *x; tY = *y;
   /* Un-swap: */
   sX = tY; sY = tX;
   /* Undo X squeeze and Y shear: */
   rX = sX/b_haoui;
   rY = sY - a_haoui*(1 - rX*rX) - c_haoui;
   *x = rX; *y = rY;
   return;
 }
}

#define a_sqout (0.5)
/* Relative amount of Y-pushing (max 0.5). */ 

void f_sqout(real* x, real* y, Opcode op)
{
 real rX, rY, sX, sY, tX, tY;
 real t, Q;
 switch(op)
 {
  case LOCORNER: *x = -1.001; *y = -1.001; return;
  case HICORNER: *x = +1.001; *y = +1.001; return;
  case DIRMAP:
   rX = *x; rY = *y;
   reduce_to_signed_unit_square(&rX,&rY);
   /* Do Y shear: */
   Q = (1-rX*rX)*(1-rY*rY);
   sX = rX;
   sY = rY + a_sqout*Q;
   /* Rotate Y: */
   tX = -sY; tY = +sX;
   *x = tX; *y = tY;
   return;
  case INVMAP:
   tX = *x; tY = *y;
   reduce_to_signed_unit_square(&tX,&tY);
   /* Un-rotate: */
   sX = +tY; sY = -tX;
   /* Undo Y shear: */
   rX = sX;
   t=a_sqout*(1 - rX*rX);
   solve_quadratic(t,-1,sY-t,  &rY,NULL);
   *x = rX; *y = rY;
   return;
 }
}

#define C_phi (-0.73736887807831990156)
#define S_phi (-0.67549029426152364229)
/* Cosine and sine of {1/phi} of a turn. */ 
/* echo 'f=(sqrt(5)-1)/2; p=4*a(1); f; p; c(f*2*p); s(f*2*p)' | bc -lq */ 

#define a_cirout (0.5)
/* Relative amount of Y-pushing (max 0.5). */ 

void f_cirout(real* x, real* y, Opcode op)
{
 real rX, rY, sX, sY, tX, tY;
 switch(op)
 {
  case LOCORNER: *x = -1.001; *y = -1.001; return;
  case HICORNER: *x = +1.001; *y = +1.001; return;
  case DIRMAP:
   rX = *x; rY = *y;
   reduce_to_unit_disk(&rX,&rY);
   /* Push up: */
   sX = rX;
   sY = rY + a_cirout*(1 - rX*rX - rY*rY);
   /* Rotate: */
   tX = + C_phi*sX + S_phi*sY;
   tY = - S_phi*sX + C_phi*sY;
   *x = tX; *y = tY;
   return;
  case INVMAP:
   tX = *x; tY = *y;
   reduce_to_unit_disk(&tX,&tY);
   /* Un-rotate: */
   sX = + C_phi*tX - S_phi*tY;
   sY = + S_phi*tX + C_phi*tY;
   /* Pull down: */
   rX = sX;
   solve_quadratic(a_cirout,-1.0,sY - a_cirout*(1-sX*sX),  &rY,NULL);
   *x = rX; *y = rY;
   return;
 }
}

#define a_cirin (0.5)
/* Relative amount of contraction towards the X axis. */ 

static int q_cirin = 1; /* TRUE to print warning for inverse. */

void f_cirin(real* x, real* y, Opcode op)
{
 real rX, rY, sX, sY, tX, tY;
 switch(op)
 {
  case LOCORNER: *x = -1.001; *y = -1.001; return;
  case HICORNER: *x = +1.001; *y = +1.001; return;
  case DIRMAP:
   rX = *x; rY = *y;
   reduce_to_unit_disk(&rX,&rY);
   /* Contract: */
   sX = rX;
   sY = rY - a_cirin*rY*(1 - rX*rX - rY*rY);
   /* Rotate: */
   tX = + C_phi*sX + S_phi*sY;
   tY = - S_phi*sX + C_phi*sY;
   *x = tX; *y = tY;
   return;
  case INVMAP:
   tX = *x; tY = *y;
   reduce_to_unit_disk(&tX,&tY);
   /* Un-rotate: */
   sX = + C_phi*tX - S_phi*tY;
   sY = + S_phi*tX + C_phi*tY;
   /* Expand: */
   /* Requires solving a cubic equation... */
   if(q_cirin) 
   {fprintf(stderr, "f_cirin: inverse not defined"); q_cirin=0;}
   return;
 }
}

#define a_contract (0.75)

void f_contract(real* x, real* y, Opcode op)
{
 switch(op)
 {
  case LOCORNER: *x = -1.001; *y = -1.001; return;
  case HICORNER: *x = +1.001; *y = +1.001; return;
  case DIRMAP:
   *x = *x * a_contract;
   *y = *y * a_contract;
   return;
  case INVMAP:
   *x = *x / a_contract;
   *y = *y / a_contract;
   return;
 }
}

void f_sqwhorl(real* x, real* y, Opcode op)
{
 real rX, rY, sX, sY;
 switch(op)
 {
  case LOCORNER: *x = -0.001; *y = -0.001; return;
  case HICORNER: *x = +1.001; *y = +1.001; return;
  case DIRMAP:
   rX = *x; rY = *y;
   reduce_to_unit_square(&rX,&rY);
   if (rY<0.5*rX)      { sX=0.5+0.5*rX;  sY=2.0*rY; }
   else if (rY<rX)     { sX=0.5+rX-rY;   sY=rX;     }
   else if (rY<2.0*rX) { sX=0.5-rY+rX;   sY=rY;     }
   else                { sX=0.5-0.5*rY;  sY=2.0*rX; }
   *x = sX; *y = sY;
   return;
  case INVMAP:
   sX = *x; sY = *y;
   reduce_to_unit_square(&sX,&sY);
   if (sX>0.5+0.5*sY)      { rX=2.0*sX-1.0; rY=0.5*sY;     }
   else if (sX>0.5)        { rX=sY;         rY=0.5+sY-sX;  }
   else if (sX>0.5-0.5*sY) { rX=sX+sY-0.5;  rY=sY;         }
   else                    { rX=0.5*sY;     rY=1.0-2.0*sX; }
   *x = rX; *y = rY;
   return;
 }
}

void f_triwhorl(real* x, real* y, Opcode op)
{
 real a,b,c, A,B,C;
 switch(op)
 {
  case LOCORNER: *x = -1.001; *y = -1.001; return;
  case HICORNER: *x = +1.001; *y = +1.001; return;
  case DIRMAP:
   cart_to_bary(*x,*y, &a,&b,&c);
   reduce_to_unit_triangle(&a,&b,&c);
   if (b > c) 
   {A = 0.5*a + b -c; B = 2*c; C = 0.5*a;}
   else
   {A = 0.5*a; B = 2*b; C = 0.5*a -b + c;}
   bary_to_cart(A,B,C, x,y);
   return;
  case INVMAP:
   cart_to_bary(*x,*y, &A,&B,&C);
   reduce_to_unit_triangle(&A,&B,&C);
   if (A > C) 
   {a = C; b = A + 0.5*B -C; c = 0.5*B;}
   else
   {a = 2*A; b = 0.5*B; c = C - A + 0.5*B;}
   bary_to_cart(a,b,c, x,y);
   return;
 }
}

void f_sqturn(real* x, real* y, Opcode op)
{
  real rX,rY,sX,sY;
 switch(op)
 {
  case LOCORNER: *x = -0.001; *y = -0.001; return;
  case HICORNER: *x = +1.001; *y = +1.001; return;
  case DIRMAP:
   rX = *x; rY = *y;
   sX = rY;
   sY = 1-rX;
   *x = sX; *y = sY;
   return;
  case INVMAP:
   sX = *x; sY = *y;
   rX = 1-sY;
   rY = sX;
   *x = rX; *y = rY;
   return;
 }
}

void f_trinop(real* x, real* y, Opcode op)
{
 real a,b,c;
 switch(op)
 {
  case LOCORNER: *x = -1.501; *y = -1.501; return;
  case HICORNER: *x = +1.501; *y = +1.501; return;
  case DIRMAP:
   cart_to_bary(*x,*y, &a,&b,&c);
   reduce_to_unit_triangle(&a,&b,&c);
   bary_to_cart(a,b,c, x,y);
   return;
  case INVMAP:
   cart_to_bary(*x,*y, &a,&b,&c);
   reduce_to_unit_triangle(&a,&b,&c);
   bary_to_cart(a,b,c, x,y);
   return;
 }
}

void f_sqnop(real* x, real* y, Opcode op)
{
 real X, Y;
 switch(op)
 {
  case LOCORNER: *x = -0.501; *y = -0.501; return;
  case HICORNER: *x = +1.501; *y = +1.501; return;
  case DIRMAP:
   X=*x; Y=*y;
   reduce_to_unit_square(&X,&Y);
   *x = X; *y = Y;
   return;
  case INVMAP:
   X=*x; Y=*y;
   reduce_to_unit_square(&X,&Y);
   *x = X; *y = Y;
   return;
 }
}

#define C_144 (-0.80901699437494742409)
#define S_144 (+0.58778525229247312918)
/* Cosine and sine of 144 degrees. */

void f_rot144(real* x, real* y, Opcode op)
{
 real X=*x, Y=*y;
 switch(op)
 {
  case LOCORNER: *x = -1.001; *y = -1.001; return;
  case HICORNER: *x = +1.001; *y = +1.001; return;
  case DIRMAP:
   *x= + C_144*X + S_144*Y;
   *y= - S_144*X + C_144*Y;
   return;
  case INVMAP:
   *x= + C_144*X - S_144*Y;
   *y= + S_144*X + C_144*Y;
   return;
 }
}

void f_rotphi(real* x, real* y, Opcode op)
{
 real X=*x, Y=*y;
 switch(op)
 {
  case LOCORNER: *x = -1.001; *y = -1.001; return;
  case HICORNER: *x = +1.001; *y = +1.001; return;
  case DIRMAP:
   *x= + C_phi*X + S_phi*Y;
   *y= - S_phi*X + C_phi*Y;
   return;
  case INVMAP:
   *x= + C_phi*X - S_phi*Y;
   *y= + S_phi*X + C_phi*Y;
   return;
 }
}

Mapfunc *func_from_name(char *name)
{
 if (strcmp(name,"henon")==0)    { return &f_henon; }
 if (strcmp(name,"haoui")==0)    { return &f_haoui; }
 if (strcmp(name,"sqout")==0)    { return &f_sqout; }
 if (strcmp(name,"cirout")==0)   { return &f_cirout; }
 if (strcmp(name,"cirin")==0)    { return &f_cirin; }
 if (strcmp(name,"contract")==0) { return &f_contract; }
 if (strcmp(name,"sqwhorl")==0)  { return &f_sqwhorl; }
 if (strcmp(name,"triwhorl")==0) { return &f_triwhorl; }
 if (strcmp(name,"sqturn")==0)   { return &f_sqturn; }
 if (strcmp(name,"rot144")==0)   { return &f_rot144; }
 if (strcmp(name,"rotphi")==0)   { return &f_rotphi; }
 if (strcmp(name,"trinop")==0)   { return &f_trinop; }
 if (strcmp(name,"sqnop")==0)    { return &f_sqnop; }
 fprintf(stderr, "unknown function name \"%s\"\n", name);
 exit(1);
 return NULL;
}

void solve_quadratic(double A, double B, double C, real *z1p, real *z2p) 
{
 double z1, z2;
 if(A == 0)
 {/* Solve linear equation: */
  z1 = -C/B;
  z2 = INFINITY;
 }
 else
 {double D,SD;
  /* Standard formula (bad for rounding...) */
  if(A<0) { A=-A; B=-B; C=-C; }
  /* Standard formula (bad for rounding...) */
  D = B*B - 4*A*C;
  if (D < 0)
  {fprintf(stderr, "quadratic without real roots %g %g %g\n",A,B,C);
   exit(1);
  }
  else
  {SD = sqrt(D);
   z1 = (-B-SD)/(2*A);
   z2 = (-B+SD)/(2*A);
  }
 }
 if(z1p != NULL) { (*z1p) = z1; }
 if(z2p != NULL) { (*z2p) = z2; }
}

#define bary_A (0.5)
  /* Negative extent of {T}. */

#define bary_H (1.5)
  /* Height of {T}. */

#define bary_W (1.73205080756887729352)
  /* Width of {T} = {sqrt(3)}. */

void cart_to_bary(real x, real y, real *a, real *b, real *c)
{
 real X, Y;
 /* Map (x,y) to (X,Y) so that {T} goes to R = ((0,0),(1,0),(0,1)): */
 Y = (y + bary_A)/bary_H;
 X = x/bary_W + 0.5*(1 - Y);
 /* Get the bary coords (*a,*b,*c) of (x,y) relative to {R}: */
 *b = X; *c = Y; *a = 1-(X+Y);
}

void bary_to_cart(real a, real b, real c, real *x, real *y)
{
 real X, Y;
 /* Get the cart coords (X,Y) of (a,b,c) relative to R = ((0,0),(1,0),(0,1)): */
 X = b; Y = c;
 /* Map (X,Y) to (*x,*y) so that {R} goes to {T}: */
 *y = bary_H*Y - bary_A;
 *x = (X - 0.5*(1 - Y))*bary_W;
}

void reduce_to_unit_disk(real *x, real *y)
{
 real X = *x, Y = *y;
 double R2 = X*X + Y*Y;
 if (R2 > 1) { X /= R2; Y /= R2; }
 assert(X*X + Y*Y < +1.00001); 
 *x = X; *y = Y;
}

void reduce_to_unit_square(real *x, real *y)
{
 real X = *x, Y = *y;
 X = X - 2*floor(X/2); Y = Y - 2*floor(Y/2);
 if (X > 1) { X = 2-X; }
 if (Y > 1) { Y = 2-Y; }
 assert((X > -0.00001) && (X < +1.00001)); 
 assert((Y > -0.00001) && (Y < +1.00001)); 
 *x = X; *y = Y;
}

void reduce_to_signed_unit_square(real *x, real *y)
{
 real X = *x, Y = *y;
 X = X - 4*floor((X+1)/4); Y = Y - 4*floor((Y+1)/4);
 if (X > 1) { X = 2-X; }
 if (Y > 1) { Y = 2-Y; }
 assert((X > -1.00001) && (X < +1.00001)); 
 assert((Y > -1.00001) && (Y < +1.00001)); 
 *x = X; *y = Y;
}

void reduce_to_unit_triangle(real *a, real *b, real *c)
{
 real A = *a, B = *b, C = *c;
 double q;
 /* Apply translations by (+3,-3,0), (+3,0,-3): */
 q = floor((B+1)/3); B = B - 3*q; A = A + 3*q;
 q = floor((C+1)/3); C = C - 3*q; A = A + 3*q;
 /* Now {B,C} are in [-1_+2], {A} is in [-3_+3]. */
 /* Apply translations by (+2,-1,-1): */
 q = floor((A+1)/2); A = A - 2*q; B = B + q; C = C + q;
 /* Now {A,B,C} are in [-1_+2]. */
 /* Fold parts with {B} or {C} above 1: */
 if (B > +1) {real d = B - 1; B = B - 2*d; A = A + d; C = C + d;}
 if (C > +1) {real d = C - 1; C = C - 2*d; B = B + d; A = A + d;}
 /* Now {A,B,C} are in [-1_1]. */
 /* Fold negative parts: */
 if (A < 00) {real d = -A; A = A + 2*d; B = B - d; C = C - d;}
 if (B < 00) {real d = -B; B = B + 2*d; A = A - d; C = C - d;}
 if (C < 00) {real d = -C; C = C + 2*d; A = A - d; B = B - d;}
 /* We must be done: */
 assert((A > -0.00001) && (A < +1.00001)); 
 assert((B > -0.00001) && (B < +1.00001)); 
 assert((C > -0.00001) && (C < +1.00001)); 
 *a = A; *b = B; *c = C;
}
