/* Last edited on 2023-01-14 12:22:17 by stolfi */

#include <memory.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <arith.h>
#include <fn.h>

/* TYPES */

typedef struct{
  float ux, uy, uz, vx, vy, vz, wx, wy, wz;
} Matrix;

/* PROTOTYPES */

void fnlinmap(AAform x, AAform y, AAform z, Matrix *m, AAform u, AAform v, AAform w);
void fnellip(real A, AAform x, real B, AAform y, real C, AAform z, AAform f);
void fnplane(real A, AAform x, real B, AAform y, real C, AAform z, real D, AAform f);

/* IMPLEMENTATIONS */

#define SQ2  1.4142135623730950488
#define SQ18 (3*SQ2)

static Matrix m1 = (Matrix){
  +2.0/3.0,  +2.0/3.0,  +1.0/3.0,
  +1.0/SQ2,  -1.0/SQ2,  00.0, 
  -1.0/SQ18, -1.0/SQ18, +4.0/SQ18
};

static Matrix m2 = (Matrix){
  +1.0/SQ18, +1.0/SQ18, +4.0/SQ18
  +1.0/SQ2,  -1.0/SQ2,  00.0, 
  -2.0/3.0,  -2.0/3.0,  +1.0/3.0,
};

char *fnname(void)
{ 
  return("fncoeur"); 
} 
  
void fndescr(FILE *f)
{
  fprintf(f, "GIVEN (x, y, z) DO\n");
  fprintf(f, "  (u1, v1, w1) = (x,y,z) M1;\n");
  fprintf(f, "  (u2, v2, w2) = (x,y,z) M2;\n");
  fprintf(f, "  e1 = (u1/a1)^2 + (v1/b1)^2 + (w1/c1)^2 - 1;\n");
  fprintf(f, "  e2 = (u2/a2)^2 + (v2/b2)^2 + (w2/c2)^2 - 1;\n");
  fprintf(f, "  p = (tanh(a3*x + b3*y + c3*z) - 1)/2; q = 1-p;\n");
  fprintf(f, "  f = p * e1 + q * e2;\n");
  fprintf(f, "RETURN (f)\n");
}

void fneval(AAform x, AAform y, AAform z, AAform f)
{
  AAform u1, v1, w1, u2, v2, w2, e1, e2, p, q;
  
  fnlinmap(x, y, z, &m1, u1, v1, w1);
  fnellip(0.5, u1, 0.5, v1, 1.0, w1, e1);
  fnlinmap(x, y, z, &m2, u2, v2, w2);
  fnellip(0.5, u2, 0.5, v2, 1.0, w2, e2);
  
  fnplane(SQ2, x, SQ2, y, -1.0, z, 0.0, p);
  aa_tanh(p, p);
  aa_trans(p, p, 1.0);
  aa_scale(p, p, 0.5);
  aa_neg(q, p);
  aa_trans(q, q, 1.0);
  aa_mul(e1, e1, p);
  aa_mul(e2, e2, q);
  aa_add(f, e1, e2);
}

void fnlinmap(AAform x, AAform y, AAform z, Matrix *m, AAform u, AAform v, AAform w)
{
  AAform t;
  aa_scale(u, x, m->ux);
  aa_scale(t, y, m->uy); aa_add(u, u, t);
  aa_scale(t, z, m->uz); aa_add(u, u, t);
  
  aa_scale(v, x, m->vx);
  aa_scale(t, y, m->vy); aa_add(v, v, t);
  aa_scale(t, z, m->vz); aa_add(v, v, t);

  aa_scale(w, x, m->wx);
  aa_scale(t, y, m->wy); aa_add(w, w, t);
  aa_scale(t, z, m->wz); aa_add(w, w, t);
}

void fnellip(real A, AAform x, real B, AAform y, real C, AAform z, AAform f)
{
  AAform s;
  aa_scale(s, x, 1.0/A); aa_sqr(f, s); ; 
  aa_scale(s, y, 1.0/B); aa_sqr(s, s); aa_add(f, f, s);
  aa_scale(s, z, 1.0/C); aa_sqr(s, s); aa_add(f, f, s);
  aa_trans(f, f, -1.0);
}

void fnplane(real A, AAform x, real B, AAform y, real C, AAform z, real D, AAform f)
{
  AAform s;
  aa_scale(f, x, A);
  aa_scale(s, y, B); aa_add(f, f, s);
  aa_scale(s, z, C); aa_add(f, f, s);
  aa_trans(f, f, D);
}
