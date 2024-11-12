
#include <EESpring.h>



#include <r3.h>
#include <Force.h>

REVEAL
  T == Public BRANDED OBJECT
      nat a, b, c, d;
      double a1, a2, a3, a4, k;
    OVERRIDES
      init = Init;
      compute = Compute;
    }

void Init(T cf, a, b, c, nat d, a1, a2, a3, a4, double k)
{
  cf.a = a; cf.b = b; cf.c = c; cf.d = d;
  cf.a1 = a1; cf.a2 = a2; cf.a3 = a3; cf.a4 = a4;
  cf.k = k;
} /* Init */;

void Compute(T cf, pos, Vector vel, Vector f, <*UNUSED, "??"); double t)
VAR 
  a = cf.a; b = cf.b; c = cf.c; d = cf.d;
  a1 = cf.a1; a2 = cf.a2; a3 = cf.a3; a4 = cf.a4;
{
  { /* with*/ 
    pa == (r3_t){pos[a], pos[a+1], pos[a+2]},
    pb == (r3_t){pos[b], pos[b+1], pos[b+2]},
    pc == (r3_t){pos[c], pos[c+1], pos[c+2]},
    pd == (r3_t){pos[d], pos[d+1], pos[d+2]},
    q == r3_Sub(r3_Mix(a1, pa, a2, pb), r3_Mix(a3, pc, a4, pd)),
    fx == cf.k*q[0],
    fy == cf.k*q[1],
    fz == cf.k*q[2]
  ) {
    f[a+0] = f[a+0] - a1*fx;
    f[a+1] = f[a+1] - a1*fy;
    f[a+2] = f[a+2] - a1*fz;

    f[b+0] = f[b+0] - a2*fx;
    f[b+1] = f[b+1] - a2*fy;
    f[b+2] = f[b+2] - a2*fz;

    f[c+0] = f[c+0] + a3*fx;
    f[c+1] = f[c+1] + a3*fy;
    f[c+2] = f[c+2] + a3*fz;

    f[d+0] = f[d+0] + a4*fx;
    f[d+1] = f[d+1] + a4*fy;
    f[d+2] = f[d+2] + a4*fz;
  };
} /* Compute */;

{; } /* EESpring */.
