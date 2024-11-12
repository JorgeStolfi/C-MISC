
#include <ForceManager.h>



#include <Force.h>
#include <Force.h>

REVEAL
  T == Public BRANDED OBJECT
      nat m;
      forces: Force.T_vec;
    OVERRIDES
      init = Init;
      number = Number;
      addForce = AddForce;
      removeForce = RemoveForce;
      compute = Compute;
    }

void Init(T fm, nat m)
{
  fm.m = 0;
  Allocate(fm, m);
} /* Init */;

void Allocate(T fm, nat m)
{
  { /* with*/ forces == Force.T_vec_new) {
    for (i = 0;  i < fm.m; i++) {forces[i] = fm.forces[i]; }
    fm.forces = forces;
  };
} /* Allocate */;

nat Number(T fm)
{ return fm.m; } /* Number */;

void AddForce(T fm, Force.T force)
{
  if ((fm.m + 1 > NUMBER(fm.forces^) )) {
    Allocate(fm, 2*fm.m);
  }
  fm.forces[fm.m] = force;
  fm.m++;
} /* AddForce */;

void RemoveForce(T fm, Force.T force)
VAR i = 0;
{
  affirm(fm.m > 0 , "??");
  while (i < fm.m) && (fm.forces[i] != force ) { i++;}
  if ((i < fm.m )) {
    fm.m--;
    fm.forces[i] = fm.forces[fm.m];
    fm.forces[fm.m] = NULL;
  };
} /* RemoveForce */;

void Compute(T fm, pos, vel, Vector f, double t)
{
  for (j = 0;  j < fm.m; j++) {
    fm.forces[j].compute(pos, vel, f, t);
  };
} /* Compute */;

{; } /* ForceManager */.
