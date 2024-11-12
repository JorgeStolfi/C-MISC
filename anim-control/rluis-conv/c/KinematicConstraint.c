
#include <KinematicConstraint.h>



#include <TimedConstraint.h>
#include <Constraint.h>

REVEAL
  T == Public BRANDED OBJECT
      nat k;
      double xi;
      double psis;
      double a, b, c, d;
    OVERRIDES
      init = Init;
      numEquations = NumEquations;
      treatEvent = TreatEvent;
      computeGrad = ComputeGrad;
      computePsi = ComputePsi;
      getGrad = GetGrad;
      dirDerivative = DirDerivative;
      getPsi = GetPsi;
      setPsi = SetPsi;
      addReactionForce = AddReactionForce;
    }

void Init(T tc, nat k, ta, double tb, p, double v, double xi)
{
  affirm(ta < tb , "??");
  NARROW(tc, TimedConstraint.T).init(ta, tb);
  tc.k = k;
  tc.xi = xi;
  tc.c = v;
  tc.d = p;
} /* Init */;

nat NumEquations(T tc)
{ 
  affirm(tc.started , "??");
  if ((tc.active )) { return 1 } else { return 0;
} /* NumEquations */;

void TreatEvent(T tc, double t, pos, Vector vel)
{
  affirm(tc.started , "??");
  if ((! tc.active )) {
    affirm(t == tc.ta , "??");
    { /* with*/ 
      k == tc.k,
      dt == tc.tb - t, dt2 == dt*dt, dt3 == dt2*dt,
      d == tc.d - pos[k],
      v == vel[k] + tc.c
    ) {
      affirm(dt > 0.0D0 , "??");
      tc.a = -2.0D0*d/dt3 + v/dt2;
      tc.b = 3.0D0*d/dt2 - (vel[k] + v)/dt;
      tc.c = vel[k];
      tc.d = pos[k];
    }
  } else {
    affirm(t == tc.tb , "??");
  }
  tc.setClock(t);
} /* TreatEvent */;

void ComputeGrad(T tc, <*UNUSED, "??"); Vector pos, <*UNUSED, "??"); double t)
{ 
  affirm(tc.started , "??");
} /* ComputeGrad */;

void ComputePsi(T tc, pos, Vector vel, double t)
{
  affirm(tc.started , "??");
  if ((tc.active )) {
    { /* with*/ 
      k == tc.k,
      xi == tc.xi,
      dt == t - tc.ta, dt2 == dt*dt, dt3 == dt2*dt,
      b2 == tc.b + tc.b,
      h == pos[k] - tc.a*dt3 - tc.b*dt2 - tc.c*dt - tc.d,
      dh == vel[k] - 3.0D0*tc.a*dt2 - b2*dt - tc.c
    ) {
      tc.psis = 6.0D0*tc.a*dt + b2 - xi*(dh + dh + xi*h);
    };
  };
} /* ComputePsi */;

void GetGrad(T tc, nat j, Vector grad)
{ 
  affirm(tc.started , "??");
  affirm(tc.active , "??");
  affirm(j == 0 , "??");
  grad[tc.k] = 1.0D0;
} /* GetGrad */;

double DirDerivative(T tc, nat j, Vector u)
{ 
  affirm(tc.started , "??");
  affirm(tc.active , "??");
  affirm(j == 0 , "??");
  return u[tc.k];
} /* DirDerivative */;

double GetPsi(T tc, nat j)
{ 
  affirm(tc.started , "??");
  affirm(tc.active , "??");
  affirm(j == 0 , "??");
  return tc.psis;
} /* GetPsi */;

void SetPsi(T tc, nat j, double psi)
{ 
  affirm(tc.started , "??");
  affirm(tc.active , "??");
  affirm(j == 0 , "??");
  tc.psis = psi;
} /* SetPsi */;

void AddReactionForce(T tc, nat j, double lambda, Vector f)
{ 
  affirm(tc.started , "??");
  affirm(tc.active , "??");
  affirm(j == 0 , "??");
  f[tc.k] = f[tc.k] - lambda ;
} /* AddReactionForce */;

{; } /* KinematicConstraint */.
