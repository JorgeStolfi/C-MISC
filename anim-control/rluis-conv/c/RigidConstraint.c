
#include <RigidConstraint.h>



#include <Constraint.h>

REVEAL
  T == Public BRANDED OBJECT
      nat n1, n2;
      double hx, hy, hz, psis;
    OVERRIDES
      init = Init;
      numEquations = NumEquations;
      computeGrad = ComputeGrad;
      computePsi = ComputePsi;
      getGrad = GetGrad;
      dirDerivative = DirDerivative;
      getPsi = GetPsi;
      setPsi = SetPsi;
      addReactionForce = AddReactionForce;
    }

void Init(T q, n1, nat n2)
{
  n1 = 3*n1; n2 = 3*n2;
  q.n1 = n1; q.n2 = n2;
} /* Init */;

nat NumEquations(<*UNUSED, "??"); T q)
{ 
  return 1;
} /* NumEquations */;

void ComputeGrad(T q, Vector pos, <*UNUSED, "??"); double t)
{
  { /* with*/ 
    n1 == q.n1, n2 == q.n2,
    x1 == pos[n1]+pos[n1], y1 == pos[n1+1]+pos[n1+1], z1 == pos[n1+2]+pos[n1+2],
    x2 == pos[n2]+pos[n2], y2 == pos[n2+1]+pos[n2+1], z2 == pos[n2+2]+pos[n2+2]
  ) {
    q.hx = x1 - x2; q.hy = y1 - y2; q.hz = z1 - z2;
  };
} /* ComputeGrad */;

void ComputePsi(T q, <*UNUSED, "??"); Vector pos, Vector vel, <*UNUSED, "??"); double t)
{
  { /* with*/ 
    n1 == q.n1, n2 == q.n2,
    vx1 == vel[n1], vy1 == vel[n1+1], vz1 == vel[n1+2],
    vx2 == vel[n2], vy2 == vel[n2+1], vz2 == vel[n2+2],
    vx12 == vx1 - vx2, vy12 == vy1 - vy2, vz12 == vz1 - vz2
  ) {
    q.psis = -2.0D0*(vx12*vx12 + vy12*vy12 + vz12*vz12);
  };
} /* ComputePsi */;

void GetGrad(T q, nat j, Vector grad)
{
  affirm(j == 0 , "??");
  { /* with*/ n1 == q.n1, n2 == q.n2 ) {
    grad[n1] = q.hx; grad[n1+1] = q.hy; grad[n1+2] = q.hz;
    grad[n2] = -q.hx; grad[n2+1] = -q.hy; grad[n2+2] = -q.hz;
  };
} /* GetGrad */;

double DirDerivative(T q, nat j, Vector u)
{
  affirm(j == 0 , "??");
  { /* with*/ n1 == q.n1, n2 == q.n2 ) {
    return q.hx*(u[n1+0] - u[n2+0]) +
           q.hy*(u[n1+1] - u[n2+1]) +
           q.hz*(u[n1+2] - u[n2+2]);
  };
} /* DirDerivative */;

double GetPsi(T q, nat j)
{
  affirm(j == 0 , "??");
  return q.psis;
} /* GetPsi */;

void SetPsi(T q, nat j, double psi)
{ 
  affirm(j == 0 , "??");
  q.psis = psi;
} /* SetPsi */;

void AddReactionForce(T q, nat j, double lambda, Vector f)
{
  affirm(j == 0 , "??");
  { /* with*/ 
    n1 == q.n1, n2 == q.n2,
    a == q.hx*lambda, b == q.hy*lambda, c == q.hz*lambda
  ) {
    f[n1] = f[n1] - a; f[n1+1] = f[n1+1] - b; f[n1+2] = f[n1+2] - c;
    f[n2] = f[n2] + a; f[n2+1] = f[n2+1] + b; f[n2+2] = f[n2+2] + c;
  };
} /* AddReactionForce */;

{; } /* RigidConstraint */.
