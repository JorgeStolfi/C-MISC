
#include <QuadricConstraint.h>



#include <Constraint.h>

REVEAL
  T == Public BRANDED OBJECT
      nat n;
      Hx, Hy, Hz,
      psis,
      double cxx, cyy, czz, cxy, cxz, cyz, cwx, cwy, cwz, cww;
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

PROCEDURE Init(T q, 
  nat n; 
  cxx, cyy, czz, cxy, cxz, cyz, cwx, cwy, cwz, cww: double
) == 
{
  q.n = 3*n;
  q.cxx = cxx; q.cyy = cyy; q.czz = czz; q.cxy = cxy; q.cxz = cxz; 
  q.cyz = cyz; q.cwx = cwx; q.cwy = cwy; q.cwz = cwz; q.cww = cww;
} /* Init */;

nat NumEquations(<*UNUSED, "??"); T q)
{ 
  return 1;
} /* NumEquations */;

void ComputeGrad(T q, Vector pos, <*UNUSED, "??"); double t)
{
  { /* with*/ 
    n == q.n,
    x == pos[n], y == pos[n+1], z == pos[n+2],
    cxx2 == q.cxx + q.cxx, 
    cyy2 == q.cyy + q.cyy, 
    czz2 == q.czz + q.czz,
    cxy == q.cxy, cxz == q.cxz, cyz == q.cyz
  ) {
    q.Hx = cxx2*x + cxy*y + cxz*z + q.cwx;
    q.Hy = cyy2*y + cxy*x + cyz*z + q.cwy;
    q.Hz = czz2*z + cxz*x + cyz*y + q.cwz;
  };
} /* ComputeGrad */;

void ComputePsi(T q, <*UNUSED, "??"); Vector pos, Vector vel, <*UNUSED, "??"); double t)
{
  { /* with*/ 
    n == q.n,
    vx == vel[n], vy == vel[n+1], vz == vel[n+2],
    cxx2 == q.cxx + q.cxx, 
    cyy2 == q.cyy + q.cyy, 
    czz2 == q.czz + q.czz,
    cxy == q.cxy, cxz == q.cxz, cyz == q.cyz,
    v1 == vx*cxx2 + vy*cxy  + vz*cxz,
    v2 == vx*cxy  + vy*cyy2 + vz*cyz,
    v3 == vx*cxz  + vy*cyz  + vz*czz2
  ) {
    q.psis = -(v1*vx + v2*vy + v3*vz);
  };
} /* ComputePsi */;

void GetGrad(T q, nat j, Vector grad)
{
  affirm(j == 0 , "??");
  grad[q.n] = q.Hx; grad[q.n+1] = q.Hy; grad[q.n+2] = q.Hz;
} /* GetGrad */;

double DirDerivative(T q, nat j, Vector u)
{
  affirm(j == 0 , "??");
  return q.Hx*u[q.n] + q.Hy*u[q.n+1] + q.Hz*u[q.n+2];
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
  { /* with*/ n == q.n ) {
    f[n+0] = f[n+0] - q.Hx*lambda;
    f[n+1] = f[n+1] - q.Hy*lambda;
    f[n+2] = f[n+2] - q.Hz*lambda;
  };
} /* AddReactionForce */;


/*--- ESPECIFIC QUADRICS ----------------------------------------------------*/

T Plane(nat n, nx, ny, double nz)
{
  { /* with*/ q == NEW) {
    q.init(n, 
      nx, ny, nz, 
      0.0D0, 0.0D0, 0.0D0, 
      0.0D0, 0.0D0, 0.0D0, 
      0.0D0
    );
    return q;
  };
} /* Plane */;

T Sphere(nat n, cx, cy, double cz)
{
  { /* with*/ q == NEW) {
    q.init(n, 
      1.0D0, 1.0D0, 1.0D0, 
      0.0D0, 0.0D0, 0.0D0, 
      -2.0D0*cx, -2.0D0*cy, -2.0D0*cz, 
      0.0D0
    );
    return q;
  };
} /* Sphere */;

T Ellipsoid(nat n, cx, cy, cz, cxx, cyy, double czz)
{
  { /* with*/ 
    p == 1.0D0/(cxx*cxx),
    q == 1.0D0/(cyy*cyy),
    r == 1.0D0/(czz*czz),
    cxz == NEW(T)
  ) {
    cxz.init(n, 
      p, q, r, 
      0.0D0, 0.0D0, 0.0D0, 
      -2.0D0*cx*p, -2.0D0*cy*q, -2.0D0*cz*r,
      0.0D0
    );
    return cxz;
  };
} /* Ellipsoid */;

T Paraboloid(nat n, cx, cy, cxx, cyy, double czz)
{
  { /* with*/ 
    p == 1.0D0/(cxx*cxx),
    q == 1.0D0/(cyy*cyy),
    r == 1.0D0/(czz*czz),
    t == NEW(T)
  ) {
    t.init(n, 
      p, q, r, 
      0.0D0, 0.0D0, 0.0D0, 
      -2.0D0*cx*p, -2.0D0*cy*q, -czz,
      0.0D0
    );
    return t;
  };
} /* Paraboloid */;

{; } /* QuadricConstraint */.
