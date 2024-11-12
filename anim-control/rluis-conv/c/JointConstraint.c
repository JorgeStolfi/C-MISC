
#include <JointConstraint.h>



#include <Constraint.h>

REVEAL
  T == Public BRANDED OBJECT
      nat n, n1, n2, n3;
      double a, b, c;
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

void Init(T q, n, n1, n2, nat n3, a, b, double c)
{
  q.n = 3*n; q.n1 = 3*n1; q.n2 = 3*n2; q.n3 = 3*n3;
  q.a = a; q.b = b; q.c = c;
} /* Init */;

nat NumEquations(<*UNUSED, "??"); T q)
{ 
  return 3;
} /* NumEquations */;

PROCEDURE ComputeGrad(<*UNUSED, "??"); T q, <*UNUSED, "??"); Vector pos,
                      <*UNUSED, "??"); double t) == 
{; } /* ComputeGrad */;
 
PROCEDURE ComputePsi(<*UNUSED, "??"); T q, <*UNUSED, "??"); pos, Vector vel,
                     <*UNUSED, "??"); double t) == 
{; } /* ComputePsi */;
 
void GetGrad(T q, nat j, Vector grad)
{
  grad[q.n +j] = 1.0D0;
  grad[q.n1+j] = -q.a;
  grad[q.n2+j] = -q.b;
  grad[q.n3+j] = -q.c;
} /* GetGrad */;

double DirDerivative(T q, nat j, Vector u)
{
  return u[q.n+j] - q.a*u[q.n1+j] - q.b*u[q.n2+j] - q.c*u[q.n3+j];
} /* DirDerivative */;

double GetPsi(<*UNUSED, "??"); T q, <*UNUSED, "??"); nat j)
{ return 0.0D0; } /* GetPsi */;

PROCEDURE SetPsi(<*UNUSED, "??"); T q, <*UNUSED, "??"); nat j, 
                 <*UNUSED, "??"); double psi) == 
{; } /* SetPsi */;

void AddReactionForce(T q, nat j, double lambda, Vector f)
{
  f[q.n +j] = f[q.n +j] -     lambda;
  f[q.n1+j] = f[q.n1+j] + q.a*lambda;
  f[q.n2+j] = f[q.n2+j] + q.b*lambda;
  f[q.n3+j] = f[q.n3+j] + q.c*lambda;
} /* AddReactionForce */;


/*--- ESPECIFIC JOINTS ------------------------------------------------------*/

T PointToPoint(n, nat n1)
{
  { /* with*/ j == NEW) {
    j.init(n, n1, 0, 0, 1.0D0, 0.0D0, 0.0D0);
    return j;
  };
} /* PointToPoint */;

T PointToLine(n, n1, nat n2, a, double b)
{
  { /* with*/ j == NEW) {
    j.init(n, n1, n2, 0, a, b, 0.0D0);
    return j;
  };
} /* PointToLine */;

T PointToFace(n, n1, n2, nat n3, a, b, double c)
{
  { /* with*/ j == NEW) {
    j.init(n, n1, n2, n3, a, b, c);
    return j;
  };
} /* PointToFace */;

{; } /* JointConstraint */.
