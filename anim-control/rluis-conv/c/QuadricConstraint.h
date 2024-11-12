
#ifndef QuadricConstraint_H
#define QuadricConstraint_H



#include <Constraint.h>

typedef
  T <: Public;
  Public == Constraint.T OBJECT METHODS
      init(nat n, cxx, cyy, czz, cxy, cxz, cyz, cwx, cwy, cwz, double cww);
        /*
          Constrains vertex "n" to lie on the quadric with equation
            cxx*x^2 + cyy*y^2 + czz*z^2 + 
            cxy*x*y + cxz*x*z + cyz*y*z +
            cwx*x   + cwy*y   + cwz*z   + cww == 0.
          Currently the value of "cww" is assumed to be such that 
          the constraint is satisfied at starting time.  So,
          the constraint actually keeps constant the value of the
          quadratic form above.  We should use constraint
          stabilization...
       */;
    }

T Plane(nat n, nx, ny, double nz);
T Sphere(nat n, cwx, cwy, double cwz);
T Ellipsoid(nat n, cwx, cwy, cwz, cxx, cyy, double czz);
T Paraboloid(nat n, cwx, cwy, cxx, cyy, double czz);
  /*
    Constrains vertex "n" to lie on the specified quadric surface,
    without friction. */
;
} /* QuadricConstraint */.

#endif
