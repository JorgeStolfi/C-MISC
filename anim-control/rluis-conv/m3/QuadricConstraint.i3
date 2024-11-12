INTERFACE QuadricConstraint;

IMPORT Constraint;

TYPE
  T <: Public;
  Public = Constraint.T OBJECT METHODS
      init(n: CARDINAL; cxx, cyy, czz, cxy, cxz, cyz, cwx, cwy, cwz, cww: LONGREAL);
        (*
          Constrains vertex "n" to lie on the quadric with equation
            cxx*x^2 + cyy*y^2 + czz*z^2 + 
            cxy*x*y + cxz*x*z + cyz*y*z +
            cwx*x   + cwy*y   + cwz*z   + cww = 0.
          Currently the value of "cww" is assumed to be such that 
          the constraint is satisfied at starting time.  So,
          the constraint actually keeps constant the value of the
          quadratic form above.  We should use constraint
          stabilization...
       *)
    END;

PROCEDURE Plane(n: CARDINAL; nx, ny, nz: LONGREAL): T;
PROCEDURE Sphere(n: CARDINAL; cwx, cwy, cwz: LONGREAL): T;
PROCEDURE Ellipsoid(n: CARDINAL; cwx, cwy, cwz, cxx, cyy, czz: LONGREAL): T;
PROCEDURE Paraboloid(n: CARDINAL; cwx, cwy, cxx, cyy, czz: LONGREAL): T;
  (*
    Constrains vertex "n" to lie on the specified quadric surface,
    without friction. *)

END QuadricConstraint.
