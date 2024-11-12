MODULE QuadricConstraint;

FROM Constraint IMPORT Vector;

REVEAL
  T = Public BRANDED OBJECT
      n: CARDINAL;
      Hx, Hy, Hz,
      psis,
      cxx, cyy, czz, cxy, cxz, cyz, cwx, cwy, cwz, cww: LONGREAL;
    OVERRIDES
      init         := Init;
      numEquations := NumEquations;
      computeGrad  := ComputeGrad;
      computePsi   := ComputePsi;
      getGrad       := GetGrad;
      dirDerivative     := DirDerivative;
      getPsi          := GetPsi;
      setPsi       := SetPsi;
      addReactionForce := AddReactionForce;
    END;

PROCEDURE Init(q: T; 
  n: CARDINAL; 
  cxx, cyy, czz, cxy, cxz, cyz, cwx, cwy, cwz, cww: LONGREAL
) =
BEGIN
  q.n := 3*n;
  q.cxx := cxx; q.cyy := cyy; q.czz := czz; q.cxy := cxy; q.cxz := cxz; 
  q.cyz := cyz; q.cwx := cwx; q.cwy := cwy; q.cwz := cwz; q.cww := cww;
END Init;

PROCEDURE NumEquations(<*UNUSED*> q: T): CARDINAL =
BEGIN 
  RETURN 1
END NumEquations;

PROCEDURE ComputeGrad(q: T; pos: Vector; <*UNUSED*> t: LONGREAL) =
BEGIN
  WITH
    n  = q.n,
    x  = pos[n], y = pos[n+1], z = pos[n+2],
    cxx2 = q.cxx + q.cxx, 
    cyy2 = q.cyy + q.cyy, 
    czz2 = q.czz + q.czz,
    cxy  = q.cxy, cxz = q.cxz, cyz = q.cyz
  DO
    q.Hx := cxx2*x + cxy*y + cxz*z + q.cwx;
    q.Hy := cyy2*y + cxy*x + cyz*z + q.cwy;
    q.Hz := czz2*z + cxz*x + cyz*y + q.cwz;
  END
END ComputeGrad;

PROCEDURE ComputePsi(q: T; <*UNUSED*> pos: Vector; vel: Vector; <*UNUSED*> t: LONGREAL) =
BEGIN
  WITH
    n  = q.n,
    vx = vel[n], vy = vel[n+1], vz = vel[n+2],
    cxx2 = q.cxx + q.cxx, 
    cyy2 = q.cyy + q.cyy, 
    czz2 = q.czz + q.czz,
    cxy  = q.cxy, cxz = q.cxz, cyz = q.cyz,
    v1 = vx*cxx2 + vy*cxy  + vz*cxz,
    v2 = vx*cxy  + vy*cyy2 + vz*cyz,
    v3 = vx*cxz  + vy*cyz  + vz*czz2
  DO
    q.psis := -(v1*vx + v2*vy + v3*vz);
  END
END ComputePsi;

PROCEDURE GetGrad(q: T; j: CARDINAL; grad: Vector) =
BEGIN
  <* ASSERT j = 0 *>
  grad[q.n] := q.Hx; grad[q.n+1] := q.Hy; grad[q.n+2] := q.Hz
END GetGrad;

PROCEDURE DirDerivative(q: T; j: CARDINAL; u: Vector): LONGREAL =
BEGIN
  <* ASSERT j = 0 *>
  RETURN q.Hx*u[q.n] + q.Hy*u[q.n+1] + q.Hz*u[q.n+2]
END DirDerivative;

PROCEDURE GetPsi(q: T; j: CARDINAL): LONGREAL =
BEGIN 
  <* ASSERT j = 0 *>
  RETURN q.psis
END GetPsi;

PROCEDURE SetPsi(q: T; j: CARDINAL; psi: LONGREAL) =
BEGIN 
  <* ASSERT j = 0 *>
  q.psis := psi
END SetPsi;

PROCEDURE AddReactionForce(q: T; j: CARDINAL; lambda: LONGREAL; f: Vector) =
BEGIN
  WITH n = q.n DO
    f[n+0] := f[n+0] - q.Hx*lambda;
    f[n+1] := f[n+1] - q.Hy*lambda;
    f[n+2] := f[n+2] - q.Hz*lambda;
  END
END AddReactionForce;


(*--- ESPECIFIC QUADRICS ----------------------------------------------------*)

PROCEDURE Plane(n: CARDINAL; nx, ny, nz: LONGREAL): T =
BEGIN
  WITH q = NEW(T) DO
    q.init(n, 
      nx, ny, nz, 
      0.0D0, 0.0D0, 0.0D0, 
      0.0D0, 0.0D0, 0.0D0, 
      0.0D0
    );
    RETURN q;
  END
END Plane;

PROCEDURE Sphere(n: CARDINAL; cx, cy, cz: LONGREAL): T =
BEGIN
  WITH q = NEW(T) DO
    q.init(n, 
      1.0D0, 1.0D0, 1.0D0, 
      0.0D0, 0.0D0, 0.0D0, 
      -2.0D0*cx, -2.0D0*cy, -2.0D0*cz, 
      0.0D0
    );
    RETURN q;
  END
END Sphere;

PROCEDURE Ellipsoid(n: CARDINAL; cx, cy, cz, cxx, cyy, czz: LONGREAL): T =
BEGIN
  WITH
    p = 1.0D0/(cxx*cxx),
    q = 1.0D0/(cyy*cyy),
    r = 1.0D0/(czz*czz),
    cxz = NEW(T)
  DO
    cxz.init(n, 
      p, q, r, 
      0.0D0, 0.0D0, 0.0D0, 
      -2.0D0*cx*p, -2.0D0*cy*q, -2.0D0*cz*r,
      0.0D0
    );
    RETURN cxz;
  END
END Ellipsoid;

PROCEDURE Paraboloid(n: CARDINAL; cx, cy, cxx, cyy, czz: LONGREAL): T =
BEGIN
  WITH
    p = 1.0D0/(cxx*cxx),
    q = 1.0D0/(cyy*cyy),
    r = 1.0D0/(czz*czz),
    t = NEW(T)
  DO
    t.init(n, 
      p, q, r, 
      0.0D0, 0.0D0, 0.0D0, 
      -2.0D0*cx*p, -2.0D0*cy*q, -czz,
      0.0D0
    );
    RETURN t;
  END
END Paraboloid;

BEGIN END QuadricConstraint.
