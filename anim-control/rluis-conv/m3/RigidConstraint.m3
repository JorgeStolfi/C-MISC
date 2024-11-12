MODULE RigidConstraint;

FROM Constraint IMPORT Vector;

REVEAL
  T = Public BRANDED OBJECT
      n1, n2: CARDINAL;
      hx, hy, hz, psis: LONGREAL;
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

PROCEDURE Init(q: T; n1, n2: CARDINAL) =
BEGIN
  n1 := 3*n1; n2 := 3*n2;
  q.n1 := n1; q.n2 := n2;
END Init;

PROCEDURE NumEquations(<*UNUSED*> q: T): CARDINAL =
BEGIN 
  RETURN 1
END NumEquations;

PROCEDURE ComputeGrad(q: T; pos: Vector; <*UNUSED*> t: LONGREAL) =
BEGIN
  WITH
    n1 = q.n1, n2 = q.n2,
    x1 = pos[n1]+pos[n1], y1 = pos[n1+1]+pos[n1+1], z1 = pos[n1+2]+pos[n1+2],
    x2 = pos[n2]+pos[n2], y2 = pos[n2+1]+pos[n2+1], z2 = pos[n2+2]+pos[n2+2]
  DO
    q.hx := x1 - x2; q.hy := y1 - y2; q.hz := z1 - z2;
  END
END ComputeGrad;

PROCEDURE ComputePsi(q: T; <*UNUSED*> pos: Vector; vel: Vector; <*UNUSED*> t: LONGREAL) =
BEGIN
  WITH
    n1 = q.n1, n2 = q.n2,
    vx1 = vel[n1], vy1 = vel[n1+1], vz1 = vel[n1+2],
    vx2 = vel[n2], vy2 = vel[n2+1], vz2 = vel[n2+2],
    vx12 = vx1 - vx2, vy12 = vy1 - vy2, vz12 = vz1 - vz2
  DO
    q.psis := -2.0D0*(vx12*vx12 + vy12*vy12 + vz12*vz12);
  END
END ComputePsi;

PROCEDURE GetGrad(q: T; j: CARDINAL; grad: Vector) =
BEGIN
  <* ASSERT j = 0 *>
  WITH n1 = q.n1, n2 = q.n2 DO
    grad[n1] :=  q.hx; grad[n1+1] :=  q.hy; grad[n1+2] :=  q.hz;
    grad[n2] := -q.hx; grad[n2+1] := -q.hy; grad[n2+2] := -q.hz;
  END
END GetGrad;

PROCEDURE DirDerivative(q: T; j: CARDINAL; u: Vector): LONGREAL =
BEGIN
  <* ASSERT j = 0 *>
  WITH n1 = q.n1, n2 = q.n2 DO
    RETURN q.hx*(u[n1+0] - u[n2+0]) +
           q.hy*(u[n1+1] - u[n2+1]) +
           q.hz*(u[n1+2] - u[n2+2])
  END
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
  <* ASSERT j = 0 *>
  WITH
    n1 = q.n1, n2 = q.n2,
    a = q.hx*lambda, b = q.hy*lambda, c = q.hz*lambda
  DO
    f[n1] := f[n1] - a; f[n1+1] := f[n1+1] - b; f[n1+2] := f[n1+2] - c;
    f[n2] := f[n2] + a; f[n2+1] := f[n2+1] + b; f[n2+2] := f[n2+2] + c;
  END
END AddReactionForce;

BEGIN END RigidConstraint.
