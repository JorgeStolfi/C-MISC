MODULE KinematicConstraint;

IMPORT TimedConstraint;
FROM Constraint IMPORT Vector;

REVEAL
  T = Public BRANDED OBJECT
      k: CARDINAL;
      xi: LONGREAL;
      psis: LONGREAL;
      a, b, c, d: LONGREAL;
    OVERRIDES
      init              := Init;
      numEquations      := NumEquations;
      treatEvent        := TreatEvent;
      computeGrad       := ComputeGrad;
      computePsi        := ComputePsi;
      getGrad            := GetGrad;
      dirDerivative          := DirDerivative;
      getPsi               := GetPsi;
      setPsi            := SetPsi;
      addReactionForce      := AddReactionForce;
    END;

PROCEDURE Init(tc: T; k: CARDINAL; ta, tb: LONGREAL; p, v: LONGREAL; xi: LONGREAL) =
BEGIN
  <* ASSERT ta < tb *>
  NARROW(tc, TimedConstraint.T).init(ta, tb);
  tc.k  := k;
  tc.xi := xi;
  tc.c := v;
  tc.d := p;
END Init;

PROCEDURE NumEquations(tc: T): CARDINAL =
BEGIN 
  <* ASSERT tc.started *>
  IF tc.active THEN RETURN 1 ELSE RETURN 0 END;
END NumEquations;

PROCEDURE TreatEvent(tc: T; t: LONGREAL; pos, vel: Vector) =
BEGIN
  <* ASSERT tc.started *>
  IF NOT tc.active THEN
    <* ASSERT t = tc.ta *>
    WITH
      k = tc.k,
      dt = tc.tb - t, dt2 = dt*dt, dt3 = dt2*dt,
      d = tc.d - pos[k],
      v = vel[k] + tc.c
    DO
      <* ASSERT dt > 0.0D0 *>
      tc.a := -2.0D0*d/dt3 + v/dt2;
      tc.b := 3.0D0*d/dt2 - (vel[k] + v)/dt;
      tc.c := vel[k];
      tc.d := pos[k];
    END
  ELSE
    <* ASSERT t = tc.tb *>
  END;
  tc.setClock(t)
END TreatEvent;

PROCEDURE ComputeGrad(tc: T; <*UNUSED*> pos: Vector; <*UNUSED*> t: LONGREAL) =
BEGIN 
  <* ASSERT tc.started *>
END ComputeGrad;

PROCEDURE ComputePsi(tc: T; pos, vel: Vector; t: LONGREAL) =
BEGIN
  <* ASSERT tc.started *>
  IF tc.active THEN
    WITH
      k = tc.k,
      xi = tc.xi,
      dt  = t - tc.ta, dt2 = dt*dt, dt3 = dt2*dt,
      b2 = tc.b + tc.b,
      h  = pos[k] - tc.a*dt3 - tc.b*dt2 - tc.c*dt - tc.d,
      dh = vel[k] - 3.0D0*tc.a*dt2 - b2*dt - tc.c
    DO
      tc.psis := 6.0D0*tc.a*dt + b2 - xi*(dh + dh + xi*h);
    END
  END
END ComputePsi;

PROCEDURE GetGrad(tc: T; j: CARDINAL; grad: Vector) =
BEGIN 
  <* ASSERT tc.started *>
  <* ASSERT tc.active *>
  <* ASSERT j = 0 *>
  grad[tc.k] := 1.0D0
END GetGrad;

PROCEDURE DirDerivative(tc: T; j: CARDINAL; u: Vector): LONGREAL =
BEGIN 
  <* ASSERT tc.started *>
  <* ASSERT tc.active *>
  <* ASSERT j = 0 *>
  RETURN u[tc.k]
END DirDerivative;

PROCEDURE GetPsi(tc: T; j: CARDINAL): LONGREAL =
BEGIN 
  <* ASSERT tc.started *>
  <* ASSERT tc.active *>
  <* ASSERT j = 0 *>
  RETURN tc.psis
END GetPsi;

PROCEDURE SetPsi(tc: T; j: CARDINAL; psi: LONGREAL) =
BEGIN 
  <* ASSERT tc.started *>
  <* ASSERT tc.active *>
  <* ASSERT j = 0 *>
  tc.psis := psi
END SetPsi;

PROCEDURE AddReactionForce(tc: T; j: CARDINAL; lambda: LONGREAL; f: Vector) =
BEGIN 
  <* ASSERT tc.started *>
  <* ASSERT tc.active *>
  <* ASSERT j = 0 *>
  f[tc.k] := f[tc.k] - lambda 
END AddReactionForce;

BEGIN END KinematicConstraint.
