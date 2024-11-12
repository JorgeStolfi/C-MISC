MODULE JointConstraint;

FROM Constraint IMPORT Vector;

REVEAL
  T = Public BRANDED OBJECT
      n, n1, n2, n3: CARDINAL;
      a, b, c: LONGREAL;
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

PROCEDURE Init(q: T; n, n1, n2, n3: CARDINAL; a, b, c: LONGREAL) =
BEGIN
  q.n := 3*n; q.n1 := 3*n1; q.n2 := 3*n2; q.n3 := 3*n3;
  q.a := a; q.b := b; q.c := c;
END Init;

PROCEDURE NumEquations(<*UNUSED*> q: T): CARDINAL =
BEGIN 
  RETURN 3
END NumEquations;

PROCEDURE ComputeGrad(<*UNUSED*> q: T; <*UNUSED*> pos: Vector;
                      <*UNUSED*> t: LONGREAL) =
BEGIN END ComputeGrad;
 
PROCEDURE ComputePsi(<*UNUSED*> q: T; <*UNUSED*> pos, vel: Vector;
                     <*UNUSED*> t: LONGREAL) =
BEGIN END ComputePsi;
 
PROCEDURE GetGrad(q: T; j: CARDINAL; grad: Vector) =
BEGIN
  grad[q.n +j] := 1.0D0;
  grad[q.n1+j] := -q.a;
  grad[q.n2+j] := -q.b;
  grad[q.n3+j] := -q.c;
END GetGrad;

PROCEDURE DirDerivative(q: T; j: CARDINAL; u: Vector): LONGREAL =
BEGIN
  RETURN u[q.n+j] - q.a*u[q.n1+j] - q.b*u[q.n2+j] - q.c*u[q.n3+j]
END DirDerivative;

PROCEDURE GetPsi(<*UNUSED*> q: T; <*UNUSED*> j: CARDINAL): LONGREAL =
BEGIN RETURN 0.0D0 END GetPsi;

PROCEDURE SetPsi(<*UNUSED*> q: T; <*UNUSED*> j: CARDINAL; 
                 <*UNUSED*> psi: LONGREAL) =
BEGIN END SetPsi;

PROCEDURE AddReactionForce(q: T; j: CARDINAL; lambda: LONGREAL; f: Vector) =
BEGIN
  f[q.n +j] := f[q.n +j] -     lambda;
  f[q.n1+j] := f[q.n1+j] + q.a*lambda;
  f[q.n2+j] := f[q.n2+j] + q.b*lambda;
  f[q.n3+j] := f[q.n3+j] + q.c*lambda;
END AddReactionForce;


(*--- ESPECIFIC JOINTS ------------------------------------------------------*)

PROCEDURE PointToPoint(n, n1: CARDINAL): T =
BEGIN
  WITH j = NEW(T) DO
    j.init(n, n1, 0, 0, 1.0D0, 0.0D0, 0.0D0);
    RETURN j;
  END
END PointToPoint;

PROCEDURE PointToLine(n, n1, n2: CARDINAL; a, b: LONGREAL): T =
BEGIN
  WITH j = NEW(T) DO
    j.init(n, n1, n2, 0, a, b, 0.0D0);
    RETURN j;
  END
END PointToLine;

PROCEDURE PointToFace(n, n1, n2, n3: CARDINAL; a, b, c: LONGREAL): T =
BEGIN
  WITH j = NEW(T) DO
    j.init(n, n1, n2, n3, a, b, c);
    RETURN j;
  END
END PointToFace;

BEGIN END JointConstraint.
