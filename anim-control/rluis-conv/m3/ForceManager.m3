MODULE ForceManager;

IMPORT Force;
FROM Force IMPORT Vector;

REVEAL
  T = Public BRANDED OBJECT
      m: CARDINAL;
      forces: REF ARRAY OF Force.T;
    OVERRIDES
      init        := Init;
      number      := Number;
      addForce    := AddForce;
      removeForce := RemoveForce;
      compute     := Compute;
    END;

PROCEDURE Init(fm: T; m: CARDINAL) =
BEGIN
  fm.m := 0;
  Allocate(fm, m);
END Init;

PROCEDURE Allocate(fm: T; m: CARDINAL) =
BEGIN
  WITH forces = NEW(REF ARRAY OF Force.T, m) DO
    FOR i := 0 TO fm.m-1 DO forces[i] := fm.forces[i] END;
    fm.forces := forces;
  END
END Allocate;

PROCEDURE Number(fm: T): CARDINAL =
BEGIN RETURN fm.m END Number;

PROCEDURE AddForce(fm: T; force: Force.T) =
BEGIN
  IF fm.m + 1 > NUMBER(fm.forces^) THEN
    Allocate(fm, 2*fm.m)
  END;
  fm.forces[fm.m] := force;
  INC(fm.m);
END AddForce;

PROCEDURE RemoveForce(fm: T; force: Force.T) =
VAR i := 0;
BEGIN
  <* ASSERT fm.m > 0 *>
  WHILE i < fm.m AND fm.forces[i] # force DO INC(i) END;
  IF i < fm.m THEN
    DEC(fm.m);
    fm.forces[i] := fm.forces[fm.m];
    fm.forces[fm.m] := NIL;
  END
END RemoveForce;

PROCEDURE Compute(fm: T; pos, vel, f: Vector; t: LONGREAL) =
BEGIN
  FOR j := 0 TO fm.m-1 DO
    fm.forces[j].compute(pos, vel, f, t)
  END
END Compute;

BEGIN END ForceManager.
