MODULE TimedConstraint;

FROM Constraint IMPORT Vector, NoEvent;

REVEAL
  T = Public BRANDED OBJECT
    OVERRIDES
      init        := Init;
      start       := Start;
      setClock    := SetClock;
      detectEvent := DetectEvent;
    END;
    
PROCEDURE Init(tc: T; ta, tb: LONGREAL) =
BEGIN
  tc.ta := ta;
  tc.tb := tb;
  tc.active := FALSE
END Init;

PROCEDURE Start(tc: T; t: LONGREAL; pos, vel: Vector) =
BEGIN 
  <* ASSERT NOT tc.started *>
  <* ASSERT NOT tc.active *>
  IF tc.ta <= t AND t < tc.tb THEN 
    (* Force constraint to begin now: *)
    tc.ta := t;
    tc.started := TRUE;
    (* Force an initial event: *)
    tc.treatEvent(t, pos, vel)
  END;
END Start;

PROCEDURE SetClock(tc: T; t: LONGREAL) =
BEGIN 
  <* ASSERT tc.started *>
  tc.active := (tc.ta <= t AND t < tc.tb)
END SetClock;

PROCEDURE DetectEvent(tc: T; 
  t0: LONGREAL; <*UNUSED*> pos0, vel0: Vector; 
  t1: LONGREAL; <*UNUSED*> pos1, vel1: Vector;
): LONGREAL =
BEGIN
  <* ASSERT tc.started *>
  IF t0 < tc.ta  AND tc.ta <= t1 THEN
    RETURN tc.ta
  ELSIF t0 < tc.tb  AND tc.tb <= t1 THEN
    RETURN tc.tb
  ELSE
    RETURN NoEvent
  END
END DetectEvent;

BEGIN END TimedConstraint.
