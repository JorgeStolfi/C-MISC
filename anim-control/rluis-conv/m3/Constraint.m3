MODULE Constraint;

REVEAL
  T = Public BRANDED OBJECT
    OVERRIDES
      start        := Start;
    END;
    
PROCEDURE Start(c: T; <*UNUSED*> t: LONGREAL; <*UNUSED*> pos, vel: Vector) =
BEGIN 
  c.started := TRUE
END Start;

BEGIN END Constraint.
