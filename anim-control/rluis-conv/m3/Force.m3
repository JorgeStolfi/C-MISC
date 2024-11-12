MODULE Force;
 
REVEAL
  T = Public BRANDED OBJECT OVERRIDES
      compute := Compute;
    END;

PROCEDURE Compute(<*UNUSED*> fo: T;
                  <*UNUSED*> pos, vel, f: Vector;
                  <*UNUSED*> t: LONGREAL) =
BEGIN END Compute;

BEGIN END Force.
