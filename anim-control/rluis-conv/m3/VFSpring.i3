INTERFACE VFSpring;

IMPORT Force;

TYPE
  T <: Public;
  Public = Force.T OBJECT METHODS
    init(a, b, c, d: CARDINAL; a1, a2, a3, k: LONGREAL);
  END;

END VFSpring.
