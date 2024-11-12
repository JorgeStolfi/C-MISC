INTERFACE EESpring;

IMPORT Force;

TYPE
  T <: Public;
  Public = Force.T OBJECT METHODS 
    init(a, b, c, d: CARDINAL; a1, a2, a3, a4, k: LONGREAL);
  END;

END EESpring.
