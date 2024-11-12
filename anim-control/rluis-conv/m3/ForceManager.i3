INTERFACE ForceManager;

IMPORT Force;

TYPE
  T <: Public;
  Public = OBJECT METHODS
      init(m: CARDINAL);
      number(): CARDINAL;
      addForce(force: Force.T);
      removeForce(force: Force.T);
      compute(pos, vel, f: Force.Vector; t: LONGREAL);
    END;

END ForceManager.
 
