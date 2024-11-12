INTERFACE Force;

IMPORT Integrator;

TYPE
  Vector = Integrator.Vector;
  T <: Public;
  Public = OBJECT METHODS
      compute(pos, vel, f: Vector; t: LONGREAL)
    END;
    
END Force.
