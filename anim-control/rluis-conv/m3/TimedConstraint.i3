INTERFACE TimedConstraint;

IMPORT Constraint;

TYPE
  T <: Public;
    (*
      A constraint that exists only between two instants "ta" and "tb".
      The "detectEvent" method will generate events at those times.  
      If the staring time lies inside [ta _ tb), the "start" method
      will force "ta" to be "t", and then it will call "treatEvent" for "ta".
    *)
      
  Public = Constraint.T OBJECT
      ta, tb: LONGREAL;
      active: (*READONLY*) BOOLEAN := FALSE;  (* TRUE between times "ta" and "tb". *)  
    METHODS
      init(ta, tb: LONGREAL);  (* Initializes "ta" and "tb", sets "active=FALSE" *)
      setClock(t: LONGREAL);   (* updates "active" for time "t"+. *)
    END;

END TimedConstraint.
