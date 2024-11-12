INTERFACE Condition;

(*
  A condition that can be used to determine whether a constraint
  is applicable or NOT.
*)

IMPORT Integrator;

CONST NoEvent = LAST(LONGREAL);

TYPE
  Vector = Integrator.Vector;
  T <: Public;
  Public = OBJECT
      started: (*READONLY*) BOOLEAN := FALSE;
    METHODS
      start(t: LONGREAL; pos, vel: Vector);
        (*
          Informs the condition that the simulation is starting at time 
          "t" in state "(pos, vel)". Must be called exactly once before 
          calling any of the other methods. *)
    
      detectEvent(
          t0: LONGREAL; pos0, vel0: Vector; 
          t1: LONGREAL; pos1, vel1: Vector;
        ): LONGREAL; 
        (*
          Returns the time of the earliest event in the range (t0 .. t1]
          associated with this constraint; or NoEvent if there is none. *)
          
      treatEvent(t: LONGREAL; pos, vel: Vector);
        (* 
          Asks the constraint to process an event that was scheduled by
          it (through "detectEvent") at time "t", and informs 
          that the state is currently "(pos,vel)". *)
    END;

END Condition.
