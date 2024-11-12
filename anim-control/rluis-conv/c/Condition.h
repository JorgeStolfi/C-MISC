
#ifndef Condition_H
#define Condition_H



/*
  A condition that can be used to determine whether a constraint
  is applicable or !.
*/

#include <Integrator.h>

CONST NoEvent == ((double.nel - 1)|?|MAX_LONGREAL);

typedef
  Vector == Integrator.Vector;
  T <: Public;
  Public == OBJECT
      started: /*READONLY*/ bool = FALSE;
    METHODS
      start(double t, pos, Vector vel);
        /*
          Informs the condition that the simulation is starting at time 
          "t" in state "(pos, vel)". Must be called exactly once before 
          calling any of the other methods. */
    
      detectEvent(
          double t0; pos0, Vector vel0, 
          double t1; pos1, Vector vel1,
        ): double; 
        /*
          Returns the time of the earliest event in the range (t0 .. t1]
          associated with this constraint; or NoEvent if there is none. */
          
      treatEvent(double t, pos, Vector vel);
        /* 
          Asks the constraint to process an event that was scheduled by
          it (through "detectEvent") at time "t", and informs 
          that the state is currently "(pos,vel)". */;
    }
;
} /* Condition */.

#endif
