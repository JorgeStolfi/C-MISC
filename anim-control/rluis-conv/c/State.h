
#ifndef State_H
#define State_H



#include <Simulator.h>

typedef
  Time == Simulator.Time;
  Position == Simulator.Position;
  Velocity == Simulator.Velocity;
  Acceleration == Simulator.Acceleration;

typedef
  T <: Public;
  Public == OBJECT
      double time;
      pos: REF Position;
      vel: REF Velocity;
      char *comment;
    METHODS
      alloc(nat nNodes): T;
        /*
          (Re)allocates position and velocity for "nNodes" nodes. */
          
      read(FILE *rd): T;
        /*
          Reads a state from "rd", (re)allocating "pos" and "vel" as needed. */
          
      write(FILE *wr);
        /*
          Writes the state TO "wr". */
          ;
    }
    
  Coords == ARRAY OF double;

T New(nat nNodes);
  /*
    Allocates a new state for a model with "nNodes" nodes. */

void Read(FILE *rd, VAR T s);
  /*
    Reads positions and velocities from "rd". */
    
void Write(FILE *wr, READONLY T s);
  /*
    Writes positions and velocities to "wr", in a format compatible
    with "Read". */
;
} /* State */.

#endif
