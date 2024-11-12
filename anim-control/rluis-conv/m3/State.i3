INTERFACE State;

IMPORT Simulator;

TYPE
  Time = Simulator.Time;
  Position = Simulator.Position;
  Velocity = Simulator.Velocity;
  Acceleration = Simulator.Acceleration;

TYPE
  T <: Public;
  Public = OBJECT
      time: LONGREAL;
      pos: REF Position;
      vel: REF Velocity;
      comment: TEXT;
    METHODS
      alloc(nNodes: CARDINAL): T;
        (*
          (Re)allocates position and velocity for "nNodes" nodes. *)
          
      read(rd: Rd.T): T;
        (*
          Reads a state from "rd", (re)allocating "pos" and "vel" as needed. *)
          
      write(wr: Wr.T);
        (*
          Writes the state TO "wr". *)
          
    END;
    
  Coords = ARRAY OF LONGREAL;

PROCEDURE New(nNodes: CARDINAL): T;
  (*
    Allocates a new state for a model with "nNodes" nodes. *)

PROCEDURE Read(rd: Rd.T; VAR s: T);
  (*
    Reads positions and velocities from "rd". *)
    
PROCEDURE Write(wr: Wr.T; READONLY s: T);
  (*
    Writes positions and velocities to "wr", in a format compatible
    with "Read". *)

END State.
