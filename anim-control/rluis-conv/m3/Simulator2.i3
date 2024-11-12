INTERFACE Simulator;

IMPORT LR3, Integrator, SystemDynamics, SystemTopology, SystemIO, 
       Constraint, Force;

TYPE
  Vector = Integrator.Vector;
  Vectors3D = SystemIO.Vectors3D;
  Vertices = REF ARRAY OF SystemTopology.Vertex;

VAR force: Vectors3D;

PROCEDURE NullForce(t: LONGREAL): Vectors3D;
PROCEDURE Init();
PROCEDURE Go(force: SystemDynamics.UserForce := NullForce);
PROCEDURE GetNumber(): CARDINAL;
PROCEDURE GetPos(i: CARDINAL): LR3.T;
PROCEDURE GetVel(i: CARDINAL): LR3.T;
PROCEDURE GetInitialPos(): Vector;
PROCEDURE GetInitialVel(): Vector;
PROCEDURE GetInitialTime(): LONGREAL;
PROCEDURE Getxi(): LONGREAL;
PROCEDURE GetVertices(): Vertices;
PROCEDURE AddForce(f: Force.T);

END Simulator.
