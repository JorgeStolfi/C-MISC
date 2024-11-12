INTERFACE Animator;

IMPORT LR3, SystemTopology, Integrator;
IMPORT SparseSymetricMatrix AS SSM, Constraints AS CT, ExternalForces AS FE;
FROM Integrator IMPORT DetectEvent, TreatEvent;

PROCEDURE Simulate(
    s: SystemTopology.T;
    t0, t1: LONGREAL;        (* Simulation starting and stopping times. *)
    stRd: Rd.T;              (* Initial state. *)
    kcRd: Rd.T;              (* Kinetic constraints, or NIL if none. *)
    efRd: Rd.T;              (* External forces, or NIL if none. *)
    mmRd: Rd.T;              (* Input factored mass matrix (NIL=compute it). *)
    fps: LONGREAL;           (* Frames per second, or 0 if none. *)
    dtMin, dtMax: LONGREAL;  (* Min and max integration step. *)
    pTol, vTol: LONGREAL;    (* Integration error tolerance per coordinate. *) 
    mmWr: Wr.T;              (* Output factored mass matrix (NIL=don't save). *)
    stWr: Wr.T;              (* Computed states. *)
    evWr: Wr.T;              (* Event trace. *)
    verbose: BOOL := FALSE; 
  );

END Animator.
