INTERFACE SystemDynamics;

IMPORT LR3, SystemTopology, Integrator;
IMPORT SparseSymetricMatrix AS SSM, Constraints AS CT, ExternalForces AS FE;
FROM Integrator IMPORT DetectEvent, TreatEvent;

TYPE
  T <: Public;
  Public = OBJECT
    METHODS

      init(
        mmRd: Rd.T := NIL;
        mmWr: Wr.T := NIL;
        verbose: BOOL := FALSE; 
      ): T;
      (*
        Computes the factored mass matrix and other tables,
        for use by the following methods. 
        
        If "mmRd" is non-NIL, its must contain a pre-computed 
        and pre-factored mass matrix.  If "mmWr" is non-NIL, 
        the matrix is written into it. *) 
        
      computeAcceleration(
          READONLY p: Position;    (* Generalized coordinates of system *)
          READONLY v: Velocity;    (* Time derivative of position. *)
          VAR f: Force;            (* IN: External forces; OUT: total forces. *)
          VAR J: ConstraintMatrix; (* Constraint gradients. *)
          VAR psi: Acceleration;   (* Constraint-relative acceleration. *)
          VAR D: DirectionMatrix;  (* Direction vectors OF reaction forces. *)
          VAR a: Acceleration;     (* OUT: Second derivatives of position. *)
          VAR h: Force;            (* OUT: Magnitude of reaction forces. *)
        );
        (*
          Computes the acceleration "a" for a state "(p,v)",
          given the external force vector "f".
          
          The acceleration is the result of "f", plus internal elastic
          and viscous forces, plus the ``reaction'' forces needed to
          satisfy kinematic constraints (such as joints between
          objects, or vertices being dragged along prescribed paths).
          
          The kinematic constraints must be presented in the form of a
          linear system of equations "J a = psi", where "J" is a
          matrix and "psi" a vector.  The "J" matrix must have one
          line for each independent one-dimensional constraint, and 
          one column for each system coordinate.
          
          The reaction force necessary to satisfy these constraints
          is computed internally by the method, as a linear combination
          of the columns of matrix "D".  The coefficients of this linear
          combination (the ``Lagrange multipliers'') are returned in "h".
          The number of columns of "D" must equal the number of rows
          of "J" (i.e. the number of one-dimensional constraints).
        *)
    END;

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

PROCEDURE         t: SystemTopology.T;    (* The system's model *)
(~)~

END SystemDynamics.
