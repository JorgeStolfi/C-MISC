INTERFACE Constraint;

IMPORT Integrator;

TYPE
  Vector = Integrator.Vector;
  T <: Public;
    (*
      An algebraic constraint on the system's state, possibly conditional
      on some algebraic predicate. *) 
      
  Public = OBJECT
      condition: Condition.T := NIL; (* Enforcement condition; "NIL" means always on. *)
      started: (*READONLY*) BOOLEAN := FALSE;
    METHODS
      activate(t: LONGREAL; pos, vel: Vector);
        (*
          Informs the constraint that it just started being enforced
          at time "t" in state "(pos, vel)".  This is necessary for 
          certain constraints of the form "keep X(q)" constant", 
          "drag vertex v so as to reach (p1,v1) at time t1", etc.
          Must be called at least once before calling any of the
          other methods below. *)
    
      linearize(t: LONGREAL := 0.0D0; pos, vel: Vector);
        (*
          Evaluates the differential form of the constraint 
            "g_j \dotprod {\ddot q} = \psi_j"
          at the given time and state. Saves the vector
          "g_j" and the term "\psi_j" internally, for 
          use by "*)
          
      computePsi(pos, vel: Vector; t: LONGREAL := 0.0D0);
        (*
          Computes the ``psi'' term for the constraint, saves it internally. *)
          
      getGrad(grad: Vector);
        (*
          Returns the gradient of the constraint, previoulsy computed
          by "computeGrad". *)
          
      dirDerivative(u: Vector): LONGREAL;
        (*
          Returns the directional derivative of the constraint
          relative to the vector "u", i.e. 
            "{{\partial H}\over{\partial q}}\dotprod u". *)
          
      getPsi(): LONGREAL;
        (*
          Returns the ``psi'' term for the constraint, previoulsy computed
          by "computePsi" or "setPsi". *)
          
      setPsi(psi: LONGREAL);
        (*
          Sets the ``psi'' term, to be returned by "getPsi". *)
      
      addReactionForce(lambda: LONGREAL; f: Vector);
        (*
          Adds to "f" the reaction force for this constraint, given
          the corresponding Lagrange multiplier "lambda". *)
    END;
    
END Constraint.
