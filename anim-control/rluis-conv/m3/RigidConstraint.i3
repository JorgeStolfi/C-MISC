INTERFACE RigidConstraint;

IMPORT Constraint;

TYPE
  T <: Public;
  Public = Constraint.T OBJECT METHODS
      init(n1, n2: CARDINAL);
        (*
          Constrains the distance between vertices n1, n2 to remain constant. *)
    END;

END RigidConstraint.
