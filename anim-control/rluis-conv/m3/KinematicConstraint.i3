INTERFACE KinematicConstraint;

IMPORT TimedConstraint;

TYPE
  T <: Public;
  Public = TimedConstraint.T OBJECT METHODS
      init(k: CARDINAL; ta, tb: LONGREAL; p, v: LONGREAL; xi: LONGREAL := 0.0D0);
        (*
          This constraint is active between times "ta" and "tb".
          It forces coordinate "k" to vary according 
          to a Bezier polynomial: from whatever value and derivative it has at
          time "ta", to value = "p" and derivative "v" at time "tb".
          Note: "k" is a generalized coordinate index, not a vertex NUMBER.
        *)
    END;

END KinematicConstraint.
