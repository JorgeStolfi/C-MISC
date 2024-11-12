INTERFACE IncompressibleConstraint;

IMPORT Constraint;

TYPE
  T <: Public;
  Public = Constraint.T OBJECT METHODS
      init(na, nb, nc, nd: CARDINAL);
        (*
          Constrains the tetrahedron with vertices "na,nb,nc,nd"
          to have constant volume. *)
          
    END;

END IncompressibleConstraint.
