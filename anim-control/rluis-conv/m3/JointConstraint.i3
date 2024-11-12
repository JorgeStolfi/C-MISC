INTERFACE JointConstraint;

IMPORT Constraint;

TYPE
  T <: Public;
  Public = Constraint.T OBJECT METHODS
      init(n, n1, n2, n3: CARDINAL; a, b, c: LONGREAL);
    END;

PROCEDURE PointToPoint(n, n1: CARDINAL): T;
  (*
    Constrains vertices n and n1 to coincide *)
    
PROCEDURE PointToLine(n, n1, n2: CARDINAL; a, b: LONGREAL): T;
  (*
    Constrains vertex n to lie along the edge (n1,n2),
    at barycentric coordinates (a,b). *)
    
PROCEDURE PointToFace(n, n1, n2, n3: CARDINAL; a, b, c: LONGREAL): T;
  (*
    Constrains vertex n to lie on the face (n1,n2,n3),
    at barycentric coordinates (a,b,c). *)

END JointConstraint.
