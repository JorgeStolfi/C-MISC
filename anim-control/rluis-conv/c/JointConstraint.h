
#ifndef JointConstraint_H
#define JointConstraint_H



#include <Constraint.h>

typedef
  T <: Public;
  Public == Constraint.T OBJECT METHODS
      init(n, n1, n2, nat n3, a, b, double c);
    }

T PointToPoint(n, nat n1);
  /*
    Constrains vertices n and n1 to coincide */
    
T PointToLine(n, n1, nat n2, a, double b);
  /*
    Constrains vertex n to lie along the edge (n1,n2),
    at barycentric coordinates (a,b). */
    
T PointToFace(n, n1, n2, nat n3, a, b, double c);
  /*
    Constrains vertex n to lie on the face (n1,n2,n3),
    at barycentric coordinates (a,b,c). */
;
} /* JointConstraint */.

#endif
