
#ifndef RigidConstraint_H
#define RigidConstraint_H



#include <Constraint.h>

typedef
  T <: Public;
  Public == Constraint.T OBJECT METHODS
      init(n1, nat n2);
        /*
          Constrains the distance between vertices n1, n2 to remain constant. */;
    }
;
} /* RigidConstraint */.

#endif
