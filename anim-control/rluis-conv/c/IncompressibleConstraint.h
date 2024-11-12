
#ifndef IncompressibleConstraint_H
#define IncompressibleConstraint_H



#include <Constraint.h>

typedef
  T <: Public;
  Public == Constraint.T OBJECT METHODS
      init(na, nb, nc, nat nd);
        /*
          Constrains the tetrahedron with vertices "na,nb,nc,nd"
          to have constant volume. */
          ;
    }
;
} /* IncompressibleConstraint */.

#endif
