
#ifndef KinematicConstraint_H
#define KinematicConstraint_H



#include <TimedConstraint.h>

typedef
  T <: Public;
  Public == TimedConstraint.T OBJECT METHODS
      init(nat k, ta, double tb, p, double v, xi: double = 0.0D0);
        /*
          This constraint is active between times "ta" and "tb".
          It forces coordinate "k" to vary according 
          to a Bezier polynomial: from whatever value and derivative it has at
          time "ta", to value == "p" and derivative "v" at time "tb".
          Note: "k" is a generalized coordinate index, not a vertex NUMBER.
        */;
    }
;
} /* KinematicConstraint */.

#endif
