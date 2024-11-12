
#ifndef TimedConstraint_H
#define TimedConstraint_H



#include <Constraint.h>

typedef
  T <: Public;
    /*
      A constraint that exists only between two instants "ta" and "tb".
      The "detectEvent" method will generate events at those times.  
      If the staring time lies inside [ta _ tb), the "start" method
      will force "ta" to be "t", and then it will call "treatEvent" for "ta".
    */
      
  Public == Constraint.T OBJECT
      double ta, tb;
      active: /*READONLY*/ bool = FALSE;  /* TRUE between times "ta" and "tb". */  
    METHODS
      init(ta, double tb);  /* Initializes "ta" and "tb", sets "activ e == FALSE" */
      setClock(double t);   /* updates "active" for time "t"+. */;
    }
;
} /* TimedConstraint */.

#endif
