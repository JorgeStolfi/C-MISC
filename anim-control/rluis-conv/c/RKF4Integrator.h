
#ifndef RKF4Integrator_H
#define RKF4Integrator_H



/* A 4th-order Runge-Kutta integrator for ordinary differential equations. */
/* Created 95/06 by R.L.Liesenfeld, C.H.Q.Forster, and J.Stolfi.           */

#include <Integrator.h>

typedef
  T <: Public;
    /*
      An adaptive Runge-Kutta integrator for ordinary differential 
      equations.
      
      The "step" method performs one step of the Runge-Kutta-Fehlberg
      method, which is theoretically accurate to fourth order.  The
      error estimate is the difference between the fourth- and
      fifth-order Runge-Kutta approximators. */
    
  Public == Integrator.T OBJECT METHODS;
    }
;
} /* RKF4Integrator */.



#endif
