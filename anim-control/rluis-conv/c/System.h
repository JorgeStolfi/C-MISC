
#ifndef System_H
#define System_H



#include <SystemDynamics.h>

typedef
  T <: Public;
  Public == SystemDynamics.T OBJECT METHODS
      init(
          nat fps; t1, dtmin, dtmax, tol, g, double spring,
          SystemTopology.FileNames *names; 
          print, verbose, collide: bool
        );
      run(char *dest, SystemDynamics.UserForce userForce);
    }
;
} /* System */.

#endif
