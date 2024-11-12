
#ifndef ForceManager_H
#define ForceManager_H



#include <Force.h>

typedef
  T <: Public;
  Public == OBJECT METHODS
      init(nat m);
      number(): nat;
      addForce(Force.T force);
      removeForce(Force.T force);
      compute(pos, vel, Force.Vector f, double t);
    }
;
} /* ForceManager */.
 

#endif
