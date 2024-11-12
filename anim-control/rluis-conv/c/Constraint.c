
#include <Constraint.h>



REVEAL
  T == Public BRANDED OBJECT
    OVERRIDES
      start = Start;
    }
    
void Start(T c, <*UNUSED, "??"); double t, <*UNUSED, "??"); pos, Vector vel)
{ 
  c.started = TRUE;
} /* Start */;

{; } /* Constraint */.
