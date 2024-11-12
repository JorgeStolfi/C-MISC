
#include <TimedConstraint.h>



#include <Constraint.h>

REVEAL
  T == Public BRANDED OBJECT
    OVERRIDES
      init = Init;
      start = Start;
      setClock = SetClock;
      detectEvent = DetectEvent;
    }
    
void Init(T tc, ta, double tb)
{
  tc.ta = ta;
  tc.tb = tb;
  tc.active = FALSE;
} /* Init */;

void Start(T tc, double t, pos, Vector vel)
{ 
  affirm(! tc.started , "??");
  affirm(! tc.active , "??");
  if ((tc.ta <= t) && (t < tc.tb )) { 
    /* Force constraint to begin now: */
    tc.ta = t;
    tc.started = TRUE;
    /* Force an initial event: */
    tc.treatEvent(t, pos, vel);
  };
} /* Start */;

void SetClock(T tc, double t)
{ 
  affirm(tc.started , "??");
  tc.active = (tc.ta <= t) && (t < tc.tb);
} /* SetClock */;

PROCEDURE DetectEvent(T tc, 
  double t0; <*UNUSED, "??"); pos0, Vector vel0, 
  double t1; <*UNUSED, "??"); pos1, Vector vel1,
): double == 
{
  affirm(tc.started , "??");
  if ((t0 < tc.ta ) && (tc.ta <= t1 )) {
    return tc.ta
  } else if ((t0 < tc.tb ) && (tc.tb <= t1 )) {
    return tc.tb
  } else {
    return NoEvent;
  };
} /* DetectEvent */;

{; } /* TimedConstraint */.
