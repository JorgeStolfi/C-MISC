
#include <EventManager.h>



typedef
  JointManager == T OBJECT
      em: T_vec;
    METHODS
      detect = Detect;
      advance = Advance;
    }
    
T Join(T_vec em)
  {
    { /* with*/ 
      t == NEW(JointManager),
      NE == em.nel,
      re == T_vec_new(NE)
    ) {
      re^ = em;
      t.em = re;
      return t;
    };
  } /* Join */;
  
PROCEDURE Detect(
    t: T
    double t0;
    Position *p0;
    Velocity *v0;
    double t1;
    Position *p1;
    Velocity *v1;
  ): Time == 
  Time *?te = ((Time.nel - 1)|?|MAX_Time);
  {
    { /* with*/ 
      em == t.em^,
      NE == em.nel
    ) {
      for (i = 0;  i < NE; i++) {
        if ((em[i] != NULL )) {
          te = min(te, em[i].detect(t0, p0, v0, t1, p1, v1);
        };
      };
    };
  } /* Detect */;

PROCEDURE Advance(
    T t;
    double t0;
    Position *p0;
    double t1;
    Position *p1;
    Velocity *v;
  ) == 
  {
    { /* with*/ 
      em == t.em^,
      NE == em.nel
    ) {
      for (i = 0;  i < NE; i++) {
        if ((em[i] != NULL )) {
          em[i].advance(t0, p0, t1, p1, v);
        };
      };
    };
  } /* Advance */;

{;
} /* EventManager */.
