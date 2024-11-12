
#include <RKF4Integrator.h>



#include <Math.h>
#include <Integrator.h>
#include <Math.h>
#include <Integrator.h>

REVEAL T == Public BRANDED OBJECT
    v: ARRAY [1..5] OF REF Velocity;
    s: REF State;
  OVERRIDES
    name = Name;
    adjustStepSize = AdjustStepSize;
    step = Step;
  }

char *Name(<*UNUSED, "??"); T ri)
  {
    return "RKF4Integrator.T";
  } /* Name */;

PROCEDURE AdjustStepSize(
    <*UNUSED, "??"); T g, 
    Time dt; 
    Coord error, tol; 
    Time dtMin, dtMax;
  ): Time == 
  CONST MaxScale == 10.0;
  CONST MinScale == 1.0/MaxScale;
  CONST Safety == 0.840896415253;
  {
    if ((error > 0.0 )) {
      { /* with*/ scale == Safety * sqrt(sqrt(tol/error)) ) {
        if ((scale <= MinScale  )) { dt = MinScale * dt
        } else if ((scale >= MaxScale  )) { dt = MaxScale * dt
        } else { dt = scale * dt;
        };
      }
    } else {
      dt = MaxScale * dt;
    }
    return min(dtMax, max(dtMin, dt));
  } /* AdjustStepSize */;

PROCEDURE Step(
    T g;
    DiffProc diff;
    Time ta;              /* Starting time */
    State *sa;    /* Starting state */
    Velocity *va; /* Derivatives at starting state */
    Time tb;              /* Final time */
    State *?sb;         /* OUT: Final state */
    Error *?er;         /* OUT: Error estimate */
  ) == 
/* Raises Abort */
  CONST
    C1a == 1.0/4.0;
    C1t == C1a;
    C2a == 3.0/32.0;
    C21 == 9.0/32.0;
    C2t == C2a + C21;
    C3a == 1932.0/2197.0;
    C31 == -7200.0/2197.0;
    C32 == 7296.0 /2197.0;
    C3t == C3a + C31 + C32;
    C4a == 439.0 /216.0;
    C41 == -8.0;
    C42 == 3680.0 /513.0;
    C43 == -845.0 /4104.0;
    C4t == C4a + C41 + C42 + C43;
    C5a == -8.0 /27.0;
    C51 == 2.0;
    C52 == -3544.0 /2565.0;
    C53 == 1859.0 /4104.0;
    C54 == -11.0 /40.0;
    C5t == C5a + C51 + C52 + C53 + C54;
    Cea == 1.0 /360.0;
    Ce2 == -128.0 /4275.0;
    Ce3 == -2197.0 /75240.0;
    Ce4 == 1.0 /50.0;
    Ce5 == 2.0 /55.0;
    Cba == 25.0 /216.0;
    Cb2 == 1408.0 /2565.0;
    Cb3 == 2197.0 /4104.0;
    Cb4 == -1.0 /5.0;
  {
    /* Allocate work areas if needed: */
    { /* with*/ n == sa.nel ) { 
      if ((g.s == NULL) || (NUMBER(g.s^) != n )) { 
        g.s = NEW(REF State, n);
        for (i = 1;  i <= 5;  i++) {
          g.v[i] = NEW(REF Velocity, n);
        };
      };
    }
  
    { /* with*/ 
      n == sa.nel,
      dt == tb - ta,
      s == g.s^,
      v1 == g.v[1]^,
      v2 == g.v[2]^,
      v3 == g.v[3]^,
      v4 == g.v[4]^,
      v5 == g.v[5]^
    ) {
      /* 
        We could do without the internal work vector "g.s" 
        by using "sb" instead.  However, if the client 
        gave us the same vector for "sa" and "sb", the 
        result would be garbage.  
        
        We could also use "v5" instead of "s" in all steps
        excet the degree 5 extrapolation, where we could use "v1". 
        But that would be too obscure...
        
        It is safer this way...
      */
      
      /* Degree 1 extrapolation: */
      for (i = 0;  i < n; i++) {
        s[i] = sa[i] + dt * C1a*va[i];
      }
      diff(ta + C1t * dt, s, v1);

      /* Degree 2 extrapolation: */
      for (i = 0;  i < n; i++) {
        s[i] = sa[i] + dt * (C2a*va[i] + C21*v1[i]);
      }
      diff(ta + C2t * dt, s, v2);

      /* Degree 3 extrapolation: */
      for (i = 0;  i < n; i++) {
        s[i] = sa[i] + dt * (C3a*va[i] + C31*v1[i] + C32*v2[i]);
      }
      diff(ta + C3t * dt, s, v3);
      
      /* Degree 4 extrapolation: */
      for (i = 0;  i < n; i++) {
        s[i] = sa[i] + dt * (C4a*va[i] + C41*v1[i] + C42*v2[i] + C43*v3[i]);
      }
      diff(ta + C4t * dt, s, v4);
      
      /* Degree 5 extrapolation: */
      for (i = 0;  i < n; i++) {
        s[i] = sa[i] + dt*(C5a*va[i] + C51*v1[i] + C52*v2[i] + C53*v3[i] + C54*v4[i]);
      }
      diff(ta + C5t * dt, s, v5);
      
      /* Final extrapolation and error estimate: */
      for (i = 0;  i < n; i++) {
        sb[i] = sa[i] + dt*(Cba*va[i] + Cb2*v2[i] + Cb3*v3[i] + Cb4*v4[i]);
        er[i] = Cea*va[i] + Ce2*v2[i] + Ce3*v3[i] + Ce4*v4[i] + Ce5*v5[i];
      }
;
    };
  } /* Step */;

{ ;
} /* RKF4Integrator */.
 
