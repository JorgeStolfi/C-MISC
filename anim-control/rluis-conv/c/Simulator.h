
#ifndef Simulator_H
#define Simulator_H



/* A second-order differential equation solver for dynamical simulations. */
/* Created 95/06 by R.L.W.Liesenfeld and J.Stolfi.                        */

/*
  This module implements a dynamical simulator,
  i.e. a numerical integrator for a system of
  second order ordinary differential equations "p'' == F(p, p',t)"
  with support for detection and treatment of client-defined
  discrete events. 

  Nomenclature: the vector "p" is called the system's {\em position vector}
  its time derivative "p'" is the {\em velocity vector} and the second
  derivative "p''" is the {\em acceleration vector}.  The pait "(p,p')"
  is the {\em state} of the system. 

  The current implementation is based on the Runge-Kutta-Fehlberg 
  fourth-order integrator "RKF4integrator.T".
*/
       
#include <Integrator.h>

typedef
  Time == Integrator.Time;
  Coord == double;
  Speed == double;
  Accel == double;
  Position == ARRAY OF Coord;
  Velocity == ARRAY OF Speed;
  Acceleration == ARRAY OF Accel;
  
PROCEDURE Simulate(
    double *?t;      /* IN: starting time; OUT: final time */
    Position *?p;      /* IN: initial position; OUT: final position */
    Velocity *?v;      /* IN: initial velocity; OUT: final velocity */
    Time dtMin, dtMax;   /* Min and max time steps */
    Coord pTol;          /* Max integration error in "p[i]" allowed per step */
    Speed vTol;          /* Max integration error in "v[i]" allowed per step */
    EvalRHSProc evalRHS;            /* Computes acceleration from time+state. */
    CheckStepProc checkStep = NULL; /* Checks and handles discrete events. */
  );
  /*
    Integrates the system, starting from time "t" and state
    "(p, v)".

    The "integrate" method will call "evalRHS" several times to
    compute the right-hand side "F(p,p',t)".  Beware that the
    states "(p,p')" submitted to "evalRHS" may lie far from the
    final integrated path.

    The "integrate" method will also call "checkStep" to detect
    possible events within each integration step, and to know when TO
    stop. See the details below. */
;
    }

typedef
  EvalRHSProc == PROCEDURE (
      double t; 
      Position *p;
      Velocity *v; 
      Acceleration *?a;
    );
    /*
      Should return in "a" the acceleration for time "t", 
      given the corresponding position "p" and velocity "v". */
  
  CheckStepProc == PROCEDURE (
      Time t0; 
      Position *p0; 
      Velocity *v0;
      Acceleration *a0;
      Time t1; 
      Position *?p1; 
      Velocity *?v1;
      Acceleration *?a1;
    ): Time;
    /*
      A procedure called by the simulator to determine the possible
      occurrence of an event in the time interval "[t0 _ t1]", given
      the integrated states and accelerations at those times.  The
      procedure should return:
      
            (0) the time "t0" itself, meaning that the integration
                should stop at time "t0" with state "(p0, v0)"; or
            
            (1) any time "te" in "[t0+dtMin _ t1)", meaning that the
                integration should be redone from "t0" to "te"; or
            
            (2) the time "t1" itself, meaning that the step should be
                taken, and the integration should then stop at time
                "t1", state "(p1, v1)"; or
                 
            (3) any time "tn" greater than "t1", meaning that the step
                should be taken, and the next step should be from "t1"
                to "tb" (or shorter).
      
      In case (1) the integration will be redone from "t0" to the
      requested time "te" and then the "checkStep" procedure will be
      called again---to inspect this truncated step, and perhaps
      refine the estimate of the event time "te".  The client must
      make sure that this loop will terminate in a reasonable amount
      of time.  In particular, in case (1) the difference "t1 - te"
      must have a positive lower bound, at least in the average sense.

      In cases (2) and (3), the "checkStep" procedure is allowed
      adjust the final state "(p1,v1)" and the acceleration "a1"
      before returning to the simulator.
      
      For instance, the procedure may adjust "p1" and "v1" so as to
      satisy any constraints or conservation laws which are implied by
      the original differential equations but could be violated in a
      long integration due to the accumulation of error.
      
      Indeed, in case (3), the "checkStep" routine *must* ensure that
      the differental equations can be integrated continuously from
      state "(p1,v1)" for some positive interval of time.  In
      particular, if "evalRHS" has a discontinuity at some time "td",
      the "checkStep" procedure must ask the simulator (through exits
      (1) and (3)) to integrate normally up to "td - epsilon", and
      then perform an integration step from "td - epsilon" to "te +
      epsilon".  When processing that last step, "checkStep" should
      handle the discontinuity, and put the system in a state from
      which integration may proceeed normally.
      
      After processing a discontinuity, "checkStep" should estimate an
      appropriate time increment "dt" for the next step---typically,
      "dtMin"---and return "t1+dt".  In the normal case, when there
      were no discontinuities, "checkStep" may just return
      "((Time.nel - 1)|?|MAX_Time)"; the simulator will automatically select the next
      time step, based on an internal estimate of the integration
      error made in the previous step. */
;
} /* Simulator */.

#endif
