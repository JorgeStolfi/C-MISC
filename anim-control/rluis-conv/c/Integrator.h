
#ifndef Integrator_H
#define Integrator_H



/*
  This interface defines the abstract data type "Integrator.T",
  an object that numerically solves ordinary differential equations
  of the form 
  |
  |    ds/dt == F(t,s)
  |
  where "s" is an unknown function from time "t" to some 
  ``state space'' "R^n", and "F" is a client-provided function.
  
  Created 95/06 by R.L.Liesenfeld, C.H.Q.Forster, and J.Stolfi
*/

EXCEPTION Abort;

typedef
  Time == double;            /* A time value */
  Coord == double;           /* A state coordinate */
  State == ARRAY OF Coord;     /* The system's state. */
  Dist == double;           /* Distance between two states */
  Error == ARRAY OF Coord;     /* Error estimate for a "State". */
  Speed == double;           /* Derivative of "Coord" relative to "Time" */
  Velocity == ARRAY OF Speed;  /* Derivative of "State" relative to "Time". */

  T == OBJECT METHODS
    /*
      A generic integrator for ordinary differential equations.
      
      This is an abstract object class that implements no methods.
      Therefore, clients should not create instances of "Integrator.T"
      directly.  Instead, they should create instances of appropriate
      concrete subclasses, such as "RKF4Integrator.T". */

      step(
          DiffProc diff;        /* Computes the right-hand side "F(t,s)" */
          Time ta;              /* Starting time */
          State *sa;    /* Starting state */
          Velocity *va; /* Derivatives at starting state */
          Time tb;              /* Final time */
          State *?sb;         /* OUT: Final state */
          Error *?eb;         /* OUT: Error estimate for "sb" */
        ); 
/* Raises Abort */
        /*
          Integrates the equation "ds/dt == F(t,s)" from time "ta"
          and state "sa" to time "tb", resulting in state "sb".
          
          The method must be given in "va" the velocity (derivative of
          state with respect to time) "F(ta,sa)", already computed.
          The method will typically call the "diff" procedure in order
          to evaluate "F(ti,si)" for additional times "ti" in the open
          range "(ta__tb)", and the corresponding states "si",
          generated internally.
          
          The method will return in "eb" an estimate of the 
          absolute integration error affecting the computed state
          "sb".  This estimate is usually computed by comparing
          two different extrapolations of "sb", and may be used
          to select the proper step size "tb - ta".
          
          The vectors "sa", "va", "sb", and "er" must be 
          memory-disjoint.
          
          The "Abort" exception is not raised by "step"; is provided
          so that "diff" may interrupt the integration at mid-step.
          
          Fixed-step integration
          ----------------------
          
          In most integrators, the "step" method assumes implicitly
          that the function "s(t)" belongs to some restricted space "S", 
          for instance polynomials of a certain small degree.
          When that assumption is not true,  the integration error
          affecting the final state "sb" may grow explosively as 
          "tb-ta" increases.  
          
          Therefore, integration over ``large'' intervals (where the
          function "s(t)" cannot be approximated by functions from
          "S") is usually performed as a sequence of smaller "step"s,
          each one starting at the time and state where the previous
          one ended.
          
          In the simplest ``fixed stepsize'' scheme, the desired range
          "[ta__tb]" is divided into "NSteps" equal intervals. 
          Here is a typical example:

            | { /* with*/ 
            |   tIni == (... initial time ...),
            |   tFin == (... final time ...),
            |   s == NEW(REF State, NCoords)^,
            |   v == NEW(REF State, NCoords)^,
            |   e == NEW(REF State, NCoords)^
            | ) {
            |   VAR t: double = tIni;
            |   {
            |     s = (... initial state ...);
            |     v = (... initial velocity ...);
            |     for (k = 1;  k <= NSteps;  k++) {
            |       { /* with*/ 
            |         r == ((double)i)/((double)NSteps),
            |         tn == (1.0d0 - r)*tIni + r*tFin
            |       ) {
            |         intg.step(Diff, t,s,v, tn,s, e);
            |         t = tn;
            |         Diff(s, v); 
            |;       }
            |;     }
            |;   }
            |; }

          For a more sophisticated scheme, with variable stepsize,
          see the comments in "adjustStepSize" below. */
         
      adjustStepSize(Time dt, error, Dist tol, dtMin, Time dtMax): Time;
        /*
          Given that a time step of size "dt" generated the specified 
          integration "error", returns the time step that should 
          make the integration error approximately "tol".
          
          The result is clipped to the range "[dtMin _ dtMax]".
          
          This method assumes that "error" is some mathematical norm
          of the estimated error vector "eb" returned by "step".  That
          is, multiplying "eb" by some positive factor should multiply
          "error" by that same factor.
          
          Here is a typical example of integration with 
          ``adaptive stepsize control'':
          
            | { /* with*/ 
            |   tIni == (... initial time ...),
            |   tFin == (... final time ...),
            |   s0 == NEW(REF State, NCoords)^,
            |   s1 == NEW(REF State, NCoords)^,
            |   v0 == NEW(REF State, NCoords)^,
            |   e == NEW(REF State, NCoords)^
            | ) {
            |   VAR t0: double = tIni; dt, double t1,
            |   {
            |     s0 = (... initial state ...);
            |     v0 = (... initial velocity ...);
            |     dt = (... initial stepsize ...);
            |     while (t0 < tb ) {
            |       t1 = min(t0 + dt, tFin);
            |       while (1) {
            |         intg.step(Diff, t0,s0,v0, t1,s1, e);
            |         { /* with*/ error == LRN.Norm) {
            |           if ((error < tol) || (dt < dtMin )) { EXIT; }
            |           dt = intg.adjustStepSize(dt, error, tol, dtMin, dtMax)
            |;         }
            |;       }
            |       t0 = t1;
            |       s0 = s1;
            |       Diff(t0, s0, v0)
            |;     }
            |;   }
            |; }

          Note that the error estimates returned by "step" may be
          unreliable if the right-hand side behaves too wildly 
          in the specified interval.  Threfore, this adaptive
          scheme should be used with caution.
          
          In particular, if "F(t,s)" has discontinuities along the
          path from "(t0,s0)" to "(t1, s1)", then it is advisable to
          break the integration into two separate steps at the
          discontinuity. */
        
      name(): char *;
        /*
          Returns the integrator's name, possibly with
          any internal parameters. */;
    }

  DiffProc == PROCEDURE (Time t, READONLY State s, VAR Velocity v); 
/* Raises Abort */
    /*
      Called by the integrator to evaluate the right-hand side of the
      differential equation. */
;
} /* Integrator */.

#endif
