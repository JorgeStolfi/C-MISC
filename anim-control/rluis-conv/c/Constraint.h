
#ifndef Constraint_H
#define Constraint_H



#include <Integrator.h>

typedef
  Vector == Integrator.Vector;
  T <: Public;
    /*
      An algebraic constraint on the system's state, possibly conditional
      on some algebraic predicate. */ 
      
  Public == OBJECT
      Condition.T condition = NULL; /* Enforcement condition; "NULL" means always on. */
      started: /*READONLY*/ bool = FALSE;
    METHODS
      activate(double t, pos, Vector vel);
        /*
          Informs the constraint that it just started being enforced
          at time "t" in state "(pos, vel)".  This is necessary for 
          certain constraints of the form "keep X(q)" constant", 
          "drag vertex v so as to reach (p1,v1) at time t1", etc.
          Must be called at least once before calling any of the
          other methods below. */
    
      linearize(t: double = 0.0D0; pos, Vector vel);
        /*
          Evaluates the differential form of the constraint 
            "g_j \dotprod {\ddot q} == \psi_j"
          at the given time and state. Saves the vector
          "g_j" and the term "\psi_j" internally, for 
          use by "*/
          
      computePsi(pos, Vector vel, t: double = 0.0D0);
        /*
          Computes the ``psi'' term for the constraint, saves it internally. */
          
      getGrad(Vector grad);
        /*
          Returns the gradient of the constraint, previoulsy computed
          by "computeGrad". */
          
      dirDerivative(Vector u): double;
        /*
          Returns the directional derivative of the constraint
          relative to the vector "u", i.e. 
            "{{\partial H}\(over){\partial q}}\dotprod u". */
          
      getPsi(): double;
        /*
          Returns the ``psi'' term for the constraint, previoulsy computed
          by "computePsi" or "setPsi". */
          
      setPsi(double psi);
        /*
          Sets the ``psi'' term, to be returned by "getPsi". */
      
      addReactionForce(double lambda, Vector f);
        /*
          Adds to "f" the reaction force for this constraint, given
          the corresponding Lagrange multiplier "lambda". */;
    }
    ;
} /* Constraint */.

#endif
