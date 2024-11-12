
#ifndef SystemDynamics_H
#define SystemDynamics_H



#include <r3.h>
#include <SystemTopology.h>
#include <Integrator.h>
#include <SparseSymetricMatrix AS SSM.h>
#include <Constraints AS CT.h>
#include <ExternalForces AS FE.h>
#include <Integrator.h>

typedef
  T <: Public;
  Public == OBJECT
    METHODS

      init(
        FILE *mmRd = NULL;
        FILE *mmWr = NULL;
        bool verbose = FALSE; 
      ): T;
      /*
        Computes the factored mass matrix and other tables,
        for use by the following methods. 
        
        If "mmRd" is non-NULL, its must contain a pre-computed 
        and pre-factored mass matrix.  If "mmWr" is non-NULL, 
        the matrix is written into it. */ 
        
      computeAcceleration(
          Position *p;    /* Generalized coordinates of system */
          Velocity *v;    /* Time derivative of position. */
          Force *f;            /* IN: External forces; OUT: total forces. */
          ConstraintMatrix *J; /* Constraint gradients. */
          Acceleration *psi;   /* Constraint-relative acceleration. */
          DirectionMatrix *D;  /* Direction vectors OF reaction forces. */
          Acceleration *a;     /* OUT: Second derivatives of position. */
          Force *h;            /* OUT: Magnitude of reaction forces. */
        );
        /*
          Computes the acceleration "a" for a state "(p,v)",
          given the external force vector "f".
          
          The acceleration is the result of "f", plus internal elastic
          and viscous forces, plus the ``reaction'' forces needed to
          satisfy kinematic constraints (such as joints between
          objects, or vertices being dragged along prescribed paths).
          
          The kinematic constraints must be presented in the form of a
          linear system of equations "J a == psi", where "J" is a
          matrix and "psi" a vector.  The "J" matrix must have one
          line for each independent one-dimensional constraint, and 
          one column for each system coordinate.
          
          The reaction force necessary to satisfy these constraints
          is computed internally by the method, as a linear combination
          of the columns of matrix "D".  The coefficients of this linear
          combination (the ``Lagrange multipliers'') are returned in "h".
          The number of columns of "D" must equal the number of rows
          of "J" (i.e. the number of one-dimensional constraints).
        */;
    }

PROCEDURE Simulate(
    SystemTopology.T s;
    double t0, t1;        /* Simulation starting and stopping times. */
    FILE *stRd;              /* Initial state. */
    FILE *kcRd;              /* Kinetic constraints, or NULL if none. */
    FILE *efRd;              /* External forces, or NULL if none. */
    FILE *mmRd;              /* Input factored mass matrix (NUL L == compute it). */
    double fps;           /* Frames per second, or 0 if none. */
    double dtMin, dtMax;  /* Min and max integration step. */
    double pTol, vTol;    /* Integration error tolerance per coordinate. */ 
    FILE *mmWr;              /* Output factored mass matrix (NUL L == don't save). */
    FILE *stWr;              /* Computed states. */
    FILE *evWr;              /* Event trace. */
    bool verbose = FALSE; 
  );

SystemTopology.T PROCEDURE         t;    /* The system's model */
(~)~
;
} /* SystemDynamics */.

#endif
