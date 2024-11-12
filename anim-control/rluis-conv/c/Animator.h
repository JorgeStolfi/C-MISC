
#ifndef Animator_H
#define Animator_H



#include <r3.h>
#include <SystemTopology.h>
#include <Integrator.h>
#include <SparseSymetricMatrix AS SSM.h>
#include <Constraints AS CT.h>
#include <ExternalForces AS FE.h>
#include <Integrator.h>

PROCEDURE Simulate(
    SystemTopology.T s,
    double t0, t1,        /* Simulation starting and stopping times. */
    FILE *stRd,              /* Initial state. */
    FILE *kcRd,              /* Kinetic constraints, or NULL if none. */
    FILE *efRd,              /* External forces, or NULL if none. */
    FILE *mmRd,              /* Input factored mass matrix (NUL L == compute it). */
    double fps,           /* Frames per second, or 0 if none. */
    double dtMin, dtMax,  /* Min and max integration step. */
    double pTol, vTol,    /* Integration error tolerance per coordinate. */ 
    FILE *mmWr,              /* Output factored mass matrix (NUL L == don't save). */
    FILE *stWr,              /* Computed states. */
    FILE *evWr,              /* Event trace. */
    bool verbose = FALSE, 
  );

#endif
