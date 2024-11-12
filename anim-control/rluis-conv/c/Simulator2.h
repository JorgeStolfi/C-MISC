
#ifndef Simulator_H
#define Simulator_H



#include <r3.h>
#include <Integrator.h>
#include <SystemDynamics.h>
#include <SystemTopology.h>
#include <SystemIO.h>
#include <
       Constraint, Force;

typedef
  Vector == Integrator.Vector;
  Vectors3D == SystemIO.Vectors3D;
  Vertices == SystemTopology.Vertex_vec;

Vectors3D *?force;

Vectors3D NullForce(double t);
void Init();
void Go(force: SystemDynamics.UserForce = NullForce);
nat GetNumber();
r3_t GetPos(nat i);
r3_t GetVel(nat i);
Vector GetInitialPos();
Vector GetInitialVel();
double GetInitialTime();
double Getxi();
Vertices GetVertices();
void AddForce(Force.T f);
;
} /* Simulator */.

#endif
