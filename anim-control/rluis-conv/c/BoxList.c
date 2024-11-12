
#include <BoxList.h>



#include <r3.h>
#include <SystemTopology.h>
#include <SystemDynamics.h>

void BuildVertexBoxes(T bx, SystemDynamics.T s, double dt)
VAR vertices = s.getVertices();
    double min, max;
{
  dt = dt/3.0D0;
  for (i = 0;  i <= ((bx.nel - 1)|?|MAX_bx);  i++) {
    { /* with*/ 
      a == vertices[i].a,
      pa == s.getPos0(a),
      va == s.getVel0(a),
      pb == s.getPos(a),
      vb == s.getVel(a),
      qa == r3_Add(pa, r3_Scale(dt, va)),
      qb == r3_Sub(pb, r3_Scale(dt, vb))
    ) {
      if ((qb[0] > pb[0] )) { min = pb[0]; max = qb[0]
      } else { min = qb[0]; max = pb[0]; }
      if ((qa[0] < min )) { min = qa[0]
      } else if ((qa[0] > max )) { max = qa[0]; }
      if ((pa[0] < min )) { min = pa[0]
      } else if ((pa[0] > max )) { max = pa[0]; }
      bx[i].min[Axis.X] = min;
      bx[i].max[Axis.X] = max;
      
      if ((qb[1] > pb[1] )) { min = pb[1]; max = qb[0]
      } else { min = qb[1]; max = pb[1]; }
      if ((qa[1] < min )) { min = qa[1]
      } else if ((qa[1] > max )) { max = qa[1]; }
      if ((pa[1] < min )) { min = pa[1]
      } else if ((pa[1] > max )) { max = pa[1]; }
      bx[i].min[Axis.Y] = min;
      bx[i].max[Axis.Y] = max;
      
      if ((qb[2] > pb[2] )) { min = pb[2]; max = qb[2]
      } else { min = qb[2]; max = pb[2]; }
      if ((qa[2] < min )) { min = qa[2]
      } else if ((qa[2] > max )) { max = qa[2]; }
      if ((pa[2] < min )) { min = pa[2]
      } else if ((pa[2] > max )) { max = pa[2]; }
      bx[i].min[Axis.Z] = min;
      bx[i].max[Axis.Z] = max;
    };
  };
} /* BuildVertexBoxes */;

void BuildEdgeBoxes(bx, T vertexBoxes, SystemTopology.T s)
VAR edges = s.edgeToVertex;
{
  for (i = 0;  i <= ((bx.nel - 1)|?|MAX_bx);  i++) {
    { /* with*/ 
      a == vertexBoxes[edges[i].a],
      b == vertexBoxes[edges[i].b]
    ) {
      bx[i].min[Axis.X] = min(a.min[Axis.X], b.min[Axis.X]);
      bx[i].max[Axis.X] = max(a.max[Axis.X], b.max[Axis.X]);
      bx[i].min[Axis.Y] = min(a.min[Axis.Y], b.min[Axis.Y]);
      bx[i].max[Axis.Y] = max(a.max[Axis.Y], b.max[Axis.Y]);
      bx[i].min[Axis.Z] = min(a.min[Axis.Z], b.min[Axis.Z]);
      bx[i].max[Axis.Z] = max(a.max[Axis.Z], b.max[Axis.Z]);
    };
  };
} /* BuildEdgeBoxes */;

void BuildFaceBoxes(bx, T edgeBoxes, SystemTopology.T s)
VAR faces = s.faceToEdge;
{
  for (i = 0;  i <= ((bx.nel - 1)|?|MAX_bx);  i++) {
    { /* with*/ 
      a == edgeBoxes[faces[i].a],
      b == edgeBoxes[faces[i].b]
    ) {
      bx[i].min[Axis.X] = min(a.min[Axis.X], b.min[Axis.X]);
      bx[i].max[Axis.X] = max(a.max[Axis.X], b.max[Axis.X]);
      bx[i].min[Axis.Y] = min(a.min[Axis.Y], b.min[Axis.Y]);
      bx[i].max[Axis.Y] = max(a.max[Axis.Y], b.max[Axis.Y]);
      bx[i].min[Axis.Z] = min(a.min[Axis.Z], b.min[Axis.Z]);
      bx[i].max[Axis.Z] = max(a.max[Axis.Z], b.max[Axis.Z]);
    };
  };
} /* BuildFaceBoxes */;

{; } /* BoxList */.
