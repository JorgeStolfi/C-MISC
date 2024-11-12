
#include <PairTable.h>



#include <SystemTopology.h>

T New(l, nat c)
{
  { /* with*/ table == NEW(T, l) ) {
    for (i = 0;  i < l; i++) {
      table[i].row = NEW(Row, c);
    }
    return table;
  };
} /* New */;

bool On(T t, nat i)
{ return t[i].on; } /* On */;

void SetOn(T t, nat i, bool on)
{ t[i].on = on; } /* SetOn */;

bool Valid(T t, i, nat j)
{ return t[i].row[j].v; } /* Valid */;

void SetValid(T t, i, nat j, bool v)
{ t[i].row[j].v = v; } /* SetValid */;

bool Overlap(T t, i, nat j)
{
  { /* with*/ flag == t[i].row[j].flag ) {
    return flag[Axis.X]) && (flag[Axis.Y]) && (flag[Axis.Z];
  };
} /* Overlap */;

nat BuildVertexFacePairs(T t, SystemTopology.T s)
VAR n = 0;
{
  for (i = 0;  i <= LAST(s.vertices^);  i++) {
    { /* with*/ 
      d == s.vertices[i].a DIV 3,
      vertexFixed == s.fixed[d]
    ) {
      t[i].on = TRUE;
      for (j = 0;  j <= LAST(s.faces^);  j++) {
        { /* with*/ 
          a == s.faces[j].a DIV 3,
          b == s.faces[j].b DIV 3,
          c == s.faces[j].c DIV 3,
          faceFixed == s.fixed[a]) && (s.fixed[b]) && (s.fixed[c],
          pair == t[i].row[j]
        ) {
          pair.v = ! (s.sameTetrahedron(a, b, c, d) OR
                         vertexFixed) && (faceFixed);
          pair.flag[Axis.X] = FALSE;
          pair.flag[Axis.Y] = FALSE;
          pair.flag[Axis.Z] = FALSE;
          if ((pair.v )) { n++;
        };
      };
    };
  }
  return n;
} /* BuildVertexFacePairs */;

nat BuildEdgeEdgePairs(T t, SystemTopology.T s)
VAR n = 0;
{
  for (i = 0;  i <= LAST(s.edges^);  i++) {
    t[i].on = TRUE;
    t[i].row[i].v = FALSE;
    { /* with*/ 
      a == s.edges[i].a DIV 3,
      b == s.edges[i].b DIV 3,
      edge1Fixed == s.fixed[a]) && (s.fixed[b]
    ) {
      for (j = i+1;  j <= LAST(s.edges^);  j++) {
        { /* with*/ 
          c == s.edges[j].a DIV 3,
          d == s.edges[j].b DIV 3,
          edge2Fixed == s.fixed[c]) && (s.fixed[d],
          pair == t[i].row[j]
        ) {
          pair.v = ! (a == c) || (a == d) || (b == c) || (b == d OR
                         s.sameTetrahedron(a, b, c, d) OR
                         edge1Fixed) && (edge2Fixed);
          pair.flag[Axis.X] = FALSE;
          pair.flag[Axis.Y] = FALSE;
          pair.flag[Axis.Z] = FALSE;
          if ((pair.v )) { n++;
        };
      };
    };
  }
  return n;
} /* BuildEdgeEdgePairs */;

{; } /* PairTable */.
