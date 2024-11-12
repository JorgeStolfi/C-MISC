
#include <IntervalList.h>



#include <SystemTopology.h>
#include <BoxList.h>
#include <PairTable.h>
#include <CandidateList.h>
#include <Candidate.h>
#include <PairTable.h>

T New(SystemTopology.T s)
{
  { /* with*/ 
    nv == 2*NUMBER(s.vertices^),
    ne == 2*NUMBER(s.edges^),
    nf == 2*NUMBER(s.faces^),
    n == nv + ne + nf,
    itv == NEW(T, n)
  ) {
    for (i = 0;  i <= nv-1 BY 2;  i++) {
      itv[i+0].object = IntervalObject.Vertex;
      itv[i+1].object = IntervalObject.Vertex;
      itv[i+0].begin = TRUE;
      itv[i+1].begin = FALSE;
      itv[i+0].i = i DIV 2;
      itv[i+1].i = itv[i].i;
    }
    
    for (i = nv;  i <= nv+ne-1 BY 2;  i++) {
      itv[i+0].object = IntervalObject.Edge;
      itv[i+1].object = IntervalObject.Edge;
      itv[i+0].begin = TRUE;
      itv[i+1].begin = FALSE;
      itv[i+0].i = (i-nv) DIV 2;
      itv[i+1].i = itv[i].i;
    }
    
    for (i = nv+ne;  i <= nv+ne+nf-1 BY 2;  i++) {
      itv[i+0].object = IntervalObject.Face;
      itv[i+1].object = IntervalObject.Face;
      itv[i+0].begin = TRUE;
      itv[i+1].begin = FALSE;
      itv[i+0].i = (i-nv-ne) DIV 2;
      itv[i+1].i = itv[i].i;
    }

    return itv;
  };
} /* New */;

PROCEDURE Process(T itv, Axis axis,
                  BoxList.T vertexBox, edgeBox, faceBox;
                  PairTable.T facePairs, edgePairs;
                  CandidateList.T *?clist) == 

  double GetValue(nat i)
  {
    { /* with*/ 
      object == itv[i].object,
      begin == itv[i].begin,
      index == itv[i].i
    ) {
      if ((object == IntervalObject.Vertex )) {
        if ((begin )) { return vertexBox[index].min[axis]
        } else { return vertexBox[index].max[axis]; }
      } else if ((object == IntervalObject.Edge )) {
        if ((begin )) { return edgeBox[index].min[axis]
        } else { return edgeBox[index].max[axis]; }
      } else {  /* object == IntervalObject.Face */
        if ((begin )) { return faceBox[index].min[axis]
        } else { return faceBox[index].max[axis];
      };
    };
  } /* GetValue */;

  void SetStatus(READONLY itvi, Interval itvj)
  VAR pairs = facePairs;
      vertexFace = TRUE;
      double lmin, lmax, cmin, cmax;
      obji = itvi.object;
      objj = itvj.object;
      l = itvi.i;
      c = itvj.i;
  {
    if ((obji == IntervalObject.Edge) && (objj == IntervalObject.Edge )) {
      pairs = edgePairs;
      vertexFace = FALSE;
      if ((l > c )) {
        l = c;
        c = itvi.i;
      }
      lmin = edgeBox[l].min[axis];
      lmax = edgeBox[l].max[axis];
      cmin = edgeBox[c].min[axis];
      cmax = edgeBox[c].max[axis];
    } else if ((obji == IntervalObject.Vertex) && (objj == IntervalObject.Face )) {
      lmin = vertexBox[l].min[axis];
      lmax = vertexBox[l].max[axis];
      cmin = faceBox[c].min[axis];
      cmax = faceBox[c].max[axis];
    } else if ((obji == IntervalObject.Face) && (objj == IntervalObject.Vertex )) {
      l = c;
      c = itvi.i;
      lmin = vertexBox[l].min[axis];
      lmax = vertexBox[l].max[axis];
      cmin = faceBox[c].min[axis];
      cmax = faceBox[c].max[axis];
    } else {
      return;
    }
    { /* with*/ overlaped == PairTable.Overlap(pairs, l, c) ) {
      pairs[l].row[c].flag[axis] = lmax > cmin) && (cmax > lmin;
      if ((??))! overlaped AND
         PairTable.On(pairs, l) AND
         PairTable.Valid(pairs, l, c) AND
         PairTable.Overlap(pairs, l, c) 
      )) {
        CandidateList.Insert(clist, (Candidate.T){vertexFace, l, c});
      };
    };
  } /* SetStatus */;

{
  for (i = 1;  i <= ((itv.nel - 1)|?|MAX_itv);  i++) {
    { /* with*/ x == GetValue) {
      VAR interval = itv[i];
          j = i-1;
      {
        while (j >= 0) && (x < GetValue) {
          SetStatus(interval, itv[j]);
          itv[j+1] = itv[j];
          j--;
        }
        if ((j < i-1 )) { itv[j+1] = interval;
      };
    };
  };
} /* Process */;

PROCEDURE ClearCandidates(VAR CandidateList.T clist,
                          facePairs, PairTable.T edgePairs) == 
VAR ptr = clist;
    CandidateList.T prev = NULL;
    PairTable.T pairs;
{
  while (ptr != NULL ) {
    { /* with*/ cd == CandidateList.GetInfo(ptr) ) {
      if ((cd.vertexFace )) { pairs = facePairs } else { pairs = edgePairs; }
      if ((??))
        PairTable.On(pairs, cd.i) AND
        PairTable.Valid(pairs, cd.i, cd.j) AND
        PairTable.Overlap(pairs, cd.i, cd.j)
      )) {
        prev = ptr
      } else {
        if ((prev != NULL )) {
          CandidateList.RemoveNext(prev)
        } else {
          clist = CandidateList.GetNext(clist);
        };
      }
      ptr = CandidateList.GetNext(ptr);
    };
  };
} /* ClearCandidates */;
  
{; } /* IntervalList */.
