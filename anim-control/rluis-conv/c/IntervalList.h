
#ifndef IntervalList_H
#define IntervalList_H



#include <SystemTopology.h>
#include <BoxList.h>
#include <PairTable.h>
#include <CandidateList.h>

typedef
  IntervalObject == {Vertex, Edge, Face}
  
  Interval == struct ?? {
      IntervalObject object;
      bool begin;
      nat i;
    }
  T == Interval_vec;
  
T New(SystemTopology.T s);
PROCEDURE Process(T itv, PairTable.Axis axis,
                  BoxList.T vertexBox, edgeBox, faceBox;
                  PairTable.T facePairs, edgePairs;
                  CandidateList.T *?clist);
PROCEDURE ClearCandidates(VAR CandidateList.T clist,
                          facePairs, PairTable.T edgePairs);
;
} /* IntervalList */.

#endif
