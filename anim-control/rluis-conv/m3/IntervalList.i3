INTERFACE IntervalList;

IMPORT SystemTopology, BoxList, PairTable, CandidateList;

TYPE
  IntervalObject = {Vertex, Edge, Face};
  
  Interval = RECORD
      object: IntervalObject;
      begin: BOOLEAN;
      i: CARDINAL;
    END;
  T = REF ARRAY OF Interval;
  
PROCEDURE New(s: SystemTopology.T): T;
PROCEDURE Process(itv: T; axis: PairTable.Axis;
                  vertexBox, edgeBox, faceBox: BoxList.T;
                  facePairs, edgePairs: PairTable.T;
                  VAR clist: CandidateList.T);
PROCEDURE ClearCandidates(VAR clist: CandidateList.T;
                          facePairs, edgePairs: PairTable.T);

END IntervalList.
