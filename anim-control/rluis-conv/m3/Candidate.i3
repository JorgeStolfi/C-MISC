INTERFACE Candidate;

IMPORT Rd, Wr;

TYPE
  T = RECORD
      vertexFace: BOOLEAN;
      i, j: CARDINAL;      (* Box indices *)
    END;
      
EXCEPTION EndOfInput;

PROCEDURE Read(rd: Rd.T; VAR c: T) RAISES {EndOfInput};
PROCEDURE Write(wr: Wr.T; c: T);

END Candidate.
