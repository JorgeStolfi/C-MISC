INTERFACE Contact;

IMPORT Rd, Wr, Force;

TYPE
  T = RECORD
      vertexFace: BOOLEAN;
      i, j: CARDINAL;
      force: Force.T;
      sign: BOOLEAN;  (* Spring direction for edge-edge collisions *)
        (* 
          sign = TRUE means remove the spring when the edge-edge determinant
          becomes positive again *)
    END;

EXCEPTION EndOfInput;

PROCEDURE Read(rd: Rd.T; VAR c: T) RAISES {EndOfInput};
PROCEDURE Write(wr: Wr.T; c: T);

END Contact.
