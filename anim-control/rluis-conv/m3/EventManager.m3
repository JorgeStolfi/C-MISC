MODULE EventManager;

TYPE
  JointManager = T OBJECT
      em: REF ARRAY OF T;
    METHODS
      detect := Detect;
      advance := Advance;
    END;
    
PROCEDURE Join(READONLY em: ARRAY OF T): T =
  BEGIN
    WITH
      t = NEW(JointManager),
      NE = NUMBER(em),
      re = NEW(REF ARRAY OF T, NE)
    DO
      re^ := em;
      t.em := re;
      RETURN t
    END
  END Join;
  
PROCEDURE Detect(
    t: T
    t0: LONGREAL;
    READONLY p0: Position;
    READONLY v0: Velocity;
    t1: LONGREAL;
    READONLY p1: Position;
    READONLY v1: Velocity;
  ): Time = 
  VAR te: Time := LAST(Time);
  BEGIN
    WITH
      em = t.em^,
      NE = NUMBER(em)
    DO
      FOR i := 0 TO NE-1 DO
        IF em[i] # NIL THEN
          te := MIN(te, em[i].detect(t0, p0, v0, t1, p1, v1)
        END
      END
    END
  END Detect;

PROCEDURE Advance(
    t: T;
    t0: LONGREAL;
    READONLY p0: Position;
    t1: LONGREAL;
    READONLY p1: Position;
    READONLY v: Velocity;
  ) =
  BEGIN
    WITH
      em = t.em^,
      NE = NUMBER(em)
    DO
      FOR i := 0 TO NE-1 DO
        IF em[i] # NIL THEN
          em[i].advance(t0, p0, t1, p1, v)
        END
      END
    END
  END Advance;

BEGIN
END EventManager.
