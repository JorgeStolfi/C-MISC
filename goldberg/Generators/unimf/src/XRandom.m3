MODULE XRandom;

(* Created June 1997 by J. Stolfi  *)
(* Based on "Random.m3" by Mark R. Brown and Bill Kalsow *)
(* See the copyright notice at the end of this file. *)

IMPORT Word, RandomReal;

CONST 
  MinSeed = LAST(INTEGER) DIV 10 + 1;
  MaxSeed = (LAST(INTEGER) DIV 10) * 9;
  (* 
    Actual seed must be in [LAST(INTEGER)/10 .. 9*LAST(INTEGER)/10]. *)

REVEAL
  T = Public BRANDED OBJECT
      i: [0..56];
      a: ARRAY [1..55] OF INTEGER;
    OVERRIDES
      init     := Init;
      integer  := Integer;
      real     := Real;
      longreal := Longreal;
      extended := Extended;
      boolean  := Boolean;
    END;

PROCEDURE Init (t: T;  seed: Seed): T =
  BEGIN
    IF (t = NIL) THEN t := NEW (T) END;
    Start (t, seed + MinSeed);
    RETURN t;
  END Init;

PROCEDURE Start (t: T;  seed: INTEGER) =
  VAR j, k: INTEGER;  i2: [1..54];
  BEGIN
    (* For an explanation of this initialization procedure see the Fortran
       program in Section 3.6 of Knuth Volume 2, second edition. *)
    <* ASSERT seed < MaxSeed *>
    <* ASSERT seed > MinSeed *>
    t.a[55] := seed;
    j := seed;
    k := 1;
    FOR i := 1 TO 54 DO
      i2 := (21 * i) MOD 55;
      t.a[i2] := k;
      k := Word.Minus (j, k); (* ignore underflow *)
      j := t.a[i2];
    END;
    FOR i := 1 TO 20 DO Next55 (t) END;
  END Start;

PROCEDURE Next55 (t: T) =
  BEGIN
    FOR j := 55 TO 32 BY  -1 DO
      (* DEC (t.a[j], t.a[j - 31]); (* ignore underflow *) *)
      t.a[j] := Word.Minus (t.a[j], t.a[j - 31]);
    END;
    FOR j := 31 TO 1 BY  -1 DO
      (* DEC (t.a[j], t.a[j + 24]); (* ignore underflow *) *)
      t.a[j] := Word.Minus (t.a[j], t.a[j + 24]);
    END;
    t.i := 56;
  END Next55;

CONST
  halfWordSize = Word.Size DIV 2;
  halfWordRadix = Word.LeftShift(1, halfWordSize);
  halfWordMask = halfWordRadix - 1;

PROCEDURE Integer (t: T;  min, max: INTEGER): INTEGER =
  VAR x, rem, range: INTEGER;
  BEGIN
    <*ASSERT min <= max *>
    LOOP
      LOOP
        DEC (t.i);
        IF t.i > 0 THEN EXIT END;
        Next55 (t);
      END;
      x := t.a[t.i];

      IF (0 < min) OR (min + LAST (INTEGER) > max) THEN
        (* the range less than LAST (INTEGER) *)
        range := max - min + 1;
        IF range < halfWordRadix THEN
          VAR xl := Word.And(x, halfWordMask);
              xh := Word.RightShift(x, halfWordSize);
              res: INTEGER;
          BEGIN
            res := xl * range;
            res := Word.RightShift(res, halfWordSize);
            res := res + xh * range;
            res := Word.RightShift(res, halfWordSize);
            RETURN min + res;
          END (* BEGIN *)
        ELSE
          x := Word.And (x, LAST (INTEGER));  (* 0 <= x <= LAST (INTEGER *)
          rem := x MOD range;                 (* 0 <= rem < range *)
          IF x - rem <= LAST (INTEGER) - range THEN
            (* got a full block:  x - rem + range <= LAST (INTEGER) *)
            RETURN min + rem;
          END
        END (* IF *)
      ELSE
        (* the range is very large, but not complete *)
        (* so, just keep trying until we find a legal value *)
        IF (min <= x) AND (x <= max) THEN RETURN x END;
      END;
    END;
  END Integer;

PROCEDURE Boolean (t: T): BOOLEAN =
  BEGIN
    LOOP
      DEC (t.i);
      IF t.i > 0 THEN EXIT END;
      Next55 (t);
    END;
    RETURN VAL (Word.And (1, t.a[t.i]), BOOLEAN);
  END Boolean;

PROCEDURE Real (t: T;  min, max: REAL): REAL =
  BEGIN
    <*ASSERT min <= max *>
    RETURN (max - min) * RandomReal.Real(t) + min;
  END Real;

PROCEDURE Longreal (t: T;  min, max: LONGREAL): LONGREAL =
  BEGIN
    <*ASSERT min <= max *>
    RETURN (max - min) * RandomReal.Longreal (t) + min;
  END Longreal;

PROCEDURE Extended (t: T;  min, max: EXTENDED): EXTENDED =
  BEGIN
    <*ASSERT min <= max *>
    RETURN (max - min) * RandomReal.Extended (t) + min;
  END Extended;

BEGIN
  <*ASSERT Word.Size MOD 2 = 0 *>
END XRandom.

(* Random.m3 was copyright (C) 1994, Digital Equipment Corp. *)
