MODULE GaussElim;

CONST EPSILON = 1.0D-8;

PROCEDURE Solve(
    VAR A: ARRAY OF ARRAY OF LONGREAL; 
    VAR x, b: ARRAY OF LONGREAL;
    n: CARDINAL;
  ) =
  VAR m: CARDINAL;
      pivot, s: LONGREAL;
  BEGIN
    <* ASSERT NUMBER(A) >= n *>
    <* ASSERT NUMBER(A[0]) >= n *>
    <* ASSERT NUMBER(b) >= n *>
    FOR k := 0 TO n-2 DO
      m := k;
      FOR i := k+1 TO n-1 DO
        IF ABS(A[i,k]) > ABS(A[m,k]) THEN m := i END
      END;
      IF m # k THEN
        WITH Am = A[m], Ak = A[k] DO
          FOR j := k+1 TO n-1 DO
            s := Am[j]; Am[j] := Ak[j]; Ak[j] := s
          END;
        END;
        s := b[m]; b[m] := b[k]; b[k] := s
      END;
      pivot := A[k,k];
      <* ASSERT pivot # 0.0D0 *>
      <* ASSERT pivot < -EPSILON OR pivot > EPSILON *>

      FOR i := k+1 TO n-1 DO
        WITH Ai = A[i], Ak = A[k] DO
          s := -Ai[k]/pivot;
          FOR j := k TO n-1 DO
            Ai[j] := Ai[j] + s*Ak[j];
          END;
        END;
        b[i] := b[i] + s*b[k];
      END
    END;

    FOR i := n-1 TO 0 BY -1 DO
      s := 0.0D0;
      FOR j := i+1 TO n-1 DO
        s := s + A[i,j]*x[j];
      END;
      x[i] := (b[i] - s)/A[i,i];
      <* ASSERT b[i] <= 0.0d0 OR b[i] >= 0.0d0 *>
      <* ASSERT A[i,i] # 0.0d0 *>
      <* ASSERT x[i] <= 0.0D0 OR x[i] >= 0.0d0 *>
    END
  END Solve;

BEGIN
END GaussElim.
