INTERFACE BezierSearch;

CONST NoRoot = LAST(LONGREAL);

PROCEDURE FindRootA(
    P0, P1, P2, P3,
    A0, A1, A2, A3,
    B0, B1, B2, B3,
    C0, C1, C2, C3: LONGREAL
  ): LONGREAL;

PROCEDURE FindRootB(
    P0, P1, P2, P3,
    A0, A1, A2, A3,
    B0, B1, B2, B3,
    C0, C1, C2, C3,
    D0, D1, D2, D3: LONGREAL;
  ): LONGREAL;

PROCEDURE FindRoot(P0, P1, P2, P3: LONGREAL; sign: BOOLEAN): LONGREAL;
  (*
    Returns a time t in [0 _ 1] where the Bezier arc changes from
    negative to positive (if sign=TRUE) or from positive to negative
    (if sign=FALSE); or NoRoot if no such t exists. *)

END BezierSearch.
