
#ifndef BezierSearch_H
#define BezierSearch_H



CONST NoRoot == ((double.nel - 1)|?|MAX_LONGREAL);

double FindRootA(
    P0, P1, P2, P3,
    A0, A1, A2, A3,
    B0, B1, B2, B3,
    C0, C1, C2, C3: double
  );

double FindRootB(
    P0, P1, P2, P3,
    A0, A1, A2, A3,
    B0, B1, B2, B3,
    C0, C1, C2, C3,
    double D0, D1, D2, D3;
  );

double FindRoot(P0, P1, P2, double P3, bool sign);
  /*
    Returns a time t in [0 _ 1] where the Bezier arc changes from
    negative to positive (if sig n == TRUE) or from positive to negative
    (if sig n == FALSE); or NoRoot if no such t exists. */

#endif
