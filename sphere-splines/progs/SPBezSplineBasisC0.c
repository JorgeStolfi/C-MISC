/* See SPBezSplineBasisC0.h */
/* Last edited on 2005-06-06 11:48:43 by stolfi */

#include <SPBezSplineBasisC0.h>

#include <SPFunction.h>
#include <SPSpline.h>
#include <SPHBezFunction.h>
#include <SPTriang.h>
#include <SPQuad.h>
#include <SPDeCasteljau.h>
#include <SPBezSplineBasis.h>
#include <SPBasic.h>
#include <vec.h>
#include <affirm.h>
#include <nat.h>
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <math.h>
  
void BuildEdgeBasisC0
  ( int deg,
    Arc e,
    Basis *bas,
    int *nb,
    char *triFile,
    Triangulation *tri
  )
  { /* Provided that {deg >= 1}, there are {deg - 1} C0 edge elements. */
    int ke;
    affirm(deg >= 1, "bad deg/cont combination");
    for (ke = 0; ke < deg-1; ke++)
      { /* The C0 edge element of index {ke} has the Bezier coefficients
          {ca_{0,ke+1,deg-ke-1}} and {cb_{0,deg-ke-1,ke+1}} set to 1,
          all other coeffs set to zero. */
        PieceDataRef_vec_t pd = MakeHBezPieces(2, deg);
        int i = 0;
        int j = ke + 1;
        int k = deg - i - j;
        fprintf(stderr, "BuildEdgeBasisC0Elem(edge = %d)\n", EdgeNum(e));
        affirm(i == 0, "inconsistent i");
        affirm((j >= 1) && (j <= deg-1), "inconsistent j");
        affirm((k >= 1) && (k <= deg-1), "inconsistent k");
        { PieceData *t0 = pd.e[0];
          BezCoeff_vec_t c0 = ((HBezFn *)t0->func)->d->c;
          t0->face = Left(e)->num;
          SetCoeffRelToEdge(tri, deg, c0, e, i, j, k, 1.0);
        }
        { PieceData *t1 = pd.e[1];
          BezCoeff_vec_t c1 = ((HBezFn *)t1->func)->d->c;
          t1->face = Left(Sym(e))->num;
          SetCoeffRelToEdge(tri, deg, c1, Sym(e), i, k, j, 1.0);
        }
        StoreBasisElem(bas, nb, triFile, tri, pd, NULL);
      }
  }

void BuildVertBasisC0
  ( int deg,
    Arc e,
    Basis *bas,
    int *nb,
    char *triFile,
    Triangulation *tri
  )
  { /* In that case, there is exactly one C0 element associated
      to the vertex {v}. The element has {m = OrgDegree(v)} pieces, one
      for each triangle incident to {v}.  In each piece, all Bezier 
      coefficients are zero, except for the coefficient associated to 
      the corner {v}, which is 1.0. */
    
    Arc a;
    int order = OrgDegree(e);
    PieceDataRef_vec_t pd = MakeHBezPieces(order, deg);
    double v[1];
    int r;

    /* The ANS construction fails for {cont == 0} and {deg < 1}: */
    affirm(deg >= 1, "invalid deg/cont combination");

    v[0] = 1.0;
    a = e;
    for (r = 0; r < order; r++)
      { PieceData *t0 = pd.e[r];
        BezCoeff_vec_t c0 = ((HBezFn *)t0->func)->d->c;
        t0->face = Left(a)->num;
        SetCoeffsVertex(a, tri, deg, c0, v, 1);
        a = Onext(a);
      }
    StoreBasisElem(bas, nb, triFile, tri, pd, &(Org(e)->pos));
  }
