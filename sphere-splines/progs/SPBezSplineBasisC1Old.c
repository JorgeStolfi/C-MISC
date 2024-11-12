/* See SPBezSplineBasisC1Old.h */
/* Last edited on 2005-06-06 11:49:39 by stolfi */

#include <SPBezSplineBasisC1Old.h>

#include <SPFunction.h>
#include <SPSpline.h>
#include <SPHBezFunction.h>
#include <SPTriang.h>
#include <SPQuad.h>
#include <SPDeCasteljau.h>
#include <SPBezSplineBasis.h>
#include <SPBezSplineBasisC1.h>
#include <SPBasic.h>
#include <vec.h>
#include <affirm.h>
#include <nat.h>
#include <bool.h>
#include <r3.h>
#include <r3x3.h>
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <math.h>

void BuildEdgeBasisD5C1Old
  ( Arc e,
    Basis *bas,
    int *nb,
    char *triFile,
    Triangulation *tri
  );
  /* Builds the basis elements associated with the edge {e},
    for degree 5 and continuity 1. */

void BuildEdgeBasisD6C1Old
  ( Arc e,
    Basis *bas,
    int *nb,
    char *triFile,
    Triangulation *tri
  );
  /* Builds the basis elements associated with the edge {e},
    for degree 6 and continuity 1. */

void BuildEdgeBasisD7C1Old
  ( Arc e,
    Basis *bas,
    int *nb,
    char *triFile,
    Triangulation *tri
  );
  /* Builds the basis elements associated with the edge {e},
    for degree 7 and continuity 1. */

void BuildEdgeBasisD8C1Old
  ( Arc e,
    Basis *bas,
    int *nb,
    char *triFile,
    Triangulation *tri
  );
  /* Builds the basis elements associated with the edge {e},
    for degree 8 and continuity 1. */

void BuildEdgeBasisD5C1Elem1
  ( Arc e,
    Basis *bas,
    int *nb,
    char *triFile,
    Triangulation *tri
  );

void BuildEdgeBasisD6C1Elem1
  ( Arc e,
    Basis *bas,
    int *nb,
    char *triFile,
    Triangulation *tri
  );

void BuildEdgeBasisD6C1Elem2
  ( Arc e,
    Basis *bas,
    int *nb,
    char *triFile,
    Triangulation *tri
  );

void BuildEdgeBasisD6C1Elem3
  ( Arc e,
    Basis *bas,
    int *nb,
    char *triFile,
    Triangulation *tri
  );

void BuildEdgeBasisD7C1Elem1
  ( Arc e,
    Basis *bas,
    int *nb,
    char *triFile,
    Triangulation *tri
  );

void BuildEdgeBasisD7C1Elem2
  ( Arc e,
    Basis *bas,
    int *nb,
    char *triFile,
    Triangulation *tri
  );

void BuildEdgeBasisD7C1Elem3
  ( Arc e,
    Basis *bas,
    int *nb,
    char *triFile,
    Triangulation *tri
  );

void BuildEdgeBasisD7C1Elem4
  ( Arc e,
    Basis *bas,
    int *nb,
    char *triFile,
    Triangulation *tri
  );

void BuildEdgeBasisD7C1Elem5
  ( Arc e,
    Basis *bas,
    int *nb,
    char *triFile,
    Triangulation *tri
  );

void BuildEdgeBasisD8C1Elem1
  ( Arc e,
    Basis *bas,
    int *nb,
    char *triFile,
    Triangulation *tri
  );

void BuildEdgeBasisD8C1Elem2
  ( Arc e,
    Basis *bas,
    int *nb,
    char *triFile,
    Triangulation *tri
  );

void BuildEdgeBasisD8C1Elem3
  ( Arc e,
    Basis *bas,
    int *nb,
    char *triFile,
    Triangulation *tri
  );

void BuildEdgeBasisD8C1Elem4
  ( Arc e,
    Basis *bas,
    int *nb,
    char *triFile,
    Triangulation *tri
  );

void BuildEdgeBasisD8C1Elem5
  ( Arc e,
    Basis *bas,
    int *nb,
    char *triFile,
    Triangulation *tri
  );

void BuildEdgeBasisD8C1Elem6
  ( Arc e,
    Basis *bas,
    int *nb,
    char *triFile,
    Triangulation *tri
  );

void BuildEdgeBasisD8C1Elem7
  ( Arc e,
    Basis *bas,
    int *nb,
    char *triFile,
    Triangulation *tri
  );

void BuildEdgeBasisC1Old
  ( int deg,
    Arc e,
    Basis *bas,
    int *nb,
    char *triFile,
    Triangulation *tri
  )
  { affirm(deg >= 5, "bad deg/cont combination");
    switch(deg)
      { case 5: BuildEdgeBasisD5C1Old(e, bas, nb, triFile, tri); break;
        case 6: BuildEdgeBasisD6C1Old(e, bas, nb, triFile, tri); break;
        case 7: BuildEdgeBasisD7C1Old(e, bas, nb, triFile, tri); break;
        case 8: BuildEdgeBasisD8C1Old(e, bas, nb, triFile, tri); break;
        default: affirm(FALSE , "unimplemented deg/cont combination");
      }
  }

void BuildVertBasisC1Old
  ( int deg,
    Arc e,
    Basis *bas,
    int *nb,
    char *triFile,
    Triangulation *tri
  )
  { /* The original ANS construction for C1 spline basis gives a set
      of {m+3} elements {bas[0..m+2]} contained in the star of {u}.
      The consruction also defines a set {VGamma(u)} of {m+3}
      `cardinal' Bezier coefficients {c[0..m+2]}, such that
      coefficient {c[j]} of element {bas[i]} is equal to {i == j}.
      Moreover, as long as {deg >= 5}, the set {VGamma(u)} is disjoint
      from the cardinal sets of all other vertices, edges, and faces
      of {tri}.

      As explained above, in each triangle {t_r} there are only 10
      coefficients {v_r[0..9]}, which may be non-zero. Actually
      {v_r[6]} and {v_r[9]} are always set to zero; and the the C0
      continuity constraints require that {v_r[0]} is the same for all
      {r}, and that {v_r[2]} is the same as {v_{r+1}[1]}, and {v_r[5]}
      is the same as {v_{r+1}[3]}, for each {r}. Also, for each {r},
      either {v_r[8]} or {v_{r+1}[7]} is a cardinal point of an edge
      element, and thus is set to zero in every vertex element.

      Therefore, we have only {4*m+1} coefficients which need to be
      computed: namely the center coeff {v_0[0]}, and also, for each
      {r}, the coeffs {v_r[1]}, {v_r[3]}, {v_r[4]}, and either
      {v_r[8]} or {v_{r+1}[7]}. Let's denote by {V(u)} this set of
      coefficients.

      The coefficients {V(u)} are still subject to three first-order
      continuity constraints (D1)--(D3) for each edge; but the {m}
      constraints (D1) that involve {v[0]} have two degrees of
      redundancy, so we are left with {4*m+1 - (3*m-2) = m+3} degrees
      of freedom.

      In the ANS construction, these degrees of freedom are realized
      by choosing {m+3} of the coefficients {V(u)} above as the
      `cardinal set' {VGamma(u)} of the vertex {u}. These coefficients
      are {v_0[0..5]}, and {v_r[5]} for {r=1..m-3}. The coefficients
      in {VGamma(u)} can be set arbitrarily, and all other
      coefficients {v_r[i]} in {V(u)\VGamma(u)} can be computed from them
      by non-singular formulas derived from the continuity
      constraints.

      In particular, after setting {v_0[0]}, {v_0[1] == v_{m-1}[2]},
      and {v_0[2]}, we use condition (D1) repeatedly to compute
      {v_r[2] == v_{r+1}[1]} for all {r = 1..m-1}. After setting
      {v_0[3] == v_{m-1}[5]}, {v_0[4]}, {v_0[5] == v_1[3]}, and
      {v_r[5] == v_{r+1}[3]} for {r = 1..m-3}, we use we use (D2)
      repeatedly forward to compute {v_r[4]} for {r=1..m-2}, then
      backwards to compute {v_{m-1}[4]}, and then once more to compute
      {v_{m-2}[5] == v_{m-1}[3]}. 

      Finally, for each {r}, we set one of the pair {v_r[8]} or
      {v_{r+1}[7]} to 0 (the one which is in the cardinal control set
      {EGamma(e_r)} of the edge {e_r}), and use (D3) to compute
      the other one. 

      In certain degenerate situations (namely, when there are two or
      more pairs of co-circular edges meeting at {u}), the constraints
      (D1)-(D3) become even more redundant. In that case we get
      additional degrees of freedom, which could be used to generate
      more basis elements beyond the {m+3} described above. We do not
      test for such conditions here, and do not generate such
      `serendipitous' elements when they exist. */
      
    auto bool_t ArcIsEdgeReference(Arc e);
      /* TRUE if the basis elements associated with the edge {e.quad}
        have used {e} (and not {Sym(e)}) as the reference arc. */

    auto bool_t VC7BelongsToEGamma(Arc e);
    auto bool_t VC8BelongsToEGamma(Arc e); 
      /* TRUE if the Bezier coefficient {v[7]} or {v[8]},
        respectively, of the corner of {Left(e)} near {Org(e)}
        is a member of the cardinal set {EGamma(e)}. */

    auto void ComputeVC7AndVC8(Arc e, double *v);
      /* Sets the coefficients {v[7]} and {v[8]}. If they belong to the 
        set Gamma, then they are zero; otherwise they are defined
        by continuity from the coefficients {v[3]} and {v[5]},
        respectively, which are assumed to be still set as in the
        triangle {Right(e)}. */

    bool_t ArcIsEdgeReference(Arc e)
      { return tri->arc.e[2*EdgeNum(e)] == e; }

    bool_t VC7BelongsToEGamma(Arc e)
      { return ArcIsEdgeReference(e); }

    bool_t VC8BelongsToEGamma(Arc e)
      { return ArcIsEdgeReference(Sym(Onext(e))); }

    void ComputeVC7AndVC8(Arc e, double *v)
      { v[6] = 0.0;
        v[9] = 0.0;
        if (VC7BelongsToEGamma(e))
          { v[7] = 0.0; }
        else
          { r3_t bar = ComputeEdgeContinuityCoeffs(Sym(e), tri);
            v[7] = bar.c[2]*v[3];
          }
        if (VC8BelongsToEGamma(e))
          { v[8] = 0.0; }
        else
          { r3_t bar = ComputeEdgeContinuityCoeffs(Onext(e), tri);
            v[8] = bar.c[1]*v[5];
          }
      }

    int m = OrgDegree(e);
    int k;

    /* The ANS construction fails for {cont == 1} and {deg < 5}: */
    affirm(deg >= 5, "invalid deg/cont combination");
    for (k = 0; k <= m + 2; k++)
      { /* Generates the {k}th element of the basis elements
          associated with vertex {Org(e)}. */
        PieceDataRef_vec_t pd = PieceDataRef_vec_new(m); /* (to be trimmed). */
        double v[10]; /* Relevant Bezier coefficients for current wedge. */
        double w[10]; /* Relevant Bezier coefficients for last wedge. */
        Arc a;
        r3_t bar;
        int r, i;
        int npieces = 0; /* Number of non-zero wedges */
        
        /* Initialize coeffs {v(0)[0..9]}: */
        for (i = 0; i < 10; i++) { v[i] = 0.0; }
        /* For elements {k=0..5}, the cardinal Gamma point is {v[k]}
          of the first wedge: */
        if (k <= 5) { v[k] = 1.0; }
        /* Compute coeffs {w[0..9] = v(m-1)[0..9]} for face {Right(e)}: */
        bar = ComputeEdgeContinuityCoeffs(e, tri);
        for (i = 0; i < 10; i++) { w[i] = 0.0; }
        w[4] = v[4]*bar.c[0] + bar.c[1]*v[1] + bar.c[2]*v[3];
        w[1] = bar.c[0]*v[2] + bar.c[1]*v[0] + bar.c[2]*v[1];
        w[0] = v[0];
        w[2] = v[1];
        w[5] = v[3];
        
        /* Set wedges {r = 0..m-3} around {Org(e)}: */
        a = e;
        for (r = 0; r < m-2; r++)
          { /* Compute {v(r)[7]} and {v(r)[8]} of wedge {r-1}: */
            ComputeVC7AndVC8(a, v);
            if (! AllCoeffsAreZero(v, 10))
              { PieceData *t0 = MakeHBezPiece(deg);
                BezCoeff_vec_t c0 = ((HBezFn *)t0->func)->d->c;
                t0->face = Left(a)->num;
                /* Set coefficients of wedge {r}: */
                SetCoeffsVertex(a, tri, deg, c0, v, 10);
                pd.e[npieces] = t0; npieces++;
              }
            /* Compute coeffs for wedge {r+1}: */
            a = Onext(a);
            bar = ComputeEdgeContinuityCoeffs(Sym(a), tri);
            { double nv2 = bar.c[0]*v[1] + bar.c[1]*v[2] + bar.c[2]*v[0];
              double nv4 = bar.c[0]*v[4] + bar.c[1]*v[5] + bar.c[2]*v[2];
              v[1] = v[2];
              v[2] = nv2;
              v[3] = v[5];
              v[4] = nv4;
              /* Cardinal point of element {k >= 6} is {v(k-5)[5]}: */
              v[5] = (k-5 == r+1 ? 1.0 : 0.0); 
            }
          }
        /* Now {Left(a)} is wedge {r = m-2}, and {v} are its coeffs. */

        /* Finish off and store the coeffs of wedge {r = m-2}: */
        /* Ensure that wedge {m-2} fits against wedge {m-1}: */
        bar = ComputeEdgeContinuityCoeffs(Sym(Onext(a)), tri);
        v[5] = (w[4] - bar.c[0]*v[4] - bar.c[2]*v[2])/bar.c[1];
        w[3] = v[5];
        ComputeVC7AndVC8(a, v);
        if (! AllCoeffsAreZero(v, 10))
          { PieceData *t0 = MakeHBezPiece(deg);
            BezCoeff_vec_t c0 = ((HBezFn *)t0->func)->d->c;
            t0->face = Left(a)->num;
            SetCoeffsVertex(a, tri, deg, c0, v, 10);
            pd.e[npieces] = t0; npieces++;
          }

        /* Finish off and store the coeffs {w[0..9]} of wedge {r = m-1}: */
        a = Onext(a);
        ComputeVC7AndVC8(a, w);
        if (! AllCoeffsAreZero(w, 10))
          { PieceData *t0 = MakeHBezPiece(deg);
            BezCoeff_vec_t c0 = ((HBezFn *)t0->func)->d->c;
            t0->face = Left(a)->num;
            SetCoeffsVertex(a, tri, deg, c0, w, 10);
            pd.e[npieces] = t0; npieces++;
          }

        /* Phew! Trim pieces vector and store the element: */
        PieceDataRef_vec_trim(&pd, npieces);
        StoreBasisElem(bas, nb, triFile, tri, pd, &(Org(e)->pos));
      }
  }

void BuildEdgeBasisD5C1Old
  ( Arc e,
    Basis *bas,
    int *nb,
    char *triFile,
    Triangulation *tri
  )
  { BuildEdgeBasisD5C1Elem1(e, bas, nb, triFile, tri); }

void BuildEdgeBasisD5C1Elem1
  ( Arc e,
    Basis *bas,
    int *nb,
    char *triFile,
    Triangulation *tri
  )
  { int deg = 5;
    PieceDataRef_vec_t pd = MakeHBezPieces(2, deg) /*NC=21*/;
    r3_t bar = ComputeEdgeContinuityCoeffs(e, tri);

    fprintf(stderr, "BuildEdgeBasisD5C1Elem1(edge = %d)\n", EdgeNum(e));
    { PieceData *t0 = pd.e[0];
      BezCoeff_vec_t c0 = ((HBezFn *)t0->func)->d->c;
      t0->face = Left(e)->num;
      SetCoeffRelToEdge(tri, deg, c0, e, 1, 2, 2, 1.0);
    }
    { PieceData *t1 = pd.e[1];
      BezCoeff_vec_t c1 = ((HBezFn *)t1->func)->d->c;
      t1->face = Left(Sym(e))->num;
      SetCoeffRelToEdge(tri, deg, c1, Sym(e), 1, 2, 2, bar.c[0]);
    }
    StoreBasisElem(bas, nb, triFile, tri, pd, NULL);
  }

void BuildEdgeBasisD6C1Old
  ( Arc e,
    Basis *bas,
    int *nb,
    char *triFile,
    Triangulation *tri
  )
  { BuildEdgeBasisD6C1Elem1(e, bas, nb, triFile, tri);
    BuildEdgeBasisD6C1Elem2(e, bas, nb, triFile, tri);
    BuildEdgeBasisD6C1Elem3(e, bas, nb, triFile, tri);
  }

void BuildEdgeBasisD6C1Elem1
  ( Arc e,
    Basis *bas,
    int *nb,
    char *triFile,
    Triangulation *tri 
  )
  { int deg = 6;
    PieceDataRef_vec_t pd = MakeHBezPieces(2, deg) /*NC=28*/;
    r3_t bar = ComputeEdgeContinuityCoeffs(e, tri);

    fprintf(stderr, "BuildEdgeBasisD6C1Elem1(edge = %d)\n", EdgeNum(e));
    { PieceData *t0 = pd.e[0];
      BezCoeff_vec_t c0 = ((HBezFn *)t0->func)->d->c;
      t0->face = Left(e)->num;
      SetCoeffRelToEdge(tri, deg, c0, e, 1,3,2, 1.0);
    }
    { PieceData *t1 = pd.e[1];
      BezCoeff_vec_t c1 = ((HBezFn *)t1->func)->d->c;
      t1->face = Left(Sym(e))->num;
      SetCoeffRelToEdge(tri, deg, c1, Sym(e), 1,2,3, bar.c[0]);
    }
    StoreBasisElem(bas, nb, triFile, tri, pd, NULL);
  }

void BuildEdgeBasisD6C1Elem2
  ( Arc e,
    Basis *bas,
    int *nb,
    char *triFile,
    Triangulation *tri
  )
  { int deg = 6;
    PieceDataRef_vec_t pd = MakeHBezPieces(2, deg) /*NC=28*/;
    r3_t bar = ComputeEdgeContinuityCoeffs(e, tri);

    fprintf(stderr, "BuildEdgeBasisD6C1Elem2(edge = %d)\n", EdgeNum(e));
    { PieceData *t0 = pd.e[0];
      BezCoeff_vec_t c0 = ((HBezFn *)t0->func)->d->c;
      t0->face = Left(e)->num;
      SetCoeffRelToEdge(tri, deg, c0, e, 1,2,3, 1.0);
    }
    { PieceData *t1 = pd.e[1];
      BezCoeff_vec_t c1 = ((HBezFn *)t1->func)->d->c;
      t1->face = Left(Sym(e))->num;
      SetCoeffRelToEdge(tri, deg, c1, Sym(e), 1,3,2, bar.c[0]);
    }
    StoreBasisElem(bas, nb, triFile, tri, pd, NULL);
  }

void BuildEdgeBasisD6C1Elem3
  ( Arc e,
    Basis *bas,
    int *nb,
    char *triFile,
    Triangulation *tri
  )
  { int deg = 6;
    PieceDataRef_vec_t pd = MakeHBezPieces(2, deg) /*NC=28*/;
    r3_t bar = ComputeEdgeContinuityCoeffs(e, tri);
    fprintf(stderr, "BuildEdgeBasisD6C1Elem3(edge = %d)\n", EdgeNum(e));
    { PieceData *t0 = pd.e[0];
      BezCoeff_vec_t c0 = ((HBezFn *)t0->func)->d->c;
      t0->face = Left(e)->num;
      SetCoeffRelToEdge(tri, deg, c0, e, 0,3,3, 1.0);
    }
    { PieceData *t1 = pd.e[1];
      BezCoeff_vec_t c1 = ((HBezFn *)t1->func)->d->c;
      t1->face = Left(Sym(e))->num;
      SetCoeffRelToEdge(tri, deg, c1, Sym(e), 0,3,3, 1.0);
      SetCoeffRelToEdge(tri, deg, c1, Sym(e), 1,2,3, bar.c[2]);
      SetCoeffRelToEdge(tri, deg, c1, Sym(e), 1,3,2, bar.c[1]);
    }
    StoreBasisElem(bas, nb, triFile, tri, pd, NULL);
  }

void BuildEdgeBasisD7C1Old
  ( Arc e,
    Basis *bas,
    int *nb,
    char *triFile,
    Triangulation *tri
  )
  { BuildEdgeBasisD7C1Elem1(e, bas, nb, triFile, tri);
    BuildEdgeBasisD7C1Elem2(e, bas, nb, triFile, tri);
    BuildEdgeBasisD7C1Elem3(e, bas, nb, triFile, tri);
    BuildEdgeBasisD7C1Elem4(e, bas, nb, triFile, tri);
    BuildEdgeBasisD7C1Elem5(e, bas, nb, triFile, tri);
  }

void BuildEdgeBasisD7C1Elem1
  ( Arc e,
    Basis *bas,
    int *nb,
    char *triFile,
    Triangulation *tri
  )
  { int deg = 7;
    r3_t bar = ComputeEdgeContinuityCoeffs(e, tri);
    PieceDataRef_vec_t pd = MakeHBezPieces(2, deg) /*NC=36*/;
    fprintf(stderr, "BuildEdgeBasisD7C1Elem1(edge = %d)\n", EdgeNum(e));
    { PieceData *t0 = pd.e[0];
      BezCoeff_vec_t c0 = ((HBezFn *)t0->func)->d->c;
      t0->face = Left(e)->num;
      SetCoeffRelToEdge(tri, deg, c0, e, 1,4,2, 1.0);
    }
    { PieceData *t1 = pd.e[1];
      BezCoeff_vec_t c1 = ((HBezFn *)t1->func)->d->c;
      t1->face = Left(Sym(e))->num;
      SetCoeffRelToEdge(tri, deg, c1, Sym(e), 1,2,4, bar.c[0]);
    }
    StoreBasisElem(bas, nb, triFile, tri, pd, NULL);
  }

void BuildEdgeBasisD7C1Elem2
  ( Arc e,
    Basis *bas,
    int *nb,
    char *triFile,
    Triangulation *tri
  )
  { int deg = 7;
    r3_t bar = ComputeEdgeContinuityCoeffs(e, tri);
    PieceDataRef_vec_t pd = MakeHBezPieces(2, deg) /*NC=36*/;
    fprintf(stderr, "BuildEdgeBasisD7C1Elem2(edge = %d)\n", EdgeNum(e));
    {
      PieceData *t0 = pd.e[0];
      BezCoeff_vec_t c0 = ((HBezFn *)t0->func)->d->c;
      t0->face = Left(e)->num;
      SetCoeffRelToEdge(tri, deg, c0, e, 1,3,3, 1.0);
    }
    { PieceData *t1 = pd.e[1];
      BezCoeff_vec_t c1 = ((HBezFn *)t1->func)->d->c;
      t1->face = Left(Sym(e))->num;
      SetCoeffRelToEdge(tri, deg, c1, Sym(e), 1,3,3, bar.c[0]);
    }
    StoreBasisElem(bas, nb, triFile, tri, pd, NULL);
  }

void BuildEdgeBasisD7C1Elem3
  ( Arc e,
    Basis *bas,
    int *nb,
    char *triFile,
    Triangulation *tri
  )
  { int deg = 7;
    r3_t bar = ComputeEdgeContinuityCoeffs(e, tri);
    PieceDataRef_vec_t pd = MakeHBezPieces(2, deg) /*NC=36*/;
    fprintf(stderr, "BuildEdgeBasisD7C1Elem3(edge = %d)\n", EdgeNum(e));
    { PieceData *t0 = pd.e[0];
      BezCoeff_vec_t c0 = ((HBezFn *)t0->func)->d->c;
      t0->face = Left(e)->num;
      SetCoeffRelToEdge(tri, deg, c0, e, 1,2,4, 1.0);
    }
    { PieceData *t1 = pd.e[1];
      BezCoeff_vec_t c1 = ((HBezFn *)t1->func)->d->c;
      t1->face = Left(Sym(e))->num;
      SetCoeffRelToEdge(tri, deg, c1, Sym(e), 1,4,2, bar.c[0]);
    }
    StoreBasisElem(bas, nb, triFile, tri, pd, NULL);
  }

void BuildEdgeBasisD7C1Elem4
  ( Arc e,
    Basis *bas,
    int *nb,
    char *triFile,
    Triangulation *tri
  )
  { int deg = 7;
    r3_t bar = ComputeEdgeContinuityCoeffs(e, tri);
    PieceDataRef_vec_t pd = MakeHBezPieces(2, deg) /*NC=36*/;
    fprintf(stderr, "BuildEdgeBasisD7C1Elem4(edge = %d)\n", EdgeNum(e));
    {
      PieceData *t0 = pd.e[0];
      BezCoeff_vec_t c0 = ((HBezFn *)t0->func)->d->c;
      t0->face = Left(e)->num;
      SetCoeffRelToEdge(tri, deg, c0, e, 0,4,3, 1.0);
    }
    { PieceData *t1 = pd.e[1];
      BezCoeff_vec_t c1 = ((HBezFn *)t1->func)->d->c;
      t1->face = Left(Sym(e))->num;
      SetCoeffRelToEdge(tri, deg, c1, Sym(e), 0,3,4, 1.0);
      SetCoeffRelToEdge(tri, deg, c1, Sym(e), 1,2,4, bar.c[2]);
      SetCoeffRelToEdge(tri, deg, c1, Sym(e), 1,3,3, bar.c[1]);
    }
    StoreBasisElem(bas, nb, triFile, tri, pd, NULL);
  }

void BuildEdgeBasisD7C1Elem5
  ( Arc e,
    Basis *bas,
    int *nb,
    char *triFile,
    Triangulation *tri
  )
  { int deg = 7;
    r3_t bar = ComputeEdgeContinuityCoeffs(e, tri);
    PieceDataRef_vec_t pd = MakeHBezPieces(2, deg) /*NC=36*/;
    fprintf(stderr, "BuildEdgeBasisD7C1Elem5(edge = %d)\n", EdgeNum(e));
    { PieceData *t0 = pd.e[0];
      BezCoeff_vec_t c0 = ((HBezFn *)t0->func)->d->c;
      t0->face = Left(e)->num;
      SetCoeffRelToEdge(tri, deg, c0, e, 0,3,4, 1.0);
    }
    { PieceData *t1 = pd.e[1];
      BezCoeff_vec_t c1 = ((HBezFn *)t1->func)->d->c;
      t1->face = Left(Sym(e))->num;
      SetCoeffRelToEdge(tri, deg, c1, Sym(e), 0,4,3, 1.0);
      SetCoeffRelToEdge(tri, deg, c1, Sym(e), 1,3,3, bar.c[2]);
      SetCoeffRelToEdge(tri, deg, c1, Sym(e), 1,4,2, bar.c[1]);
    }
    StoreBasisElem(bas, nb, triFile, tri, pd, NULL);
  }

void BuildEdgeBasisD8C1Old
  ( Arc e,
    Basis *bas,
    int *nb,
    char *triFile,
    Triangulation *tri
  )
  { BuildEdgeBasisD8C1Elem1(e, bas, nb, triFile, tri);
    BuildEdgeBasisD8C1Elem2(e, bas, nb, triFile, tri);
    BuildEdgeBasisD8C1Elem3(e, bas, nb, triFile, tri);
    BuildEdgeBasisD8C1Elem4(e, bas, nb, triFile, tri);
    BuildEdgeBasisD8C1Elem5(e, bas, nb, triFile, tri);
    BuildEdgeBasisD8C1Elem6(e, bas, nb, triFile, tri);
    BuildEdgeBasisD8C1Elem7(e, bas, nb, triFile, tri);
  }

void BuildEdgeBasisD8C1Elem1
  ( Arc e,
    Basis *bas,
    int *nb,
    char *triFile,
    Triangulation *tri
  )
  { int deg = 8;
    r3_t bar = ComputeEdgeContinuityCoeffs(e, tri);
    PieceDataRef_vec_t pd = MakeHBezPieces(2, deg) /*NC=45*/;
    fprintf(stderr, "BuildEdgeBasisD8C1Elem1(edge = %d)\n", EdgeNum(e));
    { PieceData *t0 = pd.e[0];
      BezCoeff_vec_t c0 = ((HBezFn *)t0->func)->d->c;
      t0->face = Left(e)->num;
      SetCoeffRelToEdge(tri, deg, c0, e, 1,5,2, 1.0);
    }
    { PieceData *t1 = pd.e[1];
      BezCoeff_vec_t c1 = ((HBezFn *)t1->func)->d->c;
      t1->face = Left(Sym(e))->num;
      SetCoeffRelToEdge(tri, deg, c1, Sym(e), 1,2,5, bar.c[0]);
    }
    StoreBasisElem(bas, nb, triFile, tri, pd, NULL);
  }

void BuildEdgeBasisD8C1Elem2
  ( Arc e,
    Basis *bas,
    int *nb,
    char *triFile,
    Triangulation *tri
  )
  { int deg = 8;
    r3_t bar = ComputeEdgeContinuityCoeffs(e, tri);
    PieceDataRef_vec_t pd = MakeHBezPieces(2, deg) /*NC=45*/;
    fprintf(stderr, "BuildEdgeBasisD8C1Elem2(edge = %d)\n", EdgeNum(e));
    {
      PieceData *t0 = pd.e[0];
      BezCoeff_vec_t c0 = ((HBezFn *)t0->func)->d->c;
      t0->face = Left(e)->num;
      SetCoeffRelToEdge(tri, deg, c0, e, 1,4,3, 1.0);
    }
    { PieceData *t1 = pd.e[1];
      BezCoeff_vec_t c1 = ((HBezFn *)t1->func)->d->c;
      t1->face = Left(Sym(e))->num;
      SetCoeffRelToEdge(tri, deg, c1, Sym(e), 1,3,4, bar.c[0]);
    }
    StoreBasisElem(bas, nb, triFile, tri, pd, NULL);
  }

void BuildEdgeBasisD8C1Elem3
  ( Arc e,
    Basis *bas,
    int *nb,
    char *triFile,
    Triangulation *tri
  )
  { int deg = 8;
    r3_t bar = ComputeEdgeContinuityCoeffs(e, tri);
    PieceDataRef_vec_t pd = MakeHBezPieces(2, deg) /*NC=45*/;
    fprintf(stderr, "BuildEdgeBasisD8C1Elem3(edge = %d)\n", EdgeNum(e));
    { PieceData *t0 = pd.e[0];
      BezCoeff_vec_t c0 = ((HBezFn *)t0->func)->d->c;
      t0->face = Left(e)->num;
      SetCoeffRelToEdge(tri, deg, c0, e, 1,3,4, 1.0);
    }
    { PieceData *t1 = pd.e[1];
      BezCoeff_vec_t c1 = ((HBezFn *)t1->func)->d->c;
      t1->face = Left(Sym(e))->num;
      SetCoeffRelToEdge(tri, deg, c1, Sym(e), 1,4,3, bar.c[0]);
    }
    StoreBasisElem(bas, nb, triFile, tri, pd, NULL);
  }

void BuildEdgeBasisD8C1Elem4
  ( Arc e,
    Basis *bas,
    int *nb,
    char *triFile,
    Triangulation *tri
  )
  { int deg = 8;
    r3_t bar = ComputeEdgeContinuityCoeffs(e, tri);
    PieceDataRef_vec_t pd = MakeHBezPieces(2, deg) /*NC=45*/;
    fprintf(stderr, "BuildEdgeBasisD8C1Elem4(edge = %d)\n", EdgeNum(e));
    { PieceData *t0 = pd.e[0];
      BezCoeff_vec_t c0 = ((HBezFn *)t0->func)->d->c;
      t0->face = Left(e)->num;
      SetCoeffRelToEdge(tri, deg, c0, e, 1,2,5, 1.0);
    }
    { PieceData *t1 = pd.e[1];
      BezCoeff_vec_t c1 = ((HBezFn *)t1->func)->d->c;
      t1->face = Left(Sym(e))->num;
      SetCoeffRelToEdge(tri, deg, c1, Sym(e), 1,5,2, bar.c[0]);
    }
    StoreBasisElem(bas, nb, triFile, tri, pd, NULL);
  }

void BuildEdgeBasisD8C1Elem5
  ( Arc e,
    Basis *bas,
    int *nb,
    char *triFile,
    Triangulation *tri
  )
  { int deg = 8;
    r3_t bar = ComputeEdgeContinuityCoeffs(e, tri);
    PieceDataRef_vec_t pd = MakeHBezPieces(2, deg) /*NC=45*/;
    fprintf(stderr, "BuildEdgeBasisD8C1Elem5(edge = %d)\n", EdgeNum(e));
    {
      PieceData *t0 = pd.e[0];
      BezCoeff_vec_t c0 = ((HBezFn *)t0->func)->d->c;
      t0->face = Left(e)->num;
      SetCoeffRelToEdge(tri, deg, c0, e, 0,5,3, 1.0);
    }
    { PieceData *t1 = pd.e[1];
      BezCoeff_vec_t c1 = ((HBezFn *)t1->func)->d->c;
      t1->face = Left(Sym(e))->num;
      SetCoeffRelToEdge(tri, deg, c1, Sym(e), 0,3,5, 1.0);
      SetCoeffRelToEdge(tri, deg, c1, Sym(e), 1,2,5, bar.c[2]);
      SetCoeffRelToEdge(tri, deg, c1, Sym(e), 1,3,4, bar.c[1]);
    }
    StoreBasisElem(bas, nb, triFile, tri, pd, NULL);
  }

void BuildEdgeBasisD8C1Elem6
  ( Arc e,
    Basis *bas,
    int *nb,
    char *triFile,
    Triangulation *tri
  )
  { int deg = 8;
    r3_t bar = ComputeEdgeContinuityCoeffs(e, tri);
    PieceDataRef_vec_t pd = MakeHBezPieces(2, deg) /*NC=45*/;
    fprintf(stderr, "BuildEdgeBasisD8C1Elem6(edge = %d)\n", EdgeNum(e));
    { PieceData *t0 = pd.e[0];
      BezCoeff_vec_t c0 = ((HBezFn *)t0->func)->d->c;
      t0->face = Left(e)->num;
      SetCoeffRelToEdge(tri, deg, c0, e, 0,4,4, 1.0);
    }
    { PieceData *t1 = pd.e[1];
      BezCoeff_vec_t c1 = ((HBezFn *)t1->func)->d->c;
      t1->face = Left(Sym(e))->num;
      SetCoeffRelToEdge(tri, deg, c1, Sym(e), 0,4,4, 1.0);
      SetCoeffRelToEdge(tri, deg, c1, Sym(e), 1,3,4, bar.c[2]);
      SetCoeffRelToEdge(tri, deg, c1, Sym(e), 1,4,3, bar.c[1]);
    }
    StoreBasisElem(bas, nb, triFile, tri, pd, NULL);
  }

void BuildEdgeBasisD8C1Elem7
  ( Arc e,
    Basis *bas,
    int *nb,
    char *triFile,
    Triangulation *tri
  )
  { int deg = 8;
    r3_t bar = ComputeEdgeContinuityCoeffs(e, tri);
    PieceDataRef_vec_t pd = MakeHBezPieces(2, deg) /*NC=45*/;
    fprintf(stderr, "BuildEdgeBasisD8C1Elem7(edge = %d)\n", EdgeNum(e));
    { PieceData *t0 = pd.e[0];
      BezCoeff_vec_t c0 = ((HBezFn *)t0->func)->d->c;
      t0->face = Left(e)->num;
      SetCoeffRelToEdge(tri, deg, c0, e, 0,3,5, 1.0);
    }
    { PieceData *t1 = pd.e[1];
      BezCoeff_vec_t c1 = ((HBezFn *)t1->func)->d->c;
      t1->face = Left(Sym(e))->num;
      SetCoeffRelToEdge(tri, deg, c1, Sym(e), 0,5,3, 1.0);
      SetCoeffRelToEdge(tri, deg, c1, Sym(e), 1,4,3, bar.c[2]);
      SetCoeffRelToEdge(tri, deg, c1, Sym(e), 1,5,2, bar.c[1]);
    }
    StoreBasisElem(bas, nb, triFile, tri, pd, NULL);
  }

