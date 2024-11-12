/* See SPBezSplineBasisC1New.h */
/* Last edited on 2005-06-06 11:49:52 by stolfi */

#include <SPBezSplineBasisC1New.h>

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
#include <bool.h>
#include <r3.h>
#include <r3x3.h>
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <math.h>

void ComputeEdgeC1ValElemCoeffs
  ( int kv,
    int deg,
    Arc a,
    BezCoeff_vec_t ca,
    BezCoeff_vec_t cb,
    Triangulation *tri,
    r3x3_t *frame
  );
  /* Sets the Bezier coeffs for the val-type edge element
    of degree {deg} and index {kv}, relative to the edge {e}.
    The arguments {ca} and {cb} are the Bezier coeff arrays for
    the faces {Left(a)} and {Left(b)}.  The {frame} is the 
    result of {EdgeFrame(Org(e), Dest(e))}. */

void ComputeEdgeC1DerElemCoeffs
  ( int kd,
    int deg,
    Arc e,
    BezCoeff_vec_t a,
    BezCoeff_vec_t b,
    Triangulation *tri,
    r3x3_t *frame
  );
  /* Sets the Bezier coeffs for the der-type edge element
    of degree {deg} and index {kv}, relative to the edge {e}.
    The arguments {ca} and {cb} are the Bezier coeff arrays for
    the faces {Left(a)} and {Left(b)}.  The {frame} is the 
    result of {EdgeFrame(Org(e), Dest(e))}. */

void BuildEdgeBasisC1New
  ( int deg,
    Arc e,
    Basis *bas,
    int *nb,
    char *triFile,
    Triangulation *tri
  )
  { int NBVal = deg - 5;
    int NBDer = deg - 4;
    int NB = NBVal + NBDer;
    S2Point *oe = &(Org(e)->pos);
    S2Point *de = &(Dest(e)->pos);
    r3x3_t frame = EdgeFrame(oe, de);
    int k;
    
    affirm(deg >= 5, "bad deg/cont combination");
    fprintf(stderr, "BuildEdgeBasisC1New(edge = %d)\n", EdgeNum(e));
    for (k = 0; k < NB; k++)
      { PieceDataRef_vec_t pd = MakeHBezPieces(2, deg);
        PieceData *ta = pd.e[0];
        PieceData *tb = pd.e[1];
        BezCoeff_vec_t ca = ((HBezFn *)ta->func)->d->c;
        BezCoeff_vec_t cb = ((HBezFn *)tb->func)->d->c;
        ta->face = Left(e)->num;
        tb->face = Left(Sym(e))->num;
        { int j; for (j = 0; j < ca.ne; j++) { ca.e[j] = 0.0; cb.e[j] = 0.0; } }
        if (k < NBVal)
          { int kv = k;
            ComputeEdgeC1ValElemCoeffs(kv, deg, e, ca, cb, tri, &frame);
          }
        else
          { int kd = k - NBVal;
            ComputeEdgeC1DerElemCoeffs(kd, deg, e, ca, cb, tri, &frame);
          }
        CheckPiecesForC1Continuity(tri, e, pd.e[0], pd.e[1]);
        StoreBasisElem(bas, nb, triFile, tri, pd, NULL);
      }
  }
  
void ComputeEdgeC1ValElemCoeffs
  ( int kv,
    int deg,
    Arc a,
    BezCoeff_vec_t ca,
    BezCoeff_vec_t cb,
    Triangulation *tri,
    r3x3_t *frame
  )
  {
    /* The val-type elements of index {k} for edge {a}
      has all coefficients set to zero, except 
      
        {ca_{0,deg-k,k}} (the element's cardinal value)
        {ca_{1,deg-k,k-1}, cb_{1,k-1,deg-k}}
        {ca_{1,deg-k-1,k}, cb_{1,k,deg-k-1}}
      
      The first element is set to 1, the other four are adjusted so as
      to obey the C1 constraints. The linear functionals are chosen so
      that their derivative transverse to the plane of edge {a} is zero. */
    
    Arc b = Sym(a);
    S2Point *va = &(Dest(Onext(a))->pos);
    S2Point *vb = &(Dest(Onext(b))->pos);
    int i = 0;
    int j = kv + 3;
    int k = deg - i - j;
    r3_t U = (r3_t){{frame->c[0][0], frame->c[1][0], frame->c[2][0]}}; 
    r3_t V = (r3_t){{frame->c[0][1], frame->c[1][1], frame->c[2][1]}}; 
    
    SetCoeffRelToEdge(tri, deg, ca, a,  i, j, k,  1.0);
    SetCoeffRelToEdge(tri, deg, cb, b,  i, k, j,  1.0);
                                     
    SetCoeffRelToEdge(tri, deg, ca, a,  i+1, j, k-1,  r3_dot(va, &V));
    SetCoeffRelToEdge(tri, deg, cb, b,  i+1, k-1, j,  r3_dot(vb, &V));
                                     
    SetCoeffRelToEdge(tri, deg, ca, a,  i+1, j-1, k,  r3_dot(va, &U));
    SetCoeffRelToEdge(tri, deg, cb, b,  i+1, k, j-1,  r3_dot(vb, &U));
  }

void ComputeEdgeC1DerElemCoeffs
  ( int kd,
    int deg,
    Arc a,
    BezCoeff_vec_t ca,
    BezCoeff_vec_t cb,
    Triangulation *tri,
    r3x3_t *frame
  )
  {
    /* The der-type elements of index {k} has all 
      coefficients set to zero, except 
      
        {a_{1,deg-k,k-1}} (the element's cardinal value)
        {b_{1,k-1,deg-k}}
      
      These two elements are adjusted so as to obey the C1
      constraints. The corresponding linear functional is chosen so
      that its value on the plane of edge {a} is zero, and
      its derivative transverse to that plane is 1. */
    
    Arc b = Sym(a);
    S2Point *va = &(Dest(Onext(a))->pos);
    S2Point *vb = &(Dest(Onext(b))->pos);
    int i = 1;
    int j = kd + 2;
    int k = deg - i - j;
    r3_t W = (r3_t){{frame->c[0][2], frame->c[1][2], frame->c[2][2]}}; 
    
    SetCoeffRelToEdge(tri, deg, ca, a,  i, j, k,  r3_dot(va, &W));
    SetCoeffRelToEdge(tri, deg, cb, b,  i, k, j,  r3_dot(vb, &W));
  }

r3x3_t VertexFrame(S2Point *u, S2Point *v)
  { r3_t r, s, t;
    r3x3_t m;
    int k;
    r3_t tmp;
    r = *u; 
    r3_dir(&r, &r);
    r3_decomp(v, &r, &tmp, &s);
    r3_dir(&s, &s);
    r3_cross(&r, &s, &t);
    r3_dir(&t, &t);
    for (k = 0; k < 3; k++) 
      { m.c[k][0] = r.c[k]; m.c[k][1] = s.c[k]; m.c[k][2] = t.c[k]; }
    return m;
  }

r3x3_t EdgeFrame(S2Point *u, S2Point *v)
  { r3_t r, s, t;
    r3x3_t m;
    int k;
    r3_t tmp;
    r3_decomp(u, v, &tmp, &r);
    r3_scale(1.0/r3_dot(u, &r), &r, &r);
    r3_decomp(v, u, &tmp, &s);
    r3_scale(1.0/r3_dot(v, &s), &s, &s);
    r3_cross(&r, &s, &t);
    r3_dir(&t, &t);
    for (k = 0; k < 3; k++) 
      { m.c[k][0] = r.c[k]; m.c[k][1] = s.c[k]; m.c[k][2] = t.c[k]; }
    return m;
  }

void BuildVertBasisC1New
  ( int deg,
    Arc e,
    Basis *bas,
    int *nb,
    char *triFile,
    Triangulation *tri
  )
  { /* In the new construction, the {m+3} elements associated to
      vertex {u} are divided in two subsets: three `star' elements,
      each one spanning the entire star of {u} ({m} triangles); and
      {m} `boat' elements, each of them spanning only three
      consecutive wedges (a {boat}) of the star of {u}.
      
      This construction only works when there are no two neighbors
      {p[i-1],p[i+1]} which lie on the same great circle as {u}.
      
      There is one boat element for each edge {b} out of {u}, spanning
      the triangle {Left(b)} and its neighbors {Left(a)} and
      {Left(c)}, where {a = Oprev(b)} and {c = Onext(b)}. The only
      non-zero coeffcients are {v_b[3] == v_a[5], v_b[4], v_b[5] ==
      v_c[3]}, {v_a[8]}, {v_b[7]}, {v_b[8]}, and {v_c[7]}. Coefficient
      {v_b[4]}, the cardinal coefficient of this element, is set to 1,
      and the rest is computed from the C1 continuity conditions.
      
      In the star elements, for each wedge {Left(b)} where {Org(b) = u},
      coefficients {v_b[4]} and {v_b[6]} are
      zero, for every edge {b} out of {u}. Coefficients {v_b[0]} and
      {v_b[1]} are set as follows, for each element {k = 0,1,2}:
      
        {v_b[0] = <o_b|r_k>},  (U1)
        {v_b[1] = <d_b|r_k>}   (U2)
        
      which, by C0, also imply
      
        {v_b[2] = <n_b|r_k>},   (V1)
        {v_a[1] = <p_b|r_k>}    (V2)
        
      where 
        
        {a = Oprev(b)}
        {c = Onext(b)}
        
        {d_b = Dest(b).pos},
        {p_b = Dest(Oprev(b)).pos}
        {n_b = Dest(Onext(b)).pos}
        {o_b = Org(b).pos}
        
        {r_0 = u.pos = o_b} (a unit vector),
        {r_1 = } some unit vector orthogonal to {r_0}
        {r_2 = cross(r_0, r_1)}
        
      Note that the C1 continuity conditions for coefficients
      {v_b[0..2]} are satisfied, because if we define {bn,bo,bd} by
      the equation
      
        {p_b == bn*n_b + bo*o_b + bd*d_b}, 
        
      then, for any vector {r_k}, we automatically have
      
        {<p_b|r_k> = bn*<n_b|r_k> + bo*<o_b|r_k> + bd*<d_b|r_k}  (S)
      
      Substituting (U1,U2,V1,V2) into (S), we get
      
        {v_a[1] = bn*v_b[2] + bo*v_b[0] + bd*v_b[1]}

      which is the diamond condition for those four coefficients.
      
      Once {v_b[0..2]} have been set, the C1 condition is used to
      compute coefficient {v_b[3]} from {v_b[1]}, assuming {v_b[4]}
      and {v_a[4]} are zero. Then, as in the boat elements, {v_b[7]}
      and {v_a[8]} are computed from {v_b[3]} by the C1 condition,
      taking care that the corresponding linear functional is
      orthogonal to the one used for the edge element based on those
      two coefficients. Ditto for {v_b[8]} and {v_c[7]}.
      
      The cardinal coefficients are {v_e[0]}, {v_e[1]}, and {v_e[2]},
      in the first wedge {Left(e)}. These are not strictly cardinal,
      but consider the matrix {M} where the entries in row {k} are the
      coefficients {v_e[0], v_e[1], v_e[2]} for element {k}. From
      formulas (U1,U2,V1), columns 0,1,2 of this matrix are the
      coordinates of {o_e}, {d_e}, and {n_e}, respectively, in the
      orthonormal basis {r_0,r_1,r_2}. Assuming that the triangles are
      not degenerate, this matrix has nonzero determinant; hence the
      three elements are independent among themselves, and independent
      from any other spline elements (which all have {v_e[0] == v_e[1]
      == v_e[2] == 0}). */
  
    int m = OrgDegree(e);
    S2Point *oe = &(Org(e)->pos);
    S2Point *de = &(Dest(e)->pos);
      
    int k;
    double minBeta = 1.0e-3; /* Lower limit for non-colinearity coeff. */

    auto bool_t ArcIsEdgeReference(Arc a);
      /* TRUE if the basis elements for the edge {e.quad} were built 
        using {a} (rather than {Sym(a)}) as the reference arc. */

    auto PieceData *MakePiece(Arc a, double *va);
      /* Creates a {PieceData} record for the face {Left(a)},
        with {v[0..9]} as the Bezier coefficients nearest
        to {Org(a)}> */
        
    auto double ComputeV3(Arc b, double *vb);
      /* Assumes that {vb[0..9]} are the corner coefficients for the
        face {Left(b)} near {Org(b)}.  Computes {vb[3]}, assuming that
        {vb[1]} and {vb[4]} have already been set, and {va[4]} is zero. */
        
    auto double ComputeV5(Arc b, double *vb);
      /* Assumes that {vb[0..9]} are the corner coefficients for the
        face {Left(b)} near {Org(b)}.  Computes {vb[5]}, assuming that
        {vb[2]} and {vb[4]} have already been set, and {vc[4]} is zero. */
        
    auto double ComputeV7(Arc b, double *vb);
      /* Assumes that {vb[0..9]} are the corner coefficients for the
        face {Left(b)} near {Org(b)}.  Computes {vb[7]}, assuming that
        {vb[3]} has already been set, and {vb[6]} is zero. */
    
    auto double ComputeV8(Arc b, double *vb);
      /* Assumes that {vb[0..9]} are the corner coefficients for the
        face {Left(b)} near {Org(b)}.  Computes {vb[8]}, assuming that
        {vb[5]} has already been set, and {vb[9]} is zero. */
    
    PieceData *MakePiece(Arc a, double *va)
      { PieceData *ta = MakeHBezPiece(deg);
        BezCoeff_vec_t ca = ((HBezFn *)ta->func)->d->c;
        ta->face = Left(a)->num;
        SetCoeffsVertex(a, tri, deg, ca, va, 10);
        return ta;
      }
      
    bool_t ArcIsEdgeReference(Arc a)
      { return tri->arc.e[2*EdgeNum(a)] == a; }

    double ComputeV7(Arc b, double *vb)
      { Arc c = Onext(b);
        S2Point *db = &(Dest(b)->pos);
        S2Point *dc = &(Dest(c)->pos);
        r3_t r, tmp;
        r3_decomp(oe, db, &tmp, &r);
        r3_scale(1.0/r3_dot(oe, &r), &r, &r);
        affirm(vb[6] == 0.0, "vb[6] is not zero");
        return vb[3]*r3_dot(dc, &r);
      }

    double ComputeV8(Arc b, double *vb)
      { Arc c = Onext(b);
        S2Point *db = &(Dest(b)->pos);
        S2Point *dc = &(Dest(c)->pos);
        r3_t r, tmp;
        r3_decomp(oe, dc, &tmp, &r);
        r3_scale(1.0/r3_dot(oe, &r), &r, &r);
        affirm(vb[9] == 0.0, "vb[9] is not zero");
        return vb[5]*r3_dot(db, &r);
      }
      
    double ComputeV3(Arc b, double *vb)
      { r3_t bar = ComputeEdgeContinuityCoeffs(b, tri);
        double alpha = fabs(bar.c[1]), beta = fabs(bar.c[2]);
        affirm(beta >= minBeta*(alpha + beta), "near-collinear edges"); 
        return (- bar.c[0]*vb[4] - bar.c[1]*vb[1])/bar.c[2];
      }
      
    double ComputeV5(Arc b, double *vb)
      { Arc c = Onext(b);
        r3_t bar = ComputeEdgeContinuityCoeffs(c, tri);
        double alpha = fabs(bar.c[1]), beta = fabs(bar.c[2]);
        affirm(beta >= minBeta*(alpha + beta), "near-collinear edges"); 
        return (vb[4] - bar.c[1]*vb[2])/bar.c[2];
      }

    /* The ANS construction fails for {cont == 1} and {deg < 5}: */
    affirm(deg >= 5, "invalid deg/cont combination");

    /* Generate the `boat' elements: */
    { Arc b = e;
      /* Relevant Bezier coefficients for the three wedges. */
      double va[10]; /* Previous wedge */
      double vb[10]; /* Center wedge. */
      double vc[10]; /* Next wedge. */
        
      for (k = 0; k < m; k++)
        { PieceDataRef_vec_t pd = PieceDataRef_vec_new(3);
          Arc a = Oprev(b);
          Arc c = Onext(b);
          int i;
          /* Initialize coeffs {v(0)[0..9]}: */
          for (i = 0; i < 10; i++)
            { va[i] = 0.0; vb[i] = 0.0; vc[i] = 0.0; }
          /* The cardinal Gamma point of the boat element is {vb[4]}: */
          vb[4] = 1.0;
          vb[3] = ComputeV3(b, vb); va[5] = vb[3];
          vb[5] = ComputeV5(b, vb); vc[3] = vb[5];
          /* Fix coefficients {va[8],vb[7],vb[8],vc[7]}: */
          va[8] = ComputeV8(a, va);
          vb[7] = ComputeV7(b, vb);
          vb[8] = ComputeV8(b, vb);
          vc[7] = ComputeV7(c, vc);
          pd.e[0] = MakePiece(a, va);
          pd.e[1] = MakePiece(b, vb);
          pd.e[2] = MakePiece(c, vc);
          CheckPiecesForC1Continuity(tri, b, pd.e[1], pd.e[0]);
          CheckPiecesForC1Continuity(tri, c, pd.e[2], pd.e[1]);
          StoreBasisElem(bas, nb, triFile, tri, pd, oe);
          b = Onext(b);
        }
    }
    
    /* Generate the `star' elements */
    { int k;
      r3x3_t frame = VertexFrame(oe, de);
      for (k = 0; k < 3; k++)
        { /* Generates the {k}th star element. */
          PieceDataRef_vec_t pd = PieceDataRef_vec_new(m);
          double va[10]; /* Bezier coeffs {v(r)[0..9]} for current wedge. */
          Arc a = e;
          S2Point *na = &(Dest(Onext(a))->pos);
          S2Point *da = &(Dest(a)->pos);
          r3_t Rk = (r3_t){{frame.c[0][k], frame.c[1][k], frame.c[2][k]}}; 
          int r;
          
          /* Set wedges around {Org(e)}: */
          for (r = 0; r < m; r++)
            { Arc b = Onext(a);
              /* At this point coeffs {va[0,1,3,6]} are set. */
              /* Set coeffs {va[2,4,5,7,8,9]}: */
              va[0] = r3_dot(oe, &Rk);
              va[1] = r3_dot(da, &Rk);
              va[2] = r3_dot(na, &Rk);
              va[4] = 0.0;
              va[6] = 0.0;
              va[9] = 0.0;
              va[3] = ComputeV3(a, va);
              va[5] = ComputeV5(a, va);
              va[7] = ComputeV7(a, va);
              va[8] = ComputeV8(a, va);
              pd.e[r] = MakePiece(a, va);
              
              /* Advance to next wedge: */
              a = b;
              da = na;
              na = &(Dest(Onext(a))->pos);
            }
            
          /* Everybody is conspiring to make me think I'm paranoid: */
          a = e;
          for (r = 0; r < m; r++)
            { CheckPiecesForC1Continuity(tri, a, pd.e[r], pd.e[(r + m - 1) % m]);
              a = Onext(a); 
            }

          /* Store the element: */
          StoreBasisElem(bas, nb, triFile, tri, pd, oe);
        }
    }
  }

