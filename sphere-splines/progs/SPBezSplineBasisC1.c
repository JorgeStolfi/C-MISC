/* See SPBezSplineBasisC1.h */
/* Last edited on 2005-06-06 11:49:10 by stolfi */

#include <SPBezSplineBasisC1.h>

#include <SPFunction.h>
#include <SPSpline.h>
#include <SPHBezFunction.h>
#include <SPTriang.h>
#include <SPQuad.h>
#include <SPDeCasteljau.h>
#include <SPBezSplineBasis.h>
#include <SPBezSplineBasisC1Old.h>
#include <SPBezSplineBasisC1New.h>
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
  
void BuildEdgeBasisC1
  ( int deg,
    bool_t newStyle,
    Arc e,
    Basis *bas,
    int *nb,
    char *triFile,
    Triangulation *tri
  )
  { 
    /* Provided that {deg >= 5}, and barring gegenerate cases, there
      are {2*deg - 9} C1 edge elements for each edge. Roughly
      speaking, {deg - 5} of these elements control the value of the
      function along the edge (`val-type' elements), and {deg - 4} control
      the derivative transverse to the edge (`der-type' elements).
      
      These support of all these edge elements consists of the
      two triangles {Left(e)} and {Right(e)} adjacent to {e}.
      Let {a_{i,j,k}} be the Bezier coefficients associated with the
      face {Left(e)}, labeled so that {a_{deg,0,0}} is associated with
      the corner opposite to {e}, and {a_{0,deg,0}} is associated with
      {Org(e)}. Let {b_{i,j,k}} be the Bezier coefficients of
      {Right(e)}, with {b_{deg,0,0}} associated to the corner opposite
      to {e}, and {b_{0,deg,0}} associated to {Dest(e)}.
      
      The C0 continuity conditions are equivalent to {a_{0,i,deg-i} ==
      b_{0,deg-1,i}} for {i=0..deg}. The C1 continuity conditions
      require, for each {i=0..deg-1}, that the four values
      
        {a_{0,deg-i,i}, a_{0,deg-i-1,i+1}, a_{1,deg-i-1,i}, b_{1,i,deg-i-1}}
      
      coincide with some linear finctional applied to the points
      
        {Org(e), Dest(e), Dest(Onext(e)), Dest(Oprev(e))}
        
      respectively.  */
      
    if (newStyle)
      { BuildEdgeBasisC1New(deg, e, bas, nb, triFile, tri); }
    else
      { BuildEdgeBasisC1Old(deg, e, bas, nb, triFile, tri); }
  }

void BuildVertBasisC1
  ( int deg,
    bool_t newStyle,
    Arc e,
    Basis *bas,
    int *nb,
    char *triFile,
    Triangulation *tri
  )
  { /* Provided {deg >= 5}, there are {m+3} C1 elements associated to
      the vertex {u = Org(e)}, where {m = OrgDegree(u)}.  The triangles covered
      by these elements are all incident to the vertex {u}
      
      In all pieces of these elements, all Bezier coefficients are
      zero, except for the 10 coefficients closest to the vertex {u},
      which are set to various values depending on the element and
      piece.
      
      A triangle {f} incident to {u} can be represented by the arc {b}
      such that {Org(b) == u} and {Left(b) == f}. Let the 10 relevant
      coefficients in wedge {b}, contained in face {Left(b)}, be
      denoted by {v_b[0..9]}, or just {v[0..9]} when {b} is implied
      by context. Coefficient {v[0]} lies at the corner {u = Org(b)}.
      Coefficients {v[1],v[3],v[6]} lie on the edge {b}; coefficients
      {v[2],v[5],v[9]} lie on the edge {Onext(b)}. Actually,
      coefficients {v[6]} and {v[9]} are always 0.
      
      Zero-order continuity is equivalent to ensuring that all control
      values on the periphery of the star of {u} are zero, and all
      control values along any internal edge are shared by the two
      adjacent faces. Assuming {deg >= 5}, these conditions reduce to
      requiring that, for any edge {b} out of {u},
      
        {v_b[0] == v_a[0]}
        {v_b[1] == v_a[2]}
        {v_b[3] == v_a[5]}, and
        {v_b[6] == v_a[9]},
        
      where {a = Oprev(b)}.
      
      First-order continuity is equivalent to imposing certain linear
      equations (the `diamond conditions') on control values that lie
      on the edges of the star or at most 1 step away from them. Assuming
      {deg >= 5}, these conditions reduce to the requirement that,
      for any edge {b} out of {u},
      
        {diamond(b, v_b[2], v_b[0], v_b[1], v_a[1])}   (D1)
        {diamond(b, v_b[4], v_b[1], v_b[3], v_a[4])}   (D2)
        {diamond(b, v_b[7], v_b[3], v_b[6], v_a[8])}   (D3)
        
      Here {a == Oprev(b)}, and {diamond(e, vn, vo, vd, vp)} means the
      equation {vp = bn*vn + bo*vo + bd*vd}, where {bn,bo,bd} are such
      that
      
        {Dest(Oprev(e)) == bn*Dest(Onext(e)) + bo*Org(b) + bd*Dest(b)}
      
      Alefeld, Neamtu and Schumaker have proved that these constraints
      leave {m+3} degrees of freedom --- except in degenerate
      situations when there are certain collinearities between the
      edges incident to {u}, in which case there may be a few extra 
      degrees of freedom. */
      
    if (newStyle)
      { BuildVertBasisC1New(deg, e, bas, nb, triFile, tri); }
    else
      { BuildVertBasisC1Old(deg, e, bas, nb, triFile, tri); }
  }

r3_t ComputeEdgeContinuityCoeffs(Arc e, Triangulation *tri)
  { S2Point *vp = &(Dest(Oprev(e))->pos);
    Face *f = Left(e);
    Arc fref = tri->side.e[f->num];
     r3_t bar;

    r3x3_map_col(&(f->c2b), vp, &bar);
    if (e == fref)
      { return bar; }
    else if (e == Lnext(fref))
      { return (r3_t){{bar.c[1], bar.c[2], bar.c[0]}}; }
    else if (fref == Lnext(e))
      { return (r3_t){{bar.c[2], bar.c[0], bar.c[1]}}; }
    else
      { affirm(FALSE , "face is not a triangle"); return bar; }
  }

void CheckPiecesForC1Continuity
  ( Triangulation *tri, 
    Arc a, 
    PieceData *pa, 
    PieceData *pb
  )
  {
    affirm(pa->bary && pb->bary, "non-barycentric pieces");
    { SPHBezFunction *ha = SPHBezFunction_Cast(pa->func); 
      SPHBezFunction *hb = SPHBezFunction_Cast(pb->func); 
      if ((ha != NULL) && (hb != NULL) && (ha->d->deg == hb->d->deg))
        { CheckCoeffsForC1Continuity(tri, a, ha->d->deg, ha->d->c, hb->d->c); 
          return;
        }
      else if ((ha != NULL) || (hb != NULL))
        { affirm(FALSE, "incompatible piece types"); }
      
    }
    { SPNHBezFunction *nha = SPNHBezFunction_Cast(pa->func); 
      SPNHBezFunction *nhb = SPNHBezFunction_Cast(pb->func); 
      if ((nha != NULL) && (nhb != NULL) && (nha->d->deg == nhb->d->deg))
        { /* By the theorem of Gomide and Stolfi, 1998: */
          CheckCoeffsForC1Continuity(tri, a, nha->d->deg, nha->d->c0, nhb->d->c0);
          CheckCoeffsForC1Continuity(tri, a, nha->d->deg-1, nha->d->c1, nhb->d->c1);
          return;
        }
      else if ((nha != NULL) || (nhb != NULL))
        { affirm(FALSE, "incompatible piece types"); }
    }
    affirm(FALSE, "unrecognized piece types");
  }

void CheckCoeffsForC1Continuity
  ( Triangulation *tri, 
    Arc a, 
    int deg,
    BezCoeff_vec_t ca, 
    BezCoeff_vec_t cb
  )
  { Arc b = Sym(a);
    
    auto double get_coeff(Arc e, BezCoeff_vec_t c, int i, int j, int k);
      /* Obtains from the Bezier coefficient vector {c} the entry
        whose label relative to arc {e} is {(i,j,k)}. */
    
    auto void C0_failure
      ( double ca, int ia, int ja, int ka, 
        double cb, int ib, int jb, int kb
      );

    auto void C1_failure(double cxp, double cbp, int ib, int jb, int kb);

    auto void dump_coeffs(void);
    
    double get_coeff(Arc e, BezCoeff_vec_t c, int i, int j, int k)
      { /* Note that the coefficients are numbered relative to 
          the reference edge of each face, not relative to {e}. */
        Arc fref = tri->side.e[Left(e)->num];
        BezLabel lab;
        if (e == fref)
          { lab = (BezLabel){{i,j,k}}; }
        else if (e == Lnext(fref))
          { lab = (BezLabel){{k,i,j}}; }
        else if (e == Lprev(fref))
          { lab = (BezLabel){{j,k,i}}; }
        else
          { affirm(FALSE, "inconsistent side"); }
        return c.e[SPDeCasteljau_BezLabelToIndex(lab, deg)];
      }

    void C0_failure
      ( double ca, int ia, int ja, int ka, 
        double cb, int ib, int jb, int kb
      )
      { fprintf(stderr, "C0 failure: ca(%d,%d,%d) != cb(%d,%d,%d)\n",
          ia, ja, ka, ib, jb, kb
        );
        dump_coeffs();
        affirm(FALSE, "aborted"); 
      }

    void C1_failure(double cxp, double cbp, int ib, int jb, int kb)
      { fprintf(stderr, "C1 failure: cb(%d,%d,%d) == %16.12e != %16.12e\n",
          ib, jb, kb, cbp, cxp
        );
        dump_coeffs();
        affirm(FALSE, "aborted"); 
      }

    void dump_coeffs(void)
      { int i, j, k;
        for (i = deg; i >= 0; i--)
          { for (j = deg; j > deg-i; j--) { fprintf(stderr, "    "); }
            for (j = deg-i; j >= 0; j--)
              { k = deg - i - j; 
                fprintf(stderr, " %7.4f", get_coeff(a, ca, i,j,k)); 
              }
            fprintf(stderr, "\n");
          }
        fprintf(stderr, "\n");
        for (i = 0; i <= deg; i++)
          { for (j = deg-i+1; j <=deg; j++) { fprintf(stderr, "    "); }
            for (j = 0; j <= deg-i; j++)
              { k = deg - i - j; 
                fprintf(stderr, " %7.4f", get_coeff(b, cb, i,j,k)); 
              }
            fprintf(stderr, "\n");
          }
      }
  
    r3_t bar = ComputeEdgeContinuityCoeffs(a, tri);
    int j;
    /* Check C0 continuity first, just in case: */
    for (j = 1; j <= deg; j++)
      { double cao = get_coeff(a, ca, 0, deg-j+1, j-1);
        double cad = get_coeff(a, ca, 0, deg-j, j);

        double cbo = get_coeff(b, cb, 0, j-1, deg-j+1);
        double cbd = get_coeff(b, cb, 0, j, deg-j);

        if (cao != cbo) { C0_failure(cao, 0, deg-j+1, j-1,  cbo, 0, j-1, deg-j+1); } 
        if (cad != cbd) { C0_failure(cad, 0, deg-j, j,      cbd, 0, j, deg-j); } 
      }
    /* Now check C1 continuity: */
    for (j = 1; j <= deg; j++)
      { double can = get_coeff(a, ca, 1, deg-j, j-1);
        double cao = get_coeff(a, ca, 0, deg-j+1, j-1);
        double cad = get_coeff(a, ca, 0, deg-j, j);
        
        double cbp = get_coeff(b, cb, 1, j-1, deg-j);
        double cxp = bar.c[0]*can + bar.c[1]*cao + bar.c[2]*cad;
        double mag = fabs(can) + fabs(cao) + fabs(cad) + 1.0e-30;
        if (fabs(cxp - cbp)/mag > 1.0e-13)
          { C1_failure(cxp, cbp, 1, j, deg-j-1); }
      }
  }
     
