/* See SPBezSplineBasis.h */
/* Last edited on 2005-10-27 17:08:40 by stolfi */

#include <SPBezSplineBasis.h>

#include <SPFunction.h>
#include <SPSpline.h>
#include <SPHBezFunction.h>
#include <SPNHBezFunction.h>
#include <SPTriang.h>
#include <SPQuad.h>
#include <SPDeCasteljau.h>
#include <SPBezSplineBasisC0.h>
#include <SPBezSplineBasisC1.h>
#include <SPBasic.h>

#include <vec.h>
#include <affirm.h>
#include <nat.h>

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <values.h>
#include <math.h>

Basis SPBezSplineBasis_BuildH
  ( int deg,             /* Degree of polynomials */
    int cont,            /* Continuity. */
    bool_t newStyle,       /* FALSE: old ANS, TRUE: with boat elements. */
    char *triFile,       /* Triangulation file name. */
    Triangulation *tri   /* Triangulation. */
  )
  { int nb = 0;
    int NV = tri->out.ne;
    int NT = tri->side.ne;
    int NE = tri->arc.ne / 2;
    /* The basis: */
    int NB = HSpaceDimension(deg, cont, NT);
    Basis bas = SPFunctionRef_vec_new(NB);
    /* Index tables: */
    nat_vec_t fBas = nat_vec_new(NT+1);
    nat_vec_t eBas = nat_vec_new(NE+1);
    nat_vec_t vBas = nat_vec_new(NV+1);
    /* The basis elements associated with face number {fn} are
      bas[fBas[fn]..fBas[fn+1]-1]. */
    BuildFaceBasis(&bas, &nb, triFile, tri, deg, cont, fBas);
    BuildEdgeBasis(&bas, &nb, triFile, tri, deg, cont, newStyle, eBas);
    BuildVertBasis(&bas, &nb, triFile, tri, deg, cont, newStyle, vBas);
    affirm(nb == NB, "inconsistent element count");
    affirm(bas.ne == NB, "inconsistent basis size");
    affirm(vBas.e[NV] == NB, "inconsistent index tables");
    return bas;
    
  }

Basis SPBezSplineBasis_BuildNH
  ( int deg,             /* Degree of polynomials */
    int cont,            /* Continuity. */
    bool_t newStyle,       /* FALSE: old ANS, TRUE: with boat elements. */
    char *triFile,       /* Triangulation file name. */
    Triangulation *tri   /* Triangulation. */
  )
  { 
    auto r3_t SuppNormal(Spline *fpw);
      /* Returns a vector normal to the supporting plane of {fpw}.
         May be a zero vector. */
    
    r3_t SuppNormal(Spline *fpw)
      { r4_t supp = fpw->d->supp.f;
        return (r3_t){{supp.c[1], supp.c[2], supp.c[3]}};
      }
  
    Basis hb0 = SPBezSplineBasis_BuildH(deg,   cont, newStyle, triFile, tri);
    int dim0 = hb0.ne;
    Basis hb1 = SPBezSplineBasis_BuildH(deg-1, cont, newStyle, triFile, tri);
    int dim1 = hb1.ne;
    Basis nhb = SPFunctionRef_vec_new(dim0 + dim1);
    int i, j;
    
    for (i = 0; i < hb0.ne;  i++)
      { Spline *hb0i = (Spline *)(hb0.e[i]);
        PieceDataRef_vec_t hb0pd = hb0i->d->pd;
        r3_t norm = SuppNormal(hb0i);
        Spline *nhbi;
        PieceDataRef_vec_t nhpd  = MakeNHBezPieces(hb0pd.ne, deg);
        for (j = 0; j < hb0pd.ne; j++)
          { PieceData *hb0pdj = hb0pd.e[j];
            PieceData *nhpdj = nhpd.e[j];
            nhpdj->face = hb0pdj->face;
            nhpdj->bary = hb0pdj->bary;
            affirm(hb0pdj->bary, "spline is not barycentric");
            { int k;
              SPHBezFunction *hb0pf = (HBezFn *)(hb0pdj->func);
              SPNHBezFunction *nhbpf = (NHBezFn *)(nhpdj->func);
              BezCoeff_vec_t hb0c = hb0pf->d->c;
              BezCoeff_vec_t nhc0 = nhbpf->d->c0;
              for (k = 0; k < nhc0.ne; k++) { nhc0.e[k] = hb0c.e[k]; }
            }
          }
        nhbi = SPSpline_FromPieces(triFile, tri, nhpd, &norm);
        nhb.e[i] = (SPFunction *)nhbi;
      }
      
    for (i = 0;  i < hb1.ne; i++)
      { Spline *hb1i = (Spline *)(hb1.e[i]);
        PieceDataRef_vec_t hb1pd = hb1i->d->pd;
        r3_t norm = SuppNormal(hb1i);
        Spline *nhbi;
        PieceDataRef_vec_t nhpd  = MakeNHBezPieces(hb1pd.ne, deg);
        for (j = 0; j < hb1pd.ne; j++)
          { PieceData *hb1pdj = hb1pd.e[j];
            PieceData *nhpdj = nhpd.e[j];
            nhpdj->face = hb1pdj->face;
            nhpdj->bary = hb1pdj->bary;
            affirm(hb1pdj->bary, "spline is not barycentric");
            { int k;
              SPHBezFunction *hb1pf = (HBezFn *)(hb1pdj->func);
              SPNHBezFunction *nhbpf = (NHBezFn *)(nhpdj->func);
              BezCoeff_vec_t hb1c = hb1pf->d->c;
              BezCoeff_vec_t nhc1 = nhbpf->d->c1;
              for (k = 0; k < nhc1.ne; k++) { nhc1.e[k] = hb1c.e[k]; }
            }
          }
        nhbi = SPSpline_FromPieces(triFile, tri, nhpd, &norm);
        nhb.e[dim0 + i] = (SPFunction *)nhbi;
      }
    /* !!! free(hb0.e); free(hb1.e); !!! */
    return nhb;
  }

int HSpaceDimension(int deg, int cont, int NT)
  { affirm(NT % 2 == 0, "number of triangles is not even");
    affirm(cont >= -1, "invalid continuity");
    if (deg < 0) 
      { return 0; }
    else if (cont == -1)
      { return ((deg + 1)*(deg + 2) / 2) * NT; }
    else if (cont == 0)
      { affirm(deg >= 1, "invalid deg/cont combination");
        return deg*deg * NT / 2 + 2;
      }
    else if (cont == 1) 
      { affirm(deg >= 5, "invalid deg/cont combination");
        return (deg*deg - 3*deg + 2) * NT / 2 + 6;
      }
    else
      { affirm(FALSE, "invalid continuity");
        return 0;
      }
  }

void BuildFaceBasis
  ( Basis *bas,
    int *nb,   
    char *triFile,
    Triangulation *tri,
    int deg,
    int cont,
    nat_vec_t fBas 
  )
  { Arc_vec_t side = tri->side;
    int NF = side.ne;
    int minexp = cont+1;
    int fn, i, j;

    for (fn = 0; fn < NF; fn++)
      { fBas.e[fn] = *nb;
        if (deg >= minexp)
          { Arc e = side.e[fn];
            { nat_t en0 = Edge(e)->num;
              nat_t en1 = Edge(Lnext(e))->num;
              nat_t en2 = Edge(Lprev(e))->num;
              nat_t vn0 = Org(e)->num;
              nat_t vn1 = Org(Lnext(e))->num;
              nat_t vn2 = Org(Lprev(e))->num;
              fprintf(stderr, 
                "nb = %04d face %d  edges %d %d %d  vertices %d %d %d\n",
                (*nb), fn, en0, en1, en2, vn0, vn1, vn2
              );
            }
            affirm(fn == Left(e)->num, "inconsistent face number");
            for (i = minexp; i <= deg - 2*minexp; i++)
              { for (j = minexp; j <= deg - i - minexp; j++)
                  { int k = deg - i - j;
                    BezLabel ijk = (BezLabel){{i,j,k}};
                    int indr = SPDeCasteljau_BezLabelToIndex(ijk, deg);
                    PieceDataRef_vec_t pd = MakeHBezPieces(1, deg);
                    affirm(i+j+k == deg, "wrong tot degree");
                    affirm(k >= minexp, "bad k");
                    fprintf(stderr, "  face element %04d (%d,%d,%d)\n", *nb, i, j, k);
                    { PieceData *t0 = pd.e[0];
                      BezCoeff_vec_t c0 = ((HBezFn *)t0->func)->d->c;
                      t0->face = Left(e)->num;
                      c0.e[indr] = 1.0;
                    }
                    StoreBasisElem(bas, nb, triFile, tri, pd, NULL);
                  }
              }
          }
      }
    fBas.e[NF] = *nb;
  }

void SetCoeffRelToEdge
  ( Triangulation *tri,  /* Triangulation. */
    int deg,             /* Degree of polynomial. */
    BezCoeff_vec_t c,    /* Coeffs vector for some face. */
    Arc e,               /* An edge of the face in question. */
    int i, int j, int k, /* Bezier coeff label relative to {e}. */
    double v             /* Coefficient value. */
  )
  { BezLabel ijk;
    nat_t face = Left(e)->num;
    Arc a = tri->side.e[face];
    int index;
    
    if (e == a)
      { ijk = (BezLabel){{i, j, k}}; }
    else if (Lnext(a) == e)
      { ijk = (BezLabel){{k, i, j}}; }
    else if (Lnext(e) == a)
      { ijk = (BezLabel){{j, k, i}}; }
    else
      { affirm(FALSE , "face is not a triangle"); }
      
    index = SPDeCasteljau_BezLabelToIndex(ijk, deg);
    c.e[index] = v;
  }

void BuildEdgeBasis
  ( Basis *bas,
    int *nb,
    char *triFile,
    Triangulation *tri,
    int deg,
    int cont,
    bool_t newStyle,
    nat_vec_t eBas
  )
  { Arc_vec_t arc = tri->arc;
    int NE = arc.ne/2;
    int en;

    for (en = 0; en < NE; en++)
      { eBas.e[en] = (*nb);
        if ((deg > 0) && (cont >= 0) && (deg >= 3*cont + 2))
          { Arc e = arc.e[2*en];
            { nat_t vn0 = Org(e)->num;
              nat_t vn1 = Org(Sym(e))->num;
              fprintf(stderr, "nb = %04d", (*nb));
              fprintf(stderr, " arc %d:%d", en, QuadBits(e));
              fprintf(stderr, " vertices %d %d\n", vn0, vn1);
            }
            if (cont == -1)
              { /* Nothing to do. */ }
            else if (cont == 0)
              { BuildEdgeBasisC0(deg, e, bas, nb, triFile, tri); }
            else if (cont == 1)
              { BuildEdgeBasisC1(deg, newStyle, e, bas, nb, triFile, tri); }
            else
              { affirm(FALSE , "bad continuity"); }
          }
      }
    eBas.e[NE] = (*nb);
  }

void BuildVertBasis
  ( Basis *bas,
    int *nb,
    char *triFile,
    Triangulation *tri,
    int deg,
    int cont,
    bool_t newStyle,
    nat_vec_t vBas
  )
  { Arc_vec_t out = tri->out;
    int NV = out.ne;
    int vn;

    for (vn = 0; vn < NV; vn++)
      { fprintf(stderr, "nb = %04d vertex %d\n", (*nb), vn);
        { Arc e = out.e[vn];
          vBas.e[vn] = (*nb);
          if ((cont == -1) || (deg < 0))
            { /* Nothing to do. */ }
          else if (cont == 0)
            { BuildVertBasisC0(deg, e, bas, nb, triFile, tri); }
          else if (cont == 1)
            { BuildVertBasisC1(deg, newStyle, e, bas, nb, triFile, tri); }
          else
            { affirm(FALSE, "invalid continuity"); }
        }
      }
    vBas.e[NV] = (*nb);
  }

bool_t AllCoeffsAreZero(double *v, int nv)
  { int i;
    for (i = 0; i < nv; i++) { if (v[i] != 0.0) { return FALSE; } }
    return TRUE;
  }

void SetCoeffsVertex
  ( Arc e,                 /* Reference edge for basis element. */
    Triangulation *tri,    /* Triangulation. */
    int deg,               /* Degree of polynomial. */
    BezCoeff_vec_t c,      /* Bezier coefficients. */
    double *v,             /* Values to set. */
    int nv                 /* Number of values to set. */
  )
  { int indr;
    BezLabel ijk;
    nat_t fn = Left(e)->num;
    Arc fref = tri->side.e[fn];
    int i;
    if (fref == Lnext(e))
      { for (i = 0; i < nv; i++) { c.e[i] = v[i]; } }
    else if (e == Lnext(fref))
      { for (i = 0; i < nv; i++)
          { ijk = SPDeCasteljau_IndexToBezLabel(i,deg);
            ijk = (BezLabel){{ijk.e[1], ijk.e[2], ijk.e[0]}};
            indr = SPDeCasteljau_BezLabelToIndex(ijk,deg);
            c.e[indr] = v[i];
          }
      }
    else if (fref == e)
      { for (i = 0; i < nv; i++)
          { ijk = SPDeCasteljau_IndexToBezLabel(i,deg);
            ijk = (BezLabel){{ijk.e[2], ijk.e[0], ijk.e[1]}};
            indr = SPDeCasteljau_BezLabelToIndex(ijk,deg);
            c.e[indr] = v[i];
          }
      }
  }

PieceDataRef_vec_t MakeHBezPieces(int n, int deg)
  { PieceDataRef_vec_t pd = PieceDataRef_vec_new(n);
    int i;
    for (i = 0; i < n; i++)
      { PieceData *pi = MakeHBezPiece(deg);
        pd.e[i] = pi;
      }
    return pd;
  }
  
PieceData *MakeHBezPiece(int deg)
  { int NC = SPHBezFunction_NumCoeffs(deg);
    PieceData *pi = SPSpline_PieceData_New();
    BezCoeff_vec_t c = BezCoeff_vec_new(NC);
    int k;
    pi->face = INT_MAX;
    pi->bary = TRUE;
    for (k = 0; k < NC; k++) { c.e[k] = 0.0; }
    pi->func = (SPFunction *)SPHBezFunction_FromCoeffs(deg, c);
    return pi;
  }

PieceDataRef_vec_t MakeNHBezPieces(int n, int deg)
  {
    int NC0 = SPHBezFunction_NumCoeffs(deg);
    int NC1 = SPHBezFunction_NumCoeffs(deg-1);
    PieceDataRef_vec_t pd = PieceDataRef_vec_new(n);
    int i, k;

    for (i = 0; i < n; i++)
      { PieceData *pi = SPSpline_PieceData_New();
        BezCoeff_vec_t c0 = BezCoeff_vec_new(NC0);
        BezCoeff_vec_t c1 = BezCoeff_vec_new(NC1);
        pi->face = INT_MAX;
        pi->bary = TRUE;
        for (k = 0; k < NC0; k++) { c0.e[k] = 0.0; }
        for (k = 0; k < NC1; k++) { c1.e[k] = 0.0; }
        pi->func = (SPFunction *)SPNHBezFunction_FromCoeffs(deg, c0, c1);
        pd.e[i] = pi;
      }
    return pd;
  }

void StoreBasisElem
  ( Basis *bas, 
    int *nb, 
    char *triFile, Triangulation *tri, 
    PieceDataRef_vec_t pd,
    r3_t *norm
  )
  { SPSpline *fpw = SPSpline_FromPieces(triFile, tri, pd, norm);
    SPSpline_SortPieces(fpw->d->pd);
    SPFunctionRef_vec_expand(bas, *nb);
    bas->e[*nb] = (SPFunction *)fpw;
    (*nb)++;
  }
