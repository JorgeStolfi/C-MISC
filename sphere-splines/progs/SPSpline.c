/* See SPSpline.h */
/* Last edited on 2024-11-11 05:13:40 by stolfi */

#include <SPFunction.h>
#include <SPSpline.h>
#include <SPIntegral.h>
#include <SPTriang.h>
#include <SPQuad.h>
#include <SPH3.h>

#include <affirm.h>
#include <nat.h>
#include <bool.h>
#include <filefmt.h>
#include <nget.h>
#include <fget.h>
#include <r3.h>
#include <r4.h>
/* #include <rn.h> DEBUG */
#include <r3x3.h>

#include <math.h>
#include <assert.h>
#include <limits.h>
#include <values.h>
#include <stdio.h>
#include <stdlib.h>

SPSpline_Methods *SPSpline_Methods_New(void);
SPSpline_Data *SPSpline_Data_New(void);
SPSpline *SPSpline_New(void);

SPSpline *SPSpline_FullNew(void);
  /* Allocates a new {SPSpline} object {f} (and its data record),
    with the proper {type} fields and methods record.
    The array {f->d->data} is left undefined. */

void r3x3_trmul (r3x3_t *A, r3x3_t *B, r3x3_t *R);
  /* Sets {R} to the product {A^T·B} of the tanspose of {A} and {B}. */

SPSpline *SPSpline_Cast(OBJ *f)
  { SPFunction *ff = (SPFunction *)f;
    if ((f != NULL) && isprefix(SPSpline_TypeId, ff->type))
      { return (SPSpline *)f; }
    else
      { return NULL; }
  }

bool_t SPSpline_AreCompatible(SPFunction *f, SPFunction *g, Triangulation *tri)
  { SPSpline *fpw;
    SPSpline *gpw;
    if ((fpw = SPSpline_Cast(f)) == NULL) { return FALSE; }
    if ((gpw = SPSpline_Cast(g)) == NULL) { return FALSE; }
    if (gpw->d->tri != fpw->d->tri) { return FALSE; }
    if ((tri != NULL) && (tri != gpw->d->tri)) { return FALSE; }
    return TRUE;
  }

bool_t SPSpline_IsCompatible(SPFunction *f, Triangulation *tri)
  { SPSpline *fpw;
    if ((fpw = SPSpline_Cast(f)) == NULL) { return FALSE; }
    if ((tri != NULL) && (tri != fpw->d->tri)) { return FALSE; }
    return TRUE;
  }

#define T SPSpline

double SPSpline_M_Eval(T *f, R3Point *p);
R3Gradient SPSpline_M_Grad(T *f, R3Point *p);
R3Hessian SPSpline_M_Hess(T *f, R3Point *p);
void SPSpline_M_Maple(T *f, FILE *wr);
void SPSpline_M_Add(T *f, double a, T *h);
void SPSpline_M_Scale(T *f, double a);
T *SPSpline_M_Copy(T *f);
void SPSpline_M_Write(T *f, FILE *wr);
void SPSpline_M_Free(T *f);

double clip(double x, double lo, double hi);
  /* Returns {x} if {x} lies in {[lo _ hi]}, otherwise 
    returns the closest endpoint. */

bool_t SPH3_SamePlane(SPH3_Plane *P, SPH3_Plane *Q);
  /* True iff {P} and {Q} are the same plane, minus 
    positive homogeneous scalin. */
    
bool_t CorrectSide(SPH3_Plane *Q, R3Point *p);
  /* True iff {p} lies in the non-negative halfspace of {Q}. */

bool_t DisjointCaps(SPH3_Plane *fs, SPH3_Plane *gs);
  /* TRUE if the spherical caps defined by the oriented planes
    {fs} and {gs} on the unit sphere are disjoint. */

typedef double IntegrandBothPW(PieceData *ft, PieceData *gt, S2Point *p);
  /* A procedure that returns a function value on a spherical point {p},
    given the PieceData {ft,gt} of two PW functions {f,g} for a triangle 
    (of their common triangulation) that contains {p}. */

double IntegralBothPW
  ( SPSpline *f, SPSpline *g,
    IntegrandBothPW h,
    bool_t verbose
  );
  /* Computes the integral of {H(f(p), g(p), p)} over every triangle
    {t} of the triangulation {tri == f->d->tri == g->d->tri}, using
    {SPIntegral_BySamples} with the precomputed sampling points 
    for that triangle.
    
    The procedure {h(ft, gt, p)} should compute {H(f(p), g(p), p)}
    when given the {PieceData} {ft,gt} of {f}and {g} corresponding
    to some triangle of {tri} that contains {p}.  IMPORTANT:
    {H(fp, gp, p)} must be zero when either {fp} or {gp} is zero. */

typedef double IntegrandSinglePW (PieceData *ft, S2Point *p);
  /* A procedure that returns some value on a spherical point {p},
    given the PieceData {ft} of a partial function for on a 
    triangle that contains {p}. */

double IntegralSinglePW
  ( SPSpline *f,
    IntegrandSinglePW h, 
    bool_t verbose
  );
  /* Computes the integral of {H(f(p), p)} over every triangle
    {t} of the triangulation {tri == f->d->tri}, applying
    {SPIntegral_BySamples} with its precomputed sample points.
    
    The procedure {h(ft, p)} should compute {H(f(p), p)} when given
    the {PieceData} {ft} of {f} corresponding to some triangle of
    {tri} that contains {p}. IMPORTANT: {H(fp, p)} must be zero when
    {fp} is zero. */

#define SPSpline_FileFormat "2002-11-13"

static SPSpline_Methods *SplineMths = NULL;

bool_t CorrectSide(SPH3_Plane *Q, R3Point *p)
  { double v = 
      Q->f.c[0] + 
      Q->f.c[1]*p->c[0] + 
      Q->f.c[2]*p->c[1] + 
      Q->f.c[3]*p->c[2];
    return (v > 0.0);
  }

/* OVERRIDES FOR PARENT CLASS METHODS */

double SPSpline_M_Eval(T *f, R3Point *p) 
  { affirm(isprefix(SPSpline_TypeId, f->type), "type mismatch");
    if (! CorrectSide(&(f->d->supp), p))
      { return 0.0; }
    else
      { PieceDataRef_vec_t pd = f->d->pd;
        int i = SPSpline_LocatePiece(f, p);
        if (i >= pd.ne) 
          { return 0.0; }
        else
          { PieceData *pdi = pd.e[i];
            if (pdi->bary)
              { FaceNumber tn = pdi->face;
                Face *t = Left(f->d->tri->side.e[tn]);
                r3x3_t *c2b = &(t->c2b);
                r3_t a;
                r3x3_map_col(c2b, p, &a);
                return pdi->func->m->eval(pdi->func, &a);
              }
            else
              { return pdi->func->m->eval(pdi->func, p); }
          }
      }
  }

R3Gradient SPSpline_M_Grad(T *f, R3Point *p) 
  { affirm(isprefix(SPSpline_TypeId, f->type), "type mismatch");
    if (! CorrectSide(&(f->d->supp), p))
      { return (R3Gradient){{0.0, 0.0, 0.0}}; }
    else
      { PieceDataRef_vec_t pd = f->d->pd;
        int i = SPSpline_LocatePiece(f, p);
        if (i >= pd.ne) 
          { return (R3Gradient){{0.0, 0.0, 0.0}}; }
        else
          { PieceData *pdi = pd.e[i];
            if (pdi->bary)
              { FaceNumber tn = pdi->face;
                Face *t = Left(f->d->tri->side.e[tn]);
                r3x3_t *c2b = &(t->c2b);
                r3_t a, g;
                r3x3_map_col(c2b, p, &a);
                g = pdi->func->m->grad(pdi->func, &a);
                r3x3_map_row(&g, c2b, &g);
                return g;
              }
            else
              { return pdi->func->m->grad(pdi->func, p); }
          }
      }
  }
  
R3Hessian SPSpline_M_Hess(T *f, R3Point *p) 
  { affirm(isprefix(SPSpline_TypeId, f->type), "type mismatch");
    if (! CorrectSide(&(f->d->supp), p))
      { return (R3Hessian){{0.0, 0.0, 0.0, 0.0, 0.0, 0.0}}; }
    else
      { PieceDataRef_vec_t pd = f->d->pd;
        int i = SPSpline_LocatePiece(f, p);
        if (i >= pd.ne) 
          { return (R3Hessian){{0.0, 0.0, 0.0, 0.0, 0.0, 0.0}}; }
        else
          { PieceData *pdi = pd.e[i];
            if (pdi->bary)
              { FaceNumber tn = pdi->face;
                Face *t = Left(f->d->tri->side.e[tn]);
                r3x3_t *c2b = &(t->c2b);
                r3_t a;
                R3Hessian h;
                r3x3_map_col(c2b, p, &a);
                h = pdi->func->m->hess(pdi->func, &a);
                /* Correct {h} for change of basis: */
                { r3x3_t H;
                  H.c[0][0] = h.c[0];
                  H.c[0][1] = H.c[1][0] = h.c[1];
                  H.c[1][1] = h.c[2];
                  H.c[0][2] = H.c[2][0] = h.c[3];
                  H.c[1][2] = H.c[2][1] = h.c[4];
                  H.c[2][2] = h.c[5];
                  r3x3_mul(&H, c2b, &H);
                  r3x3_trmul(c2b, &H, &H);
                  /* Should check whether it is still symmetric... */
                  h.c[0] = H.c[0][0]; 
                  h.c[1] = H.c[0][1]; 
                  h.c[2] = H.c[1][1]; 
                  h.c[3] = H.c[0][2]; 
                  h.c[4] = H.c[1][2]; 
                  h.c[5] = H.c[2][2]; 
                }
                return h;
              }
            else
              { return pdi->func->m->hess(pdi->func, p); }
          }
      }
  }

void SPSpline_M_Maple(T *f, FILE *wr)
  { Arc_vec_t out = f->d->tri->out;
    Arc_vec_t arc = f->d->tri->arc;
    int i;
    affirm(isprefix(SPSpline_TypeId, f->type), "type mismatch");
    /* Define vertex coordinate table: */
    fprintf(wr, "vertex = [\n");
    for (i = 0; i < out.ne; i++)
      { Site *u = Org(out.e[i]);
        affirm(u->num == i , "inconsistent site numbers");
        r3_gen_print(wr, &(u->pos), "%24.16e", "[", " ", "]");
        if (i != out.ne - 1) { fprintf(wr, ","); }
        fprintf(wr, "\n");
      }
    fprintf(wr, "];\n");

    /* Define face incidence table: */
    fprintf(wr, "left = [\n");
    for (i = 0; i < arc.ne; i++)
      { fprintf(wr, "%5d",  Left(arc.e[i])->num);
        if (i != arc.ne - 1) { fprintf(wr, ","); }
      }
    fprintf(wr, "];\n");

    /* Define vertex incidence table: */
    fprintf(wr, "org = [\n");
    for (i = 0; i < arc.ne;  i++)
      { fprintf(wr, "%5d",  Org(arc.e[i])->num);
        if (i != arc.ne - 1) { fprintf(wr, ","); }
      }
    fprintf(wr, "];\n");

    /* Define barycentric flag of each face: */
    fprintf(wr, "baricentric = [\n");
    for (i = 0; i < f->d->pd.ne; i++)
      { bool_t b = f->d->pd.e[i]->bary;
        fprintf(wr, "%s", (b ? "true" : "false"));
        if (i != f->d->pd.ne - 1) { fprintf(wr, ","); }
      }
    fprintf(wr, "];\n");

    /* Dump functions for each face: */
    fprintf(wr, "facef = [\n");
    for (i = 0; i < f->d->pd.ne; i++)
      { SPFunction *fi = f->d->pd.e[i]->func;
        fi->m->maple(fi, wr);
        if (i != f->d->pd.ne - 1) { fprintf(wr, ","); }
      }
    fprintf(wr, "];\n");
    fflush(wr);
  }
  
void SPSpline_M_Add(T *f, double a, T *h)
  { affirm(isprefix(SPSpline_TypeId, f->type), "type/method bug");
    affirm(isprefix(SPSpline_TypeId, h->type), "wrong operand type");
    affirm(f->d->tri == h->d->tri, "incompatible triangulations");
    { int fNP = f->d->pd.ne; /* Note that {f->d->pd.ne} will increase. */
      int hNP = h->d->pd.ne;
      int i = 0, k = fNP, j;
      for (j = 0; j < hNP; j++)
        { PieceData *hpdj = h->d->pd.e[j];
          PieceData *fpdi;
          /* Look for face {hpdj->face} among the pieces of {f}: */
          while(i < fNP)
            { fpdi = f->d->pd.e[i];
              if (fpdi->face >= hpdj->face) { break; }
              if (i > 0) { affirm(fpdi->face > f->d->pd.e[i-1]->face, "unsorted pieces"); }
              i++;
            }
          if ((i >= fNP) || (fpdi->face > hpdj->face))
            { /* New triangle, append a scaled copy of it to {f->d->pd}: */
              PieceData *fpdk = SPSpline_PieceData_New();
              PieceDataRef_vec_expand(&(f->d->pd), k);
              fpdk->face = hpdj->face;
              fpdk->bary = hpdj->bary;
              fpdk->func = hpdj->func->m->copy(hpdj->func);
              fpdk->func->m->scale(fpdk->func, a);
              f->d->pd.e[k] = fpdk; k++;
            }
          else
            { /* Old triangle, add the respective funcs: */
              affirm(hpdj->bary == fpdi->bary, "inconsistent barys");
              affirm(hpdj->face == fpdi->face, "inconsistent faces");
              fpdi->func->m->add(fpdi->func, a, hpdj->func);
            }
        }
      if (k > fNP)
        { PieceDataRef_vec_trim(&(f->d->pd), k);
          SPSpline_SortPieces(f->d->pd);
          f->d->supp = SPSpline_ComputeSupp(f->d->pd, f->d->tri->side, NULL);
        }
    }
  }
  
void SPSpline_M_Scale(T *f, double a)
  { affirm(isprefix(SPSpline_TypeId, f->type), "type/method bug");
    { int i;
      for (i = 0; i < f->d->pd.ne; i++)
        { PieceData *fpdi = f->d->pd.e[i];
          fpdi->func->m->scale(fpdi->func, a);
        }
    }
  }
  
T *SPSpline_M_Copy(T *f)
  { affirm(isprefix(SPSpline_TypeId, f->type), "type/method bug");
    { int i;
      SPSpline *g = SPSpline_FullNew();
      int NP = f->d->pd.ne;
      PieceDataRef_vec_t gpd = PieceDataRef_vec_new(NP);
      g->d->triFile = f->d->triFile;
      g->d->tri = f->d->tri;
      g->d->lastLoc = f->d->lastLoc;
      g->d->pd = gpd;
      g->d->supp = f->d->supp;
      for (i = 0; i < f->d->pd.ne; i++)
        { PieceData *fpdi = f->d->pd.e[i];
          PieceData *gpdi = SPSpline_PieceData_New();
          gpdi->face = fpdi->face;
          gpdi->bary = fpdi->bary;
          gpdi->func = fpdi->func->m->copy(fpdi->func);
          g->d->pd.e[i] = gpdi;
        }
      SPSpline_SortPieces(g->d->pd); /* Just in case */
      return g;
    }
  }
 
void SPSpline_M_Free(T *f)
  { affirm(isprefix(SPSpline_TypeId, f->type), "type/method bug");
    { PieceDataRef_vec_t pd = f->d->pd;
      int i;
      for (i = 0; i < pd.ne; i++) 
        { PieceData *pdi = pd.e[i];
          if (pdi->func != NULL) { pdi->func->m->free(pdi->func); }
          free(pdi);
        }
      free(f->d->pd.e);
      free(f->d);
      free(f);
    }
  }
   
/* CLASS-SPECIFIC METHODS */
  
void SPSpline_M_Write(T *f, FILE *wr)
  { PieceDataRef_vec_t pd = f->d->pd;
    int i;
    affirm(isprefix(SPSpline_TypeId, f->type), "type mismatch");
    filefmt_write_header(wr, "SPSpline", SPSpline_FileFormat);
    fprintf(wr, "triangulation = %s\n", f->d->triFile);
    fprintf(wr, "supp = ");
    r4_gen_print(wr, &(f->d->supp.f), "%.16g", "[", ",", "]");
    fputc('\n', wr);
    fprintf(wr, "pieces = %d", pd.ne); fputc('\n', wr);
    for (i = 0; i < pd.ne; i++) 
      { PieceData *pdi = pd.e[i];
        fprintf(wr, "face = %d\n", pdi->face);
        fprintf(wr, "barycentric = %d\n", pdi->bary);
        { SPFunction *fi = pdi->func;
          fi->m->write(fi, wr);
        }
      }
    filefmt_write_footer(wr, "SPSpline");
    fflush(wr);
  }

/* OTHER PROCS */
  
SPSpline_Methods *SPSpline_Methods_New(void)
  { void *v = malloc(sizeof(SPSpline_Methods));
    return (SPSpline_Methods *)notnull(v, "no mem for SPSpline_Methods");
  }

SPSpline_Data *SPSpline_Data_New(void)
  { void *v = malloc(sizeof(SPSpline_Data));
    return (SPSpline_Data *)notnull(v, "no mem for SPSpline_Data");
  }

SPSpline *SPSpline_New(void)
  { void *v = malloc(sizeof(SPSpline));
    return (SPSpline *)notnull(v, "no mem for SPSpline");
  }

PieceData *SPSpline_PieceData_New(void)
  { void *v = malloc(sizeof(PieceData));
    return (PieceData *)notnull(v, "no mem for PieceData");
  }

SPSpline *SPSpline_FullNew(void)
  { SPSpline *f = SPSpline_New();
    f->type = SPSpline_TypeId;
    f->d = SPSpline_Data_New();
    if (SplineMths == NULL)
      { SplineMths = SPSpline_Methods_New();
        SplineMths->fn.eval = (SPFunction_EvalMth *)&SPSpline_M_Eval;
        SplineMths->fn.grad = (SPFunction_GradMth *)&SPSpline_M_Grad;
        SplineMths->fn.hess = (SPFunction_HessMth *)&SPSpline_M_Hess;
        SplineMths->fn.maple = (SPFunction_MapleMth *)&SPSpline_M_Maple;
        /* Note: the {fn.write} method is inherited from {SPFunction}! */
        SplineMths->fn.write = (SPFunction_WriteMth *)&SPFunction_M_Write;
        SplineMths->fn.scale = (SPFunction_ScaleMth *)&SPSpline_M_Scale;
        SplineMths->fn.add = (SPFunction_AddMth *)&SPSpline_M_Add;
        SplineMths->fn.copy = (SPFunction_CopyMth *)&SPSpline_M_Copy;
        SplineMths->fn.free = (SPFunction_FreeMth *)&SPSpline_M_Free;
        
        /* Class-specific methods */
        SplineMths->write = (SPFunction_WriteMth *)&SPSpline_M_Write;
      }
    f->m = SplineMths;
    return f;
  }

SPSpline *SPSpline_Read(FILE *rd)
  { int i;
    int NP;
    SPSpline *f = SPSpline_FullNew();
    int smpOrder = SPIntegral_GetDefaultSamplingOrder();
    filefmt_read_header(rd, "SPSpline", SPSpline_FileFormat);
    f->d->triFile = nget_string(rd, "triangulation"); fget_eol(rd);
    affirm(smpOrder != 0, "must call SPIntegral_SetDefaultSamplingOrder");
    f->d->tri = SPTriang_ReadCached(f->d->triFile, smpOrder);
    f->d->lastLoc = INT_MAX;

    /* Read supporting plane: */
    nget_name_eq(rd, "supp");
    fget_skip_spaces(rd); fget_match(rd, "[");
    f->d->supp.f.c[0] = fget_double(rd);
    for (i = 1; i <= 3; i++)
      { fget_skip_spaces(rd); fget_match(rd, ",");
        f->d->supp.f.c[i] = fget_double(rd);
      }
    fget_skip_spaces(rd); fget_match(rd, "]"); fget_eol(rd);

    /* Read number of nonzero faces: */
    NP = nget_int32(rd, "pieces"); fget_eol(rd);
    { /* Read partial function data: */
      PieceDataRef_vec_t pd = PieceDataRef_vec_new(NP);
      f->d->pd = pd;
      for (i = 0; i < pd.ne; i++)
        { PieceData *pdi = SPSpline_PieceData_New();
          pdi->face = nget_int32(rd, "face"); fget_eol(rd);
          pdi->bary = nget_bool(rd, "barycentric"); fget_eol(rd);
          pdi->func = SPFunction_Read(rd);
          f->d->pd.e[i] = pdi;
        }
    }
    filefmt_read_footer(rd, "SPSpline");
    SPSpline_SortPieces(f->d->pd); /* Just in case */
    return f;
  }
  
SPSpline *SPSpline_FromPieces
  ( char *triFile,
    Triangulation *tri,
    PieceDataRef_vec_t pd,
    r3_t *norm
  )
  { SPSpline *f = SPSpline_FullNew();
    affirm(tri != NULL, "null triangulation");
    f->d->triFile = triFile;
    f->d->tri = tri;
    f->d->pd = pd;
    f->d->lastLoc = INT_MAX;
    f->d->supp = SPSpline_ComputeSupp(pd, tri->side, norm);
    SPSpline_SortPieces(f->d->pd);
    return f;
  }
  
SPSpline *SPSpline_Replicate
  ( char *triFile,
    Triangulation *tri,
    bool_t bary,
    SPFunction *func
  )
  { SPSpline *f = SPSpline_FullNew();
    affirm(tri != NULL, "null triangulation");
    f->d->triFile = triFile;
    f->d->tri = tri;
    { int NF = tri->side.ne;
      int fn;
      PieceDataRef_vec_t pd = PieceDataRef_vec_new(NF);
      for (fn = 0; fn < NF; fn++)
        { PieceData *pi = SPSpline_PieceData_New();
          pi->face = fn;
          pi->bary = bary;
          pi->func = func->m->copy(func);
          pd.e[fn] = pi;
        }
      f->d->pd = pd;
    }
    f->d->lastLoc = INT_MAX;
    f->d->supp = NoPlane;
    return f;
  }
  
Triangulation *SPSpline_CommonTriangulation(SPFunction *f, SPFunction *g)
  { SPSpline *fpw = SPSpline_Cast(f);
    SPSpline *gpw = SPSpline_Cast(g);
    if ((fpw != NULL) && (gpw != NULL))
      { return (fpw->d->tri == gpw->d->tri ? fpw->d->tri : NULL); }
    else if (fpw != NULL)
      { return fpw->d->tri; }
    else if (gpw != NULL)
      { return gpw->d->tri; }
    else
      { return NULL; }
  }

Triangulation *SPSpline_BasisTriangulation(Basis basis)
  { int i;
    Triangulation *tri = NULL;
    for (i = 0; i < basis.ne; i++)
      { SPFunction *f = basis.e[i]; 
        SPSpline *fpw = SPSpline_Cast(f);
        if (fpw != NULL)
          { Triangulation *tbi = fpw->d->tri; 
            affirm(tbi != NULL, "no triangulation");
            if (tri == NULL) 
              { /* First spline element, save its triangulation. */
                tri = tbi;
              }
            else if (tbi != tri)
              { /* Basis uses two or more distinct triangulations. */
                return NULL;
              }
            else
              { /* Same triangulation so far, OK. */ }
          }
        else
          { /* Ignore non-PW elements. */ }
      }
    return tri;
  }
  
/* INTEGRALS OF PW FUNCTIONS */
  
double SPSpline_IntegralPW
  ( SPSpline *f,
    FuncMap FMap, 
    SPFunction *w, 
    bool_t verbose
  )
  { Arc_vec_t side = f->d->tri->side;
    
    auto double H(PieceData *ft, S2Point *p);
    
    double H(PieceData *ft, S2Point *p)
      { double wval = (w == NULL ? 1.0 : w->m->eval(w, p));
        double fval;
        r3_t b; /* Barycentric coords of {p}. */

        affirm(ft != NULL, "null triangle");
        if (ft->bary)
          { Face *t = Left(side.e[ft->face]);
            r3x3_map_col(&(t->c2b), p, &b);
            fval = ft->func->m->eval(ft->func, &b);
          }
        else
          { fval = ft->func->m->eval(ft->func, p); }
        
        if (FMap.map != NULL) { fval = FMap.map(fval, p); }
        return wval * fval;
      }

    affirm(FMap.zeroPres, "FMap is not zero-preserving");
    return IntegralSinglePW(f, H, verbose);
  }

R3Point SPSpline_CentroidPW
  ( SPSpline *f,
    FuncMap FMap, 
    SPFunction *w, 
    bool_t verbose
  )
  { r3_t bar;
    int i;
    
    affirm(FMap.zeroPres, "FMap is not zero-preserving");
    for (i = 0; i < 3; i++) 
      { 
        auto double fpi(double fval, S2Point *p);
          /* Returns F(p)*p[i]: */
        double fpi(double fval, S2Point *p)
          { if (FMap.map != NULL) { fval = FMap.map(fval, p); }
            return fval * p->c[i];
          }
        FuncMap Fpi = (FuncMap){ &fpi, FMap.zeroPres, "u*p[i]" };

        bar.c[i] = SPSpline_IntegralPW(f, Fpi, w, verbose);
      }
    return bar;
  }

/* INTERNAL INTEGRATION PROCS */

double IntegralSinglePW
  ( SPSpline *f,
    IntegrandSinglePW h,
    bool_t verbose
  )
  { int fi;
    double sum = 0.0, corr = 0.0;
    Triangulation *tri = f->d->tri;
    PieceDataRef_vec_t fd = f->d->pd;
    if (verbose) 
      { fprintf(stderr, "\nIntegralSinglePW: "); }
    /* if (verbose) { fprintf(stderr, "(fN=%d)", fN); } */
    
    /* Scan face records: */
    for (fi = 0; fi < fd.ne; fi++)
      { PieceData *ft = fd.e[fi];
        FaceNumber fface = ft->face;
        Arc e = tri->side.e[fface];
        Face *ff = Left(e);
        if (verbose) { fprintf(stderr, "[%d]", fface); }
        { 
          auto double HT(S2Point *p); /* Calls {h(ft,p)}. */

          double HT(S2Point *p) { return h(ft, p); } 

          if (verbose) { fprintf(stderr, "I"); }
          SPIntegral_BySamples(HT, ff->sp, ff->wp, &sum, &corr);
        }
      }
    return sum;
  }

double IntegralBothPW
  ( SPSpline *f, SPSpline *g,
    IntegrandBothPW h,
    bool_t verbose
  )
  { int fi, gi;
    double sum = 0.0, corr = 0.0;
    if (verbose) 
      { fprintf(stderr, "\nIntegralBothPW: "); }
    affirm(f->d->tri == g->d->tri , "incompatible triangulations");
    { Triangulation *tri = f->d->tri;
      PieceDataRef_vec_t fd = f->d->pd; int fN = fd.ne;
      PieceDataRef_vec_t gd = g->d->pd; int gN = gd.ne;
      /* if (verbose) { fprintf(stderr, "(fN=%d gN=%d)", fN, gN); } */
      if (DisjointCaps(&(f->d->supp), &(g->d->supp)))
        { /* if (verbose) { fprintf(stderr, "(disjoint)"); } */
          return 0.0;
        }
      /* Just to be sure... (cost is linear anyway). */
      SPSpline_SortPieces(fd);
      SPSpline_SortPieces(gd);
      /* Merge face records: */
      fi = 0; gi = 0;
      while ((fi < fN) || (gi < gN))
        { FaceNumber fface = (fi < fN ? fd.e[fi]->face : INT_MAX);
          FaceNumber gface = (gi < gN ? gd.e[gi]->face : INT_MAX);
          PieceData *ft = (fface <= gface ? fd.e[fi] : NULL);
          PieceData *gt = (gface <= fface ? gd.e[gi] : NULL);
          if (verbose)
            { fprintf(stderr, "[%d:%d]", 
               (fface == INT_MAX ? -1 : fface), 
               (gface == INT_MAX ? -1 : gface)
              );
            }
          if ((ft == NULL) || (gt == NULL))
            { /* We can ignore this triangle, {h(ft, gt, p) == 0} in it. */
              if (verbose) { fprintf(stderr, "o"); }
            }
          else
            { 
              auto double HT(S2Point *p);
                /* Calls {h(ft,gt,p)}. */

              double HT(S2Point *p) { return h(ft, gt, p); } 

              affirm(ft->face == gt->face, "merge bug");
              { FaceNumber tface = ft->face;
                Arc e = tri->side.e[tface];
                Face *t = Left(e);

                if (verbose) { fprintf(stderr, "I"); }
                SPIntegral_BySamples(HT, t->sp, t->wp, &sum, &corr);
              }
            }
          if (ft != NULL) { fi++; }
          if (gt != NULL) { gi++; }
        }
    }
    return sum;
  }
  
/* DOT PRODUCTS OF COMPATIBLE PW FUNCTIONS */
  
double SPSpline_DotBothPW
  ( SPSpline *f,
    FuncMap FMap, 
    SPSpline *g, 
    FuncMap GMap,
    SPFunction *w, 
    bool_t verbose
  )
  { Arc_vec_t side = f->d->tri->side;
    
    auto double H(PieceData *ft, PieceData *gt, S2Point *p);
    
    double H(PieceData *ft, PieceData *gt, S2Point *p)
      { double wval = (w == NULL ? 1.0 : w->m->eval(w, p));
        double fval, gval;
        r3_t b; /* Barycentric coords of {p}. */

        affirm(ft != NULL, "null piece");
        affirm(gt != NULL, "null piece");
        affirm(ft->face == gt->face, "mismatched triangles");
        if (ft->bary || gt->bary)
          { FaceNumber tn = ft->face;
            Face *t = Left(side.e[tn]);
            r3x3_map_col(&(t->c2b), p, &b);
          }
        
        if (ft->bary)
          { fval = ft->func->m->eval(ft->func, &b); }
        else
          { fval = ft->func->m->eval(ft->func, p); }
        
        if (gt == ft)
          { gval = fval; }
        else if (gt->bary)
          { gval = gt->func->m->eval(gt->func, &b); }
        else
          { gval = gt->func->m->eval(gt->func, p); }
        
        if (FMap.map != NULL) { fval = FMap.map(fval, p); }
        if (GMap.map != NULL) { gval = GMap.map(gval, p); }
        return wval * fval * gval;
        
      }

    affirm(FMap.zeroPres, "FMap is not zero-preserving");
    affirm(GMap.zeroPres, "GMap is not zero-preserving");
    affirm (f->d->tri == g->d->tri, "incompatible triangulations");
    return IntegralBothPW(f, g, H, verbose);
  }
  
double SPSpline_VelSGrdDotBothPW
  ( SPSpline *f, 
    SPSpline *g, 
    SPFunction *w, 
    bool_t verbose
  ) 
  { Arc_vec_t side = f->d->tri->side;

    auto double H(PieceData *ft, PieceData *gt, S2Point *p);
    
    double H(PieceData *ft, PieceData *gt, S2Point *p)
      { affirm(ft != NULL, "null piece");
        affirm(gt != NULL, "null piece");
        affirm(ft->face == gt->face , "wrong triangles");
        { Face *t = Left(side.e[ft->face]);
          double wval = (w == NULL ? 1.0 : w->m->eval(w, p));
          r3_t b; /* Barycentric coords of {p}. */
          
          double dval, gval;
          
          if ((ft->bary) || (gt->bary)) { r3x3_map_col(&(t->c2b), p, &b); }

          { R3Gradient fDp;
            if (ft->bary)
              { R3Gradient fu;
                fu = ft->func->m->grad(ft->func, &b);
                r3x3_map_row(&fu, &(t->c2b), &fDp);
              }
            else
              { fDp = ft->func->m->grad(ft->func, p); }
            /* Project onto tangent plane and compute dot product: */
            r3_t ftan, tmp;
            r3_decomp(&fDp, p, &tmp, &ftan);
            dval = p->c[0]*ftan.c[1] - p->c[1]*ftan.c[0];
          }

          if (gt->bary)
            { gval = gt->func->m->eval(gt->func, &b); }
          else
            { gval = gt->func->m->eval(gt->func, p); }
            
          return wval * dval * gval;

        }
      }

    if ((f == NULL) || (g == NULL)) { return 0.0; }
    affirm (f->d->tri == g->d->tri, "incompatible triangulations");
    return IntegralBothPW(f, g, H, verbose);
  }

double SPSpline_SLapDotBothPW
  ( SPSpline *f, 
    SPSpline *g, 
    SPFunction *w, 
    bool_t verbose
  ) 
  { Arc_vec_t side = f->d->tri->side;

    auto double H(PieceData *ft, PieceData *gt, S2Point *p);
    
    double H(PieceData *ft, PieceData *gt, S2Point *p)
      { affirm(ft != NULL, "null piece");
        affirm(gt != NULL, "null piece");
        affirm(ft->face == gt->face , "wrong triangles");
        { Face *t = Left(side.e[ft->face]);
          double wval = (w == NULL ? 1.0 : w->m->eval(w, p));
          r3_t b; /* Barycentric coords of {p}. */
          R3Gradient fDp, gDp;
                
          
          if ((ft->bary) || (gt->bary)) { r3x3_map_col(&(t->c2b), p, &b); }

          if (ft->bary)
            { R3Gradient fu;
              fu = ft->func->m->grad(ft->func, &b);
              r3x3_map_row(&fu, &(t->c2b), &fDp);
            }
          else
            { fDp = ft->func->m->grad(ft->func, p); }

          if (ft == gt)
            { gDp = fDp; }
          else if (gt->bary)
            { R3Gradient gu;
              gu = gt->func->m->grad(gt->func, &b);
              r3x3_map_row(&gu, &(t->c2b), &gDp);
            }
          else
            { gDp = gt->func->m->grad(gt->func, p); }
            
          /* Project onto tangent plane and compute dot product: */
          { r3_t ftan, gtan, tmp;
            r3_decomp(&fDp, p, &tmp, &ftan);
            r3_decomp(&gDp, p, &tmp, &gtan);
            return wval * r3_dot(&ftan, &gtan);
          } 
        }
      }

    if ((f == NULL) || (g == NULL)) { return 0.0; }
    affirm (f->d->tri == g->d->tri, "incompatible triangulations");
    return IntegralBothPW(f, g, H, verbose);
  }

/* DOT PRODUCTS OF A PW FUNCTION AGAINST A GENERAL FUNCTION */
  
double SPSpline_DotSinglePW
  ( SPSpline *f,
    FuncMap FMap, 
    SPFunction *g, 
    FuncMap GMap,
    SPFunction *w, 
    bool_t verbose
  )
  { Arc_vec_t side = f->d->tri->side;
    
    auto double H(PieceData *ft, S2Point *p);
    
    double H(PieceData *ft, S2Point *p)
      { double wval = (w == NULL ? 1.0 : w->m->eval(w, p));
        double fval, gval;
        affirm(ft != NULL, "null piece");
        if (ft->bary)
          { r3_t b;
            Face *t = Left(side.e[ft->face]);
            r3x3_map_col(&(t->c2b), p, &b);
            fval = ft->func->m->eval(ft->func, &b);
          }
        else
          { fval = ft->func->m->eval(ft->func, p); }
        gval = (g == NULL ? 0.0 : g->m->eval(g, p));
        if (FMap.map != NULL) { fval = FMap.map(fval, p); }
        if (GMap.map != NULL) { gval = GMap.map(gval, p); }
        return wval * fval * gval;
      }

    affirm(FMap.zeroPres, "FMap is not zero-preserving");
    if ((g == NULL) && GMap.zeroPres) { return 0.0; }
    return IntegralSinglePW(f, H, verbose);
  }

double SPSpline_VelSGrdDotSinglePW
  ( SPSpline *f, 
    SPFunction *g, 
    SPFunction *w, 
    bool_t verbose
  ) 
  { Arc_vec_t side = f->d->tri->side;
    
    auto double H(PieceData *ft, S2Point *p);
    
    double H(PieceData *ft, S2Point *p)
      { double wval = (w == NULL ? 1.0 : w->m->eval(w, p));
        double dval, gval;
        
        /* Compute velocity dot gradient: */
        { R3Gradient fDp;
          affirm(ft != NULL, "null piece");
          if (ft->bary)
            { r3_t b;
              R3Gradient fu;
              Face *t = Left(side.e[ft->face]);
              r3x3_map_col(&(t->c2b), p, &b);
              fu = ft->func->m->grad(ft->func, &b);
              r3x3_map_row(&fu, &(t->c2b), &fDp);
            }
          else
            { fDp = ft->func->m->grad(ft->func, p); }
          /* Project gradient onto tangent plane and compute dot product: */
          r3_t ftan, tmp;
          r3_decomp(&fDp, p, &tmp, &ftan);
          dval = p->c[0]*ftan.c[1] - p->c[1]*ftan.c[0];
        }

        gval = g->m->eval(g, p);
        return wval * dval * gval;
      }

    if ((f == NULL) || (g == NULL)) { return 0.0; }
    return IntegralSinglePW(f, H, verbose);
  }

double SPSpline_SLapDotSinglePW
  ( SPSpline *f, 
    SPFunction *g, 
    SPFunction *w, 
    bool_t verbose
  ) 
  { Arc_vec_t side = f->d->tri->side;
    
    auto double H(PieceData *ft, S2Point *p);
    
    double H(PieceData *ft, S2Point *p)
      { double wval = (w == NULL ? 1.0 : w->m->eval(w, p));
        R3Gradient fDp, gDp;

        affirm(ft != NULL, "null piece");
        if (ft->bary)
          { r3_t b;
            R3Gradient fu;
            Face *t = Left(side.e[ft->face]);
            r3x3_map_col(&(t->c2b), p, &b);
            fu = ft->func->m->grad(ft->func, &b);
            r3x3_map_row(&fu, &(t->c2b), &fDp);
          }
        else
          { fDp = ft->func->m->grad(ft->func, p); }
        
        gDp = g->m->grad(g, p);
        
        /* Project onto tangent plane and compute dot product: */
        { r3_t ftan, gtan, tmp;
          r3_decomp(&fDp, p, &tmp, &ftan);
          r3_decomp(&gDp, p, &tmp, &gtan);
          return wval * r3_dot(&ftan, &gtan);
        } 
      }

    if ((f == NULL) || (g == NULL)) { return 0.0; }
    return IntegralSinglePW(f, H, verbose);
  }

double clip(double x, double lo, double hi)
  { return ((x > hi ? hi : (x < lo ? lo : x))); }
  
bool_t DisjointCaps(SPH3_Plane *fs, SPH3_Plane *gs)
  { SPH3_Plane om = Omega;
    if (SPH3_SamePlane(fs, &om) || SPH3_SamePlane(gs, &om)) { return FALSE; }
    { r3_t fval = (r3_t){{fs->f.c[1], fs->f.c[2], fs->f.c[3]}};
      r3_t gval = (r3_t){{gs->f.c[1], gs->f.c[2], gs->f.c[3]}};
      r3_t fd, gd;
      r3_dir(&fval, &fd);
      r3_dir(&gval, &gd);
      { /* Angle betwen the two normals, in [0_Pi]: */
        double thetaC = clip(r3_dot(&fd, &gd), -1.0, 1.0);
        double theta = acos(thetaC);
        /* Angular radius of each cap, in [0_Pi]: */
        double fC = clip(-(fs->f.c[0])/r3_norm(&fval), -1.0, 1.0);
        double falpha = acos(fC);
        double gC = clip(-(gs->f.c[0])/r3_norm(&gval), -1.0, 1.0);
        double galpha = acos(gC);
        /* See if they intersect: */
        return theta >= falpha + galpha;
      }
    }
  }

void SPSpline_SortPieces(PieceDataRef_vec_t pd)
  { int i, k;
    PieceData *x, *y;
    int n = pd.ne;
    for (i = 1; i < n; i++)
      { x = pd.e[i];
        k = i;
        while ((k > 0) && ((y = pd.e[k-1])->face > x->face))
          { pd.e[k] = y; k--; }
        affirm(y->face != x->face, "repeated face in spline");
        if (k < i) { pd.e[k] = x; }
      }
  }

int SPSpline_LocatePiece(SPSpline *f, S2Point *p)
  { 
    double eps = 1.0e-10;
    double minDist = HUGE_VAL;
    int minPiece = INT_MAX;
    PieceDataRef_vec_t pd = f->d->pd;
    Arc_vec_t side = f->d->tri->side;
    int i;
    
    affirm(isprefix(SPSpline_TypeId, f->type), "type mismatch");
    i = f->d->lastLoc;
    if (i < pd.ne)
      { FaceNumber tn = pd.e[i]->face;
        Face *t = Left(side.e[tn]);
        minPiece = i;
        minDist = TriDist(p, &(t->c2b));
        if (minDist <= 0.0)
          { /* fprintf(stderr, "!"); */
            return minPiece;
          }
      }
    /* fprintf(stderr, "?"); */
    for (i = 0; i < pd.ne; i++)
      { FaceNumber tn = pd.e[i]->face;
        Face *t = Left(side.e[tn]);
        double d = TriDist(p, &(t->c2b));
        if (d <= 0.0)
          { f->d->lastLoc = i;
            return i;
          }
        else if (d < minDist )
          { minDist = d; minPiece = i; }
      }
    f->d->lastLoc = minPiece;
    if (minDist < eps) { return minPiece; }
    return INT_MAX;
  }

R3Point SPSpline_SupportCentroid(PieceDataRef_vec_t pd, Arc_vec_t side)
  { r3_t b = (r3_t){{0.0, 0.0, 0.0}};
    int i, k;
    for (i = 0; i < pd.ne; i++)
      { FaceNumber tn = pd.e[i]->face;
        Arc e = side.e[tn];
        for (k = 0; k < 3; k++)
          { S2Point *p = &(Org(e)->pos);
            r3_add(p, &b, &b);
            e = Lnext(e);
          }
      }
    r3_scale(1.0/((double)3*pd.ne), &b, &b);
    return b;
  }

Sign SPSpline_PosRelPlane(PieceDataRef_vec_t pd, SPH3_Plane *Q, Arc_vec_t side)
  { Sign dclass = 0;
    int i;
    for (i = 0; i < pd.ne; i++)
      { FaceNumber tn = pd.e[i]->face;
        Arc e = side.e[tn], a = e;
        do 
          { SPH3_Point p = SPH3_FromCartesian(&(Org(a)->pos));
            Sign pclass = SPH3_Side(&p, Q); 
            if (dclass*pclass < 0) { return 0; }
            if (dclass == 0) { dclass = pclass; }
            a = Lnext(a);
          }
        while (a != e);
      }
    return dclass;
  }

S2Point SPSpline_FarthestVertex(S2Point *u, PieceDataRef_vec_t pd, Arc_vec_t side)
  { r3_t *v = u;
    double mincos = 1.0;
    int i, k;
    for (i = 0; i < pd.ne; i++)
      { FaceNumber tn = pd.e[i]->face;
        Arc e = side.e[tn];
        for (k = 0; k < 3; k++)
          { S2Point *p = &(Org(e)->pos);
            double c = r3_dot(p,u);
            if (c < mincos) { mincos = c; v = p; } 
            e = Lnext(e);
          }
      }
    return *v;
  }

SPH3_Plane SPSpline_ComputeSupp
  ( PieceDataRef_vec_t pd, 
    Arc_vec_t side, 
    r3_t *norm
  )
  { double eps = 0.00001;
    r3_t u;

    if (pd.ne == 0) { return Omega; }

    if ((norm != NULL) && (r3_L_inf_norm(norm) > 0.0))
      { u = *norm; }
    else
      { /* Compute barycenter of all vertices, normalize to sphere: */
        u = SPSpline_SupportCentroid(pd, side);
      }
      
    if (r3_L_inf_norm(&u) == 0.0) { return Omega; }
    r3_dir(&u, &u);
   
    /* Compute plane orthogonal to {norm} that contains all vertices. */
    { r3_t v = SPSpline_FarthestVertex(&u, pd, side);
      double t = r3_dot(&v, &u) - eps;
      /* Caps larger than one hemisphere are not convex:  */
      if (t < 0.0) { return Omega; }
      return (SPH3_Plane){{{-t, u.c[0], u.c[1], u.c[2]}}};
    }
  }
  
r3x3_t SPSpline_SupportFrame(PieceDataRef_vec_t pd, Arc_vec_t side)
  { r3_t u = SPSpline_SupportCentroid(pd, side);
    r3_t far;
    r3_t v, w, para;
    if (r3_L_inf_norm(&u) == 0.0) { u = (r3_t){{1.0, 1.0, 1.0}}; }
    r3_dir(&u, &u);
    far = SPSpline_FarthestVertex(&u, pd, side);
    r3_decomp(&far, &u, &para, &v);
    if (r3_norm(&v) < 1.0e-6)
      { (void)r3_throw_ortho(&u, &v); }
    else
      { r3_dir(&v, &v); }
    r3_cross(&u, &v, &w);
    r3_dir(&w, &w); /* Just in case... */
    { r3x3_t m; int i;
      for (i = 0; i < 3; i++)
        { m.c[0][i] = u.c[i]; m.c[1][i] = v.c[i]; m.c[2][i] = w.c[i]; }
      return m;
    }
  }

bool_t SPH3_SamePlane(SPH3_Plane *P, SPH3_Plane *Q)
  { int i;
    double Pm = fabs(P->f.c[0])+fabs(P->f.c[1])+fabs(P->f.c[2])+fabs(P->f.c[3]); 
    double Qm = fabs(Q->f.c[0])+fabs(Q->f.c[1])+fabs(Q->f.c[2])+fabs(Q->f.c[3]); 
    affirm((Pm > 0.0) && (Qm > 0.0), "indeterminate plane");
    for (i = 0; i < 4; i++)
      { if (Qm*P->f.c[i] != Pm*Q->f.c[i]) { return FALSE; } }
    return TRUE;
  }

/* BASIS SORTING */

void SPSpline_LexSortBasis
  ( Triangulation *tri, 
    Basis bas, 
    int ini, int lim,
    bool_t verbose
  )
  { int i;
    
    auto Sign cmp_supports(SPFunction *f, SPFunction *g);
      /* -1 if {f} should come before {g}, +1 if it should
        come after, {0} if they are equivalent. */
     
    Sign cmp_supports(SPFunction *f, SPFunction *g)
      { SPSpline *fpw = SPSpline_Cast(f);
        SPSpline *gpw = SPSpline_Cast(g);
        int fNP, gNP;
        
        if 
          ( ((fpw == NULL) || (fpw->d->tri != tri)) &&
            ((gpw == NULL) || (gpw->d->tri != tri))
          )
          { return 0; }
        else if ((fpw == NULL) || (fpw->d->tri != tri))
          { return +1; }
        else if ((gpw == NULL) || (gpw->d->tri != tri))
          { return -1; }
          
        /* Get number of triangles: */
        fNP = fpw->d->pd.ne; 
        gNP = gpw->d->pd.ne; 
          
        /* Compare supports lexically. */
        /* Assumes that pieces are sorted by increasing face number. */
        { PieceDataRef_vec_t fpd = fpw->d->pd;
          PieceDataRef_vec_t gpd = gpw->d->pd;
          int fk = fNP-1, gk = gNP-1;
          while ((fk >= 0) || (gk >= 0))
            { int fn = (fk < 0 ? -1 : fpd.e[fk]->face);
              int gn = (gk < 0 ? -1 : gpd.e[gk]->face);
              if (fn != gn) { return (fn < gn ? -1 : +1); }
              fk--; gk--;
            }
          /* Same set of triangles: */
          return 0;
        }
      }
    
    /* Heapsort of basis elements: */
    /* Make {bas} into a binary heap, with largest at the root {bas[ini]}: */
     /* fprintf(stderr, "heaping...\n"); */
    for (i = ini; i < lim; i++)
      { /* Insert {bas[i]} in heap {bas[ini..i-1]}. */
        /* Save {f = bas[i]}: */
        SPFunction *f = bas.e[i];
        /* Heap is now {bas[ini..i]} with {bas[i]} vacant. */
        /* fprintf(stderr, "  f<-[%d]", i); */
        /* Sort pieces of {f}, to simplify support comparisons: */
        { SPSpline *fpw = SPSpline_Cast(f);
          if (fpw != NULL) { SPSpline_SortPieces(fpw->d->pd); }
        }
        { /* Bubble the vacancy rootwards, until it is fit for {f}: */
          int j = i;
          while (TRUE)
            { /* Try to insert {f} into the vacancy {bas[j]}. */
              if (j <= ini) 
                { /* {bas[j]} has no father: */ break; }
              else
                { /* Get father {g = bas[k]} of the vacancy: */
                  int k = ini + (j-ini-1)/2;
                  SPFunction *g = bas.e[k];
                  if (cmp_supports(g, f) >= 0)
                    { /* {bas[j]} is fit for {f}: */ break; }
                  /* Demote father {g = bas[k]} into vacancy {bas[j]}: */
                  bas.e[j] = g;
                  /* fprintf(stderr, " [%d]<-[%d]", j, k); */
                  /* The father {bas[k]} is now vacant, try to put {f} there: */
                  j = k;
                }
            }
          /* The vacancy {bas[j]} is ok for {f}, put it there: */
          bas.e[j] = f; 
          /* fprintf(stderr, " [%d]<-f.\n", j); */
        }
      }
    
    /* fprintf(stderr, "unheaping...\n"); */
    /* Now remove elements from heap root, and store at end of {bas}: */
    for (i = lim-1; i > ini; i--)
      { int k;
        /* Save the leaf currently at {bas[i]}: */
        SPFunction *f = bas.e[i];
        /* Remove the heap root and insert it there: */
        /* fprintf(stderr, " f<-[%d]", i); */
        bas.e[i] = bas.e[ini];
        /* fprintf(stderr, " [%d]<-[%d]", i, ini); */
        /* Heap is now {bas[ini..i-1]}, with {bas[ini]} vacant. */
        /* Promote elements until we can find a place for {f}: */
        k = ini; 
        while (TRUE)
          { /* Now {bas[k]} is vacant, consider storing {f} there: */
            /* See if {bas[k]} has any chldren: */
            int j = ini + 2*(k-ini) + 1;
            if (j >= i) 
              { /* {bas[k]} has no children, OK to put {f} there: */ break; }
            else
              { /* Find largest child {g = bas[j]} of {bas[k]}: */
                SPFunction *g = bas.e[j];
                /* If {bas[j]} has a sibling, let {g} be the largest one: */
                int j2 = j + 1;
                if (j2 < i)
                  { SPFunction *g2 = bas.e[j2];
                    if (cmp_supports(g2, g) > 0) { j = j2; g = g2; }
                  }
                /* Check whether {f} can be their father: */
                if (cmp_supports(f, g) >= 0)
                  { /* Fine, can put {f} in {bas[k]}: */ break; }
                /* Promote {g}, keep going: */
                /* fprintf(stderr, " [%d]<-[%d]", k, j); */
                bas.e[k] = g; k = j;
              }
          }
        /* Now the vacancy is {bas[k]}, and it is OK to put {f} there: */
        /* fprintf(stderr, " [%d]<-f.\n", k); */
        bas.e[k] = f;
      }
    if (verbose) { SPSpline_DescribeSupports(stderr, tri, bas, ini, lim); }
    SPSpline_CheckBasisOrder(tri, bas, ini, lim); 
  }

void SPSpline_BSPSortBasis
  ( Triangulation *tri, 
    Basis bas, 
    int ini, int lim,
    bool_t verbose
  )
  { 
    auto void sort(int ini, int lim); 
      /* Sorts the elements {bas[ini..lim-1]} by recursive
         geometric split. */
      
    auto S2Point some_vertex(int ini, int lim);
      /* Returns some vertex of the supports of {bas[ini..lim-1]}. */
    
    auto S2Point farthest_vertex(S2Point *p, int ini, int lim);
      /* Returns the vertex of the supports of {bas[ini..lim-1]} that 
        is farthest from {p}. */
    
    auto SPH3_Plane bisector(S2Point *p, S2Point *q);
      /* The plane through the origin that bisects segment {p--q}. */
      
    auto void separate(int ini, int lim, SPH3_Plane *Q, int *limn, int *limp);
      /* Permutes {bas[ini..lim-1]} and finds indices {limn} and
        {limp} so that elements {bas[ini..limn-1]} have supports
        entirely contained in the negative side of {Q}; elements
        {bas[limn..limp-1]} are entirely contained in the positive
        side; and the remaining elements have supports that extend into
        both sides of {Q}. Assuming that the union of the supports of
        elements {bas[ini..lim-1]} has at least one support vertex on
        each side of {Q}, on return we will have {ini < limn < limp <=
        lim}; */
    
    void sort(int ini, int lim)
      { if (lim - ini > 1)
          { S2Point x = some_vertex(ini, lim);
            S2Point u = farthest_vertex(&x, ini, lim);
            S2Point v = farthest_vertex(&u, ini, lim);
            SPH3_Plane Q = bisector(&u, &v);
            int limn, limp;
            separate(ini, lim, &Q, &limn, &limp);
            fprintf(stderr, " [%d=%d:%d:%d]", lim-ini, limn-ini, limp-limn, lim-limp);
            /* The support can't be all on the same side of {Q}: */
            affirm(
              (limn-ini < lim-ini) || (limp-limn < lim-ini),
              "no progress in split"
            );
            sort(ini, limn);
            sort(limn, limp); 
            /* Sort the straddling elements, if they are not all clumped */
            if (limp > ini)
              { sort(limp, lim); }
            else
              { fprintf(stderr, " L[%d]", lim-ini);
                SPSpline_LexSortBasis(tri, bas, limp, lim, FALSE);
              }
          }
      }
    
    void separate(int ini, int lim, SPH3_Plane *Q, int *limn, int *limp)
      { 
        int inip = ini; /* Start of positive segment */
        int iniz = ini; /* Start of zero segment */
        int inix = ini; /* Start of unclassified segment */
        while (inix < lim)
          { /* Status is |------|+++++++|00000|???????| */
            int k = inix;
            SPFunction *f = bas.e[k];
            Sign class = 0;
            
            { SPSpline *fpw = SPSpline_Cast(bas.e[k]);
              if ((fpw != NULL) && (fpw->d->tri == tri))
                { class = SPSpline_PosRelPlane(fpw->d->pd, Q, tri->side); }
            }
            
            /* Remove {f} from the unclassified set: */
            inix ++;
            /* Now {bas[k]} is a hole at the end of the zero set. */
            if (class != 0)
              { /* Bubble the hole {ini}ward, over the zero set: */
                if (iniz != k) { bas.e[k] = bas.e[iniz]; k = iniz; }
                iniz++;
              }
            if (class < 0)
              { /* Bubble the hole {ini}ward, over the positive set: */
                if (inip < k) { bas.e[k] = bas.e[inip]; k = inip; }
                inip++;
              }
            /* Now {bas[k]} is the right place to put {f}: */
            bas.e[k] = f; 
          }
        (*limn) = inip; (*limp) = iniz;
      }
      
    SPH3_Plane bisector(S2Point *p, S2Point *q)
      { r3_t n;
        r3_sub(q, p, &n); r3_dir(&n, &n);
        return (SPH3_Plane){{{ 0.0, n.c[0], n.c[1], n.c[2] }}};
      }
      
    S2Point farthest_vertex(S2Point *p, int ini, int lim)
      { int i;
        S2Point v = *p; /* Current vertex farthest from {p} */
        double cospv = 1.0;
        for (i = ini; i < lim; i++)
          { SPSpline *fpw = SPSpline_Cast(bas.e[i]);
            if ((fpw != NULL) && (fpw->d->tri == tri))
              { PieceDataRef_vec_t pd = fpw->d->pd;
                S2Point w = SPSpline_FarthestVertex(p, pd, tri->side); 
                double cospw = r3_cos(p, &w);
                if (cospw < cospv) { v = w; cospv = cospw; }
              }
          }
        return v;
      }
      
    S2Point some_vertex(int ini, int lim)
      { int i;
        for (i = ini; i < lim; i++)
          { SPSpline *fpw = SPSpline_Cast(bas.e[i]);
            if ((fpw != NULL) && (fpw->d->tri == tri))
              { PieceDataRef_vec_t pd = fpw->d->pd;
                if (pd.ne != 0)
                  { FaceNumber fn = pd.e[0]->face;
                    return Org(tri->side.e[fn])->pos; 
                  }
              }
          }
        /* No vertices in any of those elements, just return something: */
        return (S2Point){{ 1.0, 0.0, 0.0 }};
      }
      
    /* So, let do it: */
    sort(ini, lim);
    fprintf(stderr, "\n");
    if (verbose) { SPSpline_DescribeSupports(stderr, tri, bas, ini, lim); }
    SPSpline_CheckBasisOrder(tri, bas, ini, lim); 
  }
  
void SPSpline_GDistSortBasis
  ( Triangulation *tri, 
    Basis bas, 
    int ini, int lim,
    bool_t verbose
  )
  { 
    int NF = tri->side.ne;
    
    auto void sort_elems(int ini, int lim); 
      /* Sorts the elements {bas[ini..lim-1]} by recursive
        graph-distance split. */
      
    auto void get_faces(int ini, int lim, nat_vec_t supfn, int *nsup, nat_vec_t ix);
      /* Stores in {nsup} the number of faces in the union of the
        supports of {bas[ini..lim-1]}, and in {supfn[0..nsup-1]} the
        numbers of those faces. Weird basis elements are ignored.
        Uses {ix[0..NF-1]} as a work area. */
    
    auto void bisect_faces(nat_vec_t supfn, int nsup, bool_vec_t is_low, nat_vec_t ix);
      /* Rearranges the faces {fn} listed in {supfn[0..nsup-1]} into two
        hopefully compact non-empty subsets, `low' and `high', and
        sets {is_low[fn]} accordingly. Assumes {nsup >= 2} and
        {is_sup[fn] == TRUE} iff {fn} is listed in {supfn[0..nsup-1]}.
        Uses {ix[0..NF-1]} as a work area. */
        
    auto void separate(int ini, int lim, bool_vec_t is_low, int *limn, int *limp);
      /* Permutes {bas[ini..lim-1]} and finds indices {limn} and
        {limp} so that (1) elements {bas[ini..limn-1]} have supports
        that consist entirely of `low' faces; (2) elements
        {bas[limn..limp-1]} have supports that consist entirely of
        `high' faces; and (3) the remaining elements have supports
        that contain both types of faces.  Assumes that {is_low[fn]}
        tells whether a support face is `low' or `high'.  
        
        If that the union of the supports of elements
        {bas[ini..lim-1]} uses at least one face of each type, on
        return we will have {ini < limn < limp <= lim}. Weird elements
        are sorted into group (3), as if they straddled both sets. */
      
    auto void bfs_sort_faces(nat_vec_t supfn, int nsup, nat_vec_t ix);
      /* Performs a breadth-first sort of the faces {supfn[0..nsup-1]},
        starting from the first face {supfn[0]}. 
        Uses {ix[0..NF-1]} as a work area. */

    bool_vec_t is_low = bool_vec_new(NF);
    nat_vec_t supfn = nat_vec_new(NF);
    nat_vec_t ix = nat_vec_new(NF); /* Work area for routines below. */
    
    void sort_elems(int ini, int lim)
      { int nbas = lim - ini;
        if (nbas > 1)
          { int limn, limp;
            int nsup = 0;
            get_faces(ini, lim, supfn, &nsup, ix);
            if (nsup > 1)
              { bisect_faces(supfn, nsup, is_low, ix);
                separate(ini, lim, is_low, &limn, &limp);
              }
            else
              { limn = ini; limp = ini; }
            fprintf(stderr, " [%d=%d:%d:%d]", lim-ini, limn-ini, limp-limn, lim-limp);
            if ((limn-ini < nbas) && (limp-limn < nbas) && (lim-limp < nbas))
              { /* Split was not trivial, recurse on each part: */
                sort_elems(ini,  limn);
                sort_elems(limn, limp); 
                sort_elems(limp,  lim);
              }
            else
              { /* Trivial split, fall back to lexical sort: */
                fprintf(stderr, " L[%d]", lim-ini);
                SPSpline_LexSortBasis(tri, bas, ini,  limn, FALSE);
                SPSpline_LexSortBasis(tri, bas, limn, limp, FALSE);
                SPSpline_LexSortBasis(tri, bas, limp,  lim, FALSE);
              }
          }
      }
      
    void get_faces(int ini, int lim, nat_vec_t supfn, int *nsup, nat_vec_t ix)
      { int n = 0;
        int i;
        for (i = ini; i < lim; i++)
          { SPSpline *fpw = SPSpline_Cast(bas.e[i]);
            if ((fpw != NULL) && (fpw->d->tri == tri))
              { PieceDataRef_vec_t pd = fpw->d->pd;
                int k;
                for (k = 0; k < pd.ne; k++)
                  { FaceNumber fn = pd.e[k]->face;
                    /* Face {fn} has been visited iff {supfn[ix[fn]] == fn}. */
                    int i = ix.e[fn];
                    if ((i >= n) || (supfn.e[i] != fn))
                      { ix.e[fn] = n; supfn.e[n] = fn; n++; }
                  }
              }
          }
        (*nsup) = n;
      }
    
    void bisect_faces(nat_vec_t supfn, int nsup, bool_vec_t is_low, nat_vec_t ix)
      { int iter, i;
        affirm(nsup >= 2, "cannot bisect 1 face");
        for (iter = 1; iter <= 2; iter++)
          { /* Swap first and last last face: */
            { FaceNumber fn = supfn.e[0], gn = supfn.e[nsup-1];
              supfn.e[0] = gn; supfn.e[nsup-1] = fn;
            }
            /* Sort faces by distance from first face (breadth-first search): */
            bfs_sort_faces(supfn, nsup, ix);
          }
        /* Mark first half as `low', upper half as `high': */
        { int nlow = nsup/2;
          affirm((nlow > 0) && (nlow < nsup), "invalid split");
          for (i = 0; i < nsup; i++)
            { is_low.e[supfn.e[i]] = (i < nlow); }
        }
      }
        
    void bfs_sort_faces(nat_vec_t supfn, int nsup, nat_vec_t ix)
      { int i;
        /* Mark faces {supfn[0..nsup-1]} in {ix}, by the back-index trick: */
        for (i = 0; i < nsup; i++) { ix.e[supfn.e[i]] = i; }
        /* Scan faces : */
        /* fprintf(stderr, "<%d:", nsup); */
        { int nsort = 1, nclosed = 0;
          while(nclosed < nsup)
            { FaceNumber fn = supfn.e[nclosed]; 
              Arc e, a;
              /* fprintf(stderr, "%s%d", (nclosed > 0 ? "," : ""), fn); */
              affirm(ix.e[fn] == nclosed, "inconsistent ix"); 
              /* In case the set {supfn} is disconnected: */
              if (nsort == nclosed) { nsort = nclosed + 1; }
              affirm(nclosed < nsort, "inconsistent nclosed");
              /* Enumerate sides {a} of face {fn}: */
              e = tri->side.e[fn]; a = e;
              do 
                { /* Enumerate edges {b} with {Org(b) = Org(a)}: */
                  Arc b = a;
                  do 
                    { FaceNumber gn = Left(b)->num;
                      /* Check wheter {Left(b)} is in {fsup} but unvisited: */
                      int gi = ix.e[gn];
                      if ((gi < nsup) && (gi >= nsort) && (supfn.e[gi] == gn))
                        { /* Put face {gn} in the sorted set: */
                          if (gi > nsort)
                            { FaceNumber hn = supfn.e[nsort];
                              affirm(ix.e[hn] == nsort, "inconsistent ix"); 
                              supfn.e[nsort] = gn; ix.e[gn] = nsort;
                              supfn.e[gi] = hn; ix.e[hn] = gi;
                            }
                          nsort++;
                          affirm(nsort <= nsup, "inconsistent nsort");
                        }
                      b = Onext(b);
                    }
                  while (b != a);
                  a = Lnext(a);
                } 
              while (a != e);
              nclosed++;
            }
        }
        /* fprintf(stderr, ">"); */
      }

    void separate(int ini, int lim, bool_vec_t is_low, int *limn, int *limp)
      { 
        int inip = ini; /* Start of high segment */
        int iniz = ini; /* Start of mixed segment */
        int inix = ini; /* Start of unclassified segment */
        while (inix < lim)
          { /* Status is |------|+++++++|00000|???????|, where
              {-} is all-low, {+} all-high, {0} mixed, {?} unclassified. */
            int k = inix;
            SPFunction *f = bas.e[k];
            SPSpline *fpw = SPSpline_Cast(bas.e[k]);
            bool_t has_low = FALSE, has_high = FALSE;
            if ((fpw == NULL) || (fpw->d->tri != tri))
              { has_low = TRUE; has_high = TRUE; }
            else
              { PieceDataRef_vec_t pd = fpw->d->pd;
                int i;
                for (i = 0; i < pd.ne; i++)
                  { FaceNumber fn = pd.e[i]->face;
                    bool_t fn_is_low = is_low.e[fn];
                    has_low |= fn_is_low;
                    has_high |= (! fn_is_low);
                  }
              }
            /* Remove {f} from the unclassified set: */
            inix ++;
            /* Now {bas[k]} is a hole at the end of the zero set. */
            affirm(has_low || has_high, "inconsistent low/high");
            if (! (has_low && has_high))
              { /* {f} is not mixed; bubble hole {ini}ward, over the mixed set: */
                if (iniz != k) { bas.e[k] = bas.e[iniz]; k = iniz; }
                iniz++;
              }
            if (! has_high)
              { /* {f} is all-low; bubble hole {ini}ward, over the high set: */
                if (inip < k) { bas.e[k] = bas.e[inip]; k = inip; }
                inip++;
              }
            /* Now {bas[k]} is the right place to put {f}: */
            bas.e[k] = f; 
          }
        (*limn) = inip; (*limp) = iniz;
      }
      
    /* So, let do it: */
    sort_elems(ini, lim);
    fprintf(stderr, "\n");
    free(ix.e); free(supfn.e); free(is_low.e);
    if (verbose) { SPSpline_DescribeSupports(stderr, tri, bas, ini, lim); }
    SPSpline_CheckBasisOrder(tri, bas, ini, lim); 
  }
  
bool_t SPSpline_WeirdSupport(SPFunction *f, Triangulation *tri)
  { SPSpline *fpw = SPSpline_Cast(f);
    return (fpw == NULL) || (fpw->d->tri != tri);
  }

bool_t SPSpline_SameSupport
  ( SPFunction *f, 
    SPFunction *g, 
    Triangulation *tri
  )
  { bool_t fWeird = SPSpline_WeirdSupport(f, tri);
    bool_t gWeird = SPSpline_WeirdSupport(g, tri);
    if (fWeird && gWeird) 
      { return TRUE; }
    else if (fWeird || gWeird)
      { return FALSE; }
    else
      { PieceDataRef_vec_t fpd = SPSpline_Cast(f)->d->pd;
        PieceDataRef_vec_t gpd = SPSpline_Cast(g)->d->pd;
        if (fpd.ne != gpd.ne)
          { return FALSE; }
        else 
          { /* Assumes that pieces are sorted in the same order. */
            int NP = fpd.ne, k;
            for (k = 0; k < NP; k++)
              { if (fpd.e[k]->face != gpd.e[k]->face) { return FALSE; } }
            return TRUE;
          }
      }
  }

bool_t SPSpline_SubSupport
  ( SPFunction *f, 
    SPFunction *g, 
    Triangulation *tri
  )
  { if (SPSpline_WeirdSupport(g, tri)) 
      { return TRUE; }
    else if (SPSpline_WeirdSupport(f, tri))
      { return FALSE; }
    else
      { PieceDataRef_vec_t fpd = SPSpline_Cast(f)->d->pd;
        PieceDataRef_vec_t gpd = SPSpline_Cast(g)->d->pd;
        int fNP = fpd.ne, gNP = gpd.ne;
        if (fNP > gNP)
          { return FALSE; }
        else 
          { /* Assumes that pieces are sorted by increasing face number. */
            int kf, kg = 0;
            for (kf = 0; kf < fNP; kf++)
              { FaceNumber fn = fpd.e[kf]->face, gn;
                /* Look for face {fn} in the {gp} list: */
                while ((kg < gNP) && ((gn = gpd.e[kg]->face) < fn)) { kg++; }
                if ((kg >= gNP) || (fn != gn)) { return FALSE; }
              }
            return TRUE;
          }
      }
  }

bool_t SPSpline_StrictSubSupport
  ( SPFunction *f, 
    SPFunction *g, 
    Triangulation *tri
  )
  { if (SPSpline_WeirdSupport(f, tri)) 
      { return FALSE; }
    else if (SPSpline_WeirdSupport(g, tri))
      { return TRUE; }
    else
      { PieceDataRef_vec_t fpd = SPSpline_Cast(f)->d->pd;
        PieceDataRef_vec_t gpd = SPSpline_Cast(g)->d->pd;
        int fNP = fpd.ne, gNP = gpd.ne;
        if (fNP >= gNP) 
          { return FALSE; }
        else
          { /* Assumes that pieces are sorted by increasing face number. */
            int fk = 0, gk = 0;
            while ((fk < fNP) || (gk < gNP))
              { int fn = (fk >= fNP ? INT_MAX : fpd.e[fk]->face);
                int gn = (gk >= gNP ? INT_MAX : gpd.e[gk]->face);
                if (fn < gn)
                  { /* {f} has extra stuff */ return FALSE; }
                else if (fn > gn)
                  { /* {g} has extra stuff, OK: */ gk++; }
                else
                  { /* Skip common elts: */ fk++; gk++; }
              }
            return TRUE; /* Because {fNP < gNP} */
          }
      }
  }    
 
bool_t SPSpline_NestedSupports(SPFunction *f, SPFunction *g, Triangulation *tri)
  { if (tri == NULL) { return FALSE; }
    return SPSpline_SubSupport(f,g,tri) || SPSpline_SubSupport(f,g,tri);
  }

void SPSpline_CheckBasisOrder
  ( Triangulation *tri, 
    Basis bas, 
    int ini, int lim
  )
  { int i, j;
    for (i = ini; i < lim; i++)
      { SPFunction *bi = bas.e[i];
        for (j = i+1; j < lim; j++)
          { SPFunction *bj = bas.e[j];
            affirm(! SPSpline_StrictSubSupport(bj, bi, tri), "sorting failure"); 
          }
      }
  }

void SPSpline_DescribeSupports
  ( FILE *wr,
    Triangulation *tri, 
    Basis bas, 
    int ini, int lim
  )
  { int cini = ini, clim;
    while (cini < lim)
      { SPFunction *bini = bas.e[cini];
        /* Find end of support class: */
        for (clim = cini+1; clim < lim; clim++)
          { SPFunction *blim = bas.e[clim];
            if (! SPSpline_SameSupport(bini, blim, tri)) { break; }
          }
        /* Now support class is {bas[cini..clim-1]}. */
        fprintf(wr, "  [%04d..%04d] %4d elems", cini, clim-1, clim - cini);
        SPSpline_DescribeSupport(wr, bini, tri);
        fprintf(wr, "\n");
        cini = clim;
      }
  }

void SPSpline_DescribeSupport(FILE *wr, SPFunction *f, Triangulation *tri)
  { SPSpline *fpw = SPSpline_Cast(f);
    if ((fpw == NULL) || (fpw->d->tri != tri))
      { fprintf(wr, " weird"); }
    else
      { PieceDataRef_vec_t pd = fpw->d->pd;
        int k;
        fprintf(wr, " faces (");
        for (k = 0; k < pd.ne; k++) { fprintf(wr, " %d", pd.e[k]->face); }
        fprintf(wr, " )");
        if (pd.ne > 1)
          { nat_t vCom[3], eCom[3], fCom[1];
            int nvCom, neCom, nfCom;
            FindSharedElems(&pd, tri, vCom, &nvCom, eCom, &neCom, fCom, &nfCom );
            affirm(nvCom <= 2 && neCom <= 1 && nfCom == 0, "buggy support");
            if (neCom > 0)
              { fprintf(wr, " sharing edge %d", eCom[0]); }
            else 
              { affirm(nvCom <= 1, "buggy support");
                if (nvCom > 0) { fprintf(wr, " sharing vertex %d", vCom[0]); }
              }
          }
      }
  }

void FindSharedElems
  ( PieceDataRef_vec_t *pd,
    Triangulation *tri,
    nat_t *vCom, int *nvCom, 
    nat_t *eCom, int *neCom,
    nat_t *fCom, int *nfCom
  )
  {
    auto void keep(nat_t x, nat_t *xCom, int *nxCom, int *jxCom);
      /* Used to compute intersections of vertex and edge sets.
        Looks for {x} in {xCom[jxCom..nxCom-1]}; if found,
        swaps it with {xCom[0..jxCom]} and increments {jxCom}. */

    void keep(nat_t x, nat_t *xCom, int *nxCom, int *jxCom)
      { int t;
        for(t = (*jxCom); t < (*nxCom); t++)
          { if (xCom[t] == x)
              { xCom[t] = xCom[*jxCom];
                xCom[*jxCom] = x; 
                (*jxCom)++;
              }
          }
      }

    int jvCom, jeCom, jfCom, j;

    *nvCom = 0; *neCom = 0; *nfCom = 0;
    for (j = 0; j < pd->ne; j++)
      { PieceData *pdj = pd->e[j];
        FaceNumber fn = pdj->face;
        Arc e = tri->side.e[fn];
        int k;
        /* From elems common to all previous pieces, pick those of {pd[j]}: */
        jvCom = 0; jeCom = 0; jfCom = 0;
        if (j == 0) { fCom[0] = fn; jfCom = 1; }
        for (k = 0; k < 3; k++)
          { SiteNumber vn = Org(e)->num;
            EdgeNumber en = EdgeNum(e);
            if (j == 0)
              { vCom[k] = vn; jvCom = k+1;
                eCom[k] = en; jeCom = k+1;
              }
            else
              { keep(vn, vCom, nvCom, &jvCom);
                keep(en, eCom, neCom, &jeCom);
              }
            e = Lnext(e);
          }
        *nvCom = jvCom; *neCom = jeCom; *nfCom = jfCom; 
      }
    affirm((*nfCom) <= 1, "inconsistent nfCom");
    
    if ((*nfCom) == 0) 
      { affirm((*neCom) <= 1, "inconsistent neCom"); }
    else
      { affirm((*neCom) == 3, "inconsistent neCom"); }
      
    if ((*neCom) == 3)
      { affirm((*nvCom) == 3, "inconsistent nvCom"); }
    else if ((*neCom) == 1)
      { affirm((*nvCom) == 2, "inconsistent nvCom"); }
    else
      { affirm((*nvCom) <= 1, "inconsistent nvCom"); }
  }

void SPSpline_FillVertexSupport(SPSpline *f)
  { Triangulation *tri = f->d->tri;
    PieceDataRef_vec_t *pd = &(f->d->pd);
    int k;
    SiteNumber vCom[3]; int nvCom;
    EdgeNumber eCom[3]; int neCom;
    FaceNumber fCom[3]; int nfCom;

    if (pd->ne <= 1) { /* Face element: */ return; }
    FindSharedElems(pd, tri, vCom, &nvCom, eCom, &neCom, fCom, &nfCom);
    
    /* 
    fprintf (stderr, "\n");
    fprintf (stderr, "%2d input pieces\n", oldNP);
    fprintf (stderr, "%2d shared vertices = (", nvCom);
    { int k; for (k = 0; k < nvCom; k++) { fprintf (stderr, " %2d", vCom[k]); }}
    fprintf (stderr, ")\n");
    fprintf (stderr, "%2d shared edges = (", neCom);
    { int k; for (k = 0; k < neCom; k++) { fprintf (stderr, " %2d", eCom[k]); }}
    fprintf (stderr, ")\n");
    fprintf (stderr, "%2d shared faces = (", nfCom);
    { int k; for (k = 0; k < nfCom; k++) { fprintf (stderr, " %2d", fCom[k]); }}
    fprintf (stderr, ")\n");
    */
    
    affirm(neCom <= 1, "faces share more than one edge");
    if (pd->ne == 2)
      { if (neCom > 0) { /* Edge element: */ return; } }
    else
      { affirm(neCom == 0, "three or more faces sharing an edge"); }

    /* Now there are no shared edges: */
    affirm(nvCom <= 1, "faces share no edge but more than one vertex"); 
    affirm(nvCom > 0, "support is not contained in a vertex star");

    /* Scan faces around shared vertex, adding missing ones to {f->d->pd}: */
    { SiteNumber vn = vCom[0];
      int origNP = pd->ne, totNP = origNP;
      Arc e = tri->out.e[vn], a = e;
      do 
        { FaceNumber fn = Left(a)->num;
          /* See if face {fn} was already present: */
          for (k = 0; k < origNP; k++)
            { if (pd->e[k]->face == fn) { break; }}
          if (k >= origNP) 
            { /* Face {fn} was missing, append a new piece for it: */
              PieceData *pd0 = pd->e[0];
              PieceData *pdn = SPSpline_PieceData_New();
              pdn->face = fn;
              pdn->bary = pd0->bary;
              pdn->func = pd0->func->m->copy(pd0->func);
              pdn->func->m->scale(pdn->func, 0.0);
              PieceDataRef_vec_expand(pd, totNP);
              pd->e[totNP] = pdn;
              totNP++;
            }
          a = Onext(a);
        }
      while (a != e);
      if (totNP != origNP)
        { /* fprintf (stderr, "%2d pieces -> %2d pieces\n", origNP, totNP); */
          PieceDataRef_vec_trim(pd, totNP);
          SPSpline_SortPieces(*pd);
          f->d->supp = SPSpline_ComputeSupp(*pd, tri->side, &(Org(e)->pos));
        }
    }
  }

void r3x3_trmul (r3x3_t *A, r3x3_t *B, r3x3_t *R)
  { r3x3_t RR;
    int n = 3;
    int i, j, k;
    for (i = 0; i < n; i++)
      for (j = 0; j < n; j++)
        { double s = 0.0;
          for (k = 0; k < n; k++)  s += A->c[k][i]*B->c[k][j];
          RR.c[i][j] = s;
        }
    for (i = 0; i < n; i++)
      for (j = 0; j < n; j++)
        { R->c[i][j] = RR.c[i][j]; }
  }

/* Arrays of {PieceData*}: */

vec_typeimpl(PieceDataRef_vec_t,PieceDataRef_vec,PieceDataRef);
