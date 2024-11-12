/* See SPDelaunay.h. */
/* Last edited on 2005-10-06 21:51:35 by stolfi */

#include <SPDelaunay.h>
#include <SPTriang.h>
#include <SPQuad.h>
#include <SPBasic.h>

#include <vec.h>
#include <r3.h>
#include <affirm.h>
#include <nat.h>
#include <bool.h>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define NullArc SPQuad_NullArc
#define Splice SPQuad_Splice
#define MakeEdge SPQuad_MakeEdge

void Consist(Arc e);

S2Point TetrahedronVertex(int i);
S2Point OctahedronVertex(int i); /* Slightly twisted */
S2Point IcosahedronVertex(int i);
  /* Enumerate the vertices of the regular polyhedra of unit radius. */

/* IMPLEMENTATIONS */

void Consist(Arc e)
  { Arc t = e;
    affirm(e != Onext(e), "hanging vertex - Org");
    affirm(Oprev(e) != Onext(e), "bridging vertex - Org");
    affirm(e != Dnext(e), "hanging vertex - Dest");
    affirm(Dprev(e) != Dnext(e), "bridging vertex - Dest");
    affirm(Org(e) != Dest(e), "loop edge");
    do 
      { Site *a = Org(t); 
        Site *b = Dest(t);
        Site *u = Dest(Onext(t));
        Site *v = Dest(Oprev(t));
        if (! PositiveSide(&(a->pos), &(u->pos), &(b->pos), &(v->pos)))
          { fprintf(stderr, "Consist failed\n");
            fprintf(stderr, "  a = s[%d] = ", a->num); 
            r3_print(stderr, &(a->pos)); fprintf(stderr, "\n");
            fprintf(stderr, "  b = s[%d] = ", b->num); 
            r3_print(stderr, &(b->pos)); fprintf(stderr, "\n");
            fprintf(stderr, "  u = s[%d] = ", u->num); 
            r3_print(stderr, &(u->pos)); fprintf(stderr, "\n");
            fprintf(stderr, "  v = s[%d] = ", v->num); 
            r3_print(stderr, &(v->pos)); fprintf(stderr, "\n");
            affirm(FALSE, "inconsistent Delaunay structure");
          }
        t = Onext(t);
      }
    while (t != e);
  }

Arc SPDelaunay_BuildTetrahedron(Site *a, Site *b, Site *c, Site *d)
  { Arc e1, e2, e;
    e1 = MakeEdge();
    Odata(e1) =  a;  
    Ddata(e1) =  b;
    e2 = MakeEdge();
    Odata(e2) =  b;  
    Ddata(e2) =  c;
    Splice (Sym(e1), e2);
    e = SPTriang_Connect(e2, e1);   
    if (PositiveSide(&(Org(e1)->pos), &(Org(e2)->pos), &(Dest(e2)->pos), &(d->pos)))
      { SPTriang_InsertSiteInFace(e, d); 
        Consist(e);
        Consist(Rnext(e));
        Consist(Rprev(e));
      }
    else
      { SPTriang_InsertSiteInFace(Sym(e), d);
        Consist(e);
        Consist(Lnext(e));
        Consist(Lprev(e));
      }
    return e;
  }

void SPDelaunay_Retriangulate(Arc e, int deg, Arc_vec_t *tr)
  { 
    auto Arc RemoveCorner(Arc a);
      /* Removes the corner {Org(a)}, {Dest(a)}, {Dest(Lnext(a))}
        from the face {Left(a)}, if that triangle is a Delaunay
        triangle; and returns the new edge that is 
        a side of the new triangle. Otherwise returns NullArc. */

    Arc RemoveCorner(Arc a)
      { Arc r;
        Arc b = Lnext(a);

        /* Check whether {Org(a)}, {Dest(a)}, {Dest(b)} is a 
          valid face; that is, whether all other vertices are
          on the correct side of plane: */ 
        r = Lnext(Lnext(b));
        affirm(r != a, "triangular face");
        while (PositiveSide(&(Org(r)->pos), &(Org(a)->pos), &(Dest(a)->pos), &(Dest(b)->pos)))
          { r = Lnext(r);
            if (r == a){ return SPTriang_Connect(b, a); }
          }
        return NullArc;
      }

    int nTr = 0;
    while (deg > 3)
      { Arc g = RemoveCorner(e);
        if (g != NullArc)
          { SPQuad_Arc_vec_expand(tr, nTr);
            tr->e[nTr] = g; nTr++;
            deg--; e = Sym(g);
          }
        else
          { e = Lnext(e); }
      }
    SPQuad_Arc_vec_expand(tr, nTr);
    tr->e[nTr] = e; nTr++;
  }

bool_t SamePoint(S2Point *p, S2Point *q)
  { return 
      (p->c[0] == q->c[0]) &&
      (p->c[1] == q->c[1]) &&
      (p->c[2] == q->c[2]);
  }
   
Arc SPDelaunay_InsertSite(Site *s, Arc e)
  { Arc t, res;
    Site *pri;
    /* 
      fprintf(stderr, "InsertSite(");
      r3_print(stderr, &(s->pos));
      fprintf(stderr, ")\n");
    */
    e = SPTriang_Locate(&(s->pos), e);
    if (SamePoint(&(s->pos), &(Org(e)->pos))){ return e; }
    { Site *a = Org(e);
      Site *b = Dest(e);
      Site *c = Dest(Onext(e));
      affirm(PositiveSide(&(a->pos), &(b->pos), &(c->pos), &(s->pos)), "folded triangle");
    }
    pri = Dest(e);
    SPTriang_InsertSiteInFace(e, s);
    while (1)
      { t = Dnext(e);
        if (PositiveSide(&(Org(e)->pos), &(Org(t)->pos), &(Dest(e)->pos), &(s->pos)))
          { SPTriang_SwapEdge(e); e = t; }
        else if (Org(e) == pri)
          { res = Sym(Onext(e)); break; }
        else
          { e = Lprev(Onext(e)); }
      }
    Consist(res);
    t = res;
    do { Consist(Sym(t)); t = Onext(t); } while (t != res);
    return res;
  }  

Arc SPDelaunay_BuildInc(SiteRef_vec_t site)
  { Arc e; 
    int i;
    affirm(site.ne >= 4, "too few sites");
    e = SPDelaunay_BuildTetrahedron(site.e[0], site.e[1], site.e[2], site.e[3]);
    for (i = 4; i < site.ne; i++) 
      { Site *s = site.e[i]; 
        double eps = ((double)i) * 1.0e-9;
        r3_scale(1.0 - eps, &(s->pos), &(s->pos));
        e = SPDelaunay_InsertSite(site.e[i], e);
        r3_dir(&(s->pos), &(s->pos));
      }
    return e;
  }

bool_t SPDelaunay_EdgeIsDelaunay(Arc e)
  { if (SPTriang_EdgeIsDegnerate(e))
      { return FALSE; }
    else
      { /* Edge must be convex: */
        S2Point *p = &(Org(e)->pos); 
        S2Point *q = &(Dest(e)->pos);
        S2Point *u = &(Dest(Onext(e))->pos);
        S2Point *v = &(Dest(Oprev(e))->pos);
        return PositiveSide(p, u, q, v);
      }
  }

void SPDelaunay_CheckDelaunay(Arc_vec_t arc)
  { int i;
    for (i = 0; i < arc.ne; i += 2)
      { Arc e = arc.e[i];
        if (e != NullArc)
          { affirm(SPDelaunay_EdgeIsDelaunay(e), "non-Delaunay edge"); }
      }
  } 
  
bool_t SPDelaunay_FixDelaunay(Triangulation *tri)
  { int NE = tri->arc.ne/2;
    int NF = tri->side.ne;
    
    auto bool_t swappable(Arc e);
      /* TRUE if {e} can be swapped, namely if {Left(e) != Right(e)}. */
    
    auto void push_suspect(Arc e);
      /* Marks edge {e} to be checked for swappability. */
    
    auto int num_sus(void);
      /* Number of currently suspect edges. */
    
    auto Arc pop_suspect(void);
      /* Removes an arc from the suspect list. Returns NULL if none. */
    
    auto void mark_left(Arc e);
      /* Sets {e} as the reference edge of face {Left(e)},
        and marks that face as one that needs to have its 
        coordinate matrices recomputed (because {e} was swapped). */
        
    bool_vec_t face_changed = bool_vec_new(NF);
    bool_vec_t is_suspect = bool_vec_new(NE);
    Arc_vec_t suspect = SPQuad_Arc_vec_new(NE);
    int first_sus = -1, last_sus = -1;

    void push_suspect(Arc e)
      { nat_t en = EdgeNum(e);
        if (! is_suspect.e[en])
          { is_suspect.e[en] = TRUE; 
            last_sus++; if (last_sus >= NE) { last_sus = 0; }
            suspect.e[last_sus] = e; 
            if (first_sus == -1) { first_sus = last_sus; }
          }
      }
      
    int num_sus(void)
      { if (first_sus == -1)
          { return 0; }
        else
          { return (last_sus + 1 + NE - first_sus) % NE; }
      }
    
    Arc pop_suspect(void)
      { if (first_sus == -1)
          { return NullArc; }
        else
          { Arc e = suspect.e[first_sus];
            if (first_sus == last_sus)
              { first_sus = -1; last_sus = -1; }
            else
              { first_sus = (first_sus+1) % NE; }
            suspect.e[EdgeNum(e)] = FALSE;
            return e;
          }
      }
    
    void mark_left(Arc e)
      { nat_t fn = Left(e)->num;
        tri->side.e[fn] = e;
        face_changed.e[fn] = TRUE;
      }
    
    bool_t swappable(Arc e)
      { return (Onext(e) != e) && (Dnext(e) != e); }

    bool_t failed = FALSE; /* Failed to fix all edges. */
  
    /* Initialize suspect list: */
    { int en; for (en = 0; en < NE; en++) { is_suspect.e[en] = FALSE; } }
    { int fn; for (fn = 0; fn < NF; fn++) { face_changed.e[fn] = FALSE; } }
    first_sus = -1; last_sus = -1;
    /* Mark all edges as suspect: */
    { int en; 
      for (en = 0; en < NE; en++) 
        { Arc e = tri->arc.e[2*en];
          if (e != NullArc) { push_suspect(e); }
        }
    }
    /* Check suspect edges, swap bad ones, re-suspect neighbors: */
    { Arc a; 
      int npass = 0; /* Number of re-pushed suspects since last swap. */
      while ((! failed) && ((a = pop_suspect()) != NullArc))
        { if (! SPDelaunay_EdgeIsDelaunay(a))
            { if (swappable(a))
                { Arc b = Sym(a);
                  /* Ensure {tri->out} remains consistent after swap: */
                  tri->out.e[Org(a)->num] = Onext(a);
                  tri->out.e[Org(b)->num] = Onext(b);
                  /* Swap the edge: */
                  SPTriang_SwapEdge(a);
                  /* Neighboring edges are now suspect: */
                  push_suspect(Onext(a));
                  push_suspect(Oprev(a));
                  push_suspect(Onext(b));
                  push_suspect(Oprev(b));
                  /* Schedule recomputation of matrices for {Left(a),Right(a)}: */
                  mark_left(a);
                  mark_left(b);
                  /* SPTriang_CheckTopology(tri, FALSE); */
                  fprintf (stderr, "S");
                  npass = 0;
                }
              else
                { fprintf (stderr, "!"); }
            }
          /* If still bad, keep it as suspect, hope it will go away someday: */
          if (! SPDelaunay_EdgeIsDelaunay(a))
            { push_suspect(a);
              fprintf (stderr, "#");
              npass++;
              if (npass > NE) { failed = TRUE; }
              affirm(! failed, "FixDelaunay failed");
            }
          else
            { affirm(SPDelaunay_EdgeIsDelaunay(Sym(a)), "inconsistent IsDelaunay"); }
        }
    }
    free(is_suspect.e); free(suspect.e);
    /* Recompute face data, for faces that changed: */
    { int fn;
      for (fn = 0; fn < NF; fn++) 
        { if (face_changed.e[fn])
            { Arc e = tri->side.e[fn]; 
              SPTriang_ComputeFaceMatrices(e); 
              SPTriang_ComputeFaceSamples(e, tri->smpOrder); 
            }
        }
    }
    free(face_changed.e); 
    if (failed) { return FALSE; }
    /* Check again for Delaunayhood, for some reason {failed} is not enough: */
    { int en;
      for (en = 0; en < NE; en++)
        { Arc e = tri->arc.e[2*en];
          if ((e != NullArc) && (! SPDelaunay_EdgeIsDelaunay(e)))
            { return FALSE; }
        }
    } 
    return TRUE;     
  } 

/* REGULAR POLYHEDRA */

S2Point TetrahedronVertex(int i)
  { r3_t p;
    switch(i)
      { case  0: p.c[0] = +1.0; p.c[1] = +1.0; p.c[2] = +1.0; break;
        case  1: p.c[0] = -1.0; p.c[1] = -1.0; p.c[2] = +1.0; break;
        case  2: p.c[0] = +1.0; p.c[1] = -1.0; p.c[2] = -1.0; break;
        case  3: p.c[0] = -1.0; p.c[1] = +1.0; p.c[2] = -1.0; break;
        default: affirm(FALSE , "bad vertex number");
      }
    r3_dir(&p, &p);
    return p;
  }

S2Point OctahedronVertex(int i)
  { r3_t p;
    /* Twist vertices a bit to avoid coplanar vertices during build: */
    double eps = 1.0e-6;
    switch(i)
      { case  0: p.c[0] = +1.0; p.c[1] = +eps; p.c[2] = -eps; break;
        case  1: p.c[0] = -1.0; p.c[1] = +eps; p.c[2] = -eps; break;
        case  2: p.c[0] = -eps; p.c[1] = +1.0; p.c[2] = +eps; break;
        case  3: p.c[0] = -eps; p.c[1] = -1.0; p.c[2] = +eps; break;
        case  4: p.c[0] = +eps; p.c[1] = -eps; p.c[2] = +1.0; break;
        case  5: p.c[0] = +eps; p.c[1] = -eps; p.c[2] = -1.0; break;
        default: affirm(FALSE , "bad vertex number");
      }
    r3_dir(&p, &p);
    return p;
  }
  
#define IcoA (PHI)
#define IcoB (1.0)

S2Point IcosahedronVertex(int i)
  { r3_t p;
    switch(i)
      { case  0: p.c[0] = 0.0; p.c[1] = +IcoA; p.c[2] = +IcoB; break;
        case  1: p.c[0] = 0.0; p.c[1] = +IcoA; p.c[2] = -IcoB; break;
        case  2: p.c[0] = +IcoA; p.c[1] = -IcoB; p.c[2] = 0.0; break;
        case  3: p.c[0] = -IcoA; p.c[1] = -IcoB; p.c[2] = 0.0; break;
        case  4: p.c[0] = +IcoB; p.c[1] = 0.0; p.c[2] = +IcoA; break;
        case  5: p.c[0] = +IcoB; p.c[1] = 0.0; p.c[2] = -IcoA; break;
        case  6: p.c[0] = 0.0; p.c[1] = -IcoA; p.c[2] = +IcoB; break;
        case  7: p.c[0] = 0.0; p.c[1] = -IcoA; p.c[2] = -IcoB; break;
        case  8: p.c[0] = -IcoB; p.c[1] = 0.0; p.c[2] = +IcoA; break;
        case  9: p.c[0] = -IcoB; p.c[1] = 0.0; p.c[2] = -IcoA; break;
        case 10: p.c[0] = +IcoA; p.c[1] = +IcoB; p.c[2] = 0.0; break;
        case 11: p.c[0] = -IcoA; p.c[1] = +IcoB; p.c[2] = 0.0; break;
        default: affirm(FALSE , "bad vertex number");
      }
    r3_dir(&p, &p);
    return p;
  }
  
Triangulation *SPDelaunay_RegularTetrahedron(int smpOrder)
  { SiteRef_vec_t site = SPTriang_MakeSites(4, TetrahedronVertex);
    Arc e = SPDelaunay_BuildInc(site);
    return SPTriang_FromTopology(e, smpOrder);
  }

Triangulation *SPDelaunay_RegularOctahedron(int smpOrder)
  { SiteRef_vec_t site = SPTriang_MakeSites(6, OctahedronVertex);
    Arc e = SPDelaunay_BuildInc(site);
    return SPTriang_FromTopology(e, smpOrder);
  }
 
Triangulation *SPDelaunay_RegularIcosahedron(int smpOrder)
  { SiteRef_vec_t site = SPTriang_MakeSites(12, IcosahedronVertex);
    Arc e = SPDelaunay_BuildInc(site);
    return SPTriang_FromTopology(e, smpOrder);
  }

