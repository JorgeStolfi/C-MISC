/* See SPTriangExtra.h */
/* Last edited on 2008-05-24 12:31:16 by stolfi */

#include <SPTriangExtra.h>
#include <SPTriang.h>
#include <SPIntegral.h>
#include <SPQuad.h>
#include <SPBasic.h>

#include <vec.h>
#include <r3.h>
#include <affirm.h>
#include <nat.h>

#include <values.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

Triangulation *SPTriang_Refine57Triang
  ( Triangulation *old, 
    r3_t *b5, 
    r3_t *b7
  )
  { /* Element tables of old triangulation: */
    Arc_vec_t oldOut = old->out;
    Arc_vec_t oldArc = old->arc;
    Arc_vec_t oldSide = old->side;
    /* Old element counts: */
    nat_t oldNV = oldOut.ne;
    nat_t oldNA = oldArc.ne;
    nat_t oldNF = oldSide.ne;
    /* Element counts for new triangulation: */
    nat_t NV = oldNV + 2*oldNA;
    nat_t NA = oldNA + 12*oldNA;
    nat_t NF = 7*oldNF + 2*oldNA;
    /* Element tables for new triangulation: */
    SiteRef_vec_t site = SiteRef_vec_new(NV);
      /* The vertices of {old} are also vertices of the new
        triangulation, with the original numbers {0..oldNV-1}.
        Inside each face of {old} we add three new vertices of degree 5
        and three of degree 7. If {e = Lnext^r(oldSide[i])} is any arc of {old},
        then the new 7-fold vertex in that face that is nearest to {Org(e)}
        is {v7(e) = site[oldNV + 2*(3*i + r%3)]}; and the 
        new 5-fold vertex closest to the midpoint of {e} is 
        {v5(e) = site.e[oldNV + 2*(3*i + r%3) + 1]}. */
        
    Arc_vec_t arc = SPQuad_Arc_vec_new(NA);
      /* For each old arc {e = oldArc[i]} of {old} we have a 
        corresponding new arc {d(e) = arc[i]}, with endpoints
        {v7(e) -> v7(Sym(e))}. These arcs are such that 
        {d(Sym(e)) = Sym(d(e))}; i.e., the edges of those arcs
        correspond to the edges of {old}, and have the same edge numbers.
        In addition, for each arc {e} of the old triangulation
        we have six new distinct edges (i.e. twelve new arcs).
        If {e == Lnext^r(oldSide[i])}, they are
        
          {a(e) = arc[oldNA + 12*(3*i + r%3) +  0]: v7(e) -> v5(e)} 
          {b(e) = arc[oldNA + 12*(3*i + r%3) +  2]: v5(e) -> v7(Lnext(e))}
          {c(e) = arc[oldNA + 12*(3*i + r%3) +  4]: v5(e) -> v5(Lnext(e))} 
          {f(e) = arc[oldNA + 12*(3*i + r%3) +  6]: v7(e) -> v7(Oprev(e))}
          {g(e) = arc[oldNA + 12*(3*i + r%3) +  8]: v7(e) -> v5(Sym(e))}
          {h(e) = arc[oldNA + 12*(3*i + r%3) + 10]: Org(e) -> v7(e)} */   
          
    Arc_vec_t out = SPQuad_Arc_vec_new(NV);
      /* The new canonical edges out of each old vertex {v},
        with number {v->num} in {0..oldNV-1}, is by definition
        {out[v->num] = h(oldOut[v->num])}.  For the new vertices,
        we have {out[v5(e)->num] = b(e)} and 
        {out[v7(e)->num] = a(e)} */
        
    Arc_vec_t side;
      /* The faces are found by {SPTriang_CreateTriangles} */ 
    
    auto Site *MakeNewSite(nat_t vnum, Arc e, r3_t *b);
      /* Create a new site with number {vnum} and barycentric coordinates
         {b} relative to the corners of {Left(e)}.
         Stores its address in {site.e[vnum]}. */
    
    auto void set_Onext(Arc a, Arc b);
      /* Splices so that {Onext(a) = b}. */

    auto Site *v5(Arc e, int r);
    auto Site *v7(Arc e, int r);
      /* The new 5- and 7-fold vertices in {Left(e)} closest to {Lnext^r(e)}. */
      
    auto Arc a(Arc e, int r);
    auto Arc b(Arc e, int r);
    auto Arc c(Arc e, int r);
    auto Arc d(Arc e, int r);
    auto Arc f(Arc e, int r);
    auto Arc g(Arc e, int r);
    auto Arc h(Arc e, int r);
      /* The new arcs associated with edge {Lnext^r(e)} */
      
    auto Arc spin(Arc e, int k);
      /* The arc {Lnext^k(e)}. */
      
    auto nat_t findex(Arc e);
      /* Smallest {r} such that {Lnext^r(side[fn]) == e},
        where {fn == Left(e)->num}. */
      
    auto nat_t eindex(Arc e);
      /* Smallest {r} such that {Sym^r(arc[2*en]) == e},
        where {en == EdgeNum(e)}. */
      
    Arc spin(Arc e, int k)
      { while(k > 0) { e = Lnext(e); k--; }
        return e;
      }

    nat_t findex(Arc e)
      { Arc b = oldSide.e[Left(e)->num];
        int r = 0;
        while (b != e) 
          { b = Lnext(b); r++;
            affirm(r < 3, "cannot find base arc");
          }
        return r; 
      }
    
    nat_t eindex(Arc e)
      { Arc b = oldArc.e[2*EdgeNum(e)];
        int r = 0;
        while (b != e) 
          { b = Sym(b); r++;
            affirm(r < 2, "cannot find base arc");
          }
        return r; 
      }
    
    Site *v7(Arc e, int r)
      { e = spin(e, r);
        return site.e[oldNV + 2*(3*Left(e)->num + findex(e))];
      }

    Site *v5(Arc e, int r)
      { e = spin(e, r);
        return site.e[oldNV + 2*(3*Left(e)->num + findex(e)) + 1];
      }
      
    Arc d(Arc e, int r)
      { e = spin(e, r);
        return arc.e[2*EdgeNum(e) + eindex(e)];
      }
      
    Arc a(Arc e, int r)
      { e = spin(e, r);
        return arc.e[oldNA + 12*(3*Left(e)->num + findex(e)) + 0];
      }

    Arc b(Arc e, int r)
      { e = spin(e, r);
        return arc.e[oldNA + 12*(3*Left(e)->num + findex(e)) + 2];
      }
    Arc c(Arc e, int r)
      { e = spin(e, r);
        return arc.e[oldNA + 12*(3*Left(e)->num + findex(e)) + 4];
      }
    Arc f(Arc e, int r)
      { e = spin(e, r);
        return arc.e[oldNA + 12*(3*Left(e)->num + findex(e)) + 6];
      }
    Arc g(Arc e, int r)
      { e = spin(e, r);
        return arc.e[oldNA + 12*(3*Left(e)->num + findex(e)) + 8];
      }
    Arc h(Arc e, int r)
      { e = spin(e, r);
        return arc.e[oldNA + 12*(3*Left(e)->num + findex(e)) + 10];
      }

    Site *MakeNewSite(nat_t vnum, Arc e, r3_t *b)
      { Site *v = (Site*)notnull(malloc(sizeof(Site)), "out of mem");
        S2Point *p0 = &(Org(e)->pos);
        S2Point *p1 = &(Org(Lnext(e))->pos);
        S2Point *p2 = &(Org(Lprev(e))->pos);
        double w0 = b->c[0];
        double w1 = b->c[1];
        double w2 = b->c[2];
        int k;
        v->num = vnum;
        site.e[vnum] = v;
        for (k = 0; k < 3; k++)
          { v->pos.c[k] = w0*p0->c[k] + w1*p1->c[k] + w2*p2->c[k]; }
        r3_dir(&(v->pos), &(v->pos));
        return v;
      }
    
    void set_Onext(Arc a, Arc b)
      { SPQuad_Splice(a, Oprev(b)); }

    nat_t vn, en, an, fn;
    
    /* CREATE VERTICES */
    
    /* Copy old vertices: */
    for (vn = 0; vn < oldNV; vn++)
      { site.e[vn] = Org(oldOut.e[vn]); }
    /* Create new vertices associated with old faces: */
    for (fn = 0; fn < oldNF; fn++)
      { Arc e = oldSide.e[fn];
        int r;
        for (r = 0; r < 3; r++)
          { /* Create vertices associated to arc {e} */
            nat_t v7num = oldNV + 2*(3*fn + r);
            nat_t v5num = v7num + 1;
            Site *v7e = MakeNewSite(v7num, e, b7);
            Site *v5e = MakeNewSite(v5num, e, b5);
            affirm(v7(e,0) == v7e, "inconsistent v7");
            affirm(v5(e,0) == v5e, "inconsistent v5");
            e = Lnext(e);
          }
      }

    /* CREATE EDGES, SET ORG/DEST */
    
    /* Create new edges {d(e)} associated to old edges: */
    for (en = 0; en < oldNA; en += 2)
      { Arc e = oldArc.e[en];
        Arc x = SPQuad_MakeEdge();
        EdgeNum(x) = en/2;
        arc.e[en] = x;
        arc.e[en+1] = Sym(x);
        /* Set its {Org} and {Dest} fields: */
        Odata(d(e,0)) = v7(e,0); 
        Ddata(d(e,0)) = v7(Sym(e),0);
      }
    /* Create 18 new edges per face: */
    for (fn = 0; fn < oldNF; fn++)
      { Arc e = oldSide.e[fn];
        int r;
        for (r = 0; r < 3; r++)
          { /* Create vertices associated to side {e} of face {fn}: */
            nat_t anum = oldNA + 12*(3*fn + r);
            int k;
            /* Create and store the new edges {a(e),b(e),c(e), f(e),g(e),h(e)}: */
            for (k = 0; k < 6; k++)
              { Arc x = SPQuad_MakeEdge();
                EdgeNum(x) = anum/2 + k;
                arc.e[anum + 2*k] = x;
                arc.e[anum + 2*k + 1] = Sym(x);
              }
              
            /* Set {Org(e)} and {Dest(e)} for all these new edges: */
            Odata(a(e,0)) = v7(e,0);  Ddata(a(e,0)) = v5(e,0);         
            Odata(b(e,0)) = v5(e,0);  Ddata(b(e,0)) = v7(e,1);  
            Odata(c(e,0)) = v5(e,0);  Ddata(c(e,0)) = v5(e,1);  
            Odata(f(e,0)) = v7(e,0);  Ddata(f(e,0)) = v7(Sym(e),1);  
            Odata(g(e,0)) = v7(e,0);  Ddata(g(e,0)) = v5(Sym(e),0);    
            Odata(h(e,0)) = Org(e);   Ddata(h(e,0)) = v7(e,0); 
            
            e = Lnext(e);
          }

     }

    /* SET ONEXT POINTERS */
    
    /* Set the {Onext} fields for the arcs {d(e)}: */
    for (en = 0; en < oldNA; en += 2)
      { Arc e = oldArc.e[en];
        Arc x = d(e,0);
        set_Onext(x, a(e,0));
        set_Onext(Sym(x), a(Sym(e),0));
      }
    /* Set the {Onext} fields for the other arcs: */
    for (an = 0; an < oldNA; an++)
      { Arc e = oldArc.e[an];
        { Arc x = a(e,0);
          set_Onext(x, Sym(b(e,2)));
          set_Onext(Sym(x), Sym(g(Sym(e),0)));
        }
        { Arc x = b(e,0);
          set_Onext(x, c(e,0));
          set_Onext(Sym(x), Sym(f(Sym(e),0)));
        }
        { Arc x = c(e,0);
          set_Onext(x, Sym(c(e,2)));
          set_Onext(Sym(x), Sym(a(e,1)));
        }
        { Arc x = f(e,0);
          set_Onext(x, g(e,0));
          set_Onext(Sym(x), Sym(h(Oprev(e),0)));
        }
        { Arc x = g(e,0);
          set_Onext(x, d(e,0));
          set_Onext(Sym(x), b(Sym(e),0));
        }
        { Arc x = h(e,0);
          set_Onext(x, h(Onext(e),0));
          set_Onext(Sym(x), f(e,0));
        }
      }

    /* DEFINE OUT POINTERS */
    
    /* Define new canonical edges for old vertices: */
    for (vn = 0; vn < oldNV; vn++)
      { out.e[vn] = h(oldOut.e[vn],0); }
    /* Define new canonical edges for new vertices: */
    for (an = 0; an < oldNA; an++)
      { Arc e = oldArc.e[an];
        out.e[v7(e,0)->num] = a(e,0);
        out.e[v5(e,0)->num] = b(e,0);
      }

    /* Check if all {Org} pointers are correctly set: */
    affirm(out.ne == NV, "inconsistent vertex count");
    fprintf(stderr, "vertex degrees: ");
    for (vn = 0; vn < NV; vn++)
      { Arc e = out.e[vn];
        Arc a = e;
        int order = 0;
        do 
          { affirm(Org(a)->num == vn, "inconsistent vertex nums");
            a = Onext(a);
            order++;
          } 
        while (a != e);
        fprintf(stderr, " %d", order);
      }
    fprintf(stderr, "\n");
    
    /* CREATE TRIANGLE RECORDS AND FACE TABLE */
    
    side = SPTriang_CreateFaceRecords(arc);
    affirm(side.ne == NF, "inconsistent face count");

    /* Check if all {Left} pointers are correctly set: */
    for (fn = 0; fn < NF; fn++)
      { Arc e = side.e[fn];
        Arc a = e;
        int order = 0;
        do 
          { affirm(Left(a)->num == fn, "inconsistent face nums");
            order++;
            a = Lnext(a);
          } 
        while (a != e);
        affirm(order == 3, "non-triangular faces");
      }
    
    /* Build triangle records and tables: */
    fprintf(stderr, "building tables\n");
    { /* OverlayTable finOvl = NULL; */
      void *trit = malloc(sizeof(Triangulation));
      Triangulation *tri = (Triangulation *)notnull(trit, "no mem for triangulation");
      tri->arc = arc;
      tri->out = out;
      tri->side = side;
      SPTriang_ComputeGeometryDataOfFaces(tri, old->smpOrder);
      return tri;
    }
  }

void SPTriang_Recompute57Triang
  ( Triangulation *tri, 
    Triangulation *old, 
    r3_t *b5,
    r3_t *b7
  )
  { /* Vertices of {tri} are identified, and related to vertices of
      {old}, solely on the basis of their numbers, independently on
      the topology of {tri}.

      Namely if {eOld = Lnext^r(oldSide[i])} is any arc of {old}, then
      the new 7-fold vertex in that face that is nearest to {Org(e)}
      is {v7(eOld) = site[oldNV + 2*(3*i + r%3)]}; and the new 5-fold
      vertex closest to the midpoint of {e} is {v5(eOld) =
      site.e[oldNV + 2*(3*i + r%3) + 1]}.
      
      The new positions of {v7(eOld)} and {v5(eOld)} are re-computed
      relatively to the corners of {Left(eOld)}; thus the topology
      of {old} must be the same as seen by {SPTriang_Refine57Triang}. */
  
    /* Element tables of old triangulation: */
    Arc_vec_t oldOut = old->out;
    Arc_vec_t oldArc = old->arc;
    Arc_vec_t oldSide = old->side;
    /* Old element counts: */
    nat_t oldNV = oldOut.ne;
    nat_t oldNA = oldArc.ne;
    nat_t oldNF = oldSide.ne;
    /* Element tables of new triangulation: */
    Arc_vec_t out = tri->out;
    Arc_vec_t arc = tri->arc;
    Arc_vec_t side = tri->side;
    /* Total counts for the new triangulation: */
    nat_t NV = out.ne;
    nat_t NA = arc.ne;
    nat_t NF = side.ne;
    
    auto void AdjustSite(Site *v, Arc eOld, r3_t *b);
      /* Recomputes the position of vertex {v} so that it 
        gets barycentric coordinates {(w0,w1,w2)} relative
        to the corners of the old face {Left(eOld)}. */
    
    auto nat_t eindex(Arc e);
      /* Smallest {r} such that {Sym^r(arc[2*en]) == e},
        where {en == EdgeNum(e)}. */
      
    nat_t eindex(Arc e)
      { Arc b = oldArc.e[2*EdgeNum(e)];
        int r = 0;
        while (b != e) 
          { b = Sym(b); r++;
            affirm(r < 2, "cannot find base arc");
          }
        return r; 
      }
    
    void AdjustSite(Site *v, Arc eOld, r3_t *b)
      { S2Point *p0 = &(Org(eOld)->pos);
        S2Point *p1 = &(Org(Lnext(eOld))->pos);
        S2Point *p2 = &(Org(Lprev(eOld))->pos);
        double w0 = b->c[0];
        double w1 = b->c[1];
        double w2 = b->c[2];
        int k;
        for (k = 0; k < 3; k++)
          { v->pos.c[k] = w0*p0->c[k] + w1*p1->c[k] + w2*p2->c[k]; }
        r3_dir(&(v->pos), &(v->pos));
      }
    
    nat_t vn, fn;
    
    affirm(NV == oldNV + 2*oldNA, "inconsistent vertex count");
    affirm(NA == oldNA + 12*oldNA, "inconsistent arc count");
    affirm(NF == 7*oldNF + 2*oldNA, "inconsistent face count");
    
    /* Check old vertices: */
    for (vn = 0; vn < oldNV; vn++)
      { Site *vOld = Org(oldOut.e[vn]);
        Site *vNew = Org(out.e[vn]);
        affirm(vOld == vNew, "old/new vertex mismatch");
      }
    /* Recompute new vertices associated with old faces: */
    for (fn = 0; fn < oldNF; fn++)
      { Arc eOld = oldSide.e[fn];
        int r;
        for (r = 0; r < 3; r++)
          { /* Recompute new vertices associated to arc {eOld} */
            Site *v7e = Org(out.e[oldNV + 2*(3*fn + r)]);
            Site *v5e = Org(out.e[oldNV + 2*(3*fn + r) + 1]);
            AdjustSite(v7e, eOld, b7);
            AdjustSite(v5e, eOld, b5);
            eOld = Lnext(eOld);
          }
      }
    SPTriang_ComputeGeometryDataOfFaces(tri, old->smpOrder);
  }

Triangulation *SPTriang_Refine3Triang(Triangulation *old)
  { /* Element tables of old triangulation: */
    Arc_vec_t oldOut = old->out;
    Arc_vec_t oldArc = old->arc;
    Arc_vec_t oldSide = old->side;
    /* Old element counts: */
    nat_t oldNV = oldOut.ne;
    nat_t oldNA = oldArc.ne;
    nat_t oldNF = oldSide.ne;
    /* Element counts for new triangulation: */
    nat_t NV = oldNV + oldNF;
    nat_t NA = 3*oldNA;
    nat_t NF = 3*oldNF;
    /* Element tables for new triangulation: */
    SiteRef_vec_t site = SiteRef_vec_new(NV);
      /* The old vertices of {old} are also vertices of the new
        triangulation, with the original numbers {0..oldNV-1}. Inside
        each face {i} of {old} we add one new vertex of degree 6,
        {site[oldNV + i]}. For each old arc {e}, we define  
        {v6(e) = site[oldNV + Left(e)->num]}, the new vertex 
        inserted into the left face of {e}. */
        
    Arc_vec_t arc = SPQuad_Arc_vec_new(NA);
      /* For each old arc {e = oldArc[i]} there are three new arcs:
        
          {d(e) = arc[i]}, with endpoints {v6(Sym(e)) -> v6(e)}
          {a(e) = arc[oldNA + 2*i]}, with endpoints {Org(e)} and {v6(e)}.
          {b(e) = arc[oldNA + 2*i + 1] = Sym(a(e))}.
          
        Note that {Sym(d(e)) == d(Sym(e))}, but {Sym(a(e)) != a(Sym(e))}. */
          
    Arc_vec_t out = SPQuad_Arc_vec_new(NV);
      /* The new canonical edges out of each old vertex {v},
        with number {v->num} in {0..oldNV-1}, is by definition
        {out[v->num] = a(oldOut[v->num])}.  For the new vertex 
        in face {i}, we have {out[oldNV + i] = b(oldSide[i])}. */
        
    Arc_vec_t side;
      /* The faces are found by {SPTriang_CreateTriangles} */ 
    
    auto Site *MakeNewSite(nat_t vnum, Arc e);
      /* Create a new site with number {vnum}, at the cicumcenter
        of {Left(e)}. Stores its address in {site.e[vnum]}. */
    
    auto void set_Onext(Arc a, Arc b);
      /* Splices so that {Onext(a) = b}. */

    auto Site *v6(Arc e);
      /* The new 6-fold vertices in {Left(e)}. */
      
    auto Arc a(Arc e);
    auto Arc b(Arc e);
    auto Arc d(Arc e);
      /* The new arcs associated with edge {e} */
      
    auto nat_t eindex(Arc e);
      /* Smallest {r} such that {Sym^r(arc[2*en]) == e},
        where {en == EdgeNum(e)}. */
      
    nat_t eindex(Arc e)
      { Arc b = oldArc.e[2*EdgeNum(e)];
        int r = 0;
        while (b != e) 
          { b = Sym(b); r++;
            affirm(r < 2, "cannot find base arc");
          }
        return r; 
      }
    
    Site *v6(Arc e)
      { return site.e[oldNV + Left(e)->num]; }

    Arc d(Arc e)
      { return arc.e[2*EdgeNum(e) + eindex(e)]; }
      
    Arc a(Arc e)
      { return arc.e[oldNA + 2*(2*EdgeNum(e) + eindex(e))]; }

    Arc b(Arc e)
      { return Sym(a(e)); }

    Site *MakeNewSite(nat_t vnum, Arc e)
      { Site *v = (Site*)notnull(malloc(sizeof(Site)), "out of mem");
        S2Point *p0 = &(Org(e)->pos);
        S2Point *p1 = &(Org(Lnext(e))->pos);
        S2Point *p2 = &(Org(Lprev(e))->pos);
        r3_t bar;
        v->num = vnum;
        site.e[vnum] = v;
        /* Compute barycenter of flat triangle */
        r3_add(p0, p1, &bar); r3_add(&bar, p2, &bar);
        r3_dir(&bar, &(v->pos));
        return v;
      }
    
    void set_Onext(Arc a, Arc b)
      { SPQuad_Splice(a, Oprev(b)); }

    nat_t vn, an, fn;
    
    /* CREATE VERTICES */
    
    /* Copy old vertices: */
    for (vn = 0; vn < oldNV; vn++)
      { site.e[vn] = Org(oldOut.e[vn]); }
    /* Create new vertices associated with old faces: */
    for (fn = 0; fn < oldNF; fn++)
      { Arc e = oldSide.e[fn];
        nat_t v6num = oldNV + fn;
        Site *v6e = MakeNewSite(v6num, e);
        affirm(v6(e) == v6e, "inconsistent v6");
      }

    /* CREATE EDGES, SET ORG/DEST */
    
    /* Create a new edge {d(e)} for each old *edge* {e}: */
    for (an = 0; an < oldNA; an += 2)
      { Arc e = oldArc.e[an];
        Arc x = SPQuad_MakeEdge();
        EdgeNum(x) = an/2;
        arc.e[an] = x;
        arc.e[an+1] = Sym(x);
        affirm(d(e) == x, "inconsistent d()");
        affirm(d(Sym(e)) == Sym(x), "inconsistent d(Sym())");
        Odata(d(e)) = v6(Sym(e)); 
        Ddata(d(e)) = v6(e);
      }
      
    /* Create one new edge {a(e)} for each *arc*  {e}: */
    for (an = 0; an < oldNA; an++)
      { Arc e = oldArc.e[an];
        nat_t anum = oldNA + 2*an;
        Arc x = SPQuad_MakeEdge();
        EdgeNum(x) = anum/2;
        arc.e[anum] = x;
        arc.e[anum + 1] = Sym(x);
        affirm(a(e) == x, "inconsistent a()");
        affirm(b(e) == Sym(x), "inconsistent b()");
        Odata(a(e)) = Org(e); 
        Odata(b(e)) = v6(e);         
        e = Lnext(e);
     }

    /* SET ONEXT POINTERS */
    
    for (an = 0; an < oldNA; an++)
      { Arc e = oldArc.e[an];
        /* Set the {Onext} fields for the arcs {Sym(d(e))}: */
        { Arc x = Sym(d(e));
          set_Onext(x, b(Lnext(e)));
        }
        /* Set the {Onext} fields for the arcs {a(e),b(e)}: */
        { Arc x = a(e);
          set_Onext(x, a(Onext(e)));
          set_Onext(Sym(x), Sym(d(e)));
        }
      }

    /* DEFINE THE OUT POINTERS */
    
    /* Define new canonical edges for old vertices: */
    for (vn = 0; vn < oldNV; vn++)
      { out.e[vn] = a(oldOut.e[vn]); }
    /* Define new canonical edges for new vertices: */
    for (fn = 0; fn < oldNF; fn++)
      { Arc e = oldSide.e[fn];
        out.e[v6(e)->num] = b(e);
      }
    affirm(NV = out.ne, "inconsistent vertex count");
    
    /* CREATE TRIANGLE RECORDS AND FACE TABLE */
    
    /* Build triangle records: */
    side = SPTriang_CreateFaceRecords(arc);
    affirm(side.ne == NF, "inconsistent face count");
    
    /* Build coordinate mapping tables: */
    fprintf(stderr, "building tables\n");
    { /* OverlayTable finOvl = NULL; */
      void *trit = malloc(sizeof(Triangulation));
      Triangulation *tri = (Triangulation *)notnull(trit, "no mem for triangulation");
      tri->arc = arc;
      tri->out = out;
      tri->side = side;
      SPTriang_ComputeGeometryDataOfFaces(tri, old->smpOrder);
      SPTriang_CheckTopology(tri, TRUE);
      return tri;
    }
  }

void SPTriang_CookValues(double_vec_t x, double tiny, double huge)
  { int i;
    double mid = sqrt(tiny*huge);
    double eps = tiny/mid;
    double norm = (1.0 + sqrt(1.0 + eps*eps));
    for (i = 0; i < x.ne; i++) 
      { double xi = x.e[i];
        /* Should do it more smoothly... */
        if (fabs(xi) < mid) 
          { double zi = xi/mid;
            xi = (zi + sqrt(zi*zi + eps*eps))/norm;
          }
        else
          { double zi = mid/xi;
            xi = norm/(zi + sqrt(zi*zi + eps*eps));
          }
        x.e[i] = mid*xi;
      }
  }

double SPTriang_GeomAvg(double_vec_t x)
  { double sum = 0.0;
    int i;
    for (i = 0; i < x.ne; i++) 
      { double xi = x.e[i];
        affirm(xi > 0.0, "bad x[i]");
        sum += log(xi);
      }
    return exp(sum/((double)x.ne));
  }

double SPTriang_Max(double_vec_t x)
  { double max = -INFINITY;
    int i;
    for (i = 0; i < x.ne; i++) 
      { double xi = x.e[i]; 
        if (xi > max) { max = xi; }
      }
    return max;
  }

double SPTriang_Min(double_vec_t x)
  { double min = INFINITY;
    int i;
    for (i = 0; i < x.ne; i++) 
      { double xi = x.e[i]; 
        if (xi < min) { min = xi; }
      }
    return min;
  }

double SPTriang_RatioVar(double_vec_t x)
  { double avg = SPTriang_GeomAvg(x);
    double sum = 0.0;
    int i;
    for (i = 0; i < x.ne; i++)
      { double xi = x.e[i];
        double d = (xi/avg + avg/xi)/2.0 - 1.0;
        sum += d;
      }
    return sum/((double)x.ne);
  }

double SPTriang_RatioSpr(double_vec_t x)
  { double max = -INFINITY;
    double min = INFINITY;
    int i;
    for (i = 0; i < x.ne; i++) 
      { double xi = x.e[i];
        if (xi > max) { max = xi; }
        if (xi < min) { min = xi; }
      }
    return (min/max + max/min)/2.0 - 1.0;
  }

double_vec_t SPTriang_FaceAreas(Triangulation *tri)
  { nat_t NF = tri->side.ne;
    double_vec_t area = double_vec_new(NF);
    nat_t fn;
    for (fn = 0; fn < NF; fn++)
      { Arc e = tri->side.e[fn];
        S2Point *p = &(Dest(Onext(e))->pos);
        S2Point *q = &(Org(e)->pos);
        S2Point *r = &(Dest(e)->pos);
        area.e[fn] = TriangleArea(p, q, r);
      }
    return area;
  }
    
double_vec_t SPTriang_EdgeLengths(Triangulation *tri)
  { 
    nat_t NE = tri->arc.ne/2;
    double_vec_t length = double_vec_new(NE);
    nat_t en;
    for (en = 0; en < NE; en++)
      { Arc e = tri->arc.e[2*en];
        S2Point *p = &(Org(e)->pos);
        S2Point *q = &(Dest(e)->pos);
        length.e[en] = r3_dist(p, q);
      }
    return length;
  }
  
double_vec_t SPTriang_RelAngles(Triangulation *tri)
  { 
    nat_t NA = tri->arc.ne;
    double_vec_t angle = double_vec_new(NA);
    nat_t an;
    
    auto double Angle(R3Vector *o, R3Vector *p, R3Vector *v);
      /* Angle beween linear subspaces generated
         by {(o,u)} and {(o,v)}. */
         
    double Angle(R3Vector *o, R3Vector *p, R3Vector *v)
      { R3Vector np, nv;
        r3_cross(o, p, &np);
        r3_cross(o, v, &nv);
        return acos(r3_cos(&np, &nv));
      }
    
    for (an = 0; an < NA; an++)
      { Arc e = tri->arc.e[an];
        S2Point *o = &(Org(e)->pos);
        S2Point *p = &(Dest(e)->pos);
        S2Point *q = &(Dest(Onext(e))->pos);
        double ideal = 2.0*PI/((double)OrgDegree(e));
        angle.e[an] = Angle(o, p, q)/ideal;
      }
    return angle;
  }

double_vec_t SPTriang_CrossRatios(Triangulation *tri)
  { nat_t NA = tri->arc.ne;
    double_vec_t xratio = double_vec_new(NA);
    int i;
    
    for (i = 0; i < NA; i++)
      { Arc e = tri->arc.e[i];
        Face *f = Left(e); 
        S2Point *p = &(Dest(Oprev(e))->pos);
        Arc fref = tri->side.e[f->num];
        r3_t bar;
        /* Compute bary coords of {p} relative to corners of {f}: */
        r3x3_map_col(&(f->c2b), p, &bar);
        /* Rename bary coords relative to {e}, not {fref}: */
        if (e == fref)
          { /* OK */ }
        else if (e == Lnext(fref))
          { bar = (r3_t){{bar.c[1], bar.c[2], bar.c[0]}}; }
        else if (fref == Lnext(e))
          { bar = (r3_t){{bar.c[2], bar.c[0], bar.c[1]}}; }
        else
          { affirm(FALSE , "face is not a triangle"); }
        /* Now return required cross ratio: */
        xratio.e[i] = bar.c[2]/(bar.c[1] + bar.c[2]);
      }
    return xratio;
  }

double_vec_t SPTriang_RelOppAngles(Triangulation *tri)
  { 
    nat_t NA = tri->arc.ne;
    nat_t NV = tri->out.ne;
    double_vec_t oppAngle = double_vec_new(NA);
    nat_t vn, an = 0;
    
    auto double OppAngle(R3Vector *o, R3Vector *u, R3Vector *v);
      /* External angle beween linear subspaces generated
         by {(o,u)} and {(v,o)}. */
         
    double OppAngle(R3Vector *o, R3Vector *u, R3Vector *v)
      { R3Vector np, nv;
        r3_cross(o, u, &np);
        r3_cross(v, o, &nv);
        return acos(r3_cos(&np, &nv));
      }
    
    for (vn = 0; vn < NV; vn++)
      { Arc estart = tri->out.e[vn];
        S2Point *o = &(Org(estart)->pos);
        double ideal = PI/((double)OrgDegree(estart));
        Arc e = estart;         
        Arc a = Onext(e);
        S2Point *u = &(Dest(a)->pos);
        Arc b = Onext(Onext(e));
        S2Point *v = &(Dest(b)->pos); 
        do 
          { S2Point *p = &(Dest(e)->pos);
            /* Take care that the test outcomes for {e,b} 
              is identical to the test {e,a} when {e} becomes {b}
              and {a} becomes {e}: */ 
            double du = r3_det(o, p, u);
            double dv = r3_det(o, v, p);
            affirm(r3_det(o, u, v) > 0.0, "degenerate triangle");
            /* Find consecutive arcs {a,b} out of {Org(e)}
              that bracket the opposite direction: */
            while ((a == e) || (dv <= 0.0))
              { a = b; u = v; du = -dv;
                b = Onext(b);
                v = &(Dest(b)->pos); 
                /* Take care that the test outcomes for {e,b} 
                  is identical to the test {e,a} when {e} becomes {b}
                  and {a} becomes {e}: */ 
                dv = r3_det(o, v, p);
                affirm(b != e, "degenerate vertex");
                affirm(r3_det(o, u, v) > 0.0, "degenerate triangle");
              }
            affirm(du >= 0, "degenerate vertex");
            oppAngle.e[an] = OppAngle(o, p, u) / ideal; an++;
            e = Onext(e);
          }
        while (e != estart);
      }
    affirm(an == NA, "inconsistent arc count");
    return oppAngle;
  }
