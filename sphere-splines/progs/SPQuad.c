/* See SPQuad.h */
/* Last edited on 2006-02-28 12:18:52 by stolfi */
  
/* Written by J. Stolfi on april 1993, based on an original 
  implementation by Jim Roth (DEC CADM Advanced Group, May 1986).  
  See the copyright notice at the end of this file. */ 

#include <SPQuad.h>

#include <bool.h>
#include <nat.h>
#include <affirm.h>

#include <memory.h>
#include <malloc.h>

#define Mark(e)  (((SPQuad_EdgeRec *)((e)&0xfffffffcu))->mark)

SPQuad_Arc SPQuad_MakeEdge(void)
  { SPQuad_Arc e;
    e = (SPQuad_Arc)notnull(malloc(sizeof(SPQuad_EdgeRec)), "no mem for edge");
    Onext(e) = e;
    SymDnext(e) = Sym(e);
    RotRnext(e) = Tor(e);
    TorLnext(e) = Rot(e);
    Edge(e)->num = 0;
    Odata(e) = NULL; Ldata(e) = NULL;
    Ddata(e) = NULL; Rdata(e) = NULL;
    return e;
  }

void SPQuad_DestroyEdge(SPQuad_Arc e)
  { SPQuad_Arc f = Sym(e);
    if (Onext(e) != e) SPQuad_Splice(e, Oprev(e));
    if (Onext(f) != f) SPQuad_Splice(f, Oprev(f));  
    free((char *)((e) & 0xfffffffcu));
  }

void SPQuad_Splice(SPQuad_Arc a, SPQuad_Arc b)
  { SPQuad_Arc ta, tb;
    SPQuad_Arc alpha = Rot(Onext(a));
    SPQuad_Arc beta = Rot(Onext(b));

    ta = Onext(a);
    tb = Onext(b);
    Onext(a) = tb;
    Onext(b) = ta;
    ta = Onext(alpha);
    tb = Onext(beta);
    Onext(alpha) = tb;
    Onext(beta) = ta;    
  }

SPQuad_Arc_vec_t SPQuad_RenumberEdges(SPQuad_Arc_vec_t a, nat_t start)
  {
    auto void VisitEdge(SPQuad_Arc a);
      /* If {Edge(a)} has not been visited before, renumber it. */

    auto bool_t VisitedEdge(SPQuad_Arc a);
      /* TRUE if {Edge(a)} has been visited before. */
    
    #define GUESS_NVISITED 1024
    SPQuad_Arc_vec_t visited = SPQuad_Arc_vec_new(GUESS_NVISITED);
    int nVisited = 0;
      /* {visited.e[0..nVisited-1]} contain one arc from each visited
        edge. An edge record {e} has been visited if
        {Edge(visited.e[e->num - start]) == e}. */

    bool_t VisitedEdge(SPQuad_Arc a)
      { nat_t en = Edge(a)->num;
        return 
          (en >= start) && 
          (en < start + nVisited) && 
          (Edge(visited.e[en-start]) == Edge(a));
      }

    void VisitEdge(SPQuad_Arc a)
      { if (! VisitedEdge(a))
          { /* Edge {Edge(a)} hasn't been visited yet. */
            SPQuad_Arc_vec_expand(&visited, nVisited);
            Edge(a)->num = start + nVisited;
            visited.e[nVisited] = a;
            nVisited++;
          }
      }

    int i;
    int nClosed = 0; 
      /* Edges {Edge(visited.e[0..nClosed-1])} are the edges whose
        children were already visited.  The children of {e}, by 
        definition, are {Onext(e)} and {Onext(Sym(e))}. */

    /* Put the given arcs on the visited (minus repetitions): */
    for (i = 0; i < a.ne; i++) { VisitEdge(a.e[i]); }
    /* Visit descendants of visited edges in BFS order: */
    while (nClosed < nVisited )
      { SPQuad_Arc s = visited.e[nClosed];
        VisitEdge(Onext(s));
        VisitEdge(Onext(Sym(s)));
        nClosed++;
      }
    /* Return the edges which were found: */
    SPQuad_Arc_vec_trim(&visited, nVisited);
    return visited;
  }

/* Arrays of {SPQuad_Arc}: */

vec_typeimpl(SPQuad_Arc_vec_t,SPQuad_Arc_vec,SPQuad_Arc);

/*
  Copyright notice:

  Copyright 1996,2002 Institute of Computing, Unicamp.

  Permission to use this software for any purpose is hereby granted,
  provided that any substantial copy or mechanically derived version
  of this file that is made available to other parties is accompanied
  by this copyright notice in full, and is distributed under these
  same terms.

  NOTE: this copyright notice does not claim to supersede any
  copyrights that may apply to the original DEC implementation of the
  quad-edge data structure.

  DISCLAIMER: This software is provided {as is} with no explicit or
  implicit warranty of any kind. Neither the authors nor their
  employers can be held responsible for any losses or damages that
  might be attributed to its use.

  End of copyright notice. */
