/* The quad-edge data structure (oriented surface version). */
/* Last edited on 2006-02-28 11:50:51 by stolfi */

/* The quad-edge data structure encodes the topology of a graph drawn
  on an orientable compact 2-D manifold (i.e. a borderless surface of
  finite extent) in such a way that every face is a topological disk.
  For details, see

    "Primitives for the Manipulation of General Subdivisions 
    and the Computation of Voronoi Diagrams"
    L. Guibas, J. Stolfi, ACM Transactions on Graphics, April 1985

  This implementation is specialized for oriented surfaces --- it
  does not provide the {Flip} operator.

  Based on an implementation by Jim Roth (DEC CADM Advanced Group) on
  May 1986. Modified by J. Stolfi on April 1993, and November 2002.
  See the copyright notice at the end of this file. */

#ifndef SPQuad_H
#define SPQuad_H

#include <vec.h>
#include <js.h>
#include <nat.h>

/* ARCS */

typedef nat_t SPQuad_Arc; 
  /* A directed edge, primal or dual */

vec_typedef(SPQuad_Arc_vec_t,SPQuad_Arc_vec,SPQuad_Arc);

#define SPQuad_NullArc 0

/* EDGE RECORDS */

typedef struct SPQuad_EdgeRec /* Edge record (four arcs). */
  { SPQuad_Arc next[4]; /* Topological links. */
    void *data[4];      /* Client data fields. */
    nat_t num;   /* Edge number. */
  } SPQuad_EdgeRec;

#define Edge(e) ((SPQuad_EdgeRec *)((e)&0xfffffffcu))
  /* Address of the edge record for arc {e}. */

/* EDGE ORIENTATION OPERATORS */

#define Rot(e) (SPQuad_Arc)(((e)&0xfffffffcu)+(((e)+1)&3u))
#define Sym(e) (SPQuad_Arc)(((e)&0xfffffffcu)+(((e)+2)&3u))
#define Tor(e) (SPQuad_Arc)(((e)&0xfffffffcu)+(((e)+3)&3u))
  /* The arc orientation operators: */

/* WALKING OPERATORS */

#define RotRnext(e) (Edge(e))->next[((e)+1)&3]
#define SymDnext(e) (Edge(e))->next[((e)+2)&3]
#define TorLnext(e) (Edge(e))->next[((e)+3)&3]
  /* Auxiliary ops, used below. */

#define Onext(e) (Edge(e))->next[(e)&3]
#define Rnext(e) (Tor(RotRnext(e)))
#define Dnext(e) (Sym(SymDnext(e)))
#define Lnext(e) (Rot(TorLnext(e)))
  /* The quad-edge forward-walking operators. */

#define Oprev(e) (Rot(RotRnext(e)))
#define Dprev(e) (Tor(TorLnext(e)))
#define Rprev(e) (SymDnext(e))
#define Lprev(e) (Sym(Onext(e)))
  /* The quad-edge backwards-walking operators. */

/* DATA POINTERS */

#define Odata(e) ((Edge(e))->data[(e)&3])
#define Rdata(e) ((Edge(e))->data[((e)+1)&3])
#define Ddata(e) ((Edge(e))->data[((e)+2)&3])
#define Ldata(e) ((Edge(e))->data[((e)+3)&3])
  /* The data fields of the edge record which are
    associated with the origin, right face, destination,
    and left face of arc {e}, respectively. */

/* TRAVERSAL */

SPQuad_Arc_vec_t SPQuad_RenumberEdges(SPQuad_Arc_vec_t a, nat_t start);
  /* Enumerates undirected edges reachable from the arcs {a.e[..]}, and
    stores in their {num} fields consecutive integers beginning with {start}.
    Returns an array {res} where {res.e[i]} is one primally reachable
    arc from the edge with number {start+i}.

    An arc {b} is /primally reachable/ from an arc {a} iff it can be obtained 
    from {a} by a finite number of {Sym} and {Onext} operations (no {Rot}s
    or {Tor}s).  An edge {e} is /reachable/ iff some of its arcs is reachable. */

/* ORIENTATION BITS */

#define SymBit(e) (((e)&2)>>1)  
#define RotBit(e) ((e)&1)
#define QuadBits(e) ((e)&3)
  /* Bits that distinguish the various flavors of the same edge:
      {SymBit(e)} distinguishes {e} from {Sym(e)}.
      {RotBit(e)} distinguishes {{e,Sym(e)}} from {{Rot(e),Tor(e)}}.
      {QuadBits(e) = 2*SymBit(e) + RotBit(e)} distinguishes the four arcs
        {Rot^k(e)} from each other.
    Nothing can be assumed between {QuadBits(e)} 
    and {QuadBits(Onext(e))}. */

/* EDGE AND ARC NUMBERS */

#define EdgeNum(e) (Edge(e)->num)  
#define PrimalArcNum(e) (2*(Edge(e)->num)+SymBit(e))
#define ArcNum(e) (4*(Edge(e)->num)+QuadBits(e))
  /* Numbers for edges and arcs (derived from the edge's {num} field):
      {EdgeNum(e)} is the same for all four arcs {e} with same {Edge(e)}.
      {PrimalArcNum(e)} distinguishes {e} from {Sym(e)}.
      {ArcNum(e)} distingushes all four arcs.
  */

/* MODIFICATION OPERATORS */

SPQuad_Arc SPQuad_MakeEdge(void);
  /* Create a new spherical map with a single edge, a single face,
    and two vertices. */

void SPQuad_DestroyEdge(SPQuad_Arc e);
  /* Disconnects the edge {e} from its map,
    and frees its edge record. */

void SPQuad_Splice(SPQuad_Arc a, SPQuad_Arc b);
  /* The quad-edge splice operator. */

#endif

/*
  Copyright notice:

  Copyright 1996,2002 Institute of Computing, Unicamp.

  Permission to use this software for any purpose is hereby granted,
  provided that any substantial copy or mechanically derived version
  of this file that is made available to other parties is accompanied
  by this copyright notice in full, and is distributed under these same
  terms. 

  NOTE: this copyright notice does not claim to supersede any copyrights
  that may apply to the original DEC implementation of the quad-edge
  data structure.

  DISCLAIMER: This software is provided {as is} with no explicit or
  implicit warranty of any kind.  Neither the authors nor their
  employers can be held responsible for any losses or damages
  that might be attributed to its use.

  End of copyright notice.
*/
