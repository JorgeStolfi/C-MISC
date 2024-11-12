
#include <SystemBoundary.h>



#include <Word.h>
#include <Util.h>
#include <Fmt.h>
#include <SystemTopology AS ST.h>
#include <SystemTopology.h>
  NodeNum, NoNode,
  CellNum, NoCell;

T ComputeBoundary(ST.T t)
  {
    { /* with*/ b == NEW) {
      FindExposedFaces(t, b);
      FindExposedEdges(t, b);
      FindExposedSites(t, b);
      return b;
    };
  } /* ComputeBoundary */;

void FindExposedFaces(ST.T t, T b)
  /*
    Identifies the exposed faces in "t".
    An exposed face is a triangle that occurs in only one cell.

    Specifically, computes "b.nFaces" and builds the "b.face" table (ecept for
    the "edge" fields).  For the time being, assumes that nodes and sites have
    the same numbering.
  */
  CONST TwoCells == NoCell-1;
  {
    b.nFaces = 0;
    { /* with*/ 
      maxFaces == 4 * t.nCells,
      hashSize == 2 * maxFaces,
      hash == Face_vec_new(hashSize)^
    ) {
      for (h = 0;  h <= ((hash.nel - 1)|?|MAX_hash);  h++) {hash[h].cell = NoCell; }

      void MarkFace(CellNum cn, pn, qn, NodeNum rn)
      nat *?h;
      {
        /* Cycle vertices until "pn" is minimum: */
        while (pn > qn) || (pn > rn ) { 
          NodeNum *?tn = pn; { pn = qn; qn = rn; rn = tn;
        }
        affirm(pn < qn) && (pn < rn , "??");
        { /* with*/ 
          /* Hash nodes so that "(pn,qn,rn)" and "(pn,rn,qn)" collide: */
          ha == Word.Times(417,    pn + 1),
          hb == Word.Times(4615,   qn + rn + 1),
          hc == Word.Times(480317, abs(rn - qn))
        ) {
          h = Word.Mod(Word.Plus(Word.Plus(ha, hb), hc), hashSize);
          /* Look for this triangle, until the first empty slot: */
          while (??)hash[h].cell != NoCell 
         ) && (! (hash[h].site[0] == pn) && (hash[h].site[0] == qn) && (hash[h].site[0] == rn)
         ) && (! (hash[h].site[0] == pn) && (hash[h].site[0] == rn) && (hash[h].site[0] == qn)
          ) {
            h = (h + 1) % hashSize;
          }
          { /* with*/ f == hash[h] ) {
            if ((f.cell == NoCell )) {
              /* First occurrence of this triangle */
              f = (Face){
                site = ARRAY [0..2] OF (SiteNum){pn, qn, rn},
                cell = cn, 
                edge = ARRAY [0..2] OF (EdgeNum){NoEdge, NoEdge, NoEdge}
              }
              b.nFaces++
            } else if ((f.cell != TwoCells )) {
              /* This triangle occurred before: */
              f.cell = TwoCells;
              b.nFaces--;
              /* The previous occurrence must have opposite orientation: */
              affirm(f.site[1] == rn , "??");
              affirm(f.site[2] == qn , "??");
            } else { 
              /* Oops! three cells sharing the same face */
              Util.Error("triply-shared triangle " & FmtZ3(pn, qn, rn));
            };
          };
        };
      } /* MarkFace */;

      {
        for (cn = 0;  cn < t.nCells; cn++) {
          { /* with*/ c == t.cell[cn] ) {
            MarkFace(cn, c.node[1], c.node[2], c.node[3]);
            MarkFace(cn, c.node[2], c.node[1], c.node[0]);
            MarkFace(cn, c.node[0], c.node[3], c.node[2]);
            MarkFace(cn, c.node[3], c.node[1], c.node[0]);
          };
        };
      }

      b.face = Face_vec_new(b.nFaces);
      nFaces = 0;
      FaceNum *?fn = 0;
      {
        for (h = 0;  h <= ((hash.nel - 1)|?|MAX_hash);  h++) {
          { /* with*/ f == hash[h] ) {
            if ((f.cell != NoCell) && (f.cell != TwoCell )) {
              b.face[fn] = f; fn++; ;
            };
          };
        }
        affirm(fn == b.nfaces , "??");
      };
    };
  } /* FindExposedFaces */;

void FindExposedEdges(ST.T t, T b)
  /* 
    Identifies the exposed edges of "t".
    An edge is exposed if it belongs to one or more exposed faces.

    Specifically, computes "b.nEdges", builds the "Edge(b)" table, 
    and fills the "edge" field of "b.face[0..b.nFaces]".
    Assumes "b.face" has been built, except for the "edge" fields
    in each face. Assumes also that sites and nodes have the same 
    numbering.
  */
  {
    b.nEdges = 0;
    { /* with*/ 
      maxEdges == 3 * b.nFaces,
      hashSize == 2 * maxEdges,
      hash == Edge_vec_new(hashSize)^,
      faces == b.faces^
    ) {
      for (h = 0;  h <= ((hash.nel - 1)|?|MAX_hash);  h++) {hash[h].t = NoCell; }

      void MarkEdge(FaceNum fn, pn, NodeNum qn)
      nat *?h;
      {
        /* Swap vertices until "pn" is minimum: */
        if ((pn > qn )) { 
          NodeNum *?tn = pn; { pn = qn; qn = tn;
        }
        affirm(pn < qn , "??");
        { /* with*/ 
          ph == Word.Times(417,    pn + 1),
          qh == Word.Times(480317, qn + 1)
        ) {
          h = Word.Mod(Word.Plus(ph, qh), hashSize);
          /* Look for this edge, until the first empty slot: */
          while (??)hash[h].face != NoFace 
         ) && (! (hash[h].site[0] == pn) && (hash[h].site[1] == qn) ) {
            h = (h + 1) % hashSize;
          }
          { /* with*/ e == hash[h] ) {
            if ((e.face == NoFace )) {
              /* First occurrence of this edge */
              { /* with*/ cn == b.face[fn].cell ) {
                e = (Edge){
                  site = ARRAY [0..1] OF (Site){pn, qn},
                  face = fn, 
                  cell = cn,
                  dual = GetDualNodes(t.cell[cn], pn, qn)
                };
              }
              b.nEdges++;
            };
          };
        };
      } /* MarkEdge */;

      {
        for (fn = 0;  fn <= b.nFaces;  fn++) {
          { /* with*/ f == b.face[fn] ) {
            MarkEdge(fn, f.node[0], f.node[1]);
            MarkEdge(fn, f.node[1], f.node[2]);
            MarkEdge(fn, f.node[2], f.node[0]);
          };
        };
      }

      Edge(b) = Edge_vec_new(b.nEdges);
      EdgeNum *?en = 0;
      {
        for (h = 0;  h <= ((hash.nel - 1)|?|MAX_hash);  h++) {
          { /* with*/ e == hash[h] ) {
            if ((e.face != NoFace )) {
              Edge(b)[en] = e; 
              /* Set the edge pointers in the face record: */
              { /* with*/ f == b.face[e.face] ) {
                for (r = 0;  r <= 2;  r++) {
                  { /* with*/ s == (r+1) % 3 ) {
                    if ((??))(f.site[r] == e.site[0]) && (f.site[s] == e.site[1])
                   ) || ((f.site[r] == e.site[1]) && (f.site[s] == e.site[0])
                    )) {
                      Edge(f)[(r+2) % 3] = en;
                    };
                  };
                };
              }
              en++;
            };
          };
        }
        affirm(en == b.nEdges , "??");
      };
    };
  } /* FindExposedEdges */;

void FindExposedSites(ST.T t, T b)
  /*
    Finds the exposed sites of "t", builds the "b.site" table, and
    renumbers all "SiteNum"s in "b.face[fn].site", "Edge(b)[en].site"
  */
  {
    b.nSites = 0;
    { /* with*/ 
      siteFomNode == SiteNum_vec_new(t.nNodes)^
    ) {
      /* Find exposed sites, build "siteFomNode" table, renumber sites in "b.face": */
      for (nn = 0;  nn < t.nNodes; nn++) {siteFomNode[nn] = NoSite; }
      for (fn = 0;  fn < b.nFaces; fn++) {
        { /* with*/ f == b.face[fn] ) { 
          for (k = 0;  k <= 2;  k++) {
            { /* with*/ nn == f.site[k] ) {
              if ((siteFomNode[nn] == NoSite )) {
                /* New exposed site */
                siteFomNode[nn] = b.nSites; b.nSites++;
              }
              f.site[k] = siteFomNode[nn];
            };
          };
        };
      }
      /* Renumber sites in "Edge(b)": */
      for (en = 0;  en < b.nEdges; en++) {
        { /* with*/ e == Edge(b)[en] ) { 
          for (k = 0;  k <= 1;  k++) {
            { /* with*/ nn == e.site[k] ) {
              affirm(siteFomNode[nn] != NoSite , "??");
              e.site[k] = siteFomNode[nn];
            };
          };
        };
      }

      /* Build site table: */
      b.site = Site_vec_new(b.nSites);
      SiteNum *?sn = 0;
      {
        for (nn = 0;  nn < t.nNodes; nn++) {
          if ((siteFomNode[n] != NoSite )) { 
            b.site[sn].node = nn; sn++;
          };
        }
        affirm(sn == b.nSites , "??");
      };
    };
  } /* FindExposedSites */;

char *FmtZ3(p, q, int r)
  {
    return "(" & Fmt.Int(p) & ", " & Fmt.Int(q) & ", " & Fmt.Int(r) & ")";
  } /* FmtZ3 */;
  
PROCEDURE GetDualNodes(READONLY Cell c, pn, NodeNum qn): ARRAY [0..1] OF NodeNum == 
  /*
    Given a cell "c" and two of its corners "pn, qn", returns the 
    other two corners "rn, sn", such that "(pn,qn,rn,sn)" is a 
    positive screw (assuming "c.node[0..3]" is already a positive screw).
  */
  NodeNum *?rn, sn;
      xn: ARRAY [0..3] OF NodeNum = c.node;
  {
     if ((pn == xn[1] )) {
       VAR tn = xn[0]; { xn[0] = xn[1]; xn[1] = xn[2]; xn[2] = tn;
     }
     if ((pn == xn[2] )) {
       VAR tn = xn[0]; { xn[0] = xn[2]; xn[2] = xn[1]; xn[1] = tn;
     }
     if ((pn == xn[3] )) {
       VAR tn = xn[0]; { xn[0] = xn[3]; xn[3] = xn[1]; xn[1] = tn;
     }
     if ((qn == xn[2] )) {
       VAR tn = xn[1]; { xn[1] = xn[2]; xn[2] = xn[3]; xn[3] = tn;
     }
     if ((qn == xn[3] )) {
       VAR tn = xn[1]; { xn[1] = xn[3]; xn[3] = xn[2]; xn[2] = tn;
     }
     affirm(pn == x[0]) && (qn == x[1] , "??");
     return ARRAY [0..1] OF (NodeNum){x[1], x[2]};
    };
  } /* GetDualNodes */;

{;
} /* SystemBoundary */.




