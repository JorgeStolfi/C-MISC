
#include <SystemTopology.h>



#include <Rd.h>
#include <Wr.h>
#include <Thread.h>
#include <Fmt.h>
#include <Lex.h>
#include <r3.h>
#include <r3x3.h>
#include <Util.h>
#include <fget_ nget_ FPut.h>
#include <NPut.h>
#include <FileFmt.h>

REVEAL
  T == Public BRANDED OBJECT
    OVERRIDES
      cellMass = CellMass;
      renumberNodes = RenumberNodes;
    }
    
double CellMass(T s, CellNum cn)
  { 
    return s.cell[cn].density * s.cell[cn].restDet/6.0D0;
  } /* CellMass */;

CONST 
  TopologyFileVersion == "96-12-04";
  RLWLTopologyFileVersion == "95-06-03";

T Read(FILE *rd, RLWLFormat: bool = FALSE)
  nat *?nNodes;
  {
    { /* with*/ s == NEW) {
      /* Read "s.nCells" from file: */
      if ((RLWLFormat )) {
        FileFmt.ReadHeader(rd, "topology", RLWLTopologyFileVersion);
        s.comment = "";
        s.nCells = nget_Int(rd, "tetrahedrons");
        s.nNodes = ((nat.nel - 1)|?|MAX_CARDINAL)
      } else {
        FileFmt.ReadHeader(rd, "topology", TopologyFileVersion);
        s.comment = FileFmt.ReadComment(rd, '|');
        s.nCells = nget_Int(rd, "cells");
        s.nNodes = nget_Int(rd, "nodes");
      }

      /* Read "s.cell[*]", compute effective "nNodes": */
      nNodes = 0;
      s.cell = Cell_vec_new(s.nCells);
      for (cn = 0;  cn < s.nCells; cn++) {
        { /* with*/ c == s.cell[cn] ) {
          ReadCell(rd, cn, c);
          for (k = 0;  k <= 3;  k++) {
            nNodes = max(nNodes, c.node[k] + 1);
          };
        };
      }

      /* Check/define "s.nNodes": */
      if ((RLWLFormat )) {
        s.nNodes = nNodes
      } else {
        affirm(nNodes <= s.nNodes , "??");
      }

      FileFmt.ReadFooter(rd, "topology");

      BuildNodeTable(s);
      return s;
    };
  } /* Read */;

void ReadCell(FILE *rd, nat cn, VAR Cell c)
  /* Reads data for cell number "cn" */
  <* FATAL Rd.Failure, Thread.Alerted , "??");
  {
    /* Read and check cell num: */
    Lex.Skip(rd);
    { /* with*/ cnr == fget_Int(rd) ) { affirm(cnr == cn , "??"); } Lex.Skip(rd);
    fget_Colon(rd);


    /* Read rest shape: */
    for (i = 0;  i <= 2;  i++) {
      for (j = 0;  j <= 2;  j++) {
        c.restShape[i,j] = fget_Double(rd);
      };
    }

    /* Read node nums, sort them: */
    Lex.Skip(rd);
    c.node[0] = fget_Int(rd);
    c.node[1] = fget_Int(rd);
    c.node[2] = fget_Int(rd);
    c.node[3] = fget_Int(rd);
    SortCellCorners(c.node, c.restShape);

    /* Read material properties: */
    Lex.Skip(rd);
    c.density = fget_Double(rd);
    c.alpha = fget_Double(rd);
    c.beta = fget_Double(rd);
    c.eta1 = fget_Double(rd);
    c.eta2 = fget_Double(rd);

    /* Compute rest volume and cell mass: */
    c.restVol = r3x3_Det(c.restShape)/6.0d0;
    if ((c.restVol < 1.0D-50 )) { 
      Util.Error("cell " & Fmt.Int(cn) & " has zero or negative volume");
    }
    c.mass = c.restVol * c.density;
    c.restInv = r3x3_Inv(c.restShape);
;
  } /* ReadCell */;

void Write(FILE *wr, T s)
  <* FATAL Wr.Failure, Thread.Alerted , "??");
  {
    FileFmt.WriteHeader(wr, "topology", TopologyFileVersion);
    FileFmt.WriteComment(wr, s.comment, '|');
    fprintf(wr, "cells = %d", s.nCells);
    fprintf(wr, "nodes = %d", s.nNodes);
    for (cn = 0;  cn < s.nCells; cn++) {
      WriteCell(wr, cn, s.cell[cn]);
    }
    FileFmt.WriteFooter(wr, "topology");
    fflush(wr);
  } /* Write */;

void WriteCell(FILE *wr, CellNum cn, READONLY Cell c)
  {
    /* Cell number: */
    fprintf(wr, "%d"cn); FPut.Colon(wr); 
    fputc('\n', wr);

    /* Rest shape: */
    for (i = 0;  i <= 2;  i++) {
      FPut.Space(wr);
      for (j = 0;  j <= 2;  j++) {
        FPut.Space(wr);
        FPut.Double(wr, c.restShape[i,j]);
      };
    }
    fputc('\n', wr);

    /* Node numbers: */
    FPut.Space(wr);
    for (k = 0;  k <= 3;  k++) {
      FPut.Space(wr);
      fprintf(wr, "%d"c.node[k]);
    }
    fputc('\n', wr);

    /* Material properties: */
    FPut.Space(wr, 2);
    FPut.Double(wr, c.density); FPut.Space(wr, 2);
    FPut.Double(wr, c.alpha);   FPut.Space(wr);
    FPut.Double(wr, c.beta);    FPut.Space(wr, 2);
    FPut.Double(wr, c.eta1);    FPut.Space(wr);
    FPut.Double(wr, c.eta2);    
    fputc('\n', wr);
;
  } /* WriteCell */;

void RenumberNodes(T s, nat nNodes, READONLY NodeNums map)
  {
    /* Renumber nodes in each cell: */
    for (cn = 0;  cn < s.nCells; cn++) {
      { /* with*/ c == s.cell[cn] ) { 
        for (k = 0;  k <= 3;  k++) {c.node[k] = map[c.node[k]]; }
        SortCellCorners(c.node, c.restShape);
      };
    }

    /* Recompute the node table: */
    s.nNodes = nNodes;
    BuildNodeTable(s);
  } /* RenumberNodes */;

void SortCellCorners(VAR node: ARRAY [0..3] OF NodeNum; VAR r3x3_t A)
  /* 
    Rearranges the cell corners "node[0..3]" so that "node[0]" is
    the lowest and "node[1]" is second-lowest, preserving their screw.
  */

    void CycleNodes123()
    NodeNum *?tNode; r3_t tA,
    /* Cyclically shifts node[1..3] to lower indices. */
    {
      tNode = node[1]; node[1] = node[2]; node[2] = node[3]; node[3] = tNode;
      tA = A[0]; A[0] = A[1]; A[1] = A[2]; A[2] = tA;
    } /* CycleNodes123 */;

    void SwapNodes01And23()
    NodeNum *?tNode; r3_t tA,
    /*
      Swaps node[1] with node[2] and node[3] with node[4],
      shifting the rest tetrahedron so that node[0]
      remains at the origin. */
    {
      tNode = node[0]; node[1] = node[0]; node[0] = tNode;
      tNode = node[2]; node[2] = node[3]; node[3] = tNode;
      tA = A[1];
      A[1] = r3_Sub(A[2], A[0]);
      A[2] = r3_Sub(tA, A[0]);
      A[0] = r3_Neg(A[0]);
    } /* SwapNodes01And23 */;

  {
    /* Find the lowest-numbered node, make it node[0]: */
    { /* with*/ m == min(min(node[0], node[1]), min(node[2], node[3])) ) {
      if ((m == node[1] )) { SwapNodes01And23()
      } else if ((m == node[2] )) { CycleNodes123(); SwapNodes01And23()
      } else if ((m == node[3] )) { CycleNodes123(); CycleNodes123(); SwapNodes01And23()
      } else {  affirm(m == node[0] , "??"); ;
      };
    }
    affirm(node[0] < node[1] , "??");
    affirm(node[0] < node[2] , "??");
    affirm(node[0] < node[3] , "??");
    /* Now make node[1] the second-smallest: */
    while (node[1] > node[2]) || (node[1] > node[3] ) {
      CycleNodes123(); ;
    }
    affirm(node[1] < node[2] , "??");
    affirm(node[1] < node[3] , "??");
  } /* SortCellCorners */;

void BuildNodeTable(T s)
  /*
    Builds the "s.node" table from "s.cell[*]".
    Assumes "s.nNodes" is already defined.
  */
  {
    s.node = Node_vec_new(s.nNodes);
    { /* with*/ node == s.node^ ) {
      for (nn = 0;  nn < s.nNodes; nn++) {node[nn].nCells = 0; }
      for (cn = 0;  cn < s.nCells; cn++) {
        for (k = 0;  k <= 3;  k++) {
          { /* with*/ nn == s.cell[cn].node[k], n == node[nn] ) { n.nCells++;
        } ;
      }
      for (nn = 0;  nn < s.nNodes; nn++) {
        { /* with*/ n == node[nn] ) {
          n.cell = CellNum_vec_new(n.nCells);
          n.nCells = 0;
        };
      }
      for (cn = 0;  cn < s.nCells; cn++) {
        for (k = 0;  k <= 3;  k++) {
          { /* with*/ nn == s.cell[cn].node[k], n == node[nn] ) {
            n.cell[n.nCells] = cn;
            n.nCells++;
          };
        };
      };
    };
  } /* BuildNodeTable */;

T Union(sa, T sb)
  {
    { /* with*/ 
      totCells == sa.nCells + sb.nCells,
      totNodes == sa.nNodes + sb.nNodes,

      sr == NEW(T, 
        nCells = totCells, 
        cell = Cell_vec_new(totCells),
        nNodes = totNodes, 
        node = Node_vec_new(totNodes)
      ),

      nodeShift == sa.nNodes,
      cellShift == sa.nCells
    ) {
      /* Combine cell tables: */
      { /* with*/ 
        cellA == sa.cell^,
        cellB == sb.cell^,
        cellR == sr.cell^
      ) {
        SUBARRAY(cellR, 0, sa.nCells) = cellA;
        for (cnB = 0;  cnB < sb.nCells; cnB++) {
          { /* with*/ cB == cellB[cnB], cR == cellR[cnB + cellShift] ) {
            cR = cB;
            for (k = 0;  k <= 3;  k++) {cR.node[k] = cR.node[k] + nodeShift;
          };
        };
      }

      /* Combine node tables: */
      { /* with*/ 
        nodeA == sa.node^,
        nodeB == sb.node^,
        nodeR == sr.node^
      ) {
        SUBARRAY(nodeR, 0, sa.nNodes) = nodeA;
        for (nnB = 0;  nnB < sb.nNodes; nnB++) {
          { /* with*/ nB == nodeB[nnB], nR == nodeR[nnB + nodeShift] ) {
            nR.nCells = nB.nCells;
            nR.cell = NEW(REF CellNums, nR.nCells);
            for (k = 0;  k < nR.nCells; k++) {nR.cell[k] = nR.cell[k] + cellShift;
          };
        };
      }
      return sr;
    };
  } /* Union */;

{; } /* SystemTopology */.


