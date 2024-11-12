MODULE SystemTopology;

IMPORT Rd, Wr, Thread, Fmt, Lex, LR3, LR3x3;
IMPORT Util, FGet, NGet, FPut, NPut, FileFmt;

REVEAL
  T = Public BRANDED OBJECT
    OVERRIDES
      cellMass := CellMass;
      renumberNodes := RenumberNodes;
    END;
    
PROCEDURE CellMass(s: T; cn: CellNum): LONGREAL =
  BEGIN 
    RETURN s.cell[cn].density * s.cell[cn].restDet/6.0D0
  END CellMass;

CONST 
  TopologyFileVersion =     "96-12-04";
  RLWLTopologyFileVersion = "95-06-03";

PROCEDURE Read(rd: Rd.T; RLWLFormat: BOOLEAN := FALSE): T =
  VAR nNodes: CARDINAL;
  BEGIN
    WITH s = NEW(T) DO
      (* Read "s.nCells" from file: *)
      IF RLWLFormat THEN
        FileFmt.ReadHeader(rd, "topology", RLWLTopologyFileVersion);
        s.comment := "";
        s.nCells := NGet.Int(rd, "tetrahedrons");
        s.nNodes := LAST(CARDINAL)
      ELSE
        FileFmt.ReadHeader(rd, "topology", TopologyFileVersion);
        s.comment := FileFmt.ReadComment(rd, '|');
        s.nCells := NGet.Int(rd, "cells");
        s.nNodes := NGet.Int(rd, "nodes")
      END;

      (* Read "s.cell[*]", compute effective "nNodes": *)
      nNodes := 0;
      s.cell := NEW(REF ARRAY OF Cell, s.nCells);
      FOR cn := 0 TO s.nCells-1 DO
        WITH c = s.cell[cn] DO
          ReadCell(rd, cn, c);
          FOR k := 0 TO 3 DO 
            nNodes := MAX(nNodes, c.node[k] + 1)
          END
        END
      END;

      (* Check/define "s.nNodes": *)
      IF RLWLFormat THEN
        s.nNodes := nNodes
      ELSE
        <* ASSERT nNodes <= s.nNodes *>
      END;

      FileFmt.ReadFooter(rd, "topology");

      BuildNodeTable(s);
      RETURN s
    END;
  END Read;

PROCEDURE ReadCell(rd: Rd.T; cn: CARDINAL; VAR c: Cell) =
  (* Reads data for cell number "cn" *)
  <* FATAL Rd.Failure, Thread.Alerted *>
  BEGIN
    (* Read and check cell num: *)
    Lex.Skip(rd);
    WITH cnr = FGet.Int(rd) DO <* ASSERT cnr = cn *> END; Lex.Skip(rd);
    FGet.Colon(rd);


    (* Read rest shape: *)
    FOR i := 0 TO 2 DO
      FOR j := 0 TO 2 DO 
        c.restShape[i,j] := FGet.LongReal(rd);
      END
    END;

    (* Read node nums, sort them: *)
    Lex.Skip(rd);
    c.node[0] := FGet.Int(rd);
    c.node[1] := FGet.Int(rd);
    c.node[2] := FGet.Int(rd);
    c.node[3] := FGet.Int(rd);
    SortCellCorners(c.node, c.restShape);

    (* Read material properties: *)
    Lex.Skip(rd);
    c.density := FGet.LongReal(rd);
    c.alpha   := FGet.LongReal(rd);
    c.beta    := FGet.LongReal(rd);
    c.eta1    := FGet.LongReal(rd);
    c.eta2    := FGet.LongReal(rd);

    (* Compute rest volume and cell mass: *)
    c.restVol  := LR3x3.Det(c.restShape)/6.0d0;
    IF c.restVol < 1.0D-50 THEN 
      Util.Error("cell " & Fmt.Int(cn) & " has zero or negative volume")
    END;
    c.mass := c.restVol * c.density;
    c.restInv := LR3x3.Inv(c.restShape);

  END ReadCell;

PROCEDURE Write(wr: Wr.T; s: T) =
  <* FATAL Wr.Failure, Thread.Alerted *>
  BEGIN
    FileFmt.WriteHeader(wr, "topology", TopologyFileVersion);
    FileFmt.WriteComment(wr, s.comment, '|');
    NPut.Int(wr, "cells", s.nCells);
    NPut.Int(wr, "nodes", s.nNodes);
    FOR cn := 0 TO s.nCells-1 DO
      WriteCell(wr, cn, s.cell[cn])
    END;
    FileFmt.WriteFooter(wr, "topology");
    Wr.Flush(wr)
  END Write;

PROCEDURE WriteCell(wr: Wr.T; cn: CellNum; READONLY c: Cell) =
  BEGIN
    (* Cell number: *)
    FPut.Int(wr, cn); FPut.Colon(wr); 
    FPut.EOL(wr);

    (* Rest shape: *)
    FOR i := 0 TO 2 DO 
      FPut.Space(wr);
      FOR j := 0 TO 2 DO
        FPut.Space(wr);
        FPut.LongReal(wr, c.restShape[i,j]);
      END;
    END;
    FPut.EOL(wr);

    (* Node numbers: *)
    FPut.Space(wr);
    FOR k := 0 TO 3 DO
      FPut.Space(wr);
      FPut.Int(wr, c.node[k])
    END;
    FPut.EOL(wr);

    (* Material properties: *)
    FPut.Space(wr, 2);
    FPut.LongReal(wr, c.density); FPut.Space(wr, 2);
    FPut.LongReal(wr, c.alpha);   FPut.Space(wr);
    FPut.LongReal(wr, c.beta);    FPut.Space(wr, 2);
    FPut.LongReal(wr, c.eta1);    FPut.Space(wr);
    FPut.LongReal(wr, c.eta2);    
    FPut.EOL(wr);

  END WriteCell;

PROCEDURE RenumberNodes(s: T; nNodes: CARDINAL; READONLY map: NodeNums) =
  BEGIN
    (* Renumber nodes in each cell: *)
    FOR cn := 0 TO s.nCells-1 DO 
      WITH c = s.cell[cn] DO 
        FOR k := 0 TO 3 DO c.node[k] := map[c.node[k]] END;
        SortCellCorners(c.node, c.restShape)
      END
    END;

    (* Recompute the node table: *)
    s.nNodes := nNodes;
    BuildNodeTable(s)
  END RenumberNodes;

PROCEDURE SortCellCorners(VAR node: ARRAY [0..3] OF NodeNum; VAR A: LR3x3.T) =
  (* 
    Rearranges the cell corners "node[0..3]" so that "node[0]" is
    the lowest and "node[1]" is second-lowest, preserving their screw.
  *)

    PROCEDURE CycleNodes123() =
    VAR tNode: NodeNum; tA: LR3.T;
    (* Cyclically shifts node[1..3] to lower indices. *)
    BEGIN
      tNode := node[1]; node[1] := node[2]; node[2] := node[3]; node[3] := tNode;
      tA := A[0]; A[0] := A[1]; A[1] := A[2]; A[2] := tA;
    END CycleNodes123;

    PROCEDURE SwapNodes01And23() =
    VAR tNode: NodeNum; tA: LR3.T;
    (*
      Swaps node[1] with node[2] and node[3] with node[4],
      shifting the rest tetrahedron so that node[0]
      remains at the origin. *)
    BEGIN
      tNode := node[0]; node[1] := node[0]; node[0] := tNode;
      tNode := node[2]; node[2] := node[3]; node[3] := tNode;
      tA := A[1];
      A[1] := LR3.Sub(A[2], A[0]);
      A[2] := LR3.Sub(tA, A[0]);
      A[0] := LR3.Neg(A[0]);
    END SwapNodes01And23;

  BEGIN
    (* Find the lowest-numbered node, make it node[0]: *)
    WITH m = MIN(MIN(node[0], node[1]), MIN(node[2], node[3])) DO
      IF    m = node[1] THEN SwapNodes01And23()
      ELSIF m = node[2] THEN CycleNodes123(); SwapNodes01And23()
      ELSIF m = node[3] THEN CycleNodes123(); CycleNodes123(); SwapNodes01And23()
      ELSE  <* ASSERT m = node[0] *> 
      END
    END;
    <* ASSERT node[0] < node[1] *>
    <* ASSERT node[0] < node[2] *>
    <* ASSERT node[0] < node[3] *>
    (* Now make node[1] the second-smallest: *)
    WHILE node[1] > node[2] OR node[1] > node[3] DO
      CycleNodes123(); 
    END;
    <* ASSERT node[1] < node[2] *>
    <* ASSERT node[1] < node[3] *>
  END SortCellCorners;

PROCEDURE BuildNodeTable(s: T) =
  (*
    Builds the "s.node" table from "s.cell[*]".
    Assumes "s.nNodes" is already defined.
  *)
  BEGIN
    s.node := NEW(REF ARRAY OF Node, s.nNodes);
    WITH node = s.node^ DO
      FOR nn := 0 TO s.nNodes-1 DO node[nn].nCells := 0 END;
      FOR cn := 0 TO s.nCells-1 DO
        FOR k := 0 TO 3 DO 
          WITH nn = s.cell[cn].node[k], n = node[nn] DO INC(n.nCells) END
        END; 
      END;
      FOR nn := 0 TO s.nNodes-1 DO 
        WITH n = node[nn] DO
          n.cell := NEW(REF ARRAY OF CellNum, n.nCells);
          n.nCells := 0
        END
      END;
      FOR cn := 0 TO s.nCells-1 DO
        FOR k := 0 TO 3 DO
          WITH nn = s.cell[cn].node[k], n = node[nn] DO
            n.cell[n.nCells] := cn;
            INC(n.nCells);
          END
        END
      END
    END
  END BuildNodeTable;

PROCEDURE Union(sa, sb: T): T =
  BEGIN
    WITH
      totCells = sa.nCells + sb.nCells,
      totNodes = sa.nNodes + sb.nNodes,

      sr = NEW(T, 
        nCells := totCells, 
        cell := NEW(REF ARRAY OF Cell, totCells),
        nNodes := totNodes, 
        node := NEW(REF ARRAY OF Node, totNodes)
      ),

      nodeShift = sa.nNodes,
      cellShift = sa.nCells
    DO
      (* Combine cell tables: *)
      WITH 
        cellA = sa.cell^,
        cellB = sb.cell^,
        cellR = sr.cell^
      DO
        SUBARRAY(cellR, 0, sa.nCells) := cellA;
        FOR cnB := 0 TO sb.nCells-1 DO
          WITH cB = cellB[cnB], cR = cellR[cnB + cellShift] DO
            cR := cB;
            FOR k := 0 TO 3 DO cR.node[k] := cR.node[k] + nodeShift END
          END
        END
      END;

      (* Combine node tables: *)
      WITH 
        nodeA = sa.node^,
        nodeB = sb.node^,
        nodeR = sr.node^
      DO
        SUBARRAY(nodeR, 0, sa.nNodes) := nodeA;
        FOR nnB := 0 TO sb.nNodes-1 DO
          WITH nB = nodeB[nnB], nR = nodeR[nnB + nodeShift] DO
            nR.nCells := nB.nCells;
            nR.cell := NEW(REF CellNums, nR.nCells);
            FOR k := 0 TO nR.nCells-1 DO nR.cell[k] := nR.cell[k] + cellShift END
          END
        END
      END;
      RETURN sr
    END
  END Union;

BEGIN END SystemTopology.



