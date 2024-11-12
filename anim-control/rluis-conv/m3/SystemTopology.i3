INTERFACE SystemTopology;

(*
  A "SystemTopology.T" object describes the topology and 
  material properties of a collection of elastic bodies.

  The bodies are modeled as a finite element mesh with
  tetrahedral cells. 
  
  Nomenclature:

    cell = a tetrahedral element of the mesh.
    node = a vertex of the finite element mesh.

    isolated node = a vertex that belongs to no cell.

    screw = relative orientation of four points of $\R^3$.

    positive screw = the four points define a right-handed helical path.
      
*)

IMPORT LR3x3, Rd, Wr;

TYPE
  NodeNum = CARDINAL;  NodeNums = ARRAY OF NodeNum;
  CellNum = CARDINAL;  CellNums = ARRAY OF CellNum;

CONST 
  NoNode: NodeNum = LAST(CARDINAL);
  NoCell: CellNum = LAST(CARDINAL);

TYPE
  CellNodes = ARRAY [0..3] OF NodeNum;  (* The corners of a cell *)

  Cell = RECORD
      node: CellNodes;       (* Corner nodes, semi-sorted, positive screw. *)
      density: LONGREAL;     (* Density in the rest state. *)
      alpha, beta: LONGREAL; (* Elastic moduli. *)
      eta1, eta2: LONGREAL;  (* Viscosity coefficients. *)
      restShape: LR3x3.T;    (* Maps canonical tetrahedron to rest shape. *)
      (* Computed quantities: *)
      mass: LONGREAL;        (* Cell mass = (rest density) * (rest volume). *)
      restVol: LONGREAL;     (* Rest volume. *)
      restInv: LR3x3.T;      (* Matrix inverse of "restShape" *)
    END;
    
  NodeCells = ARRAY OF CellNum;         (* The cells incident to a node *)

  Node = RECORD 
      nCells: CARDINAL;             (* Number of cells incident to this node *)
      cell: REF NodeCells;          (* Cells incident to this node *)
    END;

  T <: Public;
  Public = OBJECT
      comment: TEXT;            (* Comment text *)
      
      nCells: CARDINAL;         (* Number of cells. *)
      cell: REF ARRAY OF Cell;  (* The cells. *)

      nNodes: CARDINAL;         (* Number of distinct nodes. *)
      node: REF ARRAY OF Node;  (* The nodes *)

    METHODS

      cellMass(cell: CellNum): LONGREAL;
        (*
          Computes the mass of tetrahedron number "i". *)
          
      renumberNodes(nNodes: CARDINAL; READONLY map: NodeNums);
        (*
          Renumbers the nodes in "s".
          
          After the call, the mesh "s" will have exactly "nNodes"
          nodes.  Each original node number "nn"
          is replaced by node number "map[nn]" (which must lie in "[0..nNodes-1]").
          
          The "map" need not be injective nor surjective.  If a node
          number in "[0..nNodes-1]" does not appear in the "map", 
          a new isolated node will be created with that number.
          
          On the other hand, if "map[pn] = map[qn]", then nodes "pn"
          and "qn" will become identified.  This is OK, provided the
          result is still a valid tetrahedral mesh. *)

    END;

PROCEDURE Read(rd: Rd.T; RLWLFormat: BOOLEAN := FALSE): T;
  (*
    Parses the system's topology and material properties from "rd". 
    Reorders the vertices in each cell. Computes rest volumes and
    rest shape tensors. Builds the cell-node incidence tables.

    If "RLWLFormat = TRUE", assumes the input is in 
    Rogério's original ".tp" format (version "95-06-03"),
    with "tetrahedrons = " instead of "cells = ". *)

PROCEDURE Write(wr: Wr.T; s: T);
  (*
    Writes the system's topology and material properties to "wr",
    a format compatible with "Read". *)

PROCEDURE Union(sa, sb: T): T;
(*
  Makes a single model that is the disjoint union of "sa" and "sb".
  The nodes of "sb" are renumbered by adding "sa.nNodes" to their original
  numbers; and analogously for the cells. *)

END SystemTopology.
