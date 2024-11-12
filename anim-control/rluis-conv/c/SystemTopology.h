
#ifndef SystemTopology_H
#define SystemTopology_H



/*
  A "SystemTopology.T" object describes the topology and 
  material properties of a collection of elastic bodies.

  The bodies are modeled as a finite element mesh with
  tetrahedral cells. 
  
  Nomenclature:

    cell == a tetrahedral element of the mesh.
    node == a vertex of the finite element mesh.

    isolated node == a vertex that belongs to no cell.

    screw == relative orientation of four points of $\R^3$.

    positive screw == the four points define a right-handed helical path.
      
*/

#include <r3x3.h>
#include <Rd.h>
#include <Wr.h>

typedef
  NodeNum == nat;  NodeNums == ARRAY OF NodeNum;
  CellNum == nat;  CellNums == ARRAY OF CellNum;

CONST 
  NoNode: NodeNum == ((nat.nel - 1)|?|MAX_CARDINAL);
  NoCell: CellNum == ((nat.nel - 1)|?|MAX_CARDINAL);

typedef
  CellNodes == ARRAY [0..3] OF NodeNum;  /* The corners of a cell */

  Cell == struct ?? {
      CellNodes node;       /* Corner nodes, semi-sorted, positive screw. */
      double density;     /* Density in the rest state. */
      double alpha, beta; /* Elastic moduli. */
      double eta1, eta2;  /* Viscosity coefficients. */
      r3x3_t restShape;    /* Maps canonical tetrahedron to rest shape. */
      /* Computed quantities: */
      double mass;        /* Cell mass == (rest density) * (rest volume). */
      double restVol;     /* Rest volume. */
      r3x3_t restInv;      /* Matrix inverse of "restShape" */;
    }
    
  NodeCells == ARRAY OF CellNum;         /* The cells incident to a node */

  Node == struct ?? { 
      nat nCells;             /* Number of cells incident to this node */
      cell: REF NodeCells;          /* Cells incident to this node */;
    }

  T <: Public;
  Public == OBJECT
      char *comment;            /* Comment text */
      
      nat nCells;         /* Number of cells. */
      cell: Cell_vec;  /* The cells. */

      nat nNodes;         /* Number of distinct nodes. */
      node: Node_vec;  /* The nodes */

    METHODS

      cellMass(CellNum cell): double;
        /*
          Computes the mass of tetrahedron number "i". */
          
      renumberNodes(nat nNodes, READONLY NodeNums map);
        /*
          Renumbers the nodes in "s".
          
          After the call, the mesh "s" will have exactly "nNodes"
          nodes.  Each original node number "nn"
          is replaced by node number "map[nn]" (which must lie in "[0..nNodes-1]").
          
          The "map" need not be injective nor surjective.  If a node
          number in "[0..nNodes-1]" does not appear in the "map", 
          a new isolated node will be created with that number.
          
          On the other hand, if "map[pn] == map[qn]", then nodes "pn"
          and "qn" will become identified.  This is OK, provided the
          result is still a valid tetrahedral mesh. */
;
    }

T Read(FILE *rd, RLWLFormat: bool = FALSE);
  /*
    Parses the system's topology and material properties from "rd". 
    Reorders the vertices in each cell. Computes rest volumes and
    rest shape tensors. Builds the cell-node incidence tables.

    If "RLWLFormat == TRUE", assumes the input is in 
    Rogério's original ".tp" format (version "95-06-03"),
    with "tetrahedrons == " instead of "cells == ". */

void Write(FILE *wr, T s);
  /*
    Writes the system's topology and material properties to "wr",
    a format compatible with "Read". */

T Union(sa, T sb);
/*
  Makes a single model that is the disjoint union of "sa" and "sb".
  The nodes of "sb" are renumbered by adding "sa.nNodes" to their original
  numbers; and analogously for the cells. */
;
} /* SystemTopology */.

#endif
