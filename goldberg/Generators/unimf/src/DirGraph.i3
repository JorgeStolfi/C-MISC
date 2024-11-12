INTERFACE DirGraph;

(* Simple data structure for directed graphs *)
(* Created April 1997 by J. Stolfi  *)
(* See the copyright notice at the end of this file. *)

IMPORT Wr, XRandom, ASCII;

TYPE
    
  Count = CARDINAL;
  Counts = ARRAY OF Count;

  VertexNum = CARDINAL;  
  VertexNums = ARRAY OF VertexNum;

  Edge = ARRAY [0..1] OF CARDINAL;
  Edges = ARRAY OF Edge;
  
  EdgeNum = CARDINAL;
  EdgeNums = ARRAY OF EdgeNum;
  
  Cost = CARDINAL;
  Costs = ARRAY OF Cost;
  
  Mark = ASCII.Range;
  Marks = ARRAY OF Mark;

  T = RECORD
      NV: Count;               (* The vertices are the integers "[0..NV-1]". *)
      NE: Count;               (* The edges are numbered "[0..NE-1]" *)
      e: REF Edges;            (* The edges, in any order *)
      (* Optional parameters for max-flow problems: *)
      s: VertexNum;            (* `Source' vertex, whathever that means. *)
      t: VertexNum;            (* `Sink' vertex, whatever that means. *)
      m: REF Marks;            (* The marks on the edges, whatever that means. *)
      c: REF Costs;            (* Edge costs/caps/whatever, in the same order *)
      class: TEXT;             (* Graph class (for the "p" output line) *)
      cmt: TEXT;               (* Comment ("c" lines). *)
    END;
    
CONST
  NoMark = SET OF Mark{};

PROCEDURE Write(
    wr: Wr.T;
    READONLY g: T;
    vBase: INTEGER;
    READONLY omit: SET OF Mark := NoMark;
  );
  (*
    Writes the graph "g" to "wr". The format is
    
       c COMMENT
       c COMMENT
       ...
       p CLASS NUMVERTICES NUMEDGES
       s SOURCE
       t SINK
       a EDGE[0][0] EDGE[0][1] COST[0]
       a EDGE[1][0] EDGE[1][1] COST[1]
       a EDGE[2][0] EDGE[2][1] COST[2]
       ...
       
    All vertex numbers EDGE[i][j], SOURCE, SINK are shifted by "vBase".
    
    The COST field in the "a" lines is printed only if the
    "g.c" array is non-NIL.
  *)

PROCEDURE PermuteVertexNums(VAR g: T; rnd: XRandom.T);
  (*
    Renumbers the vertices of "g" in random order. *)

PROCEDURE ComplementVertexNums(VAR g: T);
  (*
    Renumbers all vertices in the opposite of the current order. *)

PROCEDURE ReverseEdges(VAR e: Edges);
  (*
    Reverses the direction of the edges in "e". *)

END DirGraph.
