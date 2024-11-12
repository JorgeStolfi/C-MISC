INTERFACE SystemBoundary;

(*
  A representation of the the boundary of a "SystemTopology.T" 
  model; that is, the sets of exposed faces, edges, and nodes,
  and their incidence properties.  Also procedures to obtain such 
  a representation.
  
  Nomenclature:
  
    face = a face of the mesh, belonging to one or two cells.
    edge = an edge of the mesh, shared by one or more cells.
    site = a node of the mesh (with different numbering).

    exposed face = a face that is shared by only one cell.
    exposed edge = an edge that belongs to some (two or more) exposed faces.
    exposed site = a site that belongs to some (three or more) exposed faces.

*)    
  
IMPORT SystemTopology AS ST;

TYPE
  FaceNum = CARDINAL;  EdgeNums = ARRAY OF EdgeNum;
  EdgeNum = CARDINAL;  FaceNums = ARRAY OF FaceNum;
  SiteNum = CARDINAL;  SiteNums = ARRAY OF SiteNum;
    (*
      NOTE: the number of a site in the boundary structure (its "SiteNum")
      is not the same as its number in the tetrahedral mesh (its "NodeNum"). *) 

CONST 
  NoFace: FaceNum = LAST(CARDINAL);
  NoEdge: EdgeNum = LAST(CARDINAL);
  NoSite: SiteNum = LAST(CARDINAL);

TYPE
  Face = RECORD  (* An exposed face *)
      site: ARRAY [0..2] OF SiteNum;    (* Corners, semi-sorted, ccw from outside. *)
      cell: ST.CellNum;                 (* The adjacent cell. *)
      edge: ARRAY [0..2] OF EdgeNum;    (* Side opposite to each corner. *)
    END;

  Edge = RECORD
      site: ARRAY [0..1] OF SiteNum;    (* Endpoints, sorted. *)
      face: FaceNum;                    (* Some incident face. *)
      cell: ST.CellNum;                 (* Some incident cell. *)
      dual: ARRAY [0..1] OF ST.NodeNum; (* The other two nodes of that cell. *)
      (* Note: (site[0], site[1], dual[0], dual[1]) is a positive screw. *)
    END;
    
  Site = RECORD
      node: ST.NodeNum;                 (* Node number of this site in the 3d mesh. *)
    END;
    
  T = REF RECORD

      nEdges: CARDINAL;  (* Number of exposed edges. *)
      nFaces: CARDINAL;  (* Number of exposed faces. *)
      nSites: CARDINAL;  (* Number of exposed sites. *)
  
      face: REF ARRAY OF Face;    (* The exposed faces. *)
      edge: REF ARRAY OF Edge;    (* The exposed edges. *)
      site: REF ARRAY OF Site;    (* The exposed sites. *)
      
    END;
      
PROCEDURE ComputeBoundary(t: ST.T): T;
  (*
    Identifies the boundary faces, edges and vertices of "s". *) 
    
END SystemBoundary.
