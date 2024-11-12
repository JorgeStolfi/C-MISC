
#ifndef SystemBoundary_H
#define SystemBoundary_H



/*
  A representation of the the boundary of a "SystemTopology.T" 
  model; that is, the sets of exposed faces, edges, and nodes,
  and their incidence properties.  Also procedures to obtain such 
  a representation.
  
  Nomenclature:
  
    face == a face of the mesh, belonging to one or two cells.
    edge == an edge of the mesh, shared by one or more cells.
    site == a node of the mesh (with different numbering).

    exposed face == a face that is shared by only one cell.
    exposed edge == an edge that belongs to some (two or more) exposed faces.
    exposed site == a site that belongs to some (three or more) exposed faces.

*/    
  
#include <SystemTopology AS ST.h>

typedef
  FaceNum == nat;  EdgeNums == ARRAY OF EdgeNum;
  EdgeNum == nat;  FaceNums == ARRAY OF FaceNum;
  SiteNum == nat;  SiteNums == ARRAY OF SiteNum;
    /*
      NOTE: the number of a site in the boundary structure (its "SiteNum")
      is not the same as its number in the tetrahedral mesh (its "NodeNum"). */ 

CONST 
  NoFace: FaceNum == ((nat.nel - 1)|?|MAX_CARDINAL);
  NoEdge: EdgeNum == ((nat.nel - 1)|?|MAX_CARDINAL);
  NoSite: SiteNum == ((nat.nel - 1)|?|MAX_CARDINAL);

typedef
  Face == struct ?? {  /* An exposed face */
      site: ARRAY [0..2] OF SiteNum;    /* Corners, semi-sorted, ccw from outside. */
      ST.CellNum cell;                 /* The adjacent cell. */
      edge: ARRAY [0..2] OF EdgeNum;    /* Side opposite to each corner. */;
    }

  Edge == struct ?? {
      site: ARRAY [0..1] OF SiteNum;    /* Endpoints, sorted. */
      FaceNum face;                    /* Some incident face. */
      ST.CellNum cell;                 /* Some incident cell. */
      dual: ARRAY [0..1] OF ST.NodeNum; /* The other two nodes of that cell. */
      /* Note: (site[0], site[1], dual[0], dual[1]) is a positive screw. */;
    }
    
  Site == struct ?? {
      ST.NodeNum node;                 /* Node number of this site in the 3d mesh. */;
    }
    
  T == REF  struct ?? {

      nat nEdges;  /* Number of exposed edges. */
      nat nFaces;  /* Number of exposed faces. */
      nat nSites;  /* Number of exposed sites. */
  
      face: Face_vec;    /* The exposed faces. */
      edge: Edge_vec;    /* The exposed edges. */
      site: Site_vec;    /* The exposed sites. */
      ;
    }
      
T ComputeBoundary(ST.T t);
  /*
    Identifies the boundary faces, edges and vertices of "s". */ 
    ;
} /* SystemBoundary */.

#endif
