/* SPMakeTriang -- Create spherical triangulations. */
/* Last edited on 2006-03-19 12:49:32 by stolfi */

#include <SPTriang.h>
#include <SPDelaunay.h>
#include <SPOptions.h>
#include <r3.h>
#include <SPBasic.h>
#include <bool.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* Creates triangulations for testing purposes.  The user may select
  between a few fixed geometries (a regular tetrahedron, a regular 
  icosahedron, and a deformed octahedron), or the Delaunay triangulation 
  of N random sites on the sphere. */

typedef struct Options     /* Parsed command line options */
  { char *outName;
    bool_t skewed;           /* TRUE uses a skewed distribution of points. */
    /* These are mutually exclusive: */
    int nSites;            /* Number of (random) sites */
    bool_t tetrahedron; 
    bool_t octahedron;
    bool_t icosahedron;
  } Options;

Options GetOptions(int argn, char **argc);

Triangulation *MakeTriangulation
  ( bool_t skewed,
    int nSites, 
    bool_t tetrahedron,
    bool_t octahedron,
    bool_t icosahedron
  );

SiteRef_vec_t MakeRandomSites(int n);
  /* Creates an array of random sites. */

int main(int argn, char **argc)
  {
    Options o = GetOptions(argn, argc);
    Triangulation *tri = MakeTriangulation(
        o.skewed, 
        o.nSites, 
        o.tetrahedron, 
        o.octahedron, 
        o.icosahedron
      );
    FILE *wr = open_write(txtcat(o.outName, ".tri"), TRUE);
    SPTriang_Write(wr, tri);
    fclose(wr);
    return 0;
  }
  
Triangulation *MakeTriangulation
  ( bool_t skewed,
    int nSites, 
    bool_t tetrahedron,
    bool_t octahedron,
    bool_t icosahedron
  )
  { Triangulation *tri;
    if (tetrahedron)
      { tri = SPDelaunay_RegularTetrahedron(1); }
    else if (octahedron)
      { tri = SPTriang_RegularOctahedron(1); }
    else if (icosahedron)
      { tri = SPDelaunay_RegularIcosahedron(1); }
    else
      { SiteRef_vec_t site = MakeRandomSites(nSites);
        Arc e = SPDelaunay_BuildInc(site);
        tri = SPTriang_FromTopology(e, 1);
      }
    
    if (skewed)
      { /* Modify sites by a conformal map that leaves {S^2} invariant: */
        int i; 
        double theta = PI/3.0;
        double c = cos(theta), s = sin(theta);
        for (i = 0; i < tri->out.ne; i++)
          { r3_t *p = &(Org(tri->out.e[i])->pos);
            double w = 1.0 + c*p->c[0];
            double x = c + p->c[0]; 
            double y = s*p->c[1];
            double z = s*p->c[2];
            *p = (r3_t){{ x/w, y/w, z/w }};
            /* Just to be sure... */
            r3_dir(p, p);
          }
      }
    SPTriang_ComputeGeometryDataOfFaces(tri, 1);
    return tri;
  }

SiteRef_vec_t MakeRandomSites(int nSites)
  {
    auto S2Point ThrowSite(int i);
    
    S2Point ThrowSite(int i)
      { /* Generate a random point in the unit ball: */
        r3_t p;
        do {
          p.c[0] = 2.0*drandom() - 1.0;
          p.c[1] = 2.0*drandom() - 1.0;
          p.c[2] = 2.0*drandom() - 1.0;
        } while (r3_norm_sqr(&p) > 1.0);
        r3_dir(&p, &p);
        return p;
      }
    
    srand(352563);
    srandom(352563);
    return SPTriang_MakeSites(nSites, ThrowSite);
  }


Options GetOptions(int argn, char **argc)
  {
    Options o;
    SPOptions_Parser_t *pp = SPOptions_NewParser(stderr, argn, argc);

    SPOptions_SetUsage(pp, 
      "SPMakeTriang \\\n"
      "  -outName NAME \\\n"
      "  [ -skewed ] \\\n"
      "  [ -tetrahedron | -icosahedron | -octahedron | -nSites NUM ]"
    );

    SPOptions_GetKeyword(pp, "-outName");                               
    o.outName = SPOptions_GetNext(pp);  

    o.skewed = SPOptions_TestKeyword(pp, "-skewed");

    if (SPOptions_TestKeyword(pp, "-tetrahedron"))
      { o.nSites = 4;
        o.tetrahedron = TRUE;
        o.octahedron = FALSE;
        o.icosahedron = FALSE;
      }
    else if (SPOptions_TestKeyword(pp, "-octahedron"))
      { o.nSites = 6;
        o.tetrahedron = FALSE;
        o.octahedron = TRUE;
        o.icosahedron = FALSE;
      }
    else if (SPOptions_TestKeyword(pp, "-icosahedron"))
      { o.nSites = 12;
        o.tetrahedron = FALSE;
        o.octahedron = FALSE;
        o.icosahedron = TRUE;
      }
    else
      { SPOptions_GetKeyword(pp, "-nSites");
        o.nSites = SPOptions_GetNextInt(pp, 4, 100000);
        o.tetrahedron = FALSE;
        o.octahedron = FALSE;
        o.icosahedron = FALSE;
      }
    SPOptions_Finish(pp);
    return o;
  }
