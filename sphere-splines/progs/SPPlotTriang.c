/* SPPlotTriang.c -- Plots a  triangulation  on the sphere. */
/* Last edited on 2006-03-19 12:51:21 by stolfi */

#include <SPTriang.h>
#include <SPPlot.h>
#include <SPOptions.h>
#include <SPPlotOptions.h>
#include <SPBasic.h>
#include <SPH3.h>
#include <r3.h>
#include <pswr.h>
#include <bool.h>
#include <stdio.h>
#include <math.h>

typedef struct Options /* Parsed command line options */
  { char *triName;        /* Triang file (minus ".tri" extension) or "-" */
    bool_t spherical;     /* TRUE draws spherical, FALSE draws polyhedron. */
    bool_t showBack;      /* TRUE to also plot back side. */
    char *outName;        /* Output file (minus ".ps"/".eps" extension). */
    SPPlotOptions_t plt;  /* Plot format and style options. */
  } Options;

Options GetOptions(int argn, char **argc);
Triangulation *ReadTriangulation(char *triName);

int main(int argn, char **argc)
  { Options o = GetOptions(argn, argc);
    SPPlotOptions_t *po = &(o.plt);
    Triangulation *tri = ReadTriangulation(o.triName);
    SPH3_PMap *map;
    
    /* Fix perspective parameters and compute view map: */
    SPPlot_FixView(&(po->obs), 0.0, &(po->upp), &(po->radius), 0.0, NULL, TRUE);
    map = SPPlot_PerspMap(&(po->obs), po->radius, &(po->upp));
      
    /* Open the figure stream: */
    SPPlot_Stream *fps = SPPlot_NewStream
      (po->eps, o.outName, po->paperSize, po->figSize, po->projLonLat, po->caption.ne);

    /* Set scales and start new figure: */
    double xMax = 1.0;
    double yMax = (po->projLonLat ? 0.5 : 1.0);
    pswr_new_picture(fps, -xMax,+xMax, -yMax,+yMax);
    pswr_new_picture(fps, -1.0,1.0, -1.0,1.0);

    /* Plot the triangulation: */
    { int side;
      double maxStep = (po->radius) * po->meshSize/(po->figSize/2);
      for (side = -1; side <= +1; side += 2)
        { if ((side > 0) || (o.showBack))
            { /* Choose color depending on side: */
              if (side < 0)
                { pswr_set_pen(fps, 0.0,0.7,1.0,  0.1, 0.0,0.0); }
              else
                { pswr_set_pen(fps, 0.0,0.0,0.0,  0.1, 0.0,0.0); }
              if (o.spherical)
                { SPPlot_Triangulation(fps, map, tri, maxStep, (side < 0)); }
              else
                { SPPlot_Polyhedron(fps, map, tri, (side < 0)); }
            }
        }
      if (o.spherical)
        { /* Plot sphere outline: */
          pswr_set_pen(fps, 0.0,0.0,0.0, 0.15, 0.0,0.0);
          SPPlot_Sphere(fps, map);
        }
    }
    
    /* Add caption, if it is the case: */
    if(! po->eps) { pswr_add_caption(fps, o.outName, 0.0); }
    int k; 
    for (k = 0; k < po->caption.ne; k++)
      { pswr_add_caption(fps, po->caption.e[k], 0.0); }
    pswr_close_stream(fps);
    return 0;
  }
  
Triangulation *ReadTriangulation(char *triName)
  { FILE *rd = open_read(txtcat(triName, ".tri"), TRUE);
    Triangulation *tri = SPTriang_Read(rd, 1);
    fclose(rd);
    return tri;
  }
    

Options GetOptions(int argn, char **argc)
  { Options o;
    SPOptions_Parser_t *pp = SPOptions_NewParser(stderr, argn, argc);
    SPOptions_SetUsage(pp, 
      "SPPlotTriang \\\n"
      "  -triName <name>  \\\n"
      "  -outName <name> \\\n"
      "  [ -spherical ] [ -showBack ] \\\n"
      SPPlotOptions_TriangHelp " \n"
    );

    SPOptions_GetKeyword(pp, "-triName");                           
    o.triName = SPOptions_GetNext(pp);  

    o.spherical = SPOptions_TestKeyword(pp, "-spherical");

    o.showBack = SPOptions_TestKeyword(pp, "-showBack");

    SPOptions_GetKeyword(pp, "-outName");                               
    o.outName = SPOptions_GetNext(pp);  
    
    o.plt = SPPlotOptions_TriangParse(pp);

    SPOptions_Finish(pp);
    return o;
  }

