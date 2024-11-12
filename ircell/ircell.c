/* Last edited on 2023-10-01 19:39:25 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <math.h>
#include <bool.h>
#include <pswr.h>

void ircell_draw_frame(PSStream *ps, double xCenter, double yCenter);
  /* Draws the template for the plastic frames that hold
    the cell walls together, centered at {(xCenter,yCenter)}
    in mm.  */
    
int main(int argc, char **argv)
  {
    char *prefix = "ircell_";
    FILE *wr = stdout;
    bool_t eps = TRUE;
    char *docName = NULL;
    char *paperSize = NULL;
    double xMin = -90; /* Low client {x} (mm). */
    double xMax = +90; /* High client {x} (mm). */
    double yMin = -60; /* Low client {y} (mm). */
    double yMax = +60; /* High client {y} (mm). */
    double mrg = 4.0;  /* Margin of EPS files. */
    double hCanvasSize = (xMax - xMin)*epswr_pt_per_mm + 2*mrg;
    double vCanvasSize = (yMax - yMin)*epswr_pt_per_mm + 2*mrg;
    PSStream *ps = pswr_new_stream
      ( prefix, wr, eps, docName, paperSize, FALSE, hCanvasSize, vCanvasSize );               
    pswr_new_picture(ps, xMin, xMax, yMin, yMax); 
    ircell_draw_frame(ps,-30.0,0.0);
    ircell_draw_frame(ps,+30.0,0.0);
    pswr_close_stream(ps);             
    return 0;
  }
  
void ircell_draw_frame(PSStream *ps, double xCenter, double yCenter)
  {
    /* All dimensions are in mm, */
    
    /* Circular outline: */
    double oRad = 25.0;  /* Outer radius. */
    double iRad =  3.0;  /* Center drill hole radius. */
    
    /* Slide dimensions: */
    /* double sLen = 75.25; */ /* Length. */
    double sWid = 24.85; /* Width. */
    double sThk = 1.00; /* Thickness. */
    
    /* Dimensions of enclosing triangle {T}: */
    double s30 = sin(M_PI/6.0);  /* Thirty degrees. */
    double c30 = cos(M_PI/6.0);  /* Thirty degrees. */
    double a = sThk/c30;  /* Big overhang. */
    double b = a*s30; /* Small overhang. */
    double tSide = sWid + a + 2*b; /* Side of {T}. */
    double tHeight = tSide*c30;  /* Height of {T}. */
    double tRad = 2*tHeight/3; /* Circumradius of {T}. */
    
    /* Coordinates of cell section: */
    double s120 = sin(2*M_PI/3.0);  /* Hundred and twenty degrees. */
    double c120 = cos(2*M_PI/3.0);  /* Hundred and twenty degrees. */
    int ns = 2; /* Num vertices in 1/3 of the section. */
    int nv = 3*ns; /* Total num of vertices in section. */
    double x[nv], y[nv]; /* Coords of section vertices. */
    
    x[0] = + s30*tSide - b;
    y[0] = - s30*tRad;
    
    x[1] = - s30*tSide + a;
    y[1] = - s30*tRad;
    
    for (int i = ns; i < nv; i++)
      { x[i] = + c120*x[i-ns] + s120*y[i-ns];
        y[i] = - s120*x[i-ns] + c120*y[i-ns];
      }
    for (int i = 0; i < nv; i++)
      { x[i] += xCenter; y[i] += yCenter; }
    
    /* Draw it: */
    pswr_set_pen(ps, 0.0,0.0,0.0, 0.05, 0.0,0.0);
    pswr_circle(ps, xCenter, yCenter, oRad, FALSE, TRUE);
    pswr_circle(ps, xCenter, yCenter, iRad, FALSE, TRUE);
    pswr_polygon(ps, TRUE, x, y, nv, FALSE, TRUE, TRUE);
  }
