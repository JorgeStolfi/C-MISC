/* Routines for writing PS (Encapsulated PostScript) files. */
/* Last edited on 2002-03-15 10:24:44 by stolfi */

#ifndef SOeps_H
#define SOeps_H

#include <stdio.h>

/*** ENCAPSULATED POSTSCRIPT FIGURES ***/

void eps_begin_figure (
    FILE *psfile,
    double xmin, double xmax,
    double ymin, double ymax,

    double hmin, double hmax,
    double vmin, double vmax,

    int xn, int yn
  );
  /* 
    Initializes an Encapsulated PostScript file.

    An Encapsulated PostScript file contains a single figure, with a
    BoundingBox comment, without an explicit "showpage" command;
    It is suitable for inclusion in TeX papers and such.  

    This procedure writes the preamble and bounding box, sets
    coordinate system, clip path, caption font, Postscript operators
    and constants, etc.

    Client coordinates will range over [xmin __ xmax] x [ymin __ ymax].
    Plotting area is [hmin __ hmax] x [vmin __ vmax] (in pt).
    The plotting scales dh/dx and dv/dy must be equal.

    The plotting area is divided implicitly into a grid of /xn/ by /yn/
    rectangular "cells".
  */

void eps_end_figure (FILE *psfile);
  /* Finalizes an Encapsulated PostScript file. */

/*** MULTI-PAGE POSTSCRIPT DOCUMENTS ***/

void eps_begin_document (FILE *psfile);
  /* 
    Initializes a (non-encapsulated, multi-page) PS file.
    
    A standard PostScript file contains a document with one or more pages,
    with explicit "show" commands at the end of each page, and no BoundingBox 
    comment.  It can be printed by itself.

    The procedure writes the file's preamble, with auxiliary
    Postscript operators.

    The client must call ps_begin_page and ps_end_page around each page.
    of the document.
  */

void eps_begin_page(
    FILE *psfile,
    int page,
    double xmin, double xmax,
    double ymin, double ymax,

    double hmin, double hmax,
    double vmin, double vmax,

    int xn, int yn
  );
  /* 
    Starts a new page of a PS document. 
    
    Sets coordinate system, clip path, dimension constants, etc.

    Client coordinates will range over [xmin __ xmax] x [ymin __ ymax].
    Plotting area is [hmin __ hmax] x [vmin __ vmax] (in pt).
    The plotting scales dh/dx and dv/dy must be equal.

    The plotting area is divided implicitly into a grid of /xn/ by /yn/ 
    rectangular "cells".
  */

void eps_end_page(FILE *psfile);
  /* 
    Finalizes a page: Writes page trailer line, etc. */

void eps_end_document (FILE *psfile, int npages);
  /* 
    Finalizes a multipage document.  The client must keep track 
    of the number of pages that were written to the file. */

/*** DRAWING COMMANDS ***/

void eps_begin_section (FILE *psfile, char *title);
  /* 
    Starts a new section of a plot. The title is a comment */

void eps_end_section (FILE *psfile);
  /* 
    Ends a section of a plot. */

void eps_set_pen (
    FILE *psfile,
    double r, double g, double b,
    double width,
    double dashlength,
    double dashspace
  );
  /*
    Sets pen parameters and ink color for line drawing.
    Dimensions are in *millimeters* */

void eps_draw_segment(
    FILE *psfile,
    double xa, double ya,
    double xb, double yb
  );
  /* 
    Draws segment from (xa,ya) to (xb,yb) with current pen and color */

void eps_draw_curve(
    FILE *psfile,
    double xa, double ya,
    double xb, double yb,
    double xc, double yc,
    double xd, double yd
  );
  /*
    Draws a Bezier arc with given control points, using the current pen and color. */
  
void eps_draw_rectangle(
    FILE *psfile,
    double xlo, double xhi,
    double ylo, double yhi
  );
  /* 
    Draws the outline of the given rectangle using the current pen. */

void eps_fill_rectangle(
    FILE *psfile,
    double xlo, double xhi,
    double ylo, double yhi,
    double r, double g, double b
  );
  /* Fills given rectangle with given color */

void eps_fill_and_draw_rectangle(
    FILE *psfile,
    double xlo, double xhi,
    double ylo, double yhi,
    double r, double g, double b
  );
  /* 
    Fills rectangle with given color, then
    draws its outline with current pen.
  */

void eps_fill_circle(
    FILE *psfile,
    double xc, double yc, double radius,
    double r, double g, double b
  );
  /*
    Fills the circle with given center and radius, using the given color. */
 
void eps_fill_dot(
    FILE *psfile,
    double xc, double yc, double radius,
    double r, double g, double b
  );
  /*
    Same as ps_fill_circle, but the radius is in millimeters, irrespective 
    of the current scale. */
 
void eps_draw_circle(
    FILE *psfile,
    double xc, double yc, double radius
  );
  /* 
    Draws the circle with given center and radius, using the current pen and color. */
  
void eps_draw_dot(
    FILE *psfile,
    double xc, double yc, double radius
  );
  /* 
    Same as ps_draw_circle, but the radius is in millimeters, irrespective 
    of the current scale. */
  
void eps_fill_and_draw_circle(
    FILE *psfile,
    double xc, double yc, double radius,
    double r, double g, double b
  );
  /* 
    Fills the circle with given center and radius, using the given color,
    then draws its outline, using the current pen and color. */
  
void eps_fill_and_draw_dot(
    FILE *psfile,
    double xc, double yc, double radius,
    double r, double g, double b
  );
  /* 
    Same as ps_fill_and_draw_circle, but the radius is in millimeters, irrespective 
    of the current scale. */

void eps_fill_and_draw_lune(
    FILE *psfile,
    double xc, double yc, double radius, double tilt,
    double r, double g, double b
  );
  /* 
    Fills the lune with given center, radius, and tilt, using the given color,
    then draws its outline, using the current pen and color. */
  
void eps_fill_and_draw_slice(
    FILE *psfile,
    double xc, double yc, double radius, double start, double stop,
    double r, double g, double b
  );
   /* 
    Fills the pie slice with given center, radius, and angle range,
    using the given color, then draws its outline, using the current
    pen and color. */
    
void eps_fill_polygon(
    FILE *psfile,
    double x[], double y[],
    int n,
    double r, double g, double b
  );
  /*
    Fills the polygon (x[1],y[1]),.. (x[n],y[n]) with the given color level.
  */
  
void eps_draw_polygon(
    FILE *psfile,
    double x[], double y[],
    int n
  );
  /*
    Draws the contour of the polygon (x[1],y[1]),.. (x[n],y[n]) 
    using the current pen and color.
  */
  
void eps_fill_and_draw_polygon(
    FILE *psfile,
    double x[], double y[],
    int n,
    double r, double g, double b
  );
  /*
    Fills the polygon (x[1],y[1]),.. (x[n],y[n]) with the given color level,
    then draws its contour using the current pen and color.
  */
  
void eps_fill_triangle(
    FILE *psfile,
    double xa, double ya,
    double xb, double yb,
    double xc, double yc,
    double r, double g, double b
  );
  /* Fills triangle /abc/ with given color level. */

void eps_fill_grid_cell(FILE *psfile, int xi, int yi, double r, double g, double b);
  /*
    Fills the given cell of the current cell grid with the given
    color level.
  */

void eps_draw_coord_line (FILE *psfile, char axis, double coord);
  /* 
    Draws a reference line perpendicular to the given axis 
    at the given coordinate value.
  */

void eps_draw_grid_lines(FILE *psfile);
  /* 
    Draws the cell boundaries with the current pen and color level. */

void eps_set_label_font(FILE *psfile, char *font, float size);
  /*
    Sets the name and point size of the font to be used by ps_put_label. */

void eps_put_label(
    FILE *psfile, 
    char *text, 
    double x, double y, 
    float xalign, float yalign
  );
  /*
    Prints "label" at point (x,y), using the current label font size. 
    The parameter "xalign" (resp. "yalign)" specifies which point of the string's 
    bounding box will end up at (x,y): 0.0 means the left (resp. bottom) side,
    1.0 means the right (resp. top) side.  Default is (0.5, 0.5), meaning 
    the box will be centered at (x,y). */

void eps_draw_frame (FILE *psfile);
  /* 
    Draws a frame around the plotting area. (The frame 
    will extend half a line width outside the nominal bounding box.) */

void eps_add_caption(FILE *psfile, char *txt);
  /*
    Adds a caption text below the drawing, *outside* the nominal
    bounding box. */

/*** LOW-LEVEL HACKS ***/

void eps_put_text(FILE *psfile, char *text, char *newline);
  /* 
    Writes a text string to /psfile/, in Postscript form,
    handling special chars and parentheses.
    
    Replaces any embedded '\n' by the given /newline/ string.
  */

#endif
