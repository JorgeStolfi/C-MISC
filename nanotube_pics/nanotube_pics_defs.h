// Last edited on 2019-05-13 02:31:05 by jstolfi

#define _GNU_SOURCE
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <assert.h>

#include <pswr.h>
#include <jsmath.h>
#include <jsfile.h>
#include <frgb.h>
#include <frgb_ops.h>
#include <affirm.h>

typedef struct nanotube_pics_elem_style_t
  {
    // For stage 0:
    double rwh;    // Extra widths of halos on both sides.
    // For stage 2:
    frgb_t dcol;   // Stoke color.
    frgb_t fcol;   // Fill color.
    double lw;     // Width of strokes.
    double rdot;   // Radius of dots.
    double arlen;  // Arrow length.
    double arwid;  // Arrow width.
    char *font;        // Font for text.
    double fontsize;   // Font size for text.
  } nanotube_pics_elem_style_t;
  /* Style for an element - arrow, label, atom, bond. */

typedef struct nanotube_pics_style_t
  { 
    double xtot; // Canvas width (mm).
    double ytot; // Canvas height (mm).
    
    // Size of repeating pattern of hex grid (2 atoms).
    // The midpoint of the left edge is at the midpoint of a down-slanted bond.
    double dy;
    double dx;
    
    // Length of lattice bond (distance between atom centers):
    double bdlen;
    
    // Origin of hex grid:
    double b0x, b0y;   // Coordinates of ref point of lower left cell.
    double orgx, orgy; // Coordinates of center of left (class 0) atom in lower left cell.
    
    frgb_t bgcol; //  Background color.

    // Parameters for general labels:
    nanotube_pics_elem_style_t labsty; 

    // Parameters for reference vectors and dots:
    nanotube_pics_elem_style_t rfvsty; 

    // Parameters for strip and sector boundary lines:
    nanotube_pics_elem_style_t edgsty; 

    // Parameters for basis vectors and dots:
    nanotube_pics_elem_style_t bassty; 

    // Parameters for valid types figure:
    nanotube_pics_elem_style_t vatsty; // Vectors and dots of non-isometric types.
    nanotube_pics_elem_style_t vbtsty; // Vectors and dots of enantiomers.
    nanotube_pics_elem_style_t vntsty; // Number pair labels.

    // Parameters for key variables and values:
    nanotube_pics_elem_style_t kyvsty; // Latin variables and units.
    nanotube_pics_elem_style_t kynsty; // Numbers.
    nanotube_pics_elem_style_t kyssty; // Greek variables and symbols.

    // Parameters for graphene and nanotube lattice elements:
    nanotube_pics_elem_style_t bdsty; 

    nanotube_pics_elem_style_t atsty[2]; //Atom style by class.

    // Halo parameters:
    frgb_t whcol;  //  Color of halos.

    // Shadow parameters
    double shx, shy; // Shadow isplacements.
    frgb_t shcol;    //  Color of shadows.
    double rsh;      //  Extra width of shadow, on each side.
        
  } nanotube_pics_style_t;

#define mm (72.0/25.4)
  /* One millimeter in Postscript points. */

#define px (0.1*mm)
  /* One "pixel" in Postscript points, as in the old SVG version. */

nanotube_pics_style_t *nanotube_pics_def_style(double xtot, double ytot);
  /* Defines the standard style parameters for a figure. The total
    width and heigh of the figure (in pt) will be {xtot,ytot}. */

void nanotube_pics_draw_fig_master(void);

void nanotube_pics_draw_fig_strip(int32_t n, int32_t m);

void nanotube_pics_fig_strip_choose_a1
  ( nanotube_pics_style_t *sty,
    double a0x,
    double a0y,
    double ux,
    double uy,
    double vx,
    double vy,
    double krx,
    double kry,
    double wx,
    double wy,
    double *a1xP,
    double *a1yP
  );
  /* Chooses a suitable position for the left 
    end {a1x,a1y} of the {wx,wy} vector, that avoids the margins, 
    and keeps the strip far from the key and the fundamental basis.
    
    Assumes that all components of {sty} are defined.
    Assumes that the fund basis has origin {a0x,a0y} and vectors
    {ux,uy} and {vx,vy}.  Assumes that the corner of the key nearest
    to the center of the pic is {krx,kry}.
    
    Returns the result in {*a1xP,*a1yP}. */ 
    

void nanotube_pics_draw_fig_strip_labels
  ( PSStream *ps,
    nanotube_pics_style_t *sty, 
    double a1x,
    double a1y,
    double a2x,
    double a2y,
    int32_t stage
  );
  /* Draws labels and arrows for the strip figure.
    Assumes that {(a1x,a1y)} and {(a2x,a2y)} are the corresponding
    lattice points. */

void nanotube_pics_draw_fig_master_labels
  ( PSStream *ps,
    nanotube_pics_style_t *sty, 
    double a1x,
    double a1y,
    double ux,
    double uy,
    double vx,
    double vy,
    int32_t stage
  );
  
void nanotube_pics_draw_basis_labels
  ( PSStream *ps,
    nanotube_pics_style_t *sty, 
    double a0x,
    double a0y,
    double ux,
    double uy,
    double vx,
    double vy,
    int32_t stage
  );
  /* Draws the fundamental basis vectors {(ux,uy)} and {(vx,xy)},
    with labels, based at the point {(a0x,a0y)}. */

void nanotube_pics_draw_fig_key
  ( PSStream *ps,
    nanotube_pics_style_t *sty, 
    int32_t n,
    int32_t m,
    double kxmin,
    double kxmax,
    double kymin,
    double kymax,
    int32_t stage
  );
  /* Draws the key with values of {n,m,c,alpha}.
    The top line will end at {(kx,xy)}. */

double nanotube_pics_draw_key_line
  ( PSStream *ps,
    nanotube_pics_style_t *sty, 
    double qx,
    double qy, 
    char *lab,
    double val,
    int32_t digs,
    char *unit,
    double duy,
    int32_t stage,
    nanotube_pics_elem_style_t *evsty,
    nanotube_pics_elem_style_t *ensty,
    nanotube_pics_elem_style_t *eusty
  );
  /* Draws one line of the key "{lab} = {val} {unit}".  The value will be 
    formatted with {digs} fraction digits. The label will be drawn
    with style {evsty}, the value with style {ensty}, and the unit 
    with style {eusty}.  The unit is displaced up by {duy}.
    
    Returns the estimated width of the text on the figure (in pt),
    including any leading or trailing blanks. */

void nanotube_pics_draw_lattice
  ( PSStream *ps,
    nanotube_pics_style_t *sty,
    double a1x,
    double a1y,
    double wx,
    double wy,
    int32_t stage,
    double dim
  );
  /* Draws the graphene lattice.  The reference point of the first cell
    will be at {(sty->b0x,sty->b0y)}.  The other parameters are as in
    {nanotube_pics_draw_lattice_cell}. */
  
void nanotube_pics_draw_lattice_cell
  ( PSStream *ps,
    nanotube_pics_style_t *sty,
    double x, 
    double y, 
    double a1x,
    double a1y,
    double wx,
    double wy,
    int32_t stage,
    double dim
  );
  /* Draws a lattice cell, containing two atoms at about same height:
    one of class 0 (left), one of class 1 (right). 
    
    The cell is assumed to be bounded by a rectangle of width {sty->dx} and
    height {sty->dy}, that cuts the 4 bonds that leave the two atoms
    at their midpoints. The midpoint of the leftmost bond will be at
    {(x,y)}.
 
    The parameters {a1x,a1y,wx,wy,stage,dim} are as in {nanotube_pics_draw_atom}. */
    
void nanotube_pics_draw_atom
  ( PSStream *ps,
    nanotube_pics_style_t *sty,
    double cx,
    double cy,
    double sx,
    double sy,
    double a1x,
    double a1y,
    double wx,
    double wy,
    int32_t stage,
    int32_t class,
    double dim 
  );
  /* Draws an atom of the graphene lattice and its three bonds. The
    atom will be centered at {(cx,cy)}. It will have one vertical bond
    and two diagonal bonds. Let {A,B,C} be the points on those three
    bonds that lie midway to the next atom The bounding box of {A,B,C}
    will be a rectangle of width {sx} and height {|sy|}. The vertical
    bond points up if {sy > 0}, down if {sy < 0}.
    
    The atom will be drawn only if it is inside the strip whose width
    is the segment from {(a1x,a1y)} to {(a1x,a1y)+(wx,wy)}, or on the
    edge of that strip that goes through the first point. 
    If {wx=wy=0}, the atom will be drawn anyway.
    
    If {stage} is 0, paints only the blank halo (antishadow) of the atom and bonds. 
    
    If {stage} is 1, paints only the shadow of the atom and bonds.
    
    If {stage} is 2, paints the atom and bonds, with the color of the
    the given {class}, effaced by the {dim} factor, between 0 (no dimming) 
    and 1 (fully dimmed to the background color). */

void draw_lattice_element
  ( double x, 
    double y, 
    double dx, 
    double dy,
    int32_t stage,
    double dim
  );
  /* Draws one cell of the graphene lattice.  */

void nanotube_pics_draw_strip_edges
  ( PSStream *ps,
    nanotube_pics_style_t *sty,
    double cx,
    double cy,
    double tx,
    double ty,
    int32_t stage
  );
  /* Draws two lines through {(cx,cy)} and {(cx,cy)+(tx,ty)}
    perpendicular to {(tx,ty)}, all across the figure.  
    
    The color and linewidth are defined by the style {sty}.
    The parameter {stage} is as in 
    {nanotube_pics_draw_perp_line}. */ 

// LOW-LEVEL GRAPHIC ELEMENTS

/* In all these procedures, if {stage} is 2, the line is drawn (and filled if solid)
  with the given style parameters {esty}.  If {stage} is 0, paints the object's
  white halo instead; if {stage} is 1, paints the object's shadow. */

void nanotube_pics_draw_perp_line
  ( PSStream *ps,
    nanotube_pics_style_t *sty,
    double cx,
    double cy,
    double tx,
    double ty,
    bool_t dash, 
    int32_t stage,
    double dim,
    nanotube_pics_elem_style_t *esty
  );
  /* Draws a line through {(cx,cy)} perpendicular to {(tx,ty)},
    all across the figure.  If {dash} is true the line will be dashed. */ 

void nanotube_pics_draw_ray
  ( PSStream *ps,
    nanotube_pics_style_t *sty,
    double cx,
    double cy,
    double wx,
    double wy,
    bool_t dash, 
    int32_t stage,
    double dim,
    nanotube_pics_elem_style_t *esty
  );
  /* Draws an "infinite" ray from {(cx,cy)} in the direction of {(wx,wy)}. */
  
void nanotube_pics_draw_segment
  ( PSStream *ps,
    nanotube_pics_style_t *sty,
    double ppx,
    double ppy,
    double qqx,
    double qqy,
    bool_t dash,
    bool_t arrow,
    int32_t stage,
    double dim,
    nanotube_pics_elem_style_t *esty
  );
  /* Draws a segment from {(ppx,ppy)} to {(qqx,qqy)}.  If {dash}
    is true, the line will be dashed, otherwise it will be solid.
    If {arrow} is true, also draws a arrowhead at the tip. */

void nanotube_pics_draw_dot 
  ( PSStream *ps,
    nanotube_pics_style_t *sty, 
    double cx,
    double cy, 
    int32_t stage,
    double dim,
    nanotube_pics_elem_style_t *esty
  );
  /* Draws a solid dot with center {(cx,cy)} and radius {esty->rdot}.  */
  
double nanotube_pics_draw_label
  ( PSStream *ps,
    nanotube_pics_style_t *sty, 
    double qx,
    double qy, 
    double tx,
    double ty, 
    double dst,
    char *lab,
    double xalg, 
    double yalg,
    int32_t stage,
    double dim,
    nanotube_pics_elem_style_t *esty
  );
  /* Draws the label {lab} at {(qx,qy)} with displacement
    {dst} in the direction of the vector {{tx,ty)}.
    The parameters {xalg,yalg} specify which point of the label
    will be placed at that point.  
    
    Returns the estimated width of the text on the figure (in pt),
    INCLUDING any leading or trailing blanks.  Note that the 
    text is positioned in the figure IGNORING any leading and
    trailing blanks. */

void nanotube_pics_stage_fudge
  ( nanotube_pics_style_t *sty,
    int32_t stage,
    double dim,
    nanotube_pics_elem_style_t *esty,
    double *exP,
    double *eyP,
    nanotube_pics_elem_style_t *fstyP
  );
  /* Adjusts colors and positions for a given stage of drawing.
    Assumes {esty} is the object's style when {stage} is 
    2 and {dim} is zero. Returns in {*fstyP} the style to actually use.
    Also returns in {*exP,*eyP} the displacement to apply to the object. */

PSStream *nanotube_pics_fig_open(nanotube_pics_style_t *sty, char *fname);
void nanotube_pics_fig_close(PSStream *ps);

int main(int argc, char**argv);

