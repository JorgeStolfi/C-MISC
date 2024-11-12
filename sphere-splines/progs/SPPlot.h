/* SPPlot.h -- Plots spherical functions */
/* Last edited on 2007-12-23 21:57:29 by stolfi */

#ifndef SPPlot_H
#define SPPlot_H

#include <SPTriang.h>
#include <SPQuad.h>
#include <SPFunction.h>
#include <SPBasic.h>
#include <SPH3.h>
#include <pswr.h>
#include <r3.h>
#include <bool.h>

/* POSTSCRIPT FIGURE STREAMS */

typedef PSStream SPPlot_Stream;
  /* A stream of Postscript figures (separate EPS files or a single 
    PS document). */

SPPlot_Stream *SPPlot_NewStream
  ( bool_t eps,
    char *name,
    char *paperSize, 
    double figSize,
    bool_t projLonLat,
    int nCap
  );
  /* Creates a new Postscript stream as specified by the given
    options: either a standalone document, called "{name}.ps", or a
    bunch of EPS files, called "{name}-{NNNNNN}.eps".
    The {figSize} is in millimeters. The figure will be square
    if {projLonLat==FALSE}, and with a 2:1 aspect ratio if 
    {projLonLat==TRUE}. */

double SPPlot_DefaultFigSize
  ( bool_t eps, 
    char *paperSize, 
    int nRows, 
    int nCols,
    bool_t projLonLat,
    int captionLines
  );
  /* A convenient default figure size (in mm) for the given output options.
    If {eps = FALSE}, chooses the size so as to fit the given number of rows
    and columns in the page, with 1 inch margins.  If {eps = TRUE}, ignores
    the {paperSize} and returns 150.0 mm divided by the *minimum* of 
    {nCols,nRows}. */

/* PLOT COMPONENTS */

/* In the procedures that follow, the {map} argument must be a
  pespective map that takes the observer to {[0,0,0,1]}, the
  projection plane to the plane {Z = 0}, the horizontal image axis to
  the {X} axis, and the region of interest to the interior of a
  cylinder of radius 1 centered on the {Z} axis. See {SPPlot_PerspMap}
  below. 

  However, if {map} is NULL, the projection is defined to be the 
  conversion to spherical coordinates {lon,lat,rad}, where {lon} 
  and {lat} are scaled to the range {[-1.0_+1.0]} and {[-0.5_+0.5]},
  respectively. */

void SPPlot_Sphere(SPPlot_Stream *fps, SPH3_PMap *map);
  /* Draws the outline of the unit sphere, with the line width and
    color defined by {ps_set_pen}.  If {map} is NULL, draws the
    rectangle {[-1.0_+1.0] × [-0.5_+0.5]} */
    
void SPPlot_Axis
  ( SPPlot_Stream *fps, 
    SPH3_PMap *map, 
    r3_t *dir,
    double length,
    Sign sense
  );
  /* Draws an arrow from the origin, with the specified length
    and direction.  If {map} is NULL, does nothing.
    
    If {sense == -1}, draws the arrow only if it pierces the invisible
    part of the sphere; if {sense == +1}, draws it only if it pierces
    the visible part; if {sense == 0}, always draws the arrow. */
    
void SPPlot_CoordAxes
  ( SPPlot_Stream *fps, 
    SPH3_PMap *map, 
    Sign sense
  );
  /* Draws the three Cartesian coordinate axes X, Y, Z. The {sense}
    parameter has the same meanng as in to {SPPlot_Axis}. */

void SPPlot_Dot
  ( SPPlot_Stream *fps, 
    SPH3_PMap *map,  
    S2Point *p, 
    double radius, 
    bool_t fill,
    bool_t draw,
    bool_t backSide     
  );
  /* Draws a round dot at point {p} of the sphere,
    with the given radius (in millimeters, irrespective of the
    current scale.  

    Plots only the visible part of the segment, if {backSide = FALSE};
    or only the invisible part, if {backSide = TRUE}.  If {map} is
    NULL (lon-lat plot), assumes that all points are visible. */

void SPPlot_Segment
  ( SPPlot_Stream *fps, 
    SPH3_PMap *map,  
    S2Point *p, S2Point *q, 
    double step,
    bool_t backSide     
  );
  /* Draws the spherical segment (great-circle arc) from {p} to {q},
    approximating it by straight segments at most {step} long.  

    Plots only the visible part of the segment, if {backSide = FALSE};
    or only the invisible part, if {backSide = TRUE}.  If {map} is
    NULL (lon-lat plot), assumes that all points are visible. */

void SPPlot_GreatCircle
  ( SPPlot_Stream *fps, 
    SPH3_PMap *map, 
    r3_t normal,
    double step,
    bool_t backSide
  );
  /* Draws the main circle with equation {r3_dot(p,normal)==0}. */

void SPPlot_Triangulation
  ( SPPlot_Stream *fps,    
    SPH3_PMap *map,       
    Triangulation *tri, 
    double step,          /* Maximum step size (1.0 = sphere radius). */
    bool_t backSide       /* TRUE plots the back side instead of front side. */
  );
  /* Draws all edges of the spherical triangulation {tri}.
    
    Plots only the edges of the front (visible) side of the
    triangulation, if {backSide = FALSE}; or only the edges of the
    back (invisible) side, if {backSide = TRUE}. Edges are
    approximated by line segments at most {step} long (in radians,
    i.e. relative to the sphere radius). */
    
void SPPlot_Polyhedron
  ( SPPlot_Stream *fps,
    SPH3_PMap *map,
    Triangulation *tri,
    bool_t backSide        /* TRUE plots the back side instead of front side. */
  );
  /* Writes to the given Postscript file a perspective view of the the
    spherical triangulation {tri}, projected by the given perspective
    map. Assumes that the underlying polyhedron is convex and that its
    vertices lie on the unit sphere.
    
    Plots only the edges of the front (visible) side of the
    polyhedron, if {backSide = FALSE}; or only the edges of the
    back (invisible) side, if {backSide = TRUE}. */
    
void SPPlot_LonLatGrid
  ( SPPlot_Stream *fps,
    SPH3_PMap *map, 
    int gridNLon,            /* Number of grid steps in longitude. */
    int gridNLat,            /* Number of grid steps in latitude. */
    bool_t gridDots,         /* TRUE shows lon-lat grid as face-ctr dots, FALSE as rects. */
    double dotRadius,        /* Dot radius if {gridDots}. */     
    bool_t backSide          /* TRUE plots the back side instead of front side. */
  );
  /* Draws the longitude/latitude grid with given number of cells in 
    each direction.  If {gridDots}, the grid is drawn with face-centered dots
    with the given {dotRadius} (in millimeters, irrespective of the current scale). */

void SPPlot_SupportingPlane
  ( SPPlot_Stream *fps,
    SPH3_PMap *map, 
    SPH3_Plane *supp, 
    double angStep
  );
  /* Plots the circle which is the intersection of the supporting
    plane and the sphere, approximated by line segments of at most 
    {angStep} radians. */

/* PLOTTING SPHERICAL FUNCTION VALUES AND GRADIENTS */

/* In all procedures of this section, the sphere is approximated by a
  triangular mesh, by recursively subdividing each spherical triangle
  of {tri} into an almost-regular triangular array of {N^2}
  trianglets. Each original triangle of {tri} is slightly contracted
  before being subdivided, to allow for functions which are
  discontinuous along the edges of {tri}.
  
  If {tri} is null, a regular triangulation of the sphere into octants
  is used instead.
    
  The variables {fMinObs} and {fMaxObs} are updated with the extremal
  values of {func} observed during the plot. THEY MUST BE INITIALIZED
  BY THE CLIENT before the first call with a new function. */

/* Official default number of isolines to use in {[0 _ fRange]}: */
#define DefaultIsolines 10

void SPPlot_Isolines
  ( SPPlot_Stream *fps,
    SPH3_PMap *map,          /* Perspective projection matrix. */
    ScalarField func,        /* Function to plot. */
    Triangulation *tri,      /* Reference triangulation, of NULL */
    int N,                   /* Order of mesh subdivision. */
    double fMin,             /* Don't plot isolines below this limit. */
    double fMax,             /* Don't plot isolines above this limit. */
    double fStep,            /* Isoline spacing. */
    double fStart,           /* Synchronize isolines with this level. */
    double /*I/O*/ *fMinObs, /* Min {func} value seen (caller must initialize). */
    double /*I/O*/ *fMaxObs  /* Max {func} value seen (caller must initialize). */
  );  
  /* Plots the isolines (level curves) of the spherical function {func} 
    whose levels lie in in the closed interval {[fMin _ fMax]}. The 
    isoline levels are spaced {fStep} apart, and shifted so as 
    to include value {fStart}.
  
    More precisely, the sphere is approximated by a triangular mesh,
    as said above. The function {func} is linearly interpolated from
    its values at the corners of each face of the triangular mesh, and
    the isolines of that linear interpolation are plotted as straight
    lines. The function may be discontinuous along the edges of {tri}. */
 
void SPPlot_PaintValues(
    SPPlot_Stream *fps,
    SPH3_PMap *map,          /* Perspective projection matrix. */
    ScalarField func,        /* Function to plot. */
    Triangulation *tri,      /* Reference triangulation, of NULL */
    int N,                   /* Order of mesh subdivision. */
    double fMin,             /* Nominal minimum function value. */
    Color *cMin,             /* Color to use for {fMin}. */
    bool_t clipMin,            /* TRUE omits parts below {fMin}. */
    double fMax,             /* Nominal maximum function value. */
    Color *cMax,             /* Color to use for {fMax}. */
    bool_t clipMax,            /* TRUE omits parts above {fMax}. */
    r3_t *dLight,            /* Direction towards main light source. */
    double shadow,           /* Amount of darkening by shadow. */
    double /*I/O*/ *fMinObs, /* Min {func} value seen (caller must initialize). */
    double /*I/O*/ *fMaxObs  /* Max {func} value seen (caller must initialize). */
  );  
  /* Plots the values of the spherical function {func} as color
    gradations.  
    
    Values between {fMin} and {fMax} are mapped to colors between
    {cMin} and {cMax}, by linear interpolation. Values below {fMin}
    are either left blank ({clipMin == TRUE}) or painted with {cMin}
    ({clipMin == FALSE}).  Values above {fMax} are handled in the same
    way, depending on {clipMax}.
    
    Parts of the sphere away from the direction {dLight} are darkened
    slightly, and parts facing towards it are lightened. The {shadow}
    parameter specifies the maximum relative change in the color's
    intensity, in either sense. 
    
    The sphere is approximated by a mesh of triangles, defined by
    {tri} and {N}, as explained above.  Currently,
    each triangle of the mesh is painted with a solid color,
    computed from the average of the function values at 
    the three corners. */
 
void SPPlot_Vectors
  ( SPPlot_Stream *fps,
    SPH3_PMap *map,      /* Perspective projection matrix. */
    VectorField fld,     /* Vector field defined on the sphere. */
    Triangulation *tri,  /* Reference triangulation, of NULL */
    int N,               /* Order of mesh subdivision. */
    double dotRadius,    /* Dot radius (in mm, irrespective of scale). */
    double scale         /* Scale factor for vector field. */
  );  
  /* Plots the vector field {fld}, sampled at points on the sphere.
    Each sample vector {v = fld(p)} is drawn as a dot with a whisker
    {scale*length(v)} long (in world coordinates, where the sphere
    has radius 1.0). 
    
    The sphere is approximated by a triangular mesh, defined by {tri}
    and {N}, as explained above. Then, one dot-whisker is plotted at
    the center of each mesh triangle. */

/* DRAWING IT ALL */

void SPPlot_Everything
  ( SPPlot_Stream *fps,           /* Poststcript stream, ready to plot. */
    ScalarField func,        /* The function to plot. */
    Triangulation *tri,      /* Reference triangulation or NULL */
    double relMeshSize,      /* Maximum step/triangle size (radians). */
    bool_t showTriang,       /* TRUE to plot the triangulation {tri}. */
    SPH3_Plane *supp,        /* Supporting plane, or NULL. */
    double fRange,           /* Nominal maximum absolute value of function. */
    double fStep,            /* Isoline spacing. */
    SPH3_Point *obs,         /* Observer's position. */
    SPH3_Point *upp,         /* Camera vertical reference. */
    double rad,              /* Radius of region of interest. */
    r3_t *dLight,            /* Direction towards main light source. */
    double lineWidth,        /* Nominal line width in mm. */
    int gridNLon,            /* Number of grid steps in longitude. */
    int gridNLat,            /* Number of grid steps in latitude. */
    bool_t gridDots,         /* TRUE shows lon-lat grid as face-ctr dots, FALSE as rects. */
    bool_t verbose,          /* TRUE to print messages along the way. */
    double /*I/O*/ *fMinObs, /* Min {func} value seen (caller must initialize). */
    double /*I/O*/ *fMaxObs  /* Max {func} value seen (caller must initialize). */
  );
  /* Plots `everything' about the function {func}: sphere,
    triangulation (if not NULL), values, and isolines.
    
    If {fRange > 0}, adjusts the color scale assuming that the function
    values are in the range {[-fRange _ +fRange]}. If {fRange <= 0}, uses
    {[-m _ +m]} where {m} is the actual {max(fabs(func))} (estimated
    over a large set of sample points), rounded to a nice number. */
  
/* MULTIPLE VIEWS */

void SPPlot_MultipleViews
  ( SPPlot_Stream *fps,      /* Postscript stream. */
    char *funcTag,           /* Prefix for page names. */
    ScalarField func,        /* The function to plot. */
    Triangulation *tri,      /* Triangulation to draw, or NULL */
    double relMeshSize,      /* Maximum step/triangle size (radians). */
    bool_t showTriang,       /* TRUE to plot the triangulation {tri}. */
    SPH3_Plane *supp,        /* Supporting plane, or NULL */
    double fRange,           /* Nominal maximum absolute value of function. */
    double fStep,            /* Isoline spacing. */
    SPH3_Point *obs,         /* Observer's position */
    SPH3_Point *upp,         /* Camera vertical reference */
    double rad,              /* Radius of region of interest. */
    r3_t *dLight,            /* Direction towards main light source. */
    double lineWidth,        /* Nominal line width in mm. */
    int gridNLon,            /* Number of grid steps in longitude. */
    int gridNLat,            /* Number of grid steps in latitude. */
    bool_t gridDots,         /* TRUE shows lon-lat grid as face-ctr dots, FALSE as rects. */
    int aSide,               /* First side to plot (+1=front, -1=rear, 0=none) */
    int bSide,               /* Second side to plot (ditto). */
    string_vec_t *caption,   /* Figure caption, or NULL. */
    int index,               /* Figure `index' (for "%I" expansion in {caption}). */
    double time,             /* Figure `time' (for "%T" expansion in {caption}). */
    double error,            /* Figure `error' (for "%E" expansion in {caption}). */
    double capAlign,         /* Caption alignment (0.0=left, 0.5=ctr, 1.0=right). */
    bool_t verbose,          /* TRUE to print messages along the way. */
    double /*I/O*/ *fMinObs, /* Min {func} value seen (caller must initialize). */
    double /*I/O*/ *fMaxObs  /* Max {func} value seen (caller must initialize). */
  );
  /* Plots zero, one or two views of the the function {func}.
    
    Each view is produced by {SPPlot_Everything} with the given
    parameters, except that {obs} is replaced by its
    antipode for rear views ({aSide} or {bSide} euquals -1).
    
    If {fps} is an EPS stream, calls {pswr_new_canvas} before each view.
    The {figName} parameter is set to {funcTag} plus {-f} for front,
    or {-r} for back.
    
    The observer's position is {obs} for the front view, and for the back view.
    
    If the given {caption} is not NULL, it will be printed under each
    view, after expansion (see {SPPlot_ExpandCaption}). The 
    "%N" template spec expands to {funcName}. */
  
/* PERSPECTIVE PROJECTION */

#define SPPlot_FullRadius 1.3
  /* Default radius of interest: full sphere, plus arrow heads etc. */

SPH3_PMap *SPPlot_PerspMap
  ( SPH3_Point *obs, 
    double rad,
    SPH3_Point *upp
  );
  /* Computes a perspective projection matrix, suitable for
    use as the {map} parameter in the plotting procedures above.
  
    The perspective map assumes an observer at point {obs} (which must
    be outside the unit sphere), looking towards the origin. 
    
    The resulting projective map takes the observer to the
    {Z}-infinity point {(0,0,+oo)}, and the origin to some point on the
    {Z} axis. Then perspective projection reduces to applying the map
    and discarding the {Z} coordinate.
    
    The virtual camera is tilted around the projection axis {L} (the
    line from {obs} to the origin) in such a way that point {upp} will
    project somewhere onto the positive {Y} axis. The projection is
    then scaled uniformly so that the edge of the spherical cap of {S^2}
    centered on {L} with given {rad} becomes the unit circle. */

void SPPlot_FixObserver
  ( SPH3_Point *obs, 
    double alpha, 
    SPFunction *f, 
    bool_t verbose
    );
  /* Replaces the observer position {*obs} by a suitable mixture of
    the given {*obs} and the presumed ideal observer position 
    {idealObs = SPPlot_IdealObserver(f)}.
    
    The parameter {alpha} specifies the degree of mixing: {alpha = 0}
    leaves {*obs} unchanged, {alpha = 1} replaces it by {idealObs}.
    
    However, if the given {obs} is not a valid point (i.e., if {obs =
    NoPoint = [0,0,0,0]}), then ignores it and sets {*obs = idealObs},
    as if {alpha} was 1. If {obs == NULL}, does nothing. */

SPH3_Point IdealObserver(SPFunction *f, bool_t verbose);
  /* Returns an observer position {obs} that is hopefully good for {f}.
    
    If the function {f} is non-NULL and is an {SPSpline.T} with
    a non-trivial supporting plane, then {obs} is the point at
    infinity in the direction of the plane's normal. 
    
    Otherwise, if {f} is non-NULL, and the barycenter {b} of the
    squared function {f^2(p)} is not the origin, then {obs} is the
    direction of {b}. 
    
    If all these attempts fail, returns {obs = SPPlot_DefaultObserver()}. */

SPH3_Point SPPlot_DefaultObserver(void);
  /* A default observer position: currently {[0,3,2,1]}. */

void SPPlot_FixZenith(SPH3_Point *upp, SPH3_Point *obs, bool_t verbose);
  /* Replaces {*upp}, if necessary, by a point is a valid zenith 
    refrence point for the observer {*obs}: namely, a point that 
    is not {[0,0,0,0]} and is not collinear with {*obs} and the 
    origin.
     
    If the given {*upp} is already valid, leaves it unchanged. 
    Otherwise, replaces it either {SPPlot_DefaultZenith()}, or some other
    cardinal point at infinity. If {upp == NULL}, does nothing. */

SPH3_Point SPPlot_DefaultZenith(void);
  /* A default zenith reference point: currently the {Z}-infinity
    point {[0,0,0,1]}. */

/* Maximum radius of interest less than 1 which is worth returning: */
#define MaxRadiusOfInterest 0.90

/* Minimum radius of interest that seems safe to plot: */
#define MinRadiusOfInterest 0.01

void SPPlot_FixRadiusOfInterest
  ( double *rad, 
    double alpha, 
    SPFunction *f, 
    SPH3_Point *obs,
    bool_t verbose
  );
  /* Replaces {*rad} by a suitable mixture of its given value {*rad} 
    and its presumed ideal value 
    {idealRad = SPPlot_IdealRadiusOfInterest(f,obs)}.
    
    The parameter {alpha} specifies the degree of mixing: {alpha = 0}
    leaves {*rad} unchanged, {alpha = 1} replaces it by {idealRad}.
    However, the given {*rad} is invalid (negative), ignores it and
    sets {*rad = idealRad}, as if {alpha} was 1.
    
    If {alpha != 0}, and the computed radius is bigger than
    {MaxRadiusOfInterest}, returns {FullRadius} (show the entire
    sphere, plus axes). 
    
    If {rad == NULL}, does nothing. In all other cases (including
    {alpha = 0}), makes sure that the final {*rad} is not smaller
    than {MinRadiusOfInterest}. */

double SPPlot_IdealRadiusOfInterest(SPFunction *f, SPH3_Point *obs, bool_t verbose);
  /* Returns the radius {rad} of a spherical cap {B} of {S^2} that is
    centered along the line {L} connecting {obs} to the origin, and
    contains the all the `interesting' parts of the function {f}.
    
    The cap {B} is always contained in the part of {S^2} that is
    visible from {obs}. Its radius {rad} is measured perpendicular to
    {L}. Let {D} be the flat disk whose border coincides with that of
    {B}.
    
    The radius {rad} may range from 0 (exclusive) when {B} reduces
    to the point {dir(obs)}, to 1 (inclusive), when {B} is the whole
    hemisphere of {S^2} turned towards {obs}.
    
    If the function {f} is non-NULL and is an {SPSpline.T} with
    a non-trivial supporting plane, then {B} is big enough to contain
    the cap of {S^2} determined by the positive side of that plane.
    
    Otherwise, if {f} is non-NULL, the procedure tries to choose {B}
    big enough to contains most of the integral of {f^2(p)}.
    
    If all these attempts fail, returns {rad == 1}. In particular, that
    is the case when {f} has interesting parts that extend beyond the
    horizon of {S^2} as seen from {obs}. */
    
void SPPlot_FixView
  ( SPH3_Point *obs, 
    double autoObs, 
    SPH3_Point *upp, 
    double *rad,
    double autoRad,
    SPFunction *f, 
    bool_t verbose
  );
  /* Calls {SPPlot_FixObs}, {SPPlot_FixZenith}, {SPPlot_FixRadiusOfInterest},
    with given arguments. Any of {f}, {obs}, {upp}, and {rad}
    may be null, in which case some sensible thing will happen. If {verbose}
    is true, prints the adjusted values to {stderr}. */

void SPPlot_FixRange
  ( double *fRange, 
    double *fStep,
    double alpha, 
    SPFunction *f, 
    bool_t verbose
  );
  /* Replaces {*fRange} by a suitable mixture of the given {*fRange} and
    the estimated maximum absolute value of {f(p)}, namely 
    {estRange = SPPlot_EstimatedRange(f)}.
    
    The parameter {alpha} specifies the degree of mixing: {alpha = 0}
    leaves {*fRange} unchanged, {alpha = 1} replaces it by {estRange}.
    However, if the given {*fRange} is not valid (i.e., {fRange <= 0}),
    then replaces it by {estRange}, as if {alpha} was 1.
    
    If {alpha} is not zero, the computed range {*fRange} is rounded to a
    nice value.
    
    In any case (even when {alpha = 0}), checks whether the final 
    {*fStep} is reasonable for the final {*fRange}; if not, fixes 
    it appropriately. */

double SPPlot_EstRange(SPFunction *f, bool_t verbose);
 /* Returns an estimate of the maximum of {fabs(f(p))} over the sphere. */

/* CAPTIONS */

char *SPPlot_ExpandCaptionLine
  ( char *tmp,
    int side,
    int index,
    double time,
    double error,
    double fRange,
    double fStep,
    char *name
  );
  /* Creates a copy of the template string {tmp} expanding 
    certain occurrences of "%". Specifically:

      "%T" means {time} (with 6 decimal fraction figures).
      "%E" means {error} (with 2 significant figures).
      "%S" means "f" if {side > 0}, "r" otherwise.
      "%R" means the level {fRange} of the maximum isoline.
      "%D" means the isoline spacing {fStep}.
      "%I" means the {index} (formatted as 6 digits zero-padded).
      "%N" means the {name} string.

    Other occurrences of "%" are left unchanged. */

string_vec_t SPPlot_ExpandCaption
  ( string_vec_t *tmp,
    int side,
    int index,
    double time,
    double error,
    double fRange,
    double fStep,
    char *name
  );
  /* Creates a copy of the string vector {tmp} expanding 
    each line with {SPPlot_Expand_Caption_Line}. */
    
void SPPlot_FreeCaption(string_vec_t *cap);
  /* Calls {free} on all strings {cap.e[k]}, and makes {cap} empty. */

/* MISCELLANEOUS */

S2Point SPPlot_XYZFromLatLon(double lat, double lon);
  /* Cartesian coordinates from latitude and longitude (in radians). */
  
double SPPlot_RoundToNice(double x);
  /* Rounds the number {x}, which must be positive, to a 
    nice value (namely 0.25, 0.50, or 1.00 times a power of 10). */

char *SPPlot_IterationTag(int iter);
  /* Formats {iter} as a six-digit number {NNNNNN}, except
    that {-1} is mapped to "ini" and {INT_MAX} is mapped to "fin". */

#endif
