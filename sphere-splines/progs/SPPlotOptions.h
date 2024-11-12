/* SPPlotOptions.h -- plotting options from command line. */
/* Last edited on 2005-10-02 19:40:39 by stolfi */

#ifndef SPPlotOptions_H
#define SPPlotOptions_H

#include <SPOptions.h>
#include <SPBasic.h>
#include <r3.h>
#include <SPH3.h>
#include <bool.h>
#include <values.h>

typedef struct SPPlotOptions_t
  { /* General plotting options: */
    bool_t eps;           /* TRUE to generate ".eps" instead of ".ps". */
    char *paperSize;      /* Document paper size (when {eps == FALSE}). */
    double figSize;       /* Figure size in mm. */
    double meshSize;      /* Plotting resolution in mm. */
    string_vec_t caption; /* Figure caption template. */
    /* View parameters: */
    bool_t projLonLat;    /* TRUE gives longutude/latitude plot instead of 3D. */
    SPH3_Point obs;       /* Position of observer; [0,0,0,0] if not specified. */
    double autoObs;       /* Degree of automatic adjustment for {obs}. */
    SPH3_Point upp;       /* Approximate {up} reference point, or [0,0,0,0]. */
    double radius;        /* Radius of interest, or 1.0 if whole {S^2}. */
    double autoRadius;    /* Degree of automatic adjustment for {radius}. */
    /* Parameters for PaintValues: */
    double fRange;        /* Nominal function range is {[-fRange _ fRange]}. */
    double autoRange;     /* Degree of automatic adjustment for {fRange}. */
    r3_t light;           /* Light source direction */
    double shadow;        /* Relative shadow darkening factor */
    /* Parameters for IsoLines: */
    double fStep;         /* Value difference between successive level curves */
    double lineWidth;     /* Line width in millimeters */
  } SPPlotOptions_t;  

#define SPPlotOptions_TriangHelp \
  "  [ -eps | -ps ] [ -paperSize STRING ] \\\n" \
  "  [ -caption TEXT ]... \\\n" \
  "  [ -obs OW OX OY OZ ] \\\n" \
  "  [ -up UW UX UY UZ ] \\\n" \
  "  [ -radius NUM ] \\\n" \
  "  [ -meshSize NUM ] \\\n" \
  "  [ -lineWidth MILLIMETERS ] \\\n" \
  "  [ -light DX DY DZ ] [ -shadow NUM ]"

SPPlotOptions_t SPPlotOptions_TriangParse(SPOptions_Parser_t *pp);
  /* Parses from the {pp} object a set of {SPPlotOptions}
    suitable for drawing triangulations only (without functions
    on them), as described by {SPPlotOptions_TriangHelp}. */

#define SPPlotOptions_FunctionHelp \
  "  [ -eps | -ps ] [ -paperSize STRING ] \\\n" \
  "  [ -caption TEXT ]... \\\n" \
  "  [ -obs OW OX OY OZ ] [ -autoObs NUM ] \\\n" \
  "  [ -up UW UX UY UZ ] \\\n" \
  "  [ -radius NUM ] [ -autoRadius NUM ] \\\n" \
  "  [ -fRange NUM ] [ -fStep NUM ] [ -autoRange NUM ] \\\n" \
  "  [ -figSize NUM ] [ -meshSize NUM ] \\\n" \
  "  [ -lineWidth MILLIMETERS ] \\\n" \
  "  [ -light DX DY DZ ] [ -shadow NUM ]"
    
SPPlotOptions_t SPPlotOptions_FunctionParse(SPOptions_Parser_t *pp);
  /* Parses from the {pp} object a set of {SPPlotOptions}
    suitable for drawing spherical functions, as 
    described by {SPPlotOptions_FunctionHelp}. */

#endif

