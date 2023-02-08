/* See SO2DPlot.h */
/* Last edited on 2007-12-23 22:46:42 by stolfi */ 

#include <SO2DPlot.h>

#include <SOBasic.h>
#include <SOGrid.h>
#include <SOFunction.h>
#include <SOBasic.h>

// #include <SOeps.h>

#include <dg_grid.h>

#include <vec.h>
#include <affirm.h>
#include <nat.h>
#include <pswr.h>
#include <stdio.h>
#include <math.h>

/* INTERNAL DATA TYPES */

typedef int Isoline; /* Index of an isoline. */
  /* By definition, the isoline with index 
    {k} is where the function has value {(k+0.5)*fStep}, where 
    {fStep} is the spacing between levels. */

/* Color near zero: */
static Color ZERO_PAINT_COLOR = (Color){{1.00, 1.00, 1.00}};

/* Limiting colors for color bands (brightness = 0.5): */
static Color POS_PAINT_COLOR = (Color){{1.00, 0.33, 0.00}};
static Color NEG_PAINT_COLOR = (Color){{0.00, 0.67, 1.00}};

/* Colors for isoline drawing (brightness = 0.20): */
static Color POS_DRAW_COLOR = (Color){{0.57, 0.00, 0.00}};
static Color NEG_DRAW_COLOR = (Color){{0.00, 0.17, 1.00}};

/* Scale factor from cell-relative Y (ranges over [0_1]) to plot Y. */
#define YSCALE ((double)0.70710678118654752440)

/* GLOBAL VARIABLES */

static Color *bandColor = NULL;   /* color table */
/* The band between isolines {k-1} and {k} is plotted with color 
  {bandColor[k - kPlotMin]} where {kPlotMin} is the index of the 
  first isoline to be plotted. */
  
static double isolineWidth = 0.15; /* Thickness of isolines in mm. */
/* This global variable should become unnecessary once the 
  setting of line width gets separated from the
  setting of line color in {pswr.h}. */

/* INTERNAL PROTOTYPES */

/* The procedures {SO2DPlot_ValuesInRectangle} and
  {SO2DPlot_ValuesInTriangle} plot only those isolines with indices
  in the range {kPlotMin..kPlotMax}. The points of the root cell where
  {f} lies below isoline {kPlotMin} are plotted as a single color
  band; and ditto for values above isoline {kPlotMax}.
  
  The {layer} parameter selects between color band painting ({layer == 0})
  and isoline plotting ({layer == 1}). */

void SO2DPlot_ValuesInRectangle
  ( SOPlot_File *psf,    /* Plot file. */
    SOFunction *f,       /* Function to be plotted. */
    interval_t boxR[2], /* Rectangle {C} relative to root cell. */
    Isoline kPlotMin,    /* Index of lowest isoline to plot. */
    Isoline kPlotMax,    /* Index of highest isoline to plot. */
    double fStep,        /* Value increment between isolines/bands. */
    int layer,           /* 0 = paint color bands, 1 = draw isolines. */ 
    double *fObsMin,     /* (IN/OUT) Min value seen during plot. */ 
    double *fObsMax      /* (IN/OUT) Max value seen during plot. */ 
  );
  /* Plots color bands or isolines within the rectangle boxR (in root-cell-relative
    coordinates). Also updates the values of {fObsMin} and {fObsMax},
    which must be initialized by the caller. */

void SO2DPlot_ValuesInTriangle
  ( SOPlot_File *psf,                /* Plot file. */
    double Px, double Py, double Pf, /* Coordinates and value at vertex {P}. */
    double Qx, double Qy, double Qf, /* Coordinates and value at vertex {Q}. */
    double Rx, double Ry, double Rf, /* Coordinates and value at vertex {R}. */
    Isoline kPlotMin,                /* Index of lowest isoline to plot. */
    Isoline kPlotMax,                /* Index of highest isoline to plot. */
    double fStep,                    /* Value increment between isolines/bands. */
    int layer                        /* 0 = paint color bands, 1 = draw isolines. */ 
  );
  /* Plots color bands or isolines within the triangle with corners
    {(Px,Py), (Qx,Qy), (Rx,Ry)} (in root-cell-relative coordinates).
    Function values are interpolated between the given values {Pf,Qf,Rf}
    at those corners. */

void SO2DPlot_PaintRectangle
  ( SOPlot_File *psf,
    interval_t boxR[2],
    Color *c
  );
  /* Paints the rectangle {boxR} (in root-cell-relative coordinates)
    with uniform color {c}. */

void SO2DPlot_PaintTriangle
  ( SOPlot_File *psf,
    double Ax, double Ay, 
    double Bx, double By,
    double Cx, double Cy,
    Color *c
  );
  /* Paints the triangle defined by its corners {(Ax,Ay), (Bx,By),
    and (Cx,Cy)} (in root-cell-relative coordinates) with uniform 
    color {c}. */

void SO2DPlot_DrawSegment 
  ( SOPlot_File *psf,
    double Ax, double Ay, 
    double Bx, double By,
    Color *c
  );
  /* Draws the segment with endpoints {(Ax,Ay)} and {(Bx,By)}
    (in root-cell-relative coordinates) with uniform color {c}. */

void SO2DPlot_SetColors
  ( Isoline kPlotMin, /* Index of first isoline to be drawn. */
    Isoline kPlotMax, /* Index of last isoline to be drawn. */
    double fStep      /* Value step between isolines. */
  );
  /* Initializes the color table for band painting. The band between
    isolines {-1} and {0} is painted white. If {kPlotMin < 0}, the
    band below isoline {kPlotMin} is painted with {NEG_PAINT_COLOR}.
    If {kPlotMax > 0}, the band above isoline {kPlotMax} is painted
    with {POS_PAINT_COLOR}. Intermediate isolines are painted with
    intermediate colors. */

void SO2DPlot_Subtree
  ( FILE *psf, 
    dg_tree_node_t *p, 
    dg_cell_index_t k,
    dg_rank_t rank,
    interval_t boxR[2],
    int maxDepth
  );
  /* Draws the cells of the dyadic grid rooted at node {p}, which is
    assumed to correspond to a cell {C} of index {k} and given {rank}.
    The boundary of {C} must be drawn by the client. The drawing is 
    clipped to the rectangle {boxR} (in 
    cell-root-relative coordinates).  Cells at level higher than 
    {maxDepth} are ignored, i.e. the subtree is implicitly truncated
    at level {maxDepth}. */

/* IMPLEMENTATIONS */

void SO2DPlot_Values
  ( SOPlot_File *psf,     /* Plot file. */
    SOFunction *f,        /* Function to be plotted. */
    SOGrid_Tree *tree,    /* Reference grid for adaptive domain subdivision. */ 
    interval_t boxR[2],  /* Rectangle relative to root cell. */
    double fPlotMin,      /* Minimum value for color mapping. */
    double fPlotMax,      /* Maximum value for color mapping. */
    double fStep,         /* Value increment between isolines/bands. */
    double lineWidth,     /* Line width for isoline plotting. */
    int minDepth,         /* Minimum subdivision depth for plotting. */
    int extraDepth,       /* Extra depth of subdivision for leaves of {tree}. */
    int maxDepth,         /* Maximum subdivision depth for plotting. */
    bool_t bands,           /* TRUE paints between isolines, FALSE leaves blank. */
    bool_t isolines,        /* TRUE draws the isolines, FALSE omits them. */ 
    double *fObsMin,      /* (IN/OUT) Min value seen during plot. */ 
    double *fObsMax       /* (IN/OUT) Max value seen during plot. */ 
  )
  { 
    interval_t boxC[2];   /* Coordinates of a cell {C}. */

    /* Indices of first and last isolines. */
    /* Allow for rounding errors in the parameters {fPlotMin}, {fmax}. */
    Isoline kPlotMin = (int)(ceil(fPlotMin/fStep + 0.4999999));
    Isoline kPlotMax = (int)(floor(fPlotMax/fStep + 0.5000001));
    
    int layer; /* 0 = bands, 1 = isolines; */
  
    /* Set isoline width for internal procs. */
    isolineWidth = lineWidth;

    auto void SO2DPlot_ValuesInCell
      ( dg_cell_index_t k,       /* The index of a dyadic (sub)cell. */
        dg_rank_t rank,     /* Its rank. */
        dg_tree_node_t *nd,      /* A node of {tree}, or NULL. */
        dg_rank_t lastRank  /* Rank of lowest node containing the cell. */
      );
      /* Plots color bands (if {layer == 0} or isolines (if {layer ==
        1}) of function {f}, inside the dyadic cell {C} in level
        {rank} whose index is {k}. Assumes that {nd} is the
        corresponding node in {tree}, or is NULL if such node doesn't
        exist. In either case, assumes {lastRank} is the rank of the
        lowest node of {tree} whose cell contains {C}; or is 0 if
        {tree == NULL}. */
   
    void SO2DPlot_ValuesInCell
      ( dg_cell_index_t k,       /* The index of a dyadic (sub)cell. */
        dg_rank_t rank,     /* Its rank. */
        dg_tree_node_t *nd,      /* A node of {tree}, or NULL. */
        dg_rank_t lastRank  /* Rank of lowest node containing the cell. */
      )
      {   
        int locMinDepth = lastRank + extraDepth;
        // DEBUG fprintf (stderr, "+InCell cell = %ld rank = %d", k, rank); 
        
        /* Get coordinates of plot cell relative to the root cell: */
        dg_cell_box_root_relative(2, k, boxC);
        
        // DEBUG fprintf (stderr, "  X = [%8.6f _ %8.6f]", LO(boxC[X]), HI(boxC[X]));
        // DEBUG fprintf (stderr, "  Y = [%8.6f _ %8.6f]", LO(boxC[Y]), HI(boxC[Y]));
        // DEBUG fprintf (stderr, "\n");
        
        /* Skip if cell lies outside the desired rectangle. */
        if ((LO(boxC[X]) >= HI(boxR[X])) || (HI(boxC[X]) <= LO(boxR[X]))) { return; }
        if ((LO(boxC[Y]) >= HI(boxR[Y])) || (HI(boxC[Y]) <= LO(boxR[Y]))) { return; }
        
        /* Subdivide and recurse, or plot undivided: */
        
        if ((rank < maxDepth) && ((rank < minDepth) || (rank < locMinDepth)))
          { /* Split {C} into subcells, even if it is a leaf cell. */
            dg_tree_node_t *nd0 = (nd ? LOCH(*nd): NULL);
            dg_tree_node_t *nd1 = (nd ? HICH(*nd): NULL);
            SO2DPlot_ValuesInCell(2*k,   rank+1, nd0, lastRank + (nd0 != NULL));
            SO2DPlot_ValuesInCell(2*k+1, rank+1, nd1, lastRank + (nd1 != NULL));
          }
        else 
          { /* Clip cell to given rectangle: */
            if (HI(boxC[X]) > HI(boxR[X])) { HI(boxC[X]) = HI(boxR[X]); } 
            if (LO(boxC[X]) < LO(boxR[X])) { LO(boxC[X]) = LO(boxR[X]); } 
            if (HI(boxC[Y]) > HI(boxR[Y])) { HI(boxC[Y]) = HI(boxR[Y]); }  
            if (LO(boxC[Y]) < LO(boxR[Y])) { LO(boxC[Y]) = LO(boxR[Y]); } 

            /* Plot approximate isolines in rectangle. */
            SO2DPlot_ValuesInRectangle
              ( psf, f, boxC,
                kPlotMin, kPlotMax, fStep, 
                layer, fObsMin, fObsMax
              );
          }
        // DEBUG fprintf (stderr, "-InCell cell = %ld rank = %d\n", k, rank); 
      }
     
    // DEBUG fprintf (stderr, "  minDepth = %d extraDepth = %d maxDepth = %d\n", minDepth, extraDepth, maxDepth); 
        
    if (! (isolines || bands)) { return; }

    affirm(fPlotMin < fPlotMax, "bad function range");
    affirm(kPlotMin <= kPlotMax, "no isolines in range");
    
    SO2DPlot_SetColors(kPlotMin, kPlotMax, fStep);

    if (bands)
      { layer = 0; 
        SO2DPlot_ValuesInCell (1, 0, (tree ? tree->root : NULL), 0);
      }

    if (isolines) 
      { layer = 1;
        SO2DPlot_ValuesInCell (1, 0, (tree ? tree->root : NULL), 0);
      }
    fprintf(stderr, 
      "SO2DPlot_Values: fObsMin = %26.16e, fObsMax =  %26.16e\n",
      (*fObsMin), (*fObsMax)
    );
  }

/* Number of sample points in each rectangle (either 4 or 8): */
#define NSAMPLES 8

void SO2DPlot_ValuesInRectangle
  ( SOPlot_File *psf,    /* Plot file. */
    SOFunction *f,       /* Function to be plotted. */
    interval_t boxR[2], /* Rectangle {C} relative to root cell. */
    Isoline kPlotMin,    /* Min isoline index to plot. */
    Isoline kPlotMax,    /* Max isoline index to plot. */
    double fStep,        /* Value increment between isolines/bands. */
    int layer,           /* 0 = paint color bands, 1 = draw isolines. */ 
    double *fObsMin,     /* (IN/OUT) Min value seen during plot. */ 
    double *fObsMax      /* (IN/OUT) Max value seen during plot. */ 
  )
  { 
    double p[(NSAMPLES+1)*2]; /* Sampling points, root-relative coords. */
    /* Samples {[0..NSAMPLES-1]} lie around the perimeter of {C}, */
    /* in CCW order; sample {[NSAMPLES]} is the center of the rectangle. */
    double fp[NSAMPLES+1];    /* Sampled function values: {fp[i] = f(p[i])}. */

    double *ctr = &(p[NSAMPLES*2]);
    
    double fpMin, fpMax;   /* Min and max of {fp[0..4]}. */
    Isoline kpMin, kpMax;  /* Isolines in {[fpMin _ fpMax]} are {kpMin..kmax}. */
    int i; 

    // DEBUG fprintf (stderr, "+InRectangle\n"); 

    /* Coordinates of rectangle center: */
    ctr[X] = (LO(boxR[X])+HI(boxR[X]))/2.0;
    ctr[Y] = (LO(boxR[Y])+HI(boxR[Y]))/2.0;
    
    /* Gets coordinates and function values of corners and mid-edges: */
    /* Important: {fp[0..7]} must be in counterclockwise order. */
    for (i = 0; i < NSAMPLES; i++)
      { double *pi = &(p[i*2]);
        if (NSAMPLES == 4)
          {  /* Samples are corners only. */
            pi[X] = ((i==0) || (i==3) ? LO(boxR[X]) : HI(boxR[X]));
            pi[Y] = ((i==0) || (i==1) ? LO(boxR[Y]) : HI(boxR[Y]));
          }
        else
          { /* Samples are corners and mid-edges. */
            if (i % 2 == 0)
              { /* Corner number {i/2}: */
                pi[X] = ((i==0) || (i==6) ? LO(boxR[X]) : HI(boxR[X]));
                pi[Y] = ((i==0) || (i==2) ? LO(boxR[Y]) : HI(boxR[Y]));
              }
            else
              { /* Mid-edge between corners {i} and {i+1}: */
                pi[X] = ((i==1) || (i==5) ? ctr[X] : (i==7 ? LO(boxR[X]) : HI(boxR[X])));
                pi[Y] = ((i==7) || (i==3) ? ctr[Y] : (i==1 ? LO(boxR[Y]) : HI(boxR[Y])));
              }
          }
       
        /* Pull sample point slightly inside the cell, for discntinuous splines: */
        pi[X] = 0.9999999 * pi[X] + 0.0000001 * ctr[X];
        pi[Y] = 0.9999999 * pi[Y] + 0.0000001 * ctr[Y];
        
        /* Evaluate function at sample point: */
        f->m->eval(f, pi, &fp[i]);
        // DEBUG fprintf (stderr, "    f(%8.6f,%8.6f) = %10.6f\n", pi[X], pi[Y], fp[i]); 
      }
    /* The center value: */
    f->m->eval(f, ctr, &fp[NSAMPLES]);
    // DEBUG fprintf (stderr, "    f(%8.6f,%8.6f) = %10.6f\n", ctr[X], ctr[Y], fp[NSAMPLES]); 

    /* Find min and max among those sampled values: */
    fpMin = fp[0]; fpMax = fp[0];
    for(i = 1; i < NSAMPLES+1; i++)
      { if (fp[i] < fpMin) { fpMin = fp[i]; }
        if (fp[i] > fpMax) { fpMax = fp[i]; }
      }
    
    /* Update global min and max: */
    if (fpMin < (*fObsMin)) { (*fObsMin) = fpMin; }
    if (fpMax > (*fObsMax)) { (*fObsMax) = fpMax; }
    
    /* Compute first and last isolines that enter this rectangle: */
    kpMin = (int)(ceil(fpMin/fStep - 0.4999999));
    kpMax = (int)(floor(fpMax/fStep - 0.5000001));
    
    /* Do any isolines enter the rectangle? */
    if ((kpMin > kPlotMax) || (kpMax < kPlotMin) || (kpMin > kpMax))
      { /* No. */
        if (layer == 0)
          { /* Paint whole rectangle with uniform color. */
            Color *cBand;
            int kBand = kpMin;
            if (kpMin > kPlotMax) { kBand = kPlotMax + 1; }
            if (kpMax < kPlotMin) { kBand = kPlotMin; }
            cBand = &(bandColor[kBand - kPlotMin]);
            SO2DPlot_PaintRectangle(psf, boxR, cBand);
          }
      }
    else
      { /* Paint each triangle: */
        for (i = 0; i < NSAMPLES; i++)
          { int j = (i + 1) % NSAMPLES;
            double *pi = &(p[2*i]);
            double *pj = &(p[2*j]);
            SO2DPlot_ValuesInTriangle
              ( psf, 
                ctr[X], ctr[Y], fp[NSAMPLES],
                pi[X],  pi[Y],  fp[i],
                pj[X],  pj[Y],  fp[j],
                kPlotMin, kPlotMax, fStep,
                layer
              );
          }
      }

    // DEBUG fprintf (stderr, "-InRectangle\n"); 
  }       

void SO2DPlot_ValuesInTriangle
  ( SOPlot_File *psf,   /* Plot file. */
    double Px, double Py, double Pf,
    double Qx, double Qy, double Qf,
    double Rx, double Ry, double Rf,
    Isoline kPlotMin,
    Isoline kPlotMax,
    double fStep,
    int layer
  )
  {
    double temp;
    Isoline kMin, kMax; /* Isolines that enter the triangle are {kMin..kMax}. */

    // DEBUG fprintf (stderr, "+InTriangle\n"); 

    /* Permute corners so that {Pf <= Qf <= Rf}: */
    if( Pf > Qf )
      { temp = Qx; Qx = Px; Px = temp;      
        temp = Qy; Qy = Py; Py = temp;
        temp = Qf; Qf = Pf; Pf = temp;
      }
      
    if( Qf > Rf )
      { temp = Qx; Qx = Rx; Rx = temp;
        temp = Qy; Qy = Ry; Ry = temp;
        temp = Qf; Qf = Rf; Rf = temp; 
      }
      
    if( Pf > Qf )
      { temp = Qx; Qx = Px; Px = temp;      
        temp = Qy; Qy = Py; Py = temp;
        temp = Qf; Qf = Pf; Pf = temp;
      }
    affirm( Pf <= Qf, "Pf > Qf" );
    affirm( Qf <= Rf, "Qf > Rf" ); 

    /* Compute first and last isolines that enter this rectangle: */
    kMin = (int)(ceil(Pf/fStep - 0.4999999));
    kMax = (int)(floor(Rf/fStep - 0.5000001));
    
    // DEBUG fprintf(stderr, "  Pf = %10.8f  kMin = %d  Rf = %10.8f  kMax = %d\n", Pf, kMin, Rf, kMax); 
    
    /* Do any isolines enter the triangle? */
    if ((kMin > kPlotMax) || (kMax < kPlotMin) || (kMin > kMax))
      { /* No. */
        // DEBUG fprintf(stderr, "  kPlotMin = %d  kPlotMax = %d\n", kPlotMin, kPlotMax);
        if (layer == 0)
          { /* Paint whole triangle with color: */
            Color *cBand;
            int kBand = kMin;
            if (kMin > kPlotMax) { kBand = kPlotMax + 1; }
            if (kMax < kPlotMin) { kBand = kPlotMin; }
            cBand = &(bandColor[kBand - kPlotMin]);
            SO2DPlot_PaintTriangle(psf, Px, Py, Qx, Qy, Rx, Ry, cBand);
          }
      }
    else
      { /* Break the triangle {P,Q,R} into slices at isoline levels. */
        double Lf = Pf - 1.0;       /* Function value of previous isoline. */
        double LUx = Px, LUy = Py;  /* Lower point on seg {P-R}. */
        double LVx = Px, LVy = Py;  /* Lower point on seg {P-Q} or {Q-R}. */
        bool_t Lsame = TRUE;          /* Do {LU} and {LV} coincide? */
        Isoline k;

        /* Restrict to specified isoline range: */
        if (kMin < kPlotMin) { kMin = kPlotMin; }
        if (kMax > kPlotMax) { kMax = kPlotMax; }

        for (k = kMin; k <= kMax + 1; k++)
          { /* Function value at higher isoline: */
            double Hf = (k > kMax ? Rf + 1.0 : (k + 0.5)*fStep);
            double HUx, HUy;  /* Higher point on seg {P-R}. */
            double HVx, HVy;  /* Higher point on seg {P-Q} or {Q-R}. */
            bool_t Hsame;       /* Do {HU} and {HV} coincide? */

            /* Slice index {k} is defined by: two lower points {LU =
              (LUx,LUy)} and {LV = (LVx,LVy)}; two higher points {HU =
              (HUx,HUy)} and {HV = (HVx,HVy)}; and possibly the point
              {Q}, if {Lf<Qf<Hf}. The convex region delimited by
              those points lies between isolines {k-1} and {k}, where
              the function value is {Lf} and {Hf}, respectively.
            */

            affirm(Hf >= Pf, "Hf < Pf");

            // DEBUG fprintf(stderr, "    k = %d  Hf = %10.6f\n", k, Hf);

            /* Compute points {HU,HV} on upper isoline: */
            if (Hf >= Rf)
              { /* Last band: */
                HUx = HVx = Rx; HUy = HVy = Ry;
                Hsame = TRUE;
              }
            else
              { double sU = (Hf - Pf)/(Rf - Pf), tU = 1.0 - sU;
                HUx = tU * Px + sU * Rx; 
                HUy = tU * Py + sU * Ry;
                if (Hf <= Qf)
                  { double sV = (Hf - Pf)/(Qf - Pf), tV = 1.0 - sV;
                    HVx = tV * Px + sV * Qx;
                    HVy = tV * Py + sV * Qy;
                  }
                else
                  { double sV = (Hf - Qf)/(Rf - Qf), tV = 1.0 - sV;
                    HVx = tV * Qx + sV * Rx;
                    HVy = tV * Qy + sV * Ry;
                  }
                Hsame = FALSE;
              }

            if (layer == 1)
              { /* Draw the isoline at level {k-1}, if it is non-empty. */
                if ((k > kMin) && (Lf > Pf))
                  { Color *c = (k > 0 ? &(POS_DRAW_COLOR) : &(NEG_DRAW_COLOR));
                    SO2DPlot_DrawSegment(psf, LUx, LUy, LVx, LVy, c);
                  }
              }
            else
              { /* Split the slice into trianglets, and paint them: */
                Color *c = &(bandColor[k - kPlotMin]);
                if ((Lf >= Qf) || (Hf <= Qf))
                  { /* Region is the trapezoid {LU,LV,HV,HU}. */
                    if (Lsame && Hsame)
                      { /* Nothing to do */ }
                    else if (Lsame)
                      { SO2DPlot_PaintTriangle(psf, LUx, LUy, HVx, HVy, HUx, HUy, c); }
                    else if (Hsame)
                      { SO2DPlot_PaintTriangle(psf, LUx, LUy, LVx, LVy, HVx, HVy, c); }
                    else
                      { SO2DPlot_PaintTriangle(psf, LUx, LUy, HVx, HVy, HUx, HUy, c);
                        SO2DPlot_PaintTriangle(psf, HVx, HVy, LUx, LUy, LVx, LVy, c);
                      }
                  }
                else
                  { /* Region is the convex pentagon {LU,LV,Q,HV,HU}. */
                    if (Lsame && Hsame)
                      { SO2DPlot_PaintTriangle(psf, LUx, LUy, Qx,  Qy,  HUx, HUy, c); }
                    else if (Lsame)
                      { SO2DPlot_PaintTriangle(psf, LUx, LUy, Qx,  Qy,  HVx, HVy, c);
                        SO2DPlot_PaintTriangle(psf, LUx, LUy, HVx, HVy, HUx, HUy, c);
                      }
                    else if (Hsame)
                      { SO2DPlot_PaintTriangle(psf, LUx, LUy, LVx, LVy, HVx, HVy, c);
                        SO2DPlot_PaintTriangle(psf, LVx, LVy, Qx,  Qy,  HVx, HVy, c);
                      }
                    else
                      { SO2DPlot_PaintTriangle(psf, LVx, LVy, Qx,  Qy,  LUx, LUy, c);
                        SO2DPlot_PaintTriangle(psf, LUx, LUy, Qx,  Qy,  HUx, HUy, c);
                        SO2DPlot_PaintTriangle(psf, HUx, HUy, Qx,  Qy,  HVx, HVy, c);
                      }
                  }
              }

            LUx = HUx; LUy = HUy;
            LVx = HVx; LVy = HVy;
            Lf = Hf; Lsame = Hsame;
          }
      }

    // ps_set_pen(psf, 0.0, 0.0, 0.0, 0.10, 0.0,0.0);
    // ps_draw_segment(psf, Px, Py*YSCALE, Qx, Qy*YSCALE);
    // ps_draw_segment(psf, Qx, Qy*YSCALE, Rx, Ry*YSCALE);
    // ps_draw_segment(psf, Rx, Ry*YSCALE, Px, Py*YSCALE);

    // DEBUG fprintf (stderr, "-InTriangle\n"); 

  }

void SO2DPlot_PaintRectangle
  ( SOPlot_File *psf,
    interval_t boxR[2],
    Color *c
  )
  {
    double R = c->c[0], G = c->c[1], B = c->c[2];
    if ((R >= 0) && (G >= 0) && (B >= 0))
      { double Lx = LO(boxR[X]), Ly = LO(boxR[Y]) * YSCALE;
        double Hx = HI(boxR[X]), Hy = HI(boxR[Y]) * YSCALE;
        ps_fill_rectangle(psf, Lx, Hx, Ly, Hy, R,G,B);
      }
  }
              
void SO2DPlot_PaintTriangle
  ( SOPlot_File *psf,
    double Ax, double Ay, 
    double Bx, double By,
    double Cx, double Cy,
    Color *c
  )
  {
    double R = c->c[0], G = c->c[1], B = c->c[2];
    if ((R >= 0) && (G >= 0) && (B >= 0))
      { Ay *= YSCALE;
        By *= YSCALE;
        Cy *= YSCALE;
        ps_fill_triangle(psf, Ax, Ay, Bx, By, Cx, Cy,  R,G,B);
      }
  }
              
void SO2DPlot_DrawSegment 
  ( SOPlot_File *psf,
    double Ax, double Ay, 
    double Bx, double By,
    Color *c
  )
  {
    double R = c->c[0], G = c->c[1], B = c->c[2];
    if ((R >= 0) && (G >= 0) && (B >= 0))
      { ps_set_pen(psf, R, G, B, isolineWidth, 0.0,0.0);
        ps_draw_segment(psf, Ax, Ay*YSCALE, Bx, By*YSCALE);
      }
  }

void SO2DPlot_SetColors
  ( Isoline kPlotMin, /* Index of first isoline to be drawn. */
    Isoline kPlotMax, /* Index of last isoline to be drawn. */
    double fStep      /* Value step between isolines. */
  )
  { 
    int nBands = (kPlotMax + 1) - kPlotMin + 1;
    int i;
    
    if (bandColor != NULL) { free(bandColor); }
    bandColor = (Color *)malloc(nBands*sizeof(Color)); 
    affirm(bandColor != NULL, "No memory for color table allocation");
 
    for(i = 0; i < nBands; i++)   
      { Isoline k = kPlotMin + i;  /* Index of isoline above band. */
        if (k < 0)
          { /* Values are negative. */
            double s = ((double)k)/((double)kPlotMin);
            bandColor[i] = SOPlot_ColorScale(s, &ZERO_PAINT_COLOR, &NEG_PAINT_COLOR);
          }
        else if (k > 0)
          { /* Values are positive. */
            double s = ((double)k)/((double)kPlotMax);
            bandColor[i] = SOPlot_ColorScale(s, &ZERO_PAINT_COLOR, &POS_PAINT_COLOR);
          }
        else
          { bandColor[i] = ZERO_PAINT_COLOR; }
      } 
  }

void SO2DPlot_Tree(
  SOPlot_File *psf,       /* Plot file. */
    SOGrid_Tree *tree,    /* The grid to plot. */
    interval_t boxR[2],  /* Low corner of region to plot. */
    double lineWidth,     /* Line width for grid cell boundaries. */
    int maxDepth          /* Omit cells below this depth. */
  )
  {
    ps_set_pen(psf, 0.0, 0.0, 0.0, lineWidth, 0.0, 0.0);
    ps_draw_rectangle(psf, LO(boxR[X]), HI(boxR[X]), YSCALE*LO(boxR[Y]), YSCALE*HI(boxR[Y]));
    if (tree == NULL) { return; }
    affirm(tree->d == 2, "wrong tree dimension");
    SO2DPlot_Subtree(psf, tree->root, 1, 0, boxR, maxDepth);
  }
  
void SO2DPlot_Subtree
  ( FILE *psf, 
    dg_tree_node_t *p, 
    dg_cell_index_t k,
    dg_rank_t rank,
    interval_t boxR[2],
    int maxDepth
  )
  {
    interval_t boxC[2];

    /* Ignore cells at levels below {maxDepth}: */
    if (rank >= maxDepth) { return; } 

    /* If cell is childless, there is nothing to do: */
    if ((p == NULL) || ((LOCH(*p) == NULL) && (HICH(*p) == NULL))) { return; }

    dg_cell_box_root_relative(2, k, boxC);

    /* Skip plot if cell lies outside the region of interest: */
    if ((LO(boxC[X]) >= HI(boxR[X])) || (HI(boxC[X]) <= LO(boxR[X]))) { return; }
    if ((LO(boxC[Y]) >= HI(boxR[Y])) || (HI(boxC[Y]) <= LO(boxR[Y]))) { return; }

    // ps_fill_dot(psf, (LO(boxC[X])+HI(boxC[X]))/2, (LO(boxC[Y])+HI(boxC[Y]))/2, 0.5, 0, 0, 0);

    /* Draw the interior of daughter cells first: */
    SO2DPlot_Subtree(psf, LOCH(*p), 2*k, rank + 1,   boxR, maxDepth);
    SO2DPlot_Subtree(psf, HICH(*p), 2*k+1, rank + 1, boxR, maxDepth);

    /* Now draw the boundary between the two daughter cells: */
    if ((rank % 2) == 0)
      { /* Longest axis is {X}: */
        double xMid = (LO(boxC[X])+HI(boxC[X]))/2;
        double yMin = YSCALE * (LO(boxC[Y]) < LO(boxR[Y]) ? LO(boxR[Y]) : LO(boxC[Y]));
        double yMax = YSCALE * (HI(boxC[Y]) > HI(boxR[Y]) ? HI(boxR[Y]) : HI(boxC[Y]));
        ps_draw_segment(psf, xMid, yMin, xMid, yMax);
      }
    else
      { /* Longest axis is {Y}: */
        double yMid = YSCALE * (LO(boxC[Y])+HI(boxC[Y]))/2;
        double xMin = (LO(boxC[X]) < LO(boxR[X]) ? LO(boxR[X]) : LO(boxC[X]));
        double xMax = (HI(boxC[X]) > HI(boxR[X]) ? HI(boxR[X]) : HI(boxC[X]));
        ps_draw_segment(psf, xMin, yMid, xMax, yMid);
      }
  }

SOPlot_File *SO2DPlot_FunctionFigure
  ( SOPlot_File *psf,          /* Plot file (or NULL). */
    SOFunction *f,             /* Function to plot. */
    SOGrid_Tree *tree,         /* Optional grid to plot. */
    interval_t boxR[2],       /* Low corner of region to plot. */
    double fPlotMin,           /* Nominal minimum {f} value, for color scale. */
    double fPlotMax,           /* Nominal maximum {f} value, for color scale. */
    double fStep,              /* Value step for isolines and color bands. */
    double isolineWidth,       /* Line width for isoline plotting. */
    double gridWidth,          /* Line width for grid drawing. */
    int meshDepth,             /* Depth of bisection recursion for plotting. */
    bool_t bands,                /* TRUE plots color bands. */ 
    bool_t isolines,             /* TRUE draws isolines. */ 
    bool_t grid,                 /* TRUE plots the tree */
    SOPlot_PageState *pgs,     /* Page layout and state for {psf}. */
    char *docName,             /* Document name (minus extension). */
    char *figName,             /* Figure name (minus extension). */
    double *fObsMin,           /* (IN/OUT) Min value seen during plot. */ 
    double *fObsMax            /* (IN/OUT) Max value seen during plot. */ 
  )
  { 
    int maxDepthIsoline = meshDepth;
    int minDepthIsoline = meshDepth - 6;
    int extraDepthIsoline = 2;
    int maxDepthTree = meshDepth - 2;
    psf = SOPlot_BeginFigure(
      psf, pgs, docName, figName, 
      LO(boxR[X]), HI(boxR[X]), 
      YSCALE*LO(boxR[Y]), YSCALE*HI(boxR[Y])
    );
    if (isolines || bands) 
      { SO2DPlot_Values
          ( psf, f, tree, boxR, 
            fPlotMin, fPlotMax, fStep, isolineWidth,
            minDepthIsoline, extraDepthIsoline, maxDepthIsoline,
            bands, isolines,
            fObsMin, fObsMax
          );
      }
    if (grid && (tree != NULL))
      { SO2DPlot_Tree(psf, tree, boxR, gridWidth, maxDepthTree); }
    return psf;
  }

/* Limits to avoid excessive plotting for bad choices of {fStep}: */
/* Should make them consistent with internal limits of {SO2DPlot_Values}. */
#define DefaultIsolines 10
#define MaxIsolinesInRange (3*DefaultIsolines/2)

void SO2DPlot_FixRange
  ( double *fRange, 
    double *fStep,
    double alpha, 
    SOFunction *f, 
    int verbose
  )
  { /* If no place to return result, do nothing: */
    if (fRange == NULL) { return; }
    mumble(1, "checking nominal function range...\n");
    /* If given range is invalid, ignore it: */
    if ((*fRange) <= 0.0) 
      { mumble(1, "  given range is invalid, ignoring it.\n");
        alpha = 1.0;
      }
    /* If no adjustment allowed, do nothing: */
    if (alpha != 0.0)
      { double newRange;
        { /* Estimate the maximum function value: */
          double estRange = SO2DPlot_EstRange(f, verbose);
          if (alpha == 1.0)
            { newRange = estRange; }
          else
            { double oldRange = *fRange;
              mumble(1, "  mixing with given range (alpha = %8.6f)...\n", alpha);
              newRange = alpha * estRange + (1.0 - alpha)*oldRange;
            }
        }
        affirm(newRange > 0.0, "empty range - something failed");
        mumble(1, "  rounding computed range (%10.4e) to a nice value...\n", newRange);
        { double step = SO2DPlot_RoundToNice(newRange/DefaultIsolines);
          newRange = DefaultIsolines * step;
          *fRange = newRange;
          mumble(1, "  new range is [ %10.4e _ %10.4e ].\n", -(*fRange), (*fRange));
        }
      }
    /* Fix {*fStep} if needed: */
    if (fStep == NULL) { return; }
    mumble(1, "checking isoline spacing...\n");
    if ((*fStep) <= (*fRange)/MaxIsolinesInRange)
      { (*fStep) = SO2DPlot_RoundToNice((*fRange)/DefaultIsolines);
        mumble(1, "  given spacing is too small, fixed to %10.4e.\n", *fStep);
      }
  }

double SO2DPlot_RoundToNice(double x)
  { double z = 1.0;
    /* Allow for rounding errors in the algorithm below: */
    x = 0.9999999 * x;
    /* If is too small, don't bother: */
    if (x < 1.0e-300) { return 1.0e-300; }
    /* Finds the smallest power of 10 that is no less than {x}: */
    while (z/10.0 > x) { z /= 10.0; }
    while (z <= x ) { z *= 10.0; }
    /* Finds a nice multiple of that power of 10: */
    if (x < 0.25 * z)
      { return 0.25 * z; }
    else if (x < 0.5 * z)
      { return 0.5 * z; }
    else
      { return z; }
  }

double SO2DPlot_EstRange(SOFunction *f, int verbose)
  { 
    affirm(FALSE, "not implemented yet");
    return 0.0;
    
    // SOSpline *fpw;
    // 
    // auto double func(DomPoint *p);
    // 
    // double func(DomPoint *p)
    //   { return f->m->eval(f, p); }
    //   
    // if ((fpw = SOSpline_Cast(f)) != NULL)
    //   { mumble(1, "  estimating max function value (grid sampling)...\n");
    //     return SORange_OnSOGrid(func, fpw->d->tree, 40);
    //   }
    // else
    //   { mumble(1, "  estimating max function value (root cell sampling)...\n");
    //     return SORange_OnRootCell(func, 40)
    //   }
  }
