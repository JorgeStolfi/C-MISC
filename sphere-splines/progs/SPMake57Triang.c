/* Creates a triangulation where vertices have degree 5 or 7. */
/* Last edited on 2008-05-24 12:26:16 by stolfi */

#define PROG_NAME "SPMake57Triang"

#include <SPOptions.h>
#include <r3.h>
#include <SPTriang.h>
#include <SPTriangExtra.h>
#include <SPDelaunay.h>
#include <SPOptimize.h>
#include <SPQuad.h>
#include <SPBasic.h>
#include <vec.h>
#include <affirm.h>
#include <nat.h>
#include <bool.h>
#include <stdio.h>
#include <math.h>
#include <values.h>
#include <stdlib.h>
#include <string.h>

#define PROG_USAGE \
  PROG_NAME " \\\n" \
  "  -triName NAME \\\n" \
  "  -depth NUM \\\n" \
  "  -outName NAME \\\n" \
  "  [ -shape SCALE B50 B51 B52  B70 B71 B72 ]... \\\n" \
  "  [ -swapEdges ] \\\n" \
  "  [ -optimize NITER [ -minstep NUM] [-maxStep NUM ]\\\n" \
  "     [ -length WGHT ] [ -area WGHT ]\\\n" \
  "     [ -angle WGHT ] [ -xratio WGHT ] \\\n" \
  "     [ -average | -maximum ] \\\n" \
  "  ] \\\n" \
  "  [ -verbose ] [ -plot ]\n"
    
/* 
  If {swapEdges} is true, the final triangulation (but not the
  previous ones) has its topology fixed by edge swaps so that it is a
  Delaunay diagram of the computed sites. */

typedef struct Weights  /* Weights for optimization criterion. */
  { double area;      /* Weight for face area quality metric. */
    double length;    /* Weight for edge length quality metric. */
    double angle;     /* Weight for face angle quality metric. */
    double xratio;    /* Weight for pair-angle quality metric. */
  } Weights;
  
typedef struct Options /* Parsed command line options */
  { char *triName;     /* Input triangulation (usually {icosa}). */
    int depth;         /* Number of subdivision steps. */
    char *outName;     /* Output triangulation file (minus ".tri" extension). */
    R3Point_vec_t b5;  /* Barycentric coords of new 5-fold vertices, per scale. */
    R3Point_vec_t b7;  /* Barycentric coords of new 7-fold vertices, per scale. */
    /* Optimization parameters: */
    bool_t swapEdges;  /* TRUE to force the triangulation to be Delaunay. */
    int maxIter;       /* Number of optimization iterations (0 if no opt). */
    double minStep;    /* Minimum optimization step. */
    double maxStep;    /* Maximum optimization step. */
    Weights w;         /* Optimization criterion weights. */
    bool_t average;    /* TRUE combines measures by averaging, FALSE by extremal. */
    bool_t plot;       /* TRUE to generate plot file of optimizaton progress. */
    bool_t verbose;    /* TRUE to print optimizaton trace. */
  } Options;

Options GetOptions(int argn, char **argc);

double ShapeDist
  ( R3Point_vec_t b5a, R3Point_vec_t b7a, 
    R3Point_vec_t b5b, R3Point_vec_t b7b
  );

Weights RawQualities(Triangulation *tri, Weights *w, bool_t average, bool_t print);
double Quality(Triangulation *tri, Weights *w, bool_t average);
double CombVar(double_vec_t x);
double CombSpr(double_vec_t x);
double CombAvg(double_vec_t x);
double CombMax(double_vec_t x);

void PrintMetricSummary(char *name, double_vec_t x);

void PrintOptStatus
  ( Triangulation *tri, 
    R3Point_vec_t b5best, 
    R3Point_vec_t b7best, 
    double qbest, 
    bool_t average
  );
  
Triangulation *ReadTriangulation(char *name);
  /* Reads a triangulation from file {name} plus extension ".tri". */

void WriteTriangulation(Triangulation *tri, char *name);
  /* Writes triangulation {tri} to the file {name} plus extension ".bas". */
  
void OptimizeTriang
  ( Triangulation **tri,
    int depth, 
    R3Point_vec_t b5, R3Point_vec_t b7,
    int niter, 
    double minStep, double maxStep,
    Weights *w, bool_t average,
    char *plotName,
    bool_t verbose
  );
  /* Tries to adjust the shape parameters {b5} and {b7} so as to
    improve the quality measure of the triangulation {tri[depth]}
    resulting from {depth} refinement steps from {tri[0]}, whose
    results are in {tri[1..depth]}. The quality is defined by
    {Quality(tri[depth], w, average)}. A plot of the progress is
    written to file {plotName} with extension ".plt". */

#define MAXDEPTH 6
#define MAXITER 50000
#define DEFAULTB5 (r3_t){{0.3652, 0.4581, 0.1767}}
#define DEFAULTB7 (r3_t){{0.6661, 0.2761, 0.0578}}

int main(int argn, char **argc)
  { Options o = GetOptions(argn, argc);
    Triangulation *tri[MAXDEPTH];
    int k;
    
    /* Createinitial triangulations: */
    tri[0] = ReadTriangulation(o.triName);
    for (k = 0; k < o.depth; k++)
      { tri[k+1] = SPTriang_Refine57Triang(tri[k], &(o.b5.e[k]), &(o.b7.e[k])); }
   
    /* Initialize the {b5} and {b5} variables for possible optimization: */
    R3Point_vec_t b5opt = R3Point_vec_new(o.depth);
    R3Point_vec_t b7opt = R3Point_vec_new(o.depth);
    for (k = 0; k < o.depth; k++)
      { b5opt.e[k] = o.b5.e[k]; 
        b7opt.e[k] = o.b7.e[k]; 
      }

    if (o.maxIter > 0)
      { /* Optimize triangulation: */
        fprintf(stderr, "optimizing (%d iterations)...\n", o.maxIter);
        OptimizeTriang
          ( tri, 
            o.depth, 
            b5opt, b7opt, 
            o.maxIter, o.minStep, o.maxStep,
            &(o.w), o.average, 
            (o.plot ? o.outName : NULL),
            o.verbose
          );
      }
    
    if (o.swapEdges) 
      { /* Fix edges of final triangulation, and evaluate again: */
        fprintf(stderr, "fixing edges of final triangulation...\n");
        affirm(SPDelaunay_FixDelaunay(tri[o.depth]), "FixDelaunay failed");
        { double qTot = Quality(tri[o.depth], &(o.w), o.average);
          PrintOptStatus(tri[o.depth], b5opt, b7opt, qTot, o.average);
        }
      }

    /* Write final triangulation: */
    WriteTriangulation(tri[o.depth], o.outName);
    fprintf(stderr, "done.\n");
    return 0;
  }
  
void OptimizeTriang
  ( Triangulation **tri,
    int depth, 
    R3Point_vec_t b5, R3Point_vec_t b7,
    int niter, 
    double minStep, double maxStep,
    Weights *w, bool_t average,
    char *plotName,
    bool_t verbose
  )
  { 
    auto double eval(double_vec_t x);
      /* Argument for SPOtimize_Method1. Unpacks the optimization
        argument vector {x} into the shape parameters {b5[0..depth-1]}
        and {b7[0..depth-1]}. Computes the triangulations for those
        shape parameters and returns the quality of the last stage. */
        
    auto void unpack_b5_b7(double_vec_t x, R3Point_vec_t b5, R3Point_vec_t b7);
      /* Extracts barycentric coordinates {b5[0..depth-1]} and
        {b7[0..depth-1]} from {x[0..4*depth-1]}. Specifically, {b5[k]}
        is recovered from {x[4*k]} and {x[4*k+1]}, and {b7[k]} is
        recovered from {x[4*k+2]} and {x[4*k+3]}. See {extract_bar}
        below. */
    
    auto void extract_bar(double x0, double x1, r3_t *b);
      /* Converts packed parameters {x0} and {x1} into barycentric
        coordinates {b}, ensuring these are non-negative and add to 1.0. */
    
    void unpack_b5_b7(double_vec_t x, R3Point_vec_t b5, R3Point_vec_t b7)
      { int k;
        for (k = 0; k < depth; k++)
          { r3_t *b5k = &(b5.e[k]);
            r3_t *b7k = &(b7.e[k]);
            extract_bar(x.e[4*k+0], x.e[4*k+1], b5k);
            extract_bar(x.e[4*k+2], x.e[4*k+3], b7k);
          }
      }
        
    void extract_bar(double x0, double x1, r3_t *b)
      { double b0 = fabs(x0);
        double b1 = fabs(x1);
        double eps = 0.000001;
        double b2, s;
        while (b0 > 1.0) { b0 = fabs(2.0 - b0); }
        while (b1 > 1.0) { b1 = fabs(2.0 - b1); }
        if (b0 + b1 > 1.0) { b0 = 1.0 - b0; b1 = 1.0 - b1; }
        b2 = fabs(1.0 - b0 - b1);
        s = b0 + b1 + b2 + 3.0*eps;
        (*b) = (r3_t){{(b0+eps)/s, (b1+eps)/s, (b2+eps)/s}};
      }
    
    double eval(double_vec_t x)
      { int k;
        /* Extract parameters {b5,b7} for all scales: */
        unpack_b5_b7(x, b5, b7);
        /* Recompute triangulations: */
        for (k = 0; k < depth; k++)
          { r3_t *b5k = &(b5.e[k]);
            r3_t *b7k = &(b7.e[k]);
            if (verbose) 
              { r3_gen_print(stderr, b5k, "%16.12f", "b5    = [", ",", "]\n");
                r3_gen_print(stderr, b7k, "%16.12f", "b7    = [", ",", "]\n");
              }
            SPTriang_Recompute57Triang(tri[k+1], tri[k], b5k, b7k);
          }
        
        /* Evaluate final triangulation: */
        { double qTot = Quality(tri[depth], w, average);
          if (verbose) { fprintf(stderr, "[ eval: %16.8e ]\n", qTot); }
          return qTot;
        }
      }
    
    double_vec_t x = double_vec_new(4*depth); /* Packed {b5,b7}; */
    double qbest;
    int k;
    
    /* Pack the initial guess into the {x} vector: */
    for (k = 0; k < depth; k++)
      { x.e[4*k+0] = b5.e[k].c[0]; 
        x.e[4*k+1] = b5.e[k].c[1]; 
        x.e[4*k+2] = b7.e[k].c[0]; 
        x.e[4*k+3] = b7.e[k].c[1]; 
      }
    /* Evaluate the initial guess: */
    qbest = eval(x); 
    PrintOptStatus(tri[depth], b5, b7, qbest, average);

    /* Optimize the quality by varying the shape parameters: */
    SPOptimize_Method1(eval, x, &qbest, niter, minStep, maxStep, plotName, verbose);

    /* Unpack best solution found: */
    { double qbx = eval(x); affirm(qbx == qbest, "opt fxbest buggy"); }
    fprintf(stderr, "\n");
    PrintOptStatus(tri[depth], b5, b7, qbest, average);
    free(x.e);
  }
  
void PrintOptStatus
  ( Triangulation *tri, 
    R3Point_vec_t b5best, 
    R3Point_vec_t b7best, 
    double qbest, 
    bool_t average
  )
  { Weights unitwt = (Weights){1.0, 1.0, 1.0, 1.0};
    int depth = b5best.ne;
    int k;
    for (k = 0; k < depth; k++)
      { fprintf(stderr, "  -shape %02d ", k);
        r3_gen_print(stderr, &(b5best.e[k]), "%6.4f", " ", " ", " ");
        r3_gen_print(stderr, &(b7best.e[k]), "%6.4f", " ", " ", " ");
        if (k < depth-1) { fprintf(stderr, "\\"); }
        fprintf(stderr, "\n");
      }
    (void) RawQualities(tri, &unitwt, average, TRUE); 
    fprintf(stderr, " qbest = %16.8e\n", qbest);
  }

double Quality(Triangulation *tri, Weights *w, bool_t average)
  { Weights q = RawQualities(tri, w, average, FALSE);
    q.area *= w->area;
    q.length *= w->length;
    q.angle *= w->angle;
    q.xratio *= w->xratio;
    return q.length + q.area + q.angle + q.xratio;
  }

Weights RawQualities(Triangulation *tri, Weights *w, bool_t average, bool_t print)
  { Weights q = (Weights){0,0,0,0};
    double NF = (double)tri->side.ne;
    if (print) { fprintf(stderr, "\n"); }
    if (w->area > 0.0)
      { double_vec_t x = SPTriang_FaceAreas(tri);
        double idealArea = FOURPI/NF;
        double tiny = 1.0e-5*idealArea, huge = 1.0e+5*idealArea;
        if (print) { PrintMetricSummary("area", x); }
        SPTriang_CookValues(x, tiny, huge);
        q.area = (average ? CombVar(x) : CombSpr(x));
        if (print) { fprintf(stderr, "  comb = %12.8f\n", q.area); }
        free(x.e);
      }
    if (w->length > 0.0)
      { double_vec_t x = SPTriang_EdgeLengths(tri);
        double idealLength = 4.0*sqrt(PI/NF/SQRT3);
        double tiny = 1.0e-5*idealLength, huge = 1.0e+5*idealLength;
        if (print) { PrintMetricSummary("length", x); }
        SPTriang_CookValues(x, tiny, huge);
        q.length = (average ? CombVar(x) : CombSpr(x));
        if (print) { fprintf(stderr, "  comb = %12.8f\n", q.length); }
        free(x.e);
      }
    if (w->angle > 0.0)
      { double_vec_t x = SPTriang_RelAngles(tri);
        double tiny = 1.0e-5, huge = 1.0e+5;
        if (print) { PrintMetricSummary("angle", x); }
        SPTriang_CookValues(x, tiny, huge);
        q.angle = (average ? CombVar(x) : CombSpr(x));
        if (print) { fprintf(stderr, "  comb = %12.8f\n", q.angle); }
        free(x.e);
      }
    if (w->xratio > 0.0)
      { double_vec_t x = SPTriang_CrossRatios(tri);
        double tiny = 1.0e-4, huge = 0.25e+4;
        if (print) { PrintMetricSummary("xratio", x); }
        SPTriang_CookValues(x, tiny, huge);
        { int i; for (i = 0; i < x.ne; i++) { x.e[i] = 1.0/x.e[i]; } }
        q.xratio = (average ? CombAvg(x) : CombMax(x));
        if (print) { fprintf(stderr, "  comb = %12.8f\n", q.xratio); }
        free(x.e);
      }
    if (print) { fprintf(stderr, "\n"); }
    return q;
  }
  
void PrintMetricSummary(char *name, double_vec_t x)
  { double min = INFINITY, max = -INFINITY;
    double sum = 0.0;
    int i;
    for (i = 0; i < x.ne; i++)
      { double xi = x.e[i]; 
        if (xi < min) { min = xi; }
        if (xi > max) { max = xi; }
        sum += xi;
      }
    fprintf(stderr, "  %-10s n = %6d min = %12.8f  max = %12.8f  avg = %12.8f", 
      name, x.ne, min, max, sum/((double)x.ne)
    );
  }
    
double CombVar(double_vec_t x)
  { return SPTriang_RatioVar(x); }

double CombSpr(double_vec_t x)
  { return SPTriang_RatioSpr(x); }
    
double CombAvg(double_vec_t x)
  { return SPTriang_GeomAvg(x); }
    
double CombMax(double_vec_t x)
  { return SPTriang_Max(x); }

Triangulation *ReadTriangulation(char *name)
  { FILE *rd = open_read(txtcat(name, ".tri"), TRUE);
    Triangulation *tri = SPTriang_Read(rd, 1);
    fclose(rd);
    return tri;
  }
  
void WriteTriangulation(Triangulation *tri, char *name)
  { FILE *wr = open_write(txtcat(name, ".tri"), TRUE);
    SPTriang_Write(wr, tri);
    fclose(wr);
  }


Options GetOptions(int argn, char **argc)
  {
    Options o;
    SPOptions_Parser_t *pp = SPOptions_NewParser(stderr, argn, argc);
          
    SPOptions_SetUsage(pp, PROG_USAGE);
          
    SPOptions_GetKeyword(pp, "-triName");                               
    o.triName = SPOptions_GetNext(pp);  

    SPOptions_GetKeyword(pp, "-depth");                               
    o.depth = SPOptions_GetNextInt(pp, 0, MAXDEPTH);  
       
    SPOptions_GetKeyword(pp, "-outName");                               
    o.outName = SPOptions_GetNext(pp);  

    o.swapEdges = SPOptions_TestKeyword(pp, "-swapEdges");
    
    o.maxIter = 0;
    o.w = (Weights){0,0,0,0};
    o.average = TRUE;
    if (SPOptions_TestKeyword(pp, "-optimize"))  
      { o.maxIter = SPOptions_GetNextInt(pp, 0, MAXITER);

        if (SPOptions_TestKeyword(pp, "-minStep"))  
          { o.minStep = SPOptions_GetNextDouble(pp, 0, 1.0); }
        else
          { o.minStep = 0.00005; }

        if (SPOptions_TestKeyword(pp, "-maxStep"))  
          { o.maxStep = SPOptions_GetNextDouble(pp, 0, 1.0); }
        else
          { o.maxStep = 0.02; }

        if (SPOptions_TestKeyword(pp, "-length"))  
          { o.w.length = SPOptions_GetNextDouble(pp, 0, 1.0); }
        if (SPOptions_TestKeyword(pp, "-area"))  
          { o.w.area = SPOptions_GetNextDouble(pp, 0, 1.0); }
        if (SPOptions_TestKeyword(pp, "-angle"))  
          { o.w.angle = SPOptions_GetNextDouble(pp, 0, 1.0); }
        if (SPOptions_TestKeyword(pp, "-xratio"))  
          { o.w.xratio = SPOptions_GetNextDouble(pp, 0, 1.0); }
        if (SPOptions_TestKeyword(pp, "-average"))
          { o.average = TRUE; }
        else if (SPOptions_TestKeyword(pp, "-maximum"))
          { o.average = FALSE; }
        if
          ( (o.w.area == 0.0) && (o.w.length == 0.0) && 
            (o.w.angle == 0.0) && (o.w.xratio == 0.0)
          )
          { SPOptions_Error(pp, "must specify some optimization goal"); }
      }

    { int k;
      o.b5 = R3Point_vec_new(o.depth);
      o.b7 = R3Point_vec_new(o.depth);
      for (k = 0; k < o.depth; k++)
        { o.b5.e[k] = DEFAULTB5;
          o.b7.e[k] = DEFAULTB7;
        }
      while (SPOptions_TestKeyword(pp, "-shape"))  
        { int i;
          k = SPOptions_GetNextInt(pp, 0, o.depth-1);
          for (i = 0; i < 3; i++)
            { o.b5.e[k].c[i] = SPOptions_GetNextDouble(pp, 0, 1.0); }
          for (i = 0; i < 3; i++)
            { o.b7.e[k].c[i] = SPOptions_GetNextDouble(pp, 0, 1.0); }
        }
    }
      
    o.verbose = SPOptions_TestKeyword(pp, "-verbose");
    o.plot = SPOptions_TestKeyword(pp, "-plot");

    SPOptions_Finish(pp);
    return o;
  }

