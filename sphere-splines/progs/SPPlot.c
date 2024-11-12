/* See SPPlot.h */
/* Last edited on 2024-11-11 05:13:54 by stolfi */

/* To get {asprintf}: */
#define _GNU_SOURCE

#include <SPPlot.h>

#include <SPTriang.h>
#include <SPSpline.h>
#include <SPQuad.h>
#include <SPRange.h>
#include <SPIntegral.h>
#include <SPBasic.h>
#include <vec.h>
#include <SPH3.h>
#include <pswr.h>
#include <r2.h>
#include <r3.h>
#include <r4.h>
#include <r4x4.h>
#include <r3x3.h>
#include <rmxn.h>
#include <affirm.h>
#include <nat.h>

#include <math.h>
#include <values.h>
#include <stdio.h>
#include <string.h>

#define Invisible ((Color){{-1,-1,-1}})

typedef r3_t ProjPoint;  /* Persp-mapped point {(xplot,yplot,height)}. */

/* INTERNAL PROTOTYPES */

ProjPoint Project(R3Point *p, SPH3_PMap *map);
  /* Applies the perspective projection {map} to the point {p}.
    If {map} is null, converts to longitude-latitude coordinates. */
  
R3Vector FindOrthoDir(R3Vector *d);
  /* Finds some unit vector {u} orthogonal to {d}. */

SPH3_Point SPPlot_SecondZenith(void);
  /* A second-best choice for the default zenith. */

bool_t SPPlot_IsGoodZenith(SPH3_Point *upp, SPH3_Point *obs);
  /* TRUE iff {upp} is a good zenith reference for the given {obs}. */

/* VISIBILITY 
   
  The procedures {XXXIsOK} return TRUE iff the point, segment or trianggle
  defined by the given points can be plotted.  If {map != NULL},
  they tests whether the given points lie all on the intrresting 
  part of the sphere; which is the part visible from the observer,
  if {backSide = FALSE} or the invisible part, if {backSide = TRUE}.
  
  If {map=NULL} (implying a longitude-latititude plot), the procedures
  test whether the object defined by the givcen points contains a pole
  or crosses the longitude cut line {lon = ±PI}. */

bool_t PointIsOK(S2Point *p, SPH3_PMap *map, bool_t backSide);
bool_t SegmentIsOK(S2Point *p, S2Point *q, SPH3_PMap *map, bool_t backSide);
bool_t TriangleIsOK(S2Point *p, S2Point *q, S2Point *r, SPH3_PMap *map, bool_t backSide);
  
Color InterpolateColor
  ( double f,         /* Function value */
    double fMin,      /* Minimum function value. */
    Color *cMin,      /* Color to use for {fMin}. */
    bool_t clipMin,     /* TRUE omits parts below {fMin}. */
    double fMax,      /* Maximum function VALUE. */
    Color *cMax,      /* Color to use for {fMax}. */
    bool_t clipMax      /* TRUE omits parts above {fMax}. */
  );
  /* Compute the color corresponding to {f}, by linear
    interpolation from {[fMin _ fMax]} to {[cMin _ cMax]}.
    The meaning of {clipMin} and {clipMax} is the same as
    in {SPPlot_PaintValues}. */

Color ClipColor(Color *c);
  /* Clips the color {c} to the unit cube, preserving 
    its brightness and hue. */

void PlotZerosInTriangle
  ( SPPlot_Stream *fps,
    double Px, double Py, double Pf,
    double Qx, double Qy, double Qf,
    double Rx, double Ry, double Rf
  );
  /* Given a triangle {P,Q,R} on the plane,
    and function values at its vertices, interpolates a 
    linear function through that data, and plots the set of 
    points where that function is zero. */
    
typedef void TriangleVisitScalar
  ( S2Point *P, double fP,  /* Vertex 0 (projected) and its function value. */
    S2Point *Q, double fQ,  /* Vertex 1 (projected) and its function value. */
    S2Point *R, double fR   /* Vertex 2 (projected) and its function value. */
  );
  /* A procedure that can process a triangular sphere fragment. */
  
void SPPlot_GenericPlot
  ( SPPlot_Stream *fps,
    Triangulation *tri,       /* Reference triangulation, of NULL */
    int N,                    /* Order of mesh subdivision. */
    ScalarField func,         /* Function to plot. */
    TriangleVisitScalar proc, /* Called for each sphere fragment. */
    double *fMinObs /*I/O*/,  /* Min {func} value seen (caller must initialize). */
    double *fMaxObs /*I/O*/   /* Max {func} value seen (caller must initialize). */
  );
  /* Generates the plotting mesh and calls {proc} on each piece.
    Tries to reuse values of {func} that are shared between adjacent
    triangles of the mesh.  The corners of each piece are in CCW
    order as seen from outside the sphere.
    
    If {tri} is not null, each face of {tri} is covered by a separate
    mesh that lies slightly inside the face; the meshes of adjacent faces
    do NOT share any points. 

    If {map} is NULL, produces a longitude-latitude plot. */
    
double SPPlot_RadiusOfInterestFromValues(S2Point *odr, SPFunction *f, bool_t verbose);
  /* Computes the radius of a cap centered at {odr} (at most one hemisphere)
    that contains most of {f^2(p)}. */

double SPPlot_RadiusOfInterestFromSupport(r3_t *odr, SPH3_Plane *S, bool_t verbose);
  /* Computes the radius of a cap centered at {odr} (at most one
    hemisphere) that contains as much as possible of the positive cap
    of the supporting plane plane {S}. */

SPH3_Point SPPlot_FocusFromRadiusOfInterest(SPH3_Point *obs, double rad);
  /* Computes the planar center {foc} of the spherical cap whose
     spherical center is {dir(obs)} and whose radius is {rad}. */

/* POSTSCRIPT FIGURE STREAMS */

SPPlot_Stream *SPPlot_NewStream
  ( bool_t eps,
    char *name,
    char *paperSize, 
    double figSize,
    bool_t projLonLat,
    int nCap
  )
  { /* Add caption only if there is a user caption, or it is not EPS. */
    /* Select a good figure size: */
    double figWd = figSize*(72.0/25.4);
    double figHt = (projLonLat ? figWd/2 : figWd);
    SPPlot_Stream *fps = pswr_new_stream
      ( /* name */                  name,
        /* file */                  NULL,
        /* eps */                   eps,
        /* docName */               "doc",
        /* paperSize */             paperSize,
        /* landscape */             FALSE,
        /* hPageSize, vPageSize */  figWd + 6.0, figHt + 6.0
      );
    pswr_set_canvas_layout
      ( fps,
        /* hPicSize, vPicSize */     figWd, figHt,
        /* adjustPicSize */          FALSE,
        /* hPicMargin,vPicMargin */  2.0, 2.0,
        /* captionLines */           nCap + (eps ? 0 : 1),  
        /* vCount, hCount */         0, 0  /* Let {pswr} choose it. */
      ); 
    return fps;
  }
  
double SPPlot_DefaultFigSize
  ( bool_t eps, 
    char *paperSize, 
    int nRows, 
    int nCols,
    bool_t projLonLat,
    int captionLines
  )
  { if (nRows <= 0) { nRows = 1; }
    if (nCols <= 0) { nCols = 1; }
    double hsize, vsize;
    if (eps)
      { hsize = vsize = 150.0; }
    else
      { pswr_get_paper_dimensions(paperSize, FALSE, &hsize, &vsize);
        hsize = (hsize/72.0 - 2.0)*25.0; 
        vsize = ((vsize - 10.0*nRows*captionLines)/72.0 - 2.0)*25.0;
        if (hsize < 0.0) { hsize = 10.0*nCols; }
        if (vsize < 0.0) { vsize = 10.0*nRows; }
      }
    hsize /= nCols;
    vsize /= nRows;
    /* This should be correct irrespective of the {projLonLat}: */
    return (hsize < vsize ? hsize : vsize);
  }
  
/* PLOT COMPONENTS */

void SPPlot_Sphere(SPPlot_Stream *fps, SPH3_PMap *map)
  { 
    if (map == NULL)
      { /*Plot the longitude-latitude rectangle: */
        pswr_segment(fps, -1.0, -0.5, +1.0, -0.5);
        pswr_segment(fps, +1.0, -0.5, +1.0, +0.5);
        pswr_segment(fps, +1.0, +0.5, -1.0, +0.5);
        pswr_segment(fps, -1.0, +0.5, -1.0, -0.5);
      }
    else
      { /* Plot the outline of {S^2} in perspective: */
        /* We only handle the case when the sphere projects as a circle. */
        r4x4_t P,Q;
        r3x3_t D;
        pswr_comment(fps, "SPPlot_Sphere");
        /* Compute the matrix {Q} of the 3D quadric that is the image of {S^2}. */
        P = map->inv; 
        { int i; for (i = 0; i < 4; i++) { P.c[i][0] = -P.c[i][0]; } }
        rmxn_mul_tr(4, 4, 4, &(P.c[0][0]), &(map->inv.c[0][0]), &(Q.c[0][0]));
        /* Compute the discriminant quadric {D} for {z} as a function of {(x,y)}. */
        /* The discriminant is positive inside the projection of {S^2}. */
        { int i, j; 
          for (i = 0; i < 3; i++) 
            for (j = 0; j < 3; j++) 
              { D.c[i][j] = 
                  (Q.c[3][i] + Q.c[i][3])*(Q.c[3][j] + Q.c[j][3])
                  - 4.0*Q.c[3][3]*Q.c[i][j];
              }
        }
        /* Now plot the implicit conic {p*D*p^t = 0}: */ 
        { double cA = D.c[1][1];
          double cB = D.c[1][2] + D.c[2][1];
          double cC = D.c[2][2];
          double cD = D.c[1][0] + D.c[0][1];
          double cE = D.c[2][0] + D.c[0][2];
          double cF = D.c[0][0];
          /* See what kind of conic it is: */
          if (cA*cC > 0.0)
            { /* Curve is an ellipse or circle, compute center: */
              double det = cB*cB - 4.0*cA*cC;
              double xctr = (- 2*cC*cD + cB*cE)/det;
              double yctr = (- 2*cA*cE + cB*cD)/det;
              /* Compute equation relative to center: */
              /* double rD = cD - 2*cA*xctr - cB*yctr; */
              /* double rE = cE - 2*cC*yctr - cB*xctr; */
              double rF = cA*xctr*xctr + cB*xctr*yctr + cC*yctr*yctr 
                - cD*xctr - cE*yctr + cF;
              /* Determine type of curve: */
              double m = sqrt(cA*cA + cB*cB + cC*cC);
              if ((cB/m < 1.0e-6) && (fabs(cA - cC)/m < 1.0e-6))
                { /* Curve is a circle: */
                  double r2 = -rF/cA, radius = sqrt(r2 <= 0.0 ? 0.0 : r2);
                  pswr_circle(fps, xctr, yctr, radius,  FALSE, TRUE);
                }
              else
                { /* Curve is an ellipse: */
                  fprintf(stderr, "** sphere projects as ellipse - not drawn");
                }
            }
          else
            { /* Curve is not an ellipse or circle: */
              affirm(FALSE, "sphere has funny projection");
            }
        }
      }
  }

void SPPlot_Dot
  ( SPPlot_Stream *fps, 
    SPH3_PMap *map,  
    S2Point *p, 
    double radius, 
    bool_t fill,
    bool_t draw,
    bool_t backSide     
  )
  { S2Point a;
    r3_dir(p, &a);
    if (PointIsOK(&a, map, backSide))
      { ProjPoint am = Project(&a, map);
        pswr_dot(fps, am.c[0], am.c[1], radius, fill, draw);
      }
  } 

void SPPlot_Segment
  ( SPPlot_Stream *fps, 
    SPH3_PMap *map,  
    S2Point *p, S2Point *q, 
    double step,
    bool_t backSide      
  )
  { 
    /* Plot geodesic arc step by step: */
    /* Number of steps to use: */
    int N = ceil(acos(r3_cos(p,q))/step);
    S2Point a;
    r3_dir(p, &a);
    int i;
    for (i = 1; i <= N; i++)
      {	S2Point b;
	double s = ((double)i)/((double)N);
	r3_mix(1.0 - s, p, s, q, &b);
	r3_dir(&b, &b);
        if (SegmentIsOK(&a, &b, map, backSide))
	  { ProjPoint am = Project(&a, map);
	    ProjPoint bm = Project(&b, map);
	    pswr_segment(fps, am.c[0], am.c[1], bm.c[0], bm.c[1]);
	  }
	a = b;
      }
  } 
    
void SPPlot_GreatCircle
  ( SPPlot_Stream *fps,
    SPH3_PMap *map, 
    r3_t normal,
    double step,          /* Maximum step size (1.0 = sphere radius). */
    bool_t backSide       /* Show back side instead of front side. */
  )
  { 
    pswr_comment(fps, "SPPlot_GreatCircle");
    /* Pick two orthogonal points {v,w} on the circle: */
    r3_t v, w;
    (void)r3_throw_ortho_pair(&normal, &v, &w); 
    /* Now plot the four quarter-circle arcs: */
    int i;
    for (i = 0; i < 4; i++)
      { SPPlot_Segment(fps, map, &v, &w, step, backSide);
        r3_t t = v; v = w; r3_neg(&t, &w);
      }
  }
  
void SPPlot_Triangulation
  ( SPPlot_Stream *fps,
    SPH3_PMap *map, 
    Triangulation *tri, 
    double step,          /* Maximum step size (1.0 = sphere radius). */
    bool_t backSide       /* Show back side instead of front side. */
  )
  { Arc_vec_t arc = tri->arc;
    int i;
    pswr_comment(fps, "SPPlot_Triangulation");
    for (i = 0; i < arc.ne; i += 2)
      { if (arc.e[i] != SPQuad_NullArc)
          { Arc e = arc.e[i];
            S2Point *op = &(Org(e)->pos);
            S2Point *dp = &(Dest(e)->pos); 
            SPPlot_Segment(fps, map, op, dp, step, backSide);
          }
      }
  }
  
void SPPlot_LonLatGrid
  ( SPPlot_Stream *fps,
    SPH3_PMap *map, 
    int gridNLon,            /* Number of grid steps in longitude. */
    int gridNLat,            /* Number of grid steps in latitude. */
    bool_t gridDots,         /* TRUE shows lon-lat grid as face-ctr dots, FALSE as rects. */
    double dotRadius,        /* Dot radius if {gridDots}. */     
    bool_t backSide          /* TRUE plots the back side instead of front side. */ 
  )
  { 
    if ((gridNLon == 0) || (gridNLat == 0)) { return; }
    double dLon = 2*PI/gridNLon;
    double dLat = PI/gridNLat;
    pswr_comment(fps, "SPPlot_LonLatGrid");
    int iLon, iLat;
    for (iLon = 0; iLon < gridNLon; iLon++)
      { for (iLat = 0; iLat < gridNLat; iLat++)
          { if (gridDots)
              { double lon = (iLon + 0.5)*dLon - PI;
                double lat = (iLat + 0.5)*dLat - PI/2;
                double cLat = cos(lat);
                S2Point p = (S2Point){{ cos(lon)*cLat, sin(lon)*cLat, sin(lat) }};
                SPPlot_Dot(fps, map, &p, dotRadius, TRUE, FALSE, backSide);
              }
            else
              { affirm(FALSE, "grid lines not implemented yet"); }
          }
      }
  }
  
void SPPlot_CoordAxes
  ( SPPlot_Stream *fps, 
    SPH3_PMap *map, 
    Sign sense
  )
  { R3Vector x = (R3Vector){{1.0, 0.0, 0.0}};
    R3Vector y = (R3Vector){{0.0, 1.0, 0.0}};
    R3Vector z = (R3Vector){{0.0, 0.0, 1.0}};
    SPPlot_Axis(fps, map, &x, 1.25, sense);
    SPPlot_Axis(fps, map, &y, 1.25, sense);
    SPPlot_Axis(fps, map, &z, 1.25, sense);
  }

void SPPlot_Polyhedron
  ( SPPlot_Stream *fps,
    SPH3_PMap *map,
    Triangulation *tri,
    bool_t backSide
  )
  { 
    auto void PlotEdge(Arc e);
    auto bool_t EdgeIsOK(Arc e);
    auto bool_t FaceIsOK(Arc e); 
    
    bool_t FaceIsOK(Arc e)
      { S2Point *p = &(Org(e)->pos);
        S2Point *q = &(Dest(e)->pos);
        S2Point *r = &(Org(Lprev(e))->pos);
        return TriangleIsOK(p, q, r, map, backSide);
      }
    
    bool_t EdgeIsOK(Arc e)
      { S2Point *p = &(Org(e)->pos);
        S2Point *q = &(Dest(e)->pos);
        return SegmentIsOK(p, q, map, backSide);
      }
    
    void PlotEdge(Arc e)
      { ProjPoint op = Project(&(Org(e)->pos), map); 
        ProjPoint dp = Project(&(Dest(e)->pos), map); 
        pswr_segment(fps, op.c[0], op.c[1], dp.c[0], dp.c[1]);
      }

    int i;
    Arc_vec_t arc = tri->arc;

    if (map == NULL) { return; }
    pswr_comment(fps, "SPPlot_Polyhedron");
    
    /* Plot edges: */
    for (i = 0;  i < arc.ne; i += 2)
      { Arc e = arc.e[i];
        if ((e != SPQuad_NullArc) && (! EdgeIsOK(e)))
          { PlotEdge(e); }
      }
  }

void SPPlot_SupportingPlane
  ( SPPlot_Stream *fps,
    SPH3_PMap *map, 
    SPH3_Plane *supp, 
    double angStep
  )
  { R3Point ctr = (R3Point){{supp->f.c[1], supp->f.c[2], supp->f.c[3]}};
    double rad;
    R3Vector u = (R3Vector){{0.0, 0.0, 0.0}}, v;
    int i;
    
    pswr_comment(fps, "SPPlot_SupportingPlane");
    /* Normalize length of {ctr} so that it is the circle center on the plane: */
    { double r2 = r3_norm_sqr(&ctr);
      double W = supp->f.c[0];
      if ((r2 == 0.0) || (W*W >= r2)) { /* Cap is empty: */ return; }
      r3_scale(-W/r2, &ctr, &ctr);
    }
    /* Compute radius {rad} of flat circle: */
    { double rad2 = (1.0 - r3_norm_sqr(&ctr));
      rad = (rad2 <= 0.0 ? 0.0 : sqrt(rad2));
    }
    
    /* Pick ortho vectors {u,v} orthogonal to {ctr} with length {rad}: */
    { (void)r3_throw_ortho_pair(&ctr, &u, &v);
      /* Normalize vectors {u,v} to length {rad}: */
      r3_scale(rad, &u, &u);
      r3_scale(rad, &v, &v);
    }
    /* Now plot the circle.  Watch out for small or zero {rad}. */
    { R3Point w, p, q;
      int N = ceil(TWOPI/angStep);
      r3_add(&ctr, &u, &p);
      for (i = 1; i <= N; i++)
        { double theta = TWOPI*((double)i)/((double)(N));
          double ct = cos(theta), st = sin(theta);
          r3_mix(ct, &u, st, &v, &w); 
          r3_add(&ctr, &w, &q);
          SPPlot_Segment(fps, map, &p, &q, 1.0, FALSE);
          p = q;
        }
    }
  }

void SPPlot_GenericPlot
  ( SPPlot_Stream *fps,
    Triangulation *tri,       /* Reference triangulation, of NULL */
    int N,                    /* Order of mesh subdivision. */
    ScalarField func,         /* Function to plot. */
    TriangleVisitScalar proc, /* Called for each sphere fragment. */
    double *fMinObs /*I/O*/,  /* Min {func} value seen (caller must initialize). */
    double *fMaxObs /*I/O*/   /* Max {func} value seen (caller must initialize). */
  )
  { 
    auto double callFunc(S2Point *p);
      /* Calls {func(p)} and updates {fMaxObs, fMinObs} */

    auto void processTriangle(S2Point *P, S2Point *Q, S2Point *R);
      /* Subdivides the given triangle into {N^2} trianglets and calls {proc}
        on them. */

    double callFunc(S2Point *p)
      { double f = func(p);
        if (f > (*fMaxObs)) { (*fMaxObs) = f; }
        if (f < (*fMinObs)) { (*fMinObs) = f; }
        return f;
      }

    S2Point_vec_t S = S2Point_vec_new(N + 1); /* Saved points. */
    double_vec_t fS = double_vec_new(N + 1);  /* Their function values. */
    
    void processTriangle(S2Point *P, S2Point *Q, S2Point *R)
      { double fN = (double)N;
        int i;
        for (i = 0; i <= N; i++)
          { int j;
            for (j = 0; j <= N-i; j++) 
              { int k = N - i - j;
                r3_t T = (r3_t)
                  {{(P->c[0]*i + Q->c[0]*j + R->c[0]*k)/fN,
                    (P->c[1]*i + Q->c[1]*j + R->c[1]*k)/fN,
                    (P->c[2]*i + Q->c[2]*j + R->c[2]*k)/fN
                  }};
                double fT;
                r3_dir(&T, &T); 
                fT = callFunc(&T);
                if (i > 0)
                  { /* Plot triangles using points from previous row: */
                    /* Be sure to preserve orientation rel. to {P,Q,R}: */
                    r3_t *V = &(S.e[j]);   double fV = fS.e[j];  
                    r3_t *W = &(S.e[j+1]); double fW = fS.e[j+1];
                    if (j > 0)
                      { r3_t *U = &(S.e[j-1]); double fU = fS.e[j-1];
                        proc(V, fV, U, fU, &T, fT);
                      }
                    proc(W, fW, V, fV, &T, fT);
                  }
                /* Save point and function value: */
                S.e[j] = T; fS.e[j] = fT;
              }
          }
      }

    if (tri != NULL)
      { /* Process each part of the given triangulation */
        int NT = tri->side.ne;
        double epsilon = 0.0001; /* Relative perturbation to avoid edges */
        int i, j;
        /* Plot the function inside each triangle: */
        for (i = 0; i < NT; i++)
          { Arc e = tri->side.e[i];
            S2Point P = Dest(Onext(e))->pos;
            S2Point Q = Org(e)->pos;
            S2Point R = Dest(e)->pos;
            
            /* Perturb points so that they lie slightly inside the triangle: */
            for (j = 0; j < 3; j++)
              { double Bj = (P.c[j] + Q.c[j] + R.c[j])/3.0;
                P.c[j] = (1.0-epsilon)*P.c[j] + epsilon*Bj;
                Q.c[j] = (1.0-epsilon)*Q.c[j] + epsilon*Bj;
                R.c[j] = (1.0-epsilon)*R.c[j] + epsilon*Bj;
              }
            r3_dir(&P, &P);
            r3_dir(&Q, &Q); 
            r3_dir(&R, &R);

            processTriangle(&P, &Q, &R);
          }
      }
    else
      { /* Subdivide the regular octahedral triangulation. */
        int x, y, z;
        for (x = -1; x <= +1; x += 2)
          { for (y = -1; y <= +1; y += 2)
              { for (z = -1; z <= +1; z += 2)
                  { S2Point P = (S2Point){{(double)x, 0.0, 0.0}};
                    S2Point Q = (S2Point){{0.0, (double)y, 0.0}};
                    S2Point R = (S2Point){{0.0, 0.0, (double)z}};
                    /* Take care to call in CCW order: */
                    if (x*y*z > 0)
                      { processTriangle(&P, &Q, &R); }
                    else                          
                      { processTriangle(&R, &Q, &P); }
                  } 
              }
          }
      }
  }

#define AxisNBarbs 20
#define AxisHeadLength 0.150
#define AxisHeadRadius 0.015

void SPPlot_Axis
  ( SPPlot_Stream *fps, 
    SPH3_PMap *map, 
    r3_t *dir,
    double length,
    Sign sense
  )
  { R3Point c;
    
    if (map == NULL) { return; }
    pswr_comment(fps, "SPPlot_Axis");
    /* Forward part of arrow: */
    c = *dir; r3_dir(&c, &c); /* Point where arrow exits sphere. */
    if ((sense == 0) || PointIsOK(&c, map, sense < 0))
      { /* Draw forward part of arrow: */
        R3Point b, t;
        R3Vector u, v;
        ProjPoint cm, tm;
        int i;
        cm = Project(&c, map);
        r3_scale(length, dir, &t); /* Forward endpoint of arrow. */
        tm = Project(&t, map);
        /* Draw arrow shaft: */
        pswr_segment(fps, cm.c[0], cm.c[1], tm.c[0], tm.c[1]);
        r3_scale(length - AxisHeadLength, dir, &b); /* Start of arrowhead. */
        /* Pick two orthogonal directions, orthogonal to {dir}: */
        u = FindOrthoDir(dir);
        r3_cross(dir, &u, &v);
        /* Generate arrowhead as a cone of segments: */  
        for (i = 1; i <= AxisNBarbs; i++)
          { double alpha = TWOPI*((double)i)/((double)AxisNBarbs);
            double c = cos(alpha);
            double s = sin(alpha);
            R3Vector uv;
            R3Point p;
            ProjPoint pm;
            r3_mix(c*AxisHeadRadius, &u, s*AxisHeadRadius, &v, &uv);
            r3_add(&b, &uv, &p);
            pm = Project(&p, map);
            pswr_segment(fps, tm.c[0], tm.c[1], pm.c[0], pm.c[1]);
          }
      }
    /* Backward part of arrow: */
    r3_neg(dir, &c); r3_dir(&c, &c); /* Point where arrow enters sphere. */
    if ((sense == 0) || PointIsOK(&c, map, sense < 0))
      { /* Draw backward part of arrow */
        R3Point t; 
        ProjPoint cm, tm;
        r3_scale(-length, dir, &t); /* Backward endpoint of arrow. */
        cm = Project(&c, map);
        tm = Project(&t, map);
        pswr_segment(fps, cm.c[0], cm.c[1], tm.c[0], tm.c[1]);
      }
  }

/* Limits to avoid excessive plotting for bad choices of {fStep}: */
#define MaxIsolinesInTriangle (2*DefaultIsolines)
#define MaxIsolinesInRange (2*DefaultIsolines)

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
  )
  {
    bool_t complained = FALSE;

    auto void TriangleIsolines
      ( S2Point *P, double fP,
        S2Point *Q, double fQ,
        S2Point *R, double fR 
      );
      /* Plots isolines of the affine scalar field in a flat triangle,
        given the corner points and their function values. */
    
    void TriangleIsolines
      ( S2Point *P, double fP,
        S2Point *Q, double fQ,
        S2Point *R, double fR 
      )
      { 
        if (! TriangleIsOK(P, Q, R, map, FALSE)) { return; }
        
        ProjPoint Pm = Project(P, map);
        ProjPoint Qm = Project(Q, map);
        ProjPoint Rm = Project(R, map);
        double levMin, levMax;
        int kMin, kMax, k;
        
        if (fMin > fMax) { return; }
        
        /* Find function range in triangle (assume linear interpolation): */
        levMin = (fP < fQ ? fP : fQ);
        levMin = (levMin < fR ? levMin : fR);
        
        levMax = (fP > fQ ? fP : fQ);
        levMax = (levMax > fR ? levMax : fR);
        
        /* Don't plot isolines outside range {[fMin _ fMax]}: */
        if (levMin < fMin) { levMin = fMin; }
        if (levMax > fMax) { levMax = fMax; }
        
        if (levMin <= levMax)
          { /* Sanity check to avoid excessive plotting: */
            if ((levMax - levMin) > fStep * MaxIsolinesInTriangle)
              { if (! complained) 
                  { fprintf(stderr, "** too many isolines in triangle");
                    fprintf(stderr, " levMin = %g  levMax = %g", levMin, levMax);
                    fprintf(stderr, " fStep = %g\n", fStep);
                    complained = TRUE;
                  }
                levMin = 0.0; levMax = 0.0;
              }
            /* Find number range of isolines: */
            kMin = ceil((levMin - fStart)/fStep);
            kMax = floor((levMax - fStart)/fStep);
            /* Plot the isolines: */
            for (k = kMin; k <= kMax; k++)
              { double fLevel = fStart + ((double)k) * fStep;
                PlotZerosInTriangle(
                  fps,
                  Pm.c[0], Pm.c[1], fP - fLevel,
                  Qm.c[0], Qm.c[1], fQ - fLevel,
                  Rm.c[0], Rm.c[1], fR - fLevel
                );
              }
          }
      } 
      
    pswr_comment(fps, "SPPlot_Isolines");
    SPPlot_GenericPlot
      ( fps, tri, N,
        func, TriangleIsolines, fMinObs, fMaxObs
      );
  }
  
void SPPlot_PaintValues(
    SPPlot_Stream *fps,
    SPH3_PMap *map,          /* Perspective projection matrix. */
    ScalarField func,        /* Function to plot. */
    Triangulation *tri,      /* Reference triangulation, of NULL */
    int N,                   /* Order of mesh subdivision. */
    double fMin,             /* Nominal minimum function value. */
    Color *cMin,             /* Color to use for {fMin}. */
    bool_t clipMin,          /* TRUE omits parts below {fMin}. */
    double fMax,             /* Nominal maximum function value. */
    Color *cMax,             /* Color to use for {fMax}. */
    bool_t clipMax,          /* TRUE omits parts above {fMax}. */
    r3_t *dLight,            /* Direction towards main light source. */
    double shadow,           /* Amount of darkening by shadow. */
    double /*I/O*/ *fMinObs, /* Min {func} value seen (caller must initialize). */
    double /*I/O*/ *fMaxObs  /* Max {func} value seen (caller must initialize). */
  )
  {
    auto void TrianglePaint
      ( S2Point *P, double fP,
        S2Point *Q, double fQ,
        S2Point *R, double fR
      );
      /* Paints color-mapped values of an affine scalar field in a 
        flat triangle, given the corner points and their 
        function values. */
    
    void TrianglePaint
      ( S2Point *P, double fP,
        S2Point *Q, double fQ,
        S2Point *R, double fR
      )
      { 
        if (! TriangleIsOK(P, Q, R, map, FALSE)) { return; }
        
        ProjPoint Pm = Project(P, map);
        ProjPoint Qm = Project(Q, map);
        ProjPoint Rm = Project(R, map);
        double fMed;
        Color tc;
      
        fMed = (fP + fQ + fR) / 3.0;
        tc = InterpolateColor(fMed, fMin, cMin, clipMin, fMax, cMax, clipMax);
        if (tc.c[0] < 0.0)
          { /* invisible color */ return; }
        else
          { S2Point dMed; 
            double illum;
            Color ac;
            r3_add(P, Q, &dMed);
            r3_add(R, &dMed, &dMed);
            r3_dir(&dMed, &dMed);
            illum = 1.0 + shadow * r3_cos(&dMed, dLight);
            ac = (Color){{illum*tc.c[0], illum*tc.c[1], illum*tc.c[2]}};
            ac = ClipColor(&ac);
            pswr_set_fill_color(fps, ac.c[0], ac.c[1], ac.c[2]);
            pswr_triangle
              ( fps, 
                Pm.c[0], Pm.c[1], 
                Qm.c[0], Qm.c[1], 
                Rm.c[0], Rm.c[1], 
                TRUE, FALSE
              );
          }
      }
  
    pswr_comment(fps, "SPPlot_PaintValues");
    SPPlot_GenericPlot
      ( fps, tri, N,
        func, TrianglePaint, fMinObs, fMaxObs
      );
  }
  
void SPPlot_Vectors
  ( SPPlot_Stream *fps,
    SPH3_PMap *map,      /* Perspective projection matrix. */
    VectorField fld,     /* Vector field defined on the sphere. */
    Triangulation *tri,  /* Reference triangulation, of NULL */
    int N,               /* Order of point grid inside each mesh triangle. */
    double dotRadius,    /* Dot radius (in mm, irrespective of scale). */
    double scale         /* Scale factor for vector field. */
  )
  {
    auto void ProcessTriangle
      ( S2Point *P, double fP,
        S2Point *Q, double fQ,
        S2Point *R, double fR
      );
      /* Plots the vector field at the center of the given triangle.
        The parameters {fP}, {fQ}, {fR} are ignored. */
    
    auto double DummyFunc(S2Point *p);
      /* A dummy function to evaluate at mesh corners */
  
    void ProcessTriangle
      ( S2Point *P, double fP,
        S2Point *Q, double fQ,
        S2Point *R, double fR
      )
      { r3_t C; int k;
        for (k = 0; k < 3; k++)
          { C.c[k] = (P->c[k] + Q->c[k] + R->c[k])/3.0; }
        r3_dir(&C, &C);
        if (! PointIsOK(&C, map, FALSE)) { return; }
          { S2Point T; /* Tip of whisker. */
            ProjPoint Cm, Tm;
            R3Vector v = fld(&C);
            r3_scale(scale, &v, &T);
            r3_add(&C, &T, &T);
            Cm = Project(&C, map);
            Tm = Project(&T, map);
            pswr_set_fill_color(fps, 0,0,0);
            pswr_dot(fps, Cm.c[0], Cm.c[1], dotRadius,  TRUE, FALSE);
            pswr_segment(fps, Cm.c[0], Cm.c[1], Tm.c[0], Tm.c[1]);
          }
      }
      
    double DummyFunc(S2Point *T)
      { return 0.0; }
    
    double fMinObs = INFINITY;
    double fMaxObs = -INFINITY;
    
    pswr_comment(fps, "SPPlot_Vectors");
    SPPlot_GenericPlot
      ( fps, tri, N,
        DummyFunc, ProcessTriangle, &fMinObs, &fMaxObs
      );
  }

Color InterpolateColor
  ( double f,         /* Function value */
    double fMin,      /* Minimum function value. */
    Color *cMin,      /* Color to use for {fMin}. */
    bool_t clipMin,     /* TRUE omits parts below {fMin}. */
    double fMax,      /* Maximum abs function value. */
    Color *cMax,      /* Color to use for {fMax}. */
    bool_t clipMax      /* TRUE omits parts above {fMax}. */
  )
  { if (f < fMin)
      { if (clipMin) { return Invisible; } else { return *cMin; } }
    else if (f > fMax)
      { if (clipMax) { return Invisible; } else { return *cMax; } }
    else
      { double s = (f - fMin)/(fMax - fMin);
        double t = 1.0 - s;
        return (Color)
          {{
            s*cMax->c[0] + t*cMin->c[0],
            s*cMax->c[1] + t*cMin->c[1],
            s*cMax->c[2] + t*cMin->c[2]
          }};
      }
  }
  
Color ClipColor(Color *c)
  { 
    double m = r3_L_inf_norm(c);
    if (m <= 1.0)
      { return *c; } 
    else
      { Color a = (Color){{c->c[0]/m, c->c[1]/m, c->c[2]/m}};
        double yc = 0.3*c->c[0] + 0.6*c->c[1] + 0.1*c->c[2];
        double ya = yc/m;
        double s = (yc < 1.0 ? yc : 1.0) - ya;
        return (Color)
          {{
            a.c[0] + s*(1.0 - a.c[0]),
            a.c[1] + s*(1.0 - a.c[1]),
            a.c[2] + s*(1.0 - a.c[2])
          }};
      }
  }

void PlotZerosInTriangle
  ( SPPlot_Stream *fps,
    double Px, double Py, double Pf,
    double Qx, double Qy, double Qf,
    double Rx, double Ry, double Rf
  )
  { if (((Pf == 0.0) + (Qf == 0.0) + (Rf == 0.0)) >= 2 )
      { if ((Pf == 0.0) && (Qf == 0.0) && (Rf == 0.0))
          { /* Ignore triangle. */ return;
            /* pswr_fill_triangle(fps, Px, Py, Qx, Qy, Rx, Ry, 0,0,0); */
          }
        if ((Pf == 0.0) && (Qf == 0.0))
          { pswr_segment(fps, Px, Py, Qx, Qy); }
        if ((Qf == 0.0) && (Rf == 0.0)) 
          { pswr_segment(fps, Qx, Qy, Rx, Ry); }
        if ((Pf == 0.0) && (Rf == 0.0))
          { pswr_segment(fps, Px, Py, Rx, Ry); }
        return;
      }
    if ((Pf >= 0.0) && (Qf >= 0.0) && (Rf >= 0.0))
      { /* Function is entirely non-negative */ return; }
    if ((Pf <= 0.0) && (Qf <= 0.0) && (Rf <= 0.0))
      { /* Function is entirely non-positive */ return; }
    /* Now there is at least one positive corner and one negative corner */
    /* Permute corners so that {Rf <= Pf <= Qf}: */
    if (Rf > Pf)
      { double temp;
        temp = Px; Px = Rx; Rx = temp;
        temp = Py; Py = Ry; Ry = temp;
        temp = Pf; Pf = Rf; Rf = temp;
      }
    if (Pf > Qf)
      { double temp;
        temp= Px; Px = Qx; Qx = temp;
        temp= Py; Py = Qy; Qy = temp;
        temp= Pf; Pf = Qf; Qf = temp;
      }
    if (Rf > Pf)
      { double temp;
        temp = Px; Px = Rx; Rx = temp;
        temp = Py; Py = Ry; Ry = temp;
        temp = Pf; Pf = Rf; Rf = temp;
      }

    affirm(Rf <= Pf , "sort error(Rf,Pf)");
    affirm(Pf <= Qf , "sort error(Pf,Qf)");
    affirm(Rf < 0.0 , "sign error(Rf)");
    affirm(Qf > 0.0 , "sign error(Qf)");
    /* Now swap q,r so that {sign(Pf) != sign(Rf)}: */
    if (Pf < 0.0 )
      { double temp;
        temp = Qx; Qx = Rx; Rx = temp;
        temp = Qy; Qy = Ry; Ry = temp;
        temp = Qf; Qf = Rf; Rf = temp;
      }
    /*now {Qf != 0.0, Rf != 0.0} */
    /* also {(Pf == 0.0) || ((sign(Pf) != sign(Qf)) && (sign(Pf) != sign(Rf))}. */
    affirm(Qf != 0.0 , "zero error(Qf)");
    affirm(Rf != 0.0 , "zero error(Rf)");
    affirm((Qf > 0.0) != (Rf > 0.0) , "rel sign error(Qf,Rf)");
    if (Pf != 0.0) 
      { affirm(((Pf > 0.0) == (Qf > 0.0)) && ((Pf > 0.0) != (Rf > 0.0)), "cond");
        { double SPR = Rf/(Rf - Pf);
          double SRP = Pf/(Pf - Rf);
          double PRx = Px*SPR + Rx*SRP;
          double PRy = Py*SPR + Ry*SRP;
          double SQR = Rf/(Rf - Qf);
          double SRQ = Qf/(Qf - Rf);
          double QRx = Qx*SQR + Rx*SRQ;
          double QRy = Qy*SQR + Ry*SRQ;
          pswr_segment(fps, PRx, PRy, QRx, QRy);
        }
      }
    else
      { double SQR = Rf/(Rf - Qf);
        double SRQ = Qf/(Qf - Rf);
        double QRx = Qx*SQR + Rx*SRQ;
        double QRy = Qy*SQR + Ry*SRQ;
      
        affirm((Qf > 0.0) != (Rf > 0.0) , "rel sign error (Qf,Rf)");
        pswr_segment(fps, Px, Py, QRx, QRy);
      }
  }

ProjPoint Project(R3Point *p, SPH3_PMap *map)
  { 
    if (map == NULL)
      { /* Longitude-latitude conversion: */
        double r = hypot(p->c[0], p->c[1]);
        double R = hypot(r, p->c[2]);
        if (r > 0)
          { double lon = atan2(p->c[1], p->c[0])/PI;
            double lat = atan2(p->c[2], r)/PI;
            return (ProjPoint){{lon, lat, R}};
	  }
	else if (p->c[2] > 0)
	  { return (ProjPoint){{0.0, +1.0, R}}; }
        else if (p->c[2] < 0)
	  { return (ProjPoint){{0.0, -1.0, R}}; }
        else 
	  { return (ProjPoint){{0.0, 00.0, R}}; }
      }
    else
      { 
	/* Perspective projection: */
    	SPH3_Point q = (SPH3_Point){{{1.0, p->c[0], p->c[1], p->c[2]}}};
    	SPH3_Point r = SPH3_MapPoint(&q, map);
    	double w = r.h.c[0];
    	return (ProjPoint){{r.h.c[1]/w, r.h.c[2]/w, r.h.c[3]/w}};
      }
  }

bool_t PointIsOK(S2Point *p, SPH3_PMap *map, bool_t backSide)
  { 
    if (map == NULL) { return TRUE; }
    SPH3_Point ph = SPH3_FromCartesian(p);
    /* The observer position is row 3 of {map->inv}: */
    r4x4_t *m = &(map->inv);
    SPH3_Point obs = (SPH3_Point){{{m->c[3][0], m->c[3][1], m->c[3][2], m->c[3][3]}}};
    /* The horizon plane as seen from that observer: */
    SPH3_Plane H = (SPH3_Plane){{{-obs.h.c[0], obs.h.c[1], obs.h.c[2], obs.h.c[3]}}};
    return (r4_dot(&(ph.h), &(H.f)) >= 0.0) != backSide;
  }            

bool_t SegmentIsOK(S2Point *p, S2Point *q, SPH3_PMap *map, bool_t backSide)
  {
    if (map != NULL)
      { if (! PointIsOK(p, map, backSide)) { return FALSE; }
        if (! PointIsOK(q, map, backSide)) { return FALSE; }
      }
    else
      { /* Check for pole-crossing or longitude fold-over: */
        if ((p->c[1] == 0) && (q->c[1] == 0))
          { /* Segment lies on the {Y=0} plane. Check X of both points: */
            if (p->c[0] <= 0) { return FALSE; }
            if (q->c[0] <= 0) { return FALSE; }
          }
        else if (p->c[1]*q->c[1] <= 0)
	  { /* Segment crosses the {Y=0} plane. Check X of crossing: */
            double mz = p->c[0]*q->c[1] - p->c[1]*q->c[0];
            if (p->c[1] < q->c[1]) { mz = -mz; }
            if (mz >= 0) { return FALSE; }
	  }
      }
    return TRUE;
  }

bool_t TriangleIsOK(S2Point *p, S2Point *q, S2Point *r, SPH3_PMap *map, bool_t backSide)
  { 
    if (map != NULL)
      { /* Persp-plot - triangle is OK if obs and org are separated by plane of triang: */
        ProjPoint Pm = Project(p, map);
        ProjPoint Qm = Project(q, map);
        ProjPoint Rm = Project(r, map);
	SPH3_Point Ph = SPH3_FromCartesian(&Pm);
	SPH3_Point Qh = SPH3_FromCartesian(&Qm);
	SPH3_Point Rh = SPH3_FromCartesian(&Rm);
	SPH3_Plane H = SPH3_PlaneFromThreePoints(&Ph, &Qh, &Rh);
	SPH3_Point pObs = (SPH3_Point){{{0.0, 0.0, 0.0, 1.0}}}; /* Projected obs. */
	SPH3_Point pCtr = (SPH3_Point){{{1.0, 0.0, 0.0, 0.0}}}; /* Proj. pt in sphere. */
	Sign sCtr = SPH3_Side(&pCtr, &H);
	Sign sObs = SPH3_Side(&pObs, &H);
	return (sCtr != sObs) != backSide;
      }
    else
      { /* Lon/lat plot - triangle is bad iff any side crosses the lon cutline: */
        if (! SegmentIsOK(p, q, map, backSide)) { return FALSE; }
	if (! SegmentIsOK(q, r, map, backSide)) { return FALSE; }
	if (! SegmentIsOK(r, p, map, backSide)) { return FALSE; }
	return TRUE;
      }
  }

void SPPlot_Everything
  ( SPPlot_Stream *fps,           /* Poststcript stream, ready to plot. */
    ScalarField func,        /* The function to plot. */
    Triangulation *tri,      /* Triangulation to draw, or NULL */
    double relMeshSize,      /* Maximum step/triangle size (radians). */
    bool_t showTriang,         /* TRUE to plot the triangulation {tri}. */
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
    bool_t verbose,            /* TRUE to print messages along the way. */
    double /*I/O*/ *fMinObs, /* Min {func} value seen (caller must initialize). */
    double /*I/O*/ *fMaxObs  /* Max {func} value seen (caller must initialize). */
  )
  { Color FMaxColor = (Color){{1.0, 0.5, 0.4}};
    Color FMedColor = (Color){{0.9, 0.9, 0.9}};
    Color FMinColor = (Color){{0.1, 0.9, 0.9}};

    /* Triangle area and side (assuming uniform triangulation): */
    double triArea = FOURPI/((double)(tri == NULL ? 8 : tri->side.ne));
    double triDiam = sqrt(4.0*triArea/SQRT3);
    /* Mesh order needed to reach the given {relMeshSize}: */
    int N = ceil(triDiam/relMeshSize);

    SPH3_PMap *map = SPPlot_PerspMap(obs, rad, upp);
    affirm(fRange > 0.0, "invalid fRange");
    affirm(fStep > 0.0, "invalid fStep");

    mumble("color scale range [ %12.4e _ %12.4e ]\n", -fRange, fRange);
    mumble("isoline spacing %12.4e\n", fStep);
    double shadow = (map == NULL ? 0.0 : 0.1);
    
    if (map != NULL)
      { mumble("Back axes...\n");;
        pswr_set_pen(fps, 0.0,0.0,0.0, 2.0 * lineWidth, 0,0);
        SPPlot_CoordAxes(fps, map, -1);
      }

    mumble("Negative shading...\n");;
    SPPlot_PaintValues
      ( fps, map,
        func, 
        tri, N,
        /*fMin*/ -fRange, 
        /*cMin*/ &FMinColor,
        /*clipMin*/ FALSE,
        /*fMax*/ 0.0,
        /*cMax*/ &FMedColor,
        /*clipMax*/ TRUE,
        /*dLight*/ dLight,
        /*shadow*/ shadow,
        fMinObs, 
        fMaxObs
      ); 
    mumble("Positive shading...\n");;
    SPPlot_PaintValues
      ( fps, map,
        func, 
        tri, N,
        /*fMin*/ 0.0,
        /*cMin*/ &FMedColor,
        /*clipMin*/ TRUE,
        /*fMax*/ fRange,
        /*cMax*/ &FMaxColor,
        /*clipMax*/ FALSE,
        /*dLight*/ dLight,
        /*shadow*/ shadow,
        fMinObs,
        fMaxObs
      ); 

    if (showTriang && (tri != NULL))
      { mumble("Triangulation...\n");;
        pswr_set_pen(fps, 0.0,0.0,0.0, lineWidth, 0,0);
        SPPlot_Triangulation(fps, map, tri, relMeshSize, FALSE);
      }

    if (supp != NULL)
      { mumble("Supporting plane...\n");;
        pswr_set_pen(fps, 0.0,0.7,0.0, lineWidth, 0,0);
        SPPlot_SupportingPlane(fps, map, supp, relMeshSize);
      }

    mumble("Negative isolines...\n");;
    pswr_set_pen(fps, 0.0,0.0,1.0, 0.75 * lineWidth, 0,0);
    SPPlot_Isolines
      ( fps, map,
        func, 
        tri, N,
        /* fMin */ -fRange,
        /* fMax */ 0.0,
        /*fStep*/ fStep,
        /*fStart*/ -0.5 * fStep,
        fMinObs,
        fMaxObs
      ); 
    mumble("Positive isolines...\n");;
    pswr_set_pen(fps, 0.5,0.0,0.0, 0.75 * lineWidth, 0,0);
    SPPlot_Isolines
      ( fps, map,
        func, 
        tri, N,
        /* fMin */ 0.0,
        /* fMax */ fRange,
        /*fStep*/ fStep,
        /*fStart*/ 0.5 * fStep,
        fMinObs,
        fMaxObs
      ); 

    if ((gridNLon != 0) && (gridNLat != 0))
      { mumble("Lon-Lat grid...\n");
        pswr_set_pen(fps, 0.5,0.0,0.0, 0.75 * lineWidth, 0,0);
        pswr_set_fill_color(fps, 0.0,0.0,0.0);
        SPPlot_LonLatGrid(fps, map, gridNLon, gridNLat, gridDots, 2.0*lineWidth, FALSE); 
      }

    mumble("Sphere outline...\n");
    pswr_set_pen(fps, 0.0,0.0,0.0, lineWidth, 0,0);
    SPPlot_Sphere(fps, map);

    mumble("Front axes...\n");
    pswr_set_pen(fps, 0.0,0.0,0.0, 2.0 * lineWidth, 0,0);
    SPPlot_CoordAxes(fps, map, +1);
  }

void SPPlot_MultipleViews
  ( SPPlot_Stream *fps,      /* Postscript stream. */
    char *funcTag,           /* Prefix for for canvas names. */
    ScalarField func,        /* The function to plot. */
    Triangulation *tri,      /* Triangulation to draw, or NULL */
    double relMeshSize,      /* Maximum step/triangle size (radians). */
    bool_t showTriang,       /* TRUE to plot the triangulation {tri}. */
    SPH3_Plane *supp,        /* Supporting plane, or NULL */
    double fRange,           /* Nominal maximum absolute value of function. */
    double fStep,            /* Isoline spacing. */
    SPH3_Point *obs,         /* Observer's position, or NULL. */
    SPH3_Point *upp,         /* Camera vertical reference */
    double rad,              /* Radius of region of interest. */
    r3_t *dLight,            /* Direction towards main light source. */
    double lineWidth,        /* Nominal line width in mm. */
    int gridNLon,            /* Number of grid steps in longitude. */
    int gridNLat,            /* Number of grid steps in latitude. */
    bool_t gridDots,         /* TRUE shows lon-lat grid as face-ctr dots, FALSE as rects. */
    int aSide,               /* First view to plot (+1=front, -1=rear, 0=none) */
    int bSide,               /* Second view to plot (ditto). */
    string_vec_t *caption,   /* Figure caption, or NULL. */
    int index,               /* Figure `index' (for "%I" expansion in {caption}). */
    double time,             /* Figure `time' (for "%T" expansion in {caption}). */
    double error,            /* Figure `error' (for "%E" expansion in {caption}). */
    double capAlign,         /* Caption alignment (0.0=left, 0.5=ctr, 1.0=right). */
    bool_t verbose,          /* TRUE to print messages along the way. */
    double /*I/O*/ *fMinObs, /* Min {func} value seen (caller must initialize). */
    double /*I/O*/ *fMaxObs  /* Max {func} value seen (caller must initialize). */
  )
  { 
    mumble("color scale range [ %12.4e _ %12.4e ]\n", -fRange, fRange);
    mumble("isoline spacing %12.4e\n", fStep);
    bool_t projLonLat = (obs == NULL);

    int is;
    for (is = 0; is < 2; is++)
      { int side = (is == 0 ? aSide : bSide);
        if (side != 0)
          { /* Get observer for this view: */
            SPH3_Point obslocal, *obsx;
            if (projLonLat) 
              { obsx = NULL; }
            else
              { obsx = &obslocal; (*obsx) = (*obs); 
                if (side < 0) 
                  { /* Flip observer across origin: */
                    int k;
                    for (k = 1; k <= 3; k++) { obsx->h.c[k] = - obsx->h.c[k]; }
                  }
	      }
            /* Start new figure: */
            if (fps->eps)
              { /* Start new canvas (hence a new file): */
                char *figTag = txtcat(funcTag, (side == 1 ? "-f" : "-r"));
                pswr_new_canvas(fps, figTag);
                free(figTag);
              }
            double xMax = 1.0;
            double yMax = (projLonLat ? 0.5 : 1.0);
            pswr_new_picture(fps, -xMax,+xMax, -yMax,+yMax);
            /* Plot view: */
            SPPlot_Everything(
              fps,
              func, 
              tri, relMeshSize, showTriang,
              supp,
              fRange, fStep,
              obsx, upp,
              rad,
              dLight, 
              lineWidth,
              gridNLon, gridNLat, gridDots,
              verbose,
              fMinObs, fMaxObs
            );
            /* Add caption, etc: */
            if (caption != NULL)
              { string_vec_t xcap = 
                  SPPlot_ExpandCaption
                    ( caption, side, index, time, error, fRange, fStep, funcTag );
                int k; 
                for (k = 0; k < xcap.ne; k++)
                  { pswr_add_caption(fps, xcap.e[k], capAlign); }
                SPPlot_FreeCaption(&xcap);
              }
          }
      }
  }

/* PERSPECTIVE PROJECTION */

SPH3_PMap *SPPlot_PerspMap
  ( SPH3_Point *obs, 
    double rad,
    SPH3_Point *upp
  )
  { 
    if (obs == NULL) { return NULL; }
    SPH3_Point foc = SPPlot_FocusFromRadiusOfInterest(obs, rad);
    SPH3_PMap *map = notnull(malloc(sizeof(SPH3_PMap)), " no mem");
    (*map) = SPH3_PerspMap(obs, &foc, upp);
    /* Now combine {map} with a uniform scale of {1/rad}: */
    { int i, j;
      r4x4_t *md = &(map->dir);
      r4x4_t *mi = &(map->inv);
      for (i = 0; i < 4; i++)
        for (j = 1; j < 4; j++)
          { md->c[i][j] /= rad; mi->c[j][i] *= rad; }
    }
    return map;
  }

void SPPlot_FixObserver
  ( SPH3_Point *obs, 
    double alpha, 
    SPFunction *f, 
    bool_t verbose
  )
  { /* If no place to return result, do nothing: */
    if (obs == NULL) { return; }
    mumble("checking observer...\n");
    /* If given {obs} is invalid, ignore it: */
    if (r4_norm(&(obs->h)) == 0.0) 
      { mumble("  no observer given, ignoring it.\n");
        alpha = 1.0;
      }
    /* If no adjustment allowed, do nothing: */
    if (alpha == 0.0) { return; }
    { /* Define the ideal observer {idealObs}: */
      SPH3_Point idealObs = IdealObserver(f, verbose);
      if (alpha == 1.0)
        { *obs = idealObs; return; }
      else
        { SPH3_Point oldObs = *obs;
          mumble("  mixing with given observer (alpha = %8.6f)...\n", alpha);
          /* Normalize homog coords, to make {alpha} meaningful: */
          r4_dir(&(idealObs.h), &(idealObs.h));
          r4_dir(&(oldObs.h), &(oldObs.h));
          r4_mix(alpha, &(idealObs.h), 1.0-alpha, &(oldObs.h), &(obs->h));
          return;
        }
    }
  }
  
SPH3_Point IdealObserver(SPFunction *f, bool_t verbose)
  { SPSpline *fpw;
    mumble("  looking for ideal observer position...\n");
    if (f != NULL) 
      { if ((fpw = SPSpline_Cast(f)) != NULL)
          { /* It is an {SPSpline}, try its supporting plane: */
            SPH3_Plane *S = &(fpw->d->supp);
            R3Vector dir = (R3Vector){{S->f.c[1], S->f.c[2], S->f.c[3]}};
            if (r3_norm(&dir) > 0.0)
              { mumble("    using supporting plane normal.\n");
                return (SPH3_Point){{{0.0, dir.c[0], dir.c[1], dir.c[2]}}};
              }
            else
              { mumble("    supporting plane is undefined.\n"); }
          }
        /* Try the barycenter of the function values squared: */
        mumble("    computing barycenter of function squared...\n");
        { 
          auto double f2(double fv, S2Point *p);
          
          double f2(double fv, S2Point *p)
            { return fv*fv; }
          
          FuncMap F2 = (FuncMap){&f2, TRUE, "u^2"};
          R3Vector fbar = SPFunction_Centroid(f, F2, NULL, NULL);
          if (r3_norm(&fbar) > 0.0)
            { mumble("    using barycenter of function squared.\n");
              r3_dir(&fbar, &fbar);
              return (SPH3_Point){{{0.0, fbar.c[0], fbar.c[1], fbar.c[2]}}};
            }
          else
            { mumble("    function is too uniform.\n"); }
        }
      }
    /* No luck, pick a default observer: */
    mumble("    using default observer position.\n");
    return SPPlot_DefaultObserver();
  }

SPH3_Point SPPlot_DefaultObserver(void)
  { R3Vector d = (R3Vector){{3.0, 2.0, 1.0}};
    r3_dir(&d, &d);
    return (SPH3_Point){{{0.0, d.c[0], d.c[1], d.c[2]}}};
  }

SPH3_Point SPPlot_DefaultZenith(void)
  { return (SPH3_Point){{{0.0, 0.0, 0.0, 1.0}}}; }

SPH3_Point SPPlot_SecondZenith(void)
  { return (SPH3_Point){{{0.0, 0.0, 1.0, 0.0}}}; }

bool_t SPPlot_IsGoodZenith(SPH3_Point *upp, SPH3_Point *obs)
  { SPH3_Point org = Origin; 
    if ((upp == NULL) || (r4_norm(&(upp->h)) == 0.0)) { return FALSE; }
    if (upp->h.c[0] < 0.0) { return FALSE; }
    if ((obs == NULL) || (r4_norm(&(obs->h)) == 0.0))
      { return TRUE; }
    else
      { R3Vector u = SPH3_Dir(&org, upp);
        R3Vector v = SPH3_Dir(&org, obs);
        if (fabs(r3_dot(&u, &v)) >= 0.9999) { return FALSE; }
      }
    return TRUE;
  }

void SPPlot_FixZenith(SPH3_Point *upp, SPH3_Point *obs, bool_t verbose)
  { SPH3_Point r;
    /* If no place to put result, do nothing: */
    if (upp == NULL) { return; }
    mumble("checking zenith reference...\n");
    if (SPPlot_IsGoodZenith(upp, obs)) { return; }
    r = SPPlot_DefaultZenith(); 
    if (SPPlot_IsGoodZenith(&r, obs)) { *upp = r; return; }
    r = SPPlot_SecondZenith(); 
    if (SPPlot_IsGoodZenith(&r, obs)) { *upp = r; return; }
    affirm(FALSE, "zenith failure");
  }

void SPPlot_FixRadiusOfInterest
  ( double *rad, 
    double alpha, 
    SPFunction *f, 
    SPH3_Point *obs, 
    bool_t verbose
  )
  { /* If no place to put result, do nothing: */
    if (rad == NULL) { return; }
    mumble("checking radius of interest...\n");
    /* If the given radius is invalid, ignore it: */
    if ((*rad) < 0.0) 
      { mumble("  given radius (%.5f) is invalid, ignoring it.\n", (*rad));
        alpha = 1.0;
      }
    /* Combine radius and ideal radius: */
    if (alpha != 0.0)
      { double idealRad = SPPlot_IdealRadiusOfInterest(f, obs, verbose);
        mumble("  mixing ideal radius (%.5f) with given radius (alpha = %8.6f)...\n", idealRad, alpha);
        *rad = alpha*idealRad + (1.0 - alpha)*(*rad);
        if ((*rad) > MaxRadiusOfInterest)
          { mumble
              ( "  computed radius (%.5f) is too big, using %.3f.\n", 
                (*rad), SPPlot_FullRadius
              );
            (*rad) = SPPlot_FullRadius;
            
          }
      }
    /* Make sure that the computed {rad} is not too small: */
    if ((*rad) < MinRadiusOfInterest)
      { mumble("  computed radius (%.5f) is too small, using %.5f.\n", (*rad), MinRadiusOfInterest);
        (*rad) = MinRadiusOfInterest;
      }
  }

double SPPlot_IdealRadiusOfInterest(SPFunction *f, SPH3_Point *obs, bool_t verbose)
  { SPSpline *fpw;
    mumble("  computing ideal radius of interest...\n");
    if ((obs == NULL) || (r4_norm(&(obs->h)) == 0.0))
      { mumble("    observer is undefined: ideal radius is 1.\n");
        return 1.0;
      }
    affirm(obs->h.c[0] >= 0, "observer beyond infinity");
    /* Get angular radius of interest from function {f}: */
    { S2Point odr;
      double rvis, rint;
      /* Compute direction of observer {odr} and visible cap radius {rvis}: */
      { R3Point obsc = (r3_t){{obs->h.c[1], obs->h.c[2], obs->h.c[3]}};
        double obsm = r3_norm(&obsc);
        double cvis = obs->h.c[0]/obsm;
        r3_scale(1/obsm, &obsc, &odr);
        if (cvis < -1.0) { cvis = -1.0; }
        affirm(cvis > -0.00001, "observer inside sphere?");
        if (cvis >  0.0) { cvis =  0.0; }
        rvis = sqrt(1.0 - cvis*cvis);
      }
      /* Compute angular width {rint} of interesting cap: */ 
      if (f == NULL)
        { mumble("    function is NULL, interesting part is whole sphere.\n");
          rint = 1.0;
        } 
      else if ((fpw = SPSpline_Cast(f)) != NULL)
        { rint = SPPlot_RadiusOfInterestFromSupport(&odr, &(fpw->d->supp), verbose); }
      else
        { rint = SPPlot_RadiusOfInterestFromValues(&odr, f, verbose); }
      /* Clip interesting cap radius {rint} to visible cap radius {rvis}: */
      if (rint >= rvis)
        { mumble("    region of iterest covers the entire visible part of S^2.\n");
          return rvis;
        }
      else
        { return rint; }
    }
  }
  
double SPPlot_RadiusOfInterestFromValues(r3_t *odr, SPFunction *f, bool_t verbose)
  { 
    auto double f2(double fv, S2Point *p);
      /* Function value squared */
    double f2(double fv, S2Point *p) { return fv*fv; }
    FuncMap F2 = (FuncMap){f2, TRUE, "u^2"};

    auto double f2d(double fv, S2Point *p);
      /* Function value squared, times {odr}-coordinate: */
    double f2d(double fv, S2Point *p)
      { double dv = r3_dot(p, odr);
        return fv*fv*dv;
      }
    FuncMap F2d = (FuncMap){f2d, TRUE, "u^2*dv"};

    auto double f2d2(double fv, S2Point *p);
      /* Function value squared, times {odr}-coordinate squared: */
    double f2d2(double fv, S2Point *p)
      { double dv = r3_dot(p, odr);
        return fv*fv*dv*dv;
      }
    FuncMap F2d2 = (FuncMap){f2d2, TRUE, "u^2*dv^2"};

    mumble("    computing moments of f^2(p)...\n");
    { 
      double If2 = SPFunction_Integral(f, F2, NULL, NULL, FALSE);
      double If2d = SPFunction_Integral(f, F2d, NULL, NULL, FALSE);
      double If2d2 = SPFunction_Integral(f, F2d2, NULL, NULL, FALSE);
      if (If2 <= 0.0)
        { mumble("      function is zero.\n");
          return 1.0;
        }
      else
        { double avg = If2d/If2;
          double dev = sqrt(fabs(If2d2/If2 - avg*avg));
          double cmin = avg - 2*dev;
          mumble("      projection avg = %6.3f dev = %6.3f.\n", avg, dev);
          if (cmin <= 0.0) 
            { mumble("      function is too spread out.\n");
              return 1.0;
            }
          else if (cmin >= 1.0)
            { return 0.0; }
          else
            { return sqrt(1.0 - cmin*cmin); }
        }
    }
  }

double SPPlot_RadiusOfInterestFromSupport(r3_t *odr, SPH3_Plane *S, bool_t verbose)
  { R3Vector cS = (R3Vector){{S->f.c[1], S->f.c[2], S->f.c[3]}};
    double mS = r3_norm(&cS);
    mumble("    computing radius of interest from supporting plane...\n");
    if (mS <= 0.0)
      { fprintf(stderr, "      supporting plane is undefined.\n");
        return 1.0; 
      }
    else 
      { double dS = -S->f.c[0]/mS; /* {dS} = signed dist of {S} from origin. */
        if (dS < 0.0)
          { fprintf(stderr, "      support cap is bigger than hemisphere.\n");
            return 1.0;
          }
        else
          { S2Point tip; 
            double tsAng, otAng, osAng;
            fprintf(stderr, "       expanding radius to include support cap.\n");
            tsAng = (dS >= 1.0 ? 0.0 : acos(dS)); /* Angular rad of support. */
            r3_scale(1.0/mS, &cS, &tip); /* {tip} is spher. ctr. of support cap. */
            { double otCos = r3_dot(&tip, odr);
              if (otCos < -1.0) { otCos = -1.0; }
              if (otCos > 1.0) { otCos = 1.0; }
              otAng = acos(otCos); /* Angle between {odr} and {tip}: */
            }
            osAng = otAng + tsAng; /* Angular radius of cap of interest. */
            if (osAng > PI/2)
              { mumble("      region of interest bigger that one hemisphere.\n");
                return 1.0;
              }
            else
              { return sin(osAng); }
          }
      }
  }

void SPPlot_FixView
  ( SPH3_Point *obs, 
    double alphaObs, 
    SPH3_Point *upp, 
    double *rad,
    double alphaRad,
    SPFunction *f, 
    bool_t verbose
  )
  {
    if (obs != NULL)
      { SPPlot_FixObserver(obs, alphaObs, f, verbose);
        if (verbose)
          { r4_gen_print(stderr, &(obs->h), "%6.3f", "using obs = [ ", ", ", " ]\n"); }
      }
    if (upp != NULL)
      { SPPlot_FixZenith(upp, obs, verbose);
        if (verbose)
          { r4_gen_print(stderr, &(upp->h), "%6.3f", "using upp = [ ", ", ", " ]\n"); }
      }
    if (rad != NULL) 
      { SPPlot_FixRadiusOfInterest(rad, alphaRad, f, obs, verbose); 
        mumble("using radius = %8.4f\n", *rad);
      }
  }

SPH3_Point SPPlot_FocusFromRadiusOfInterest(SPH3_Point *obs, double rad)
  { 
    if (rad >= 1.0)
      { return Origin; }
    else
      { S2Point odr;
        /* Compute direction of observer {odr}: */
        { R3Point obsc = (r3_t){{obs->h.c[1], obs->h.c[2], obs->h.c[3]}};
          r3_dir(&obsc, &odr);
        }
        if (rad <= 0.0)
          { return SPH3_FromCartesian(&odr); }
        else
          { double cint = sqrt(1.0 - rad*rad);
            r3_scale(cint, &odr, &odr);
            return SPH3_FromCartesian(&odr);
          }
      }
  }

void SPPlot_FixRange
  ( double *fRange, 
    double *fStep,
    double alpha, 
    SPFunction *f, 
    bool_t verbose
  )
  { /* If no place to return result, do nothing: */
    if (fRange == NULL) { return; }
    mumble("checking nominal function range...\n");
    /* If given range is invalid, ignore it: */
    if ((*fRange) <= 0.0) 
      { mumble("  given range (%.5f) is invalid, ignoring it.\n", (*fRange));
        alpha = 1.0;
      }
    /* If no adjustment allowed, do nothing: */
    if ((alpha != 0.0) && (f != NULL))
      { double newRange;
        { /* Estimate the maximum function value: */
          double estRange = SPPlot_EstRange(f, verbose);
          if (alpha == 1.0)
            { newRange = estRange; }
          else
            { double oldRange = *fRange;
              mumble("  mixing ideal range (%.5f) with given range (alpha = %8.6f)...\n", estRange, alpha);
              newRange = alpha * estRange + (1.0 - alpha)*oldRange;
            }
        }
        affirm(newRange > 0.0, "empty range - something failed");
        mumble("  rounding computed range (%10.4e) to a nice value...\n", newRange);
        { double step = SPPlot_RoundToNice(newRange/DefaultIsolines);
          newRange = DefaultIsolines * step;
          *fRange = newRange;
          mumble("  new range is [ %10.4e _ %10.4e ].\n", -(*fRange), (*fRange));
        }
      }
    /* Fix {*fStep} if needed: */
    if (fStep == NULL) { return; }
    mumble("checking isoline spacing...\n");
    if ((*fStep) <= (*fRange)/MaxIsolinesInRange)
      { (*fStep) = SPPlot_RoundToNice((*fRange)/DefaultIsolines);
        mumble("  given spacing is too small, fixed to %10.4e.\n", *fStep);
      }
  }

double SPPlot_EstRange(SPFunction *f, bool_t verbose)
  { SPSpline *fpw;
  
    auto double func(S2Point *p);
    
    double func(S2Point *p)
      { return f->m->eval(f, p); }

    if ((fpw = SPSpline_Cast(f)) != NULL)
      { mumble("  estimating max function value (triangulation sampling)...\n");
        return SPRange_OnTriangulation(func, fpw->d->tri, 40);
      }
    else
      { mumble("  estimating max function value (sphere sampling)...\n");
        return SPRange_OnSphere(func, 40);
      }
  }

/* MISCELLANEOUS */

R3Vector FindOrthoDir(R3Vector *d)
  { int iMax, iMed, iMin;
    R3Vector u;
    /* Find the max, med, and min coordinates of {d}: */
    if (fabs(d->c[0]) >= fabs(d->c[1]))
      { iMax = 0; iMed = 1; } 
    else
      { iMax = 1; iMed = 0; }
    if (fabs(d->c[2]) <= fabs(d->c[iMed]))
      { iMin = 2; }
    else if (fabs(d->c[2]) <= fabs(d->c[iMax]))
      { iMin = iMed; iMed = 2; }
    else
      { iMin = iMed; iMed = iMax; iMax = 2; }
    /* Build result: */
    u.c[iMax] = -d->c[iMed]; u.c[iMed] = d->c[iMax]; u.c[iMin] = 0.0;
    r3_dir(&u, &u);
    return u;
  }

S2Point SPPlot_XYZFromLatLon(double lat, double lon)
  { double cLat = cos(lat); double sLat = sin(lat);
    double cLon = cos(lon); double sLon = sin(lon);
    return (S2Point){{cLat*cLon, cLat*sLon, sLat}};
  }

double SPPlot_RoundToNice(double x)
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

char *SPPlot_IterationTag(int iter)
  { if (iter == -1) 
      { return txtcat("ini",""); }
    else if (iter == INT_MAX)
      { return txtcat("fin",""); }
    else
      { affirm(iter >= 0, "bad iteration number");
        return fmt_int(iter,6);
      }
  }

string_vec_t SPPlot_ExpandCaption
  ( string_vec_t *tmp,
    int side,
    int index,
    double time,
    double error,
    double fRange,
    double fStep,
    char *name
  )
  { int nCap = tmp->ne;
    string_vec_t cap = string_vec_new(nCap);
    int k;
    for (k = 0; k < nCap; k++)
      { 
        cap.e[k] = 
          SPPlot_ExpandCaptionLine
           ( tmp->e[k], side, index, time, error, fRange, fStep, name );
      }
    return cap;
  }

char *SPPlot_ExpandCaptionLine
  ( char *tmp,
    int side,
    int index,
    double time,
    double error,
    double fRange,
    double fStep,
    char *name
  )
  {
    int nt = strlen(tmp);
    char_vec_t cap = char_vec_new(nt);
    int it = 0, ic = 0;
    while (it < nt)
      { char c = tmp[it]; it++;
        if ((c == '%') && (it < nt)) 
          { char d = tmp[it]; it++;
            if (d == '%')
              { char_vec_expand(&cap, ic);
                cap.e[ic] = '%'; ic++;
              }
            else 
              { char *sub = NULL;
                if (d == 'T')
                  { asprintf(&sub, "%+7.4f", time); }
                else if (d == 'E')
                  { asprintf(&sub, "%+8.1e", error); }
                else if (d == 'S')
                  { asprintf(&sub, "%s", (side > 0 ? "f" : "r")); }
                else if (d == 'I')
                  { asprintf(&sub, "%06d", index); }
                else if (d == 'N')
                  { asprintf(&sub, "%s", name); }
                else if (d == 'R')
                  { asprintf(&sub, "%10.3e", fRange); }
                else if (d == 'D')
                  { asprintf(&sub, "%10.3e", fStep); }
                else
                  { asprintf(&sub, "%%%c", d); }
                int ns = strlen(sub);
                char_vec_expand(&cap, ic + ns - 1);
                int is;
                for (is = 0; is < ns; is++) { cap.e[ic] = sub[is]; ic++; }
                free(sub);
              }
          }
        else
          { char_vec_expand(&cap, ic);
            cap.e[ic] = c; ic++;
          }
      }
    char_vec_expand(&cap, ic);
    cap.e[ic] = '\000'; ic++;
    char_vec_trim(&cap, ic);
    return cap.e;
  }
  
void SPPlot_FreeCaption(string_vec_t *cap)
  {
    int nCap = cap->ne;
    int k;
    for (k = 0; k < nCap; k++)
      { 
        if (cap->e[k] != NULL) { free(cap->e[k]); cap->e[k] = NULL; }
      }
    string_vec_trim(cap, 0);
  }
