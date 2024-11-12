/* Checks whether the {grad} and {eval} methods are consistent. */
/* Last edited on 2024-11-11 05:14:12 by stolfi */

#define PROG_NAME "SPGradTest"

#include <SPFunction.h>
#include <SPSpline.h>
#include <SPTriang.h>
#include <SPPlot.h>
#include <SPOptions.h>
#include <SPPlotOptions.h>
#include <SPApprox.h>
#include <SPIntegral.h>
#include <SPBasic.h>
#include <affirm.h>
#include <nat.h>
#include <r3.h>
#include <r4.h>
#include <SPH3.h>
#include <math.h>
#include <stdio.h>
#include <values.h>

typedef struct Options
  { char *sfnName;        /* Function file (minus ".sfn" extension) or "-". */
    double delta;         /* Size of perturbation for numerical gradient. */
    char *outName;        /* Prefix for output file names. */
    int smpOrder;         /* Sampling order for test points. */
    /* Plotting options: */
    SPPlotOptions_t plt;  /* Plot format and style options. */
  } Options;
    
typedef struct GradInfo 
  { S2Point p;       /* A point on the sphere. */
    double fp;       /* Function values at {p}. */
    R3Gradient cg;   /* The Cartesian gradient at {p}. */
    r3_t dp;         /* A perturbation for {p} (not necessarily tangent). */
    r3_t q;          /* The perturbed point {p + dp}. */
    double fq;       /* Function values at {{q}. */
    double cgdp;     /* Dot product of {cg} and {dp}. */
    double error;    /* Difference between {fq} and {fp+cgdp}. */
  } GradInfo;

typedef struct SGrdInfo 
  { S2Point p;       /* A point on the sphere. */
    S2Gradient sg;   /* The spherical gradient at {p}. */
    double cosgp;    /* Cosine of {sg} and {p}. */
  } SGrdInfo;

Options GetOptions(int argn, char **argc);

void ComputeGradInfos(SPFunction *fn, S2Point *p, double delta, bool_t tg, GradInfo *GI, SGrdInfo *SI);
  /* Tests the gradient method of {fn} with displacements along three directions at {p},
    returns the data in {GI[0..2]}. If {tg} is true, the three directions are tangent
    to the sphere, at 120 degrees; otherwise they are directed along the coordinate 
    axes. */
   
void CheckGradError(SPFunction *fn, Options *o, r3_t *cgMaxErr, double *sgMaxCos);

void PlotGradError(SPFunction *fn, Options *o, r3_t cgMaxErr, double sgMaxCos);

GradInfo TestGrad(SPFunction *fn, S2Point *p, double delta, r3_t *dir);
SGrdInfo TestSGrd(SPFunction *fn, S2Point *p);

void PrintGradInfo(FILE *wr, GradInfo *GI);
void PrintSGrdInfo(FILE *wr, SGrdInfo *SI);

Triangulation *GetTriangulation(SPFunction *f);
SPH3_Plane GetSupportingPlane(SPFunction *f);

int main(int argn, char **argc)
  { Options o = GetOptions(argn, argc);
    SPIntegral_SetDefaultSamplingOrder(o.smpOrder);
    SPFunction *fn = SPApprox_ReadFunction(txtcat(o.sfnName, ".sfn"));
    r3_t cgErrMax;
    double sgAngMax;
    CheckGradError(fn, &o, &cgErrMax, &sgAngMax);
    PlotGradError(fn, &o, cgErrMax, sgAngMax);
    return 0;
  }

void ComputeGradInfos(SPFunction *fn, S2Point *p, double delta, bool_t tg, GradInfo *GI, SGrdInfo *SI)
  { 
    r3_t v[3];
    int k;
    if (tg)
      { double cos60 = 0.5;
        double sin60 = sqrt(3.0)/2.0;
        /* Generate two orthogonal vectors in the tangent space at {p}: */
        r3_t u[2];
        (void)r3_throw_ortho_pair(p, &(u[0]), &(u[1]));
        /* Generate three equally spaced vectors in the tangent space {p}: */
        v[0] = u[0]; 
        r3_mix(-cos60, &(u[0]), +sin60, &(u[1]), &(v[1]));
        r3_mix(-cos60, &(u[0]), -sin60, &(u[1]), &(v[2]));
      }
    else
      { for (k = 0; k < 3; k++)
          { v[k] = (r3_t){{ 0.0, 0.0, 0.0 }}; v[k].c[k] = 1.0; }
      }
    /* Test the gradient with those displacements: */
    for (k = 0; k < 3; k++) 
      { GI[k] = TestGrad(fn, p, delta, &(v[k])); }
    /* Test the spherical gradient: */
    (*SI) = TestSGrd(fn, p);
  }

void CheckGradError(SPFunction *fn, Options *o, r3_t *cgMaxErr, double *sgMaxCos)
  { int i, k;

    int cgNumErr = 0;         /* Number of Cartesian gradient tests along each axis. */
    r3_t cgTotErr2;           /* Sum of error squared in Cartesian gradient, in each coord. */
    GradInfo cgMaxErrGI[3];   /* Point of maximum error, in each coord. */
    
    int sgNumCos = 0;         /* Total number of spherical gradient flatness tests. */
    double sgTotCos2;         /* Sum of {cos(p,grad)^2}. */
    SGrdInfo sgMaxCosSI;      /* Point of maximum angular deviation from tangent space. */

    S2Point_vec_t sp = S2Point_vec_new(0); /* Sample points. */
    
    /* Generate sample points: */
    { double_vec_t wp = double_vec_new(0); /* Sample weights for integration. */
      int ns = 0;
      SPIntegral_SampleOctants(o->smpOrder, &sp, &wp, &ns);
      S2Point_vec_trim(&sp, ns);
      free(wp.e);
    }
    
    /* Evaluate gradient at sample points: */
    GradInfo GI[3]; 
    cgTotErr2 = (r3_t){{ 0.0, 0.0, 0.0 }};
    (*cgMaxErr) = (r3_t){{ -1.0, -1.0, -1.0 }};
    SGrdInfo SI;
    sgTotCos2 = 0;
    (*sgMaxCos) = -1.0;
    for (i = 0; i < sp.ne; i++)
      { r3_t p = sp.e[i];
        
        /* Check Cartesian gradient with axial steps: */
        ComputeGradInfos(fn, &p, o->delta, FALSE, GI, &SI);
        for (k = 0; k < 3; k++)
          { double error = GI[k].error;
            cgTotErr2.c[k] += error*error;
            if (fabs(error) > (*cgMaxErr).c[k])
              { (*cgMaxErr).c[k] = fabs(error); cgMaxErrGI[k] = GI[k]; }
          }
        cgNumErr++;
        
        /* Check flatness of spherical gradient with axial steps: */
        double cosgp = SI.cosgp;
        sgTotCos2 += cosgp*cosgp;
        if (fabs(cosgp) > (*sgMaxCos))
          { (*sgMaxCos) = fabs(cosgp); sgMaxCosSI = SI; }
        sgNumCos++;
      }
    r3_t cgAvgErr;
    for (k = 0; k < 3; k++) 
      { cgAvgErr.c[k] = sqrt(cgTotErr2.c[k]/((double)cgNumErr));
        fprintf(stderr, "gradient error on axis %d:", k);
        fprintf(stderr, " avg = %14.7e max = %14.7e\n", cgAvgErr.c[k], (*cgMaxErr).c[k]);
        fprintf(stderr, "\n");
        PrintGradInfo(stderr, &(cgMaxErrGI[k]));
        fprintf(stderr, "\n");
      }   
    
    double sgAvgCos = sqrt(sgTotCos2/((double)sgNumCos));
    fprintf(stderr, "cos(gradient-point angle):");
    fprintf(stderr, " avg = %14.7e max = %14.7e\n", sgAvgCos, (*sgMaxCos));
    fprintf(stderr, "\n");
    PrintSGrdInfo(stderr, &sgMaxCosSI);
  }
  
void PrintGradInfo(FILE *wr, GradInfo *GI)
  {
    fprintf(wr, "  p =    "); r3_print(wr, &(GI->p)); 
    fprintf(wr, "  f(p)  = %16.12f\n", GI->fp);
    fprintf(wr, "  q =    "); r3_print(wr, &(GI->q)); 
    fprintf(wr, "  f(q)  = %16.12f\n", GI->fq);
    fprintf(wr, "  dp =   "); r3_print(wr, &(GI->dp)); 
    fprintf(wr, "\n");
    fprintf(wr, "  grad = "); r3_print(wr, &(GI->cg)); 
    fprintf(wr, "\n");
    fprintf(wr, "  df(p)   = %16.12f\n", GI->fq - GI->fp);
    fprintf(wr, "  grad*dp = %16.12f\n", GI->cgdp);
    fprintf(wr, "\n");
    fprintf(wr, "  error = %14.7e\n", GI->error);
  }

void PrintSGrdInfo(FILE *wr, SGrdInfo *SI)
  {
    fprintf(wr, "  p =    "); r3_print(wr, &(SI->p)); 
    fprintf(wr, "\n");
    fprintf(wr, "  SGrd = "); r3_print(wr, &(SI->sg)); 
    fprintf(wr, "\n");
    fprintf(wr, "  cosgp = %14.7e\n", SI->cosgp);
  }

void PlotGradError(SPFunction *fn, Options *o, r3_t cgMaxErr, double sgMaxCos) 
  { 
    
    SPPlotOptions_t *po = &(o->plt);
    Triangulation *tri = GetTriangulation(fn);
    SPH3_Plane fnSupp = GetSupportingPlane(fn);
    double eMinObs, eMaxObs;
    
    /* Fix observer and zenith if needed: */
    SPPlot_FixView(&(po->obs), 0.0, &(po->upp), NULL, 0.0, NULL, TRUE);

    /* Open the figure stream: */
    SPPlot_Stream *fps = SPPlot_NewStream
      ( po->eps, o->outName, po->paperSize, po->figSize, po->projLonLat, po->caption.ne );

    /* Compute mesh size in sphere units: */
    double relMeshSize = po->meshSize/(po->figSize/2);

    r4_gen_print(stderr, &(fnSupp.f), "%6.3f", "support = [ ", ", ", " ]\n");
    eMinObs =   MAXDOUBLE;
    eMaxObs = - MAXDOUBLE;

    if (po->radius < 1.0) { fprintf(stderr, "** \"-radius\" ignored\n"); }

    GradInfo GI[3]; SGrdInfo SI;
    int k;
    for (k = 0; k < 4; k++)
      {
        /* Choose the isoline spacing: */
        double rawMax = (k < 3 ? cgMaxErr.c[k] : sgMaxCos);
        fprintf(stderr, "max observed value = %16.10e\n", rawMax);
        double isoStep = SPPlot_RoundToNice(rawMax/DefaultIsolines);
        double isoMax = isoStep * DefaultIsolines;
        fprintf(stderr, "plot isoMax = %16.10e isoStep = %16.10e\n", isoMax, isoStep);

        auto double evalGradError(S2Point *p);

        double evalGradError(S2Point *p)
          { ComputeGradInfos(fn, p, o->delta, FALSE, GI, &SI);
            if (k < 3) 
              { return GI[k].error; }
            else
              { return SI.cosgp; }
          }

        pswr_new_canvas(fps, NULL);
        
        SPPlot_MultipleViews
          ( /* {fps} */          fps, 
            /* {funcTag} */      o->sfnName,
            /* {func} */         evalGradError,
            /* {tri} */          tri, 
            /* {relMeshSize} */  relMeshSize,
            /* {showTriang} */   (tri != NULL),
            /* {supp} */         &fnSupp,
            /* {fRange} */       isoMax,
            /* {fStep} */        isoStep,
            /* {obs} */          &(po->obs),
            /* {upp} */          &(po->upp),
            /* {rad} */          po->radius,
            /* {dLight} */       &(po->light),
            /* {lineWidth} */    po->lineWidth,
            /* {gridNLon} */     0,
            /* {gridNLat} */     0,
            /* {gridDots} */     FALSE,
            /* {aSide} */        +1, 
            /* {bSide} */        -1, 
            /* {caption} */      &(po->caption), 
            /* {index} */        0, 
            /* {time} */         0.0,
            /* {error} */        rawMax,
            /* {capAlign} */     0.5,
            /* {verbose} */      TRUE,
            /* {fMinObs} */      &eMinObs,
            /* {fMaxObs} */      &eMaxObs
          );
      }

    pswr_close_stream(fps);
    fprintf(stderr, "observed error extrema:");
    fprintf(stderr, " min = %14.7e max = %14.7e\n", eMinObs, eMaxObs);
    fprintf(stderr, "\n");
  }
  
GradInfo TestGrad(SPFunction *fn, S2Point *p, double delta, r3_t *dir) 
  { GradInfo GI;
    GI.p = (*p);
    GI.fp = fn->m->eval(fn, &(GI.p));
    GI.cg = fn->m->grad(fn, p);
    
    r3_scale(delta, dir, &(GI.dp));
    r3_add(&(GI.p), &(GI.dp), &(GI.q));
    GI.fq = fn->m->eval(fn, &(GI.q));
    r3_sub(&(GI.q), &(GI.p), &(GI.dp));
    GI.cgdp = r3_dot(&(GI.cg), &(GI.dp));
    GI.error = (GI.fq - (GI.fp + GI.cgdp))/delta;
    return GI;
  }
  
SGrdInfo TestSGrd(SPFunction *fn, S2Point *p) 
  { SGrdInfo SI;
    SI.p = (*p);
    SI.sg = SPFunction_SGrd(fn, p);
    SI.cosgp = r3_dot(&(SI.sg), &(SI.p))/(r3_norm(&(SI.sg)) + 1.0e-100);
    return SI;
  }


Triangulation *GetTriangulation(SPFunction *f)
  { SPSpline *fpw;
    if ((fpw = SPSpline_Cast(f)) != NULL) 
      { return fpw->d->tri; }
    else
      { return NULL; }
  }

SPH3_Plane GetSupportingPlane(SPFunction *f)
  { SPSpline *fpw;
    if ((fpw = SPSpline_Cast(f)) != NULL)
      { return fpw->d->supp; }
    else
      { fprintf(stderr, "** warning - no supporting plane - type = %s\n", f->type);
        return Omega;
      }
  }

Options GetOptions(int argn, char **argc)
  { Options o;
    SPOptions_Parser_t *pp = SPOptions_NewParser(stderr, argn, argc);

    SPOptions_SetUsage(pp, 
      PROG_NAME " \\\n"
      "  -sfnName NAME \\\n"
      "  -delta NUM [ -smpOrder NUM ] \\\n"
      "  -outName NAME \\\n"
      SPPlotOptions_FunctionHelp " \n"
    );

    if (SPOptions_TestKeyword(pp, "-sfnName"))
      { o.sfnName = SPOptions_GetNext(pp); }
    else
      { o.sfnName = "-"; }  

    if (SPOptions_TestKeyword(pp, "-smpOrder"))
      { o.smpOrder = SPOptions_GetNextInt(pp, 5, 1000); }
    else
      { o.smpOrder = 100; }

    SPOptions_GetKeyword(pp, "-delta");
    o.delta = SPOptions_GetNextDouble(pp, 0.0, 1.0);

    SPOptions_GetKeyword(pp, "-outName");                               
    o.outName = SPOptions_GetNext(pp);  

    o.plt = SPPlotOptions_FunctionParse(pp);
    
    SPOptions_Finish(pp);
    return o;
  }
  
