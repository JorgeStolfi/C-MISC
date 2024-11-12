/* Checks whether the {grad} and {eval} methods are consistent. */
/* Last edited on 2024-11-11 05:14:06 by stolfi */

#define PROG_NAME "SPHessTest"

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
    double delta;         /* Size of perturbation for numerical hessian. */
    char *outName;        /* Prefix for output file names. */
    int smpOrder;         /* Sampling order for test points. */
    /* Plotting options: */
    SPPlotOptions_t plt;  /* Plot format and style options. */
  } Options;
    
typedef struct HessInfo 
  { S2Point p;         /* A point on the sphere. */
    double fp;         /* Function value at {p}. */
    R3Gradient g;      /* The {R^3} gradient at {p}. */
    R3Hessian H;       /* The {R^3} Hessian at {p}. */
    r3_t dp;           /* A tangent perturbation for {p}. */
    r3_t q;            /* The perturbed point {p + du + dv}. */
    double fq;         /* Function value at {q}. */
    double Tq;         /* Quadratic Taylor form at {p} evaluated at {q}. */
    double error;      /* Difference between {Tq} and {fq}. */
  } HessInfo;

Options GetOptions(int argn, char **argc);

void ComputeHessInfos(SPFunction *fn, S2Point *p, double delta, bool_t tg, HessInfo *HI);
  /* Tests the Hessian method of {fn} with six displacements of at
    {p}, returns the data in {HI[0..5]}. If {tg} is TRUE, the
    displacements are {delta}-long vectors tangential to the sphere,
    else they are all 1- and 2-axis {delta} steps. */
   
void CheckHessError(SPFunction *fn, Options *o, r6_t *maxErr);
void PlotHessError(SPFunction *fn, Options *o, r6_t maxErr);

HessInfo TestHessian(SPFunction *fn, S2Point *p, double delta, r3_t *dir);
void PrintHessInfo(FILE *wr, HessInfo *HI);
Triangulation *GetTriangulation(SPFunction *f);
SPH3_Plane GetSupportingPlane(SPFunction *f);

int main(int argn, char **argc)
  { Options o = GetOptions(argn, argc);
    SPIntegral_SetDefaultSamplingOrder(o.smpOrder);
    SPFunction *fn = SPApprox_ReadFunction(txtcat(o.sfnName, ".sfn"));
    r6_t maxErr;
    CheckHessError(fn, &o, &maxErr);
    PlotHessError(fn, &o, maxErr);
    return 0;
  }

void ComputeHessInfos(SPFunction *fn, S2Point *p, double delta, bool_t tg, HessInfo *HI)
  { 
    int k;
    r3_t v[6];
    if (tg)
      {
        /* Generate two orthogonal vectors in the tangent space at {p}: */
        r3_t u[2];
        (void)r3_throw_ortho_pair(p, &(u[0]), &(u[1]));
        /* Generate six equally spaced vectors in the tangent space {p}: */
        double cos30 = sqrt(3.0)/2.0;
        double sin30 = 0.5;
        v[0] = u[0]; 
        r3_mix(+cos30, &(u[0]), +sin30, &(u[1]), &(v[1]));
        r3_mix(+sin30, &(u[0]), +cos30, &(u[1]), &(v[2]));
        v[3] = u[1];
        r3_mix(-sin30, &(u[0]), +cos30, &(u[1]), &(v[4]));
        r3_mix(-cos30, &(u[0]), +sin30, &(u[1]), &(v[5]));
      }
    else
      { int i, j;
        k = 0;
        for (j = 0; j < 3; j++)
          for (i = 0; i <= j; i++)
            { 
              v[k] = (r3_t){{ 0.0, 0.0, 0.0 }};
              v[k].c[i] = 1.0; 
              v[k].c[j] = 1.0;
              k++;
            }
      }
    /* Test the Hessian with those displacements: */
    for (k = 0; k < 6; k++) 
      { HI[k] = TestHessian(fn, p, delta, &(v[k])); }
  }

void CheckHessError(SPFunction *fn, Options *o, r6_t *maxErr)
  { int i, k;
    r6_t totErr2;           /* Sum of error squared in each Hessian element. */
    int numErr = 0;         /* Number of errors summed in each Hessian element. */
    HessInfo maxErrHI[6];   /* Point of maximum error in each Hessian element. */
    S2Point_vec_t sp = S2Point_vec_new(0); /* Sample points. */
    
    /* Generate sample points: */
    { double_vec_t wp = double_vec_new(0); /* Sample weights for integration. */
      int ns = 0;
      SPIntegral_SampleOctants(o->smpOrder, &sp, &wp, &ns);
      S2Point_vec_trim(&sp, ns);
      free(wp.e);
    }
    
    /* Get value, gradient and Hessian at sample points and directions: */
    (*maxErr) = (r6_t){{ -1.0, 1.0, -1.0, -1.0, -1.0, -1.0 }};
    totErr2 = (r6_t){{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }};
    HessInfo HI[6];
    for (i = 0; i < sp.ne; i++)
      { r3_t p = sp.e[i];
        ComputeHessInfos(fn, &p, o->delta, FALSE, HI);
        for (k = 0; k < 6; k++)
          { double error = HI[k].error;
            totErr2.c[k] += error*error;
            if (fabs(error) > (*maxErr).c[k])
              { (*maxErr).c[k] = fabs(error); maxErrHI[k] = HI[k]; }
          }
        numErr++;
      }
    
    for (k = 0; k < 6; k++)
      { double avgErr = sqrt(totErr2.c[k]/((double)numErr));
        fprintf(stderr, "hessian component %d:", k);
        fprintf(stderr, " avg error = %14.7e\n", avgErr);
        fprintf(stderr, " max error = %14.7e\n", (*maxErr).c[k]);
        fprintf(stderr, "\n");
        PrintHessInfo(stderr, maxErrHI);
      }
  }
  
void PrintHessInfo(FILE *wr, HessInfo *HI)
  {
    fprintf(wr, "  p =    "); r3_print(wr, &(HI->p)); 
    fprintf(wr, "  f(p)  = %16.12f\n", HI->fp);
    
    fprintf(wr, "  grad = "); r3_print(wr, &(HI->g)); 
    fprintf(wr, "\n");
    
    fprintf(wr, "  hess = "); r6_print(wr, &(HI->H)); 
    fprintf(wr, "\n");
    
    fprintf(wr, "  dp =   "); r3_print(wr, &(HI->dp)); 
    fprintf(wr, "\n");
    
    fprintf(wr, "  q =    "); r3_print(wr, &(HI->q)); 
    fprintf(wr, "  f(q)  = %16.12f\n", HI->fq);
    fprintf(wr, "\n");
    
    fprintf(wr, "  Tq = %16.12f\n", HI->Tq);
    fprintf(wr, "\n");
    fprintf(wr, "  error = %14.7e\n", HI->error);
  }

void PlotHessError(SPFunction *fn, Options *o, r6_t maxErr) 
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
    eMinObs =   INFINITY;
    eMaxObs = - INFINITY;

    if (po->radius < 1.0) { fprintf(stderr, "** \"-radius\" ignored\n"); }

    HessInfo HI[6];
    int k;
    for (k = 0; k < 6; k++)
      {
        /* Choose the isoline spacing: */
        double rawMax = maxErr.c[k];
        fprintf(stderr, "max observed value = %16.10e\n", rawMax);
        double isoStep = SPPlot_RoundToNice(rawMax/DefaultIsolines);
        double isoMax = isoStep * DefaultIsolines;
        fprintf(stderr, "plot isoMax = %16.10e isoStep = %16.10e\n", isoMax, isoStep);


        auto double evalHessError(S2Point *p);

        double evalHessError(S2Point *p)
          { ComputeHessInfos(fn, p, o->delta, FALSE, HI);
            return HI[k].error; 
          }

        pswr_new_canvas(fps, NULL);

        SPPlot_MultipleViews
          ( /* {fps} */          fps, 
            /* {funcTag} */      o->sfnName,
            /* {func} */         evalHessError,
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
  
HessInfo TestHessian(SPFunction *fn, S2Point *p, double delta, r3_t *dir) 
  { HessInfo HI;
    HI.p = *p;
    HI.fp = fn->m->eval(fn, &(HI.p));
    HI.g = fn->m->grad(fn, p);
    HI.H = fn->m->hess(fn, p);
    
    /* Generate the displaced point and evaluate {f} there: */
    r3_scale(delta, dir, &(HI.dp));
    r3_add(&(HI.p), &(HI.dp), &(HI.q));
    HI.fq = fn->m->eval(fn, &(HI.q));
    
    /* Evaluate the quadratic form: */
    HI.Tq = 
      HI.fp + 
      r3_dot(&(HI.g), &(HI.dp)) + 
      HI.H.c[0]*HI.dp.c[0]*HI.dp.c[0]/2 + 
      HI.H.c[1]*HI.dp.c[0]*HI.dp.c[1] +
      HI.H.c[2]*HI.dp.c[1]*HI.dp.c[1]/2 + 
      HI.H.c[3]*HI.dp.c[0]*HI.dp.c[2] + 
      HI.H.c[4]*HI.dp.c[1]*HI.dp.c[2] +
      HI.H.c[5]*HI.dp.c[2]*HI.dp.c[2]/2;
    
    HI.error = (HI.fq - HI.Tq)/r3_norm_sqr(&HI.dp);
    return HI;
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
  
