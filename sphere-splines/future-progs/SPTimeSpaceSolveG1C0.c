/* Solves a space-time differential equation on S^2 with g=3, c=1 time basis. */
/* Last edited on 2023-10-15 03:25:54 by stolfi */

#define PROG_NAME "SPTimeSpaceSolveG1C0"

/*
  This program solves a separable time-space differential equation
  on the sphere
  
    { D(f)(e,t) == TMap(f(e,t), e,t) }              (1)
    
  where {f} is the function to be determined, which depends on time
  {t} and a point {e} on the sphere; {D} is a linear operator that
  depends on {f} and its partial derivatives; and {TMap} is a given
  function. 
  
  An example is the heat diffusion equation on a rigidly rotating
  sphere, where the left-hand side is:
  
    { TD(f)(p) - F·SLap(f)(p) + R·V(e)¤SGrd(f)(p) - L·f(p) } (2)
  
  where 
  
    {p} stands for the pair {(e,t)}; 
    
    {TD} denotes the partial derivative with respect to time;
    
    {F} is the diffusion coefficient;
    {SLap} is the spherical Laplace-Beltrami operator;
    
    {R} is the angular speed;
    {V} is the velocity field on the sphere for unit angular speed; 
    {SGrd} is the spherical gradient operator;
    {¤} is scalar product in the sphere's tangent space; 
    
    {L} is the coefficient of spontaneous growth/decay.
  
  Here the right-hand side {TMap(f(p), p)} is the
  external heat input/loss.
  
  PETROV-GALERKIN APPROACH
  
  We look for an approximation {h} to the solution {f}, consisting of
  a linear combination of /finite space-time elements/,
  
    { h(p) = SUM { c[iju] · phi[iju](p) } } 
   
  where each {c[iju]} is a real coefficient (to be determined), and
  each finite element {phi[iju]} is a known function of bounded support
  in time and space. (The element index is written {iju} to emphasize
  that the element grid stretches both in space {i} and time {j},
  possibly with several elements {u} for each time step.)
  
  After Petrov and Galerkin, we look for the coefficients {c[iju]} that anihilate
  scalar products
  
    { <D(h) - H | theta[rsv]> }
    
  where {H} is the function {p -> TMap(h(p), p)}, and {theta} is
  another family of finite space-time elements, the /gauge functions/,
  which are `dense enough' in a sense that will become clear later on.
  Like the basis elements {phi[iju]}, each gauge function {theta[rsv]} is
  identified by a spatial (locus) index {r} and a temporal (epoch)
  index {s}, and there are {q} functions for each locus and epoch,
  identified by an index {v} in {0..q-1}.
  
  The Petrov-Galerkin criterion reduces to solving the  non-linear system
  
    { M c == b }        (3)
    
  where {M} is a fixed matrix
  
    { M[rsv,iju] == < D(phi[iju]) | theta[rsv] > }
        
  and {b} is a time-varying vector

    { b[rsv] == <H | theta[rsv]> }.  (4)
           
  Note that this system is not trivial since the right-hand side {b}
  depends on the approximate solution {h}, and therefore on the
  coefficient vector {c}, usually in a non-linear way. So (3) must be
  solved iteratively.
  
  Specifically, we consider a basis {phi} of the tensor-product type,
  derived from a basis {sigma[0..m-1]} of spatial elements and a basis
  {tau[0..n-1,0..q-1]} of temporal elements. The spatial basis {sigma}
  consists of spherical polynomial splines, homogeneous or general,
  defined on some irregular triangulation {T}. The temporal basis
  {tau} consists of univariate splines, with elements centered on {n}
  equally-spaced /epochs/, with {q} elements per epoch.
  
  Namely, the index {iju} is decomposed into a triple {(i,j,u)} where
  {i} in {0..m-1} selects a spatial element, {j} in {0..n-1} indicates
  an epoch, and {u} in {0..q-1} selects a temporal element for epoch
  {j}. We then define
  
    { phi[iju](p) = sigma[i](e) · tau[j,u](t) }  (5)
    
  Similarly, for the gauge functions {theta[rsv]} we use tensor
  products {rho[i] · pi[ju]} of some spatial (spherical) basis
  {rho[0..m-1]} and a temporal (univariate) basis {pi[0..n-1,0..q-1]}.
  
  One advantage of tensor-style bases is that the scalar-product
  integrals can be separated, so that
  
    { <phi[iju] | theta[rsv]> = <sigma[i] | rho[r]> · <tau[j,u] | pi[s,v]> }
  
  PROGRESSIVE SOLUTION
  
  System (3) is usually too big to solve at once. Fortunately the
  boundary conditions are of the initial-value type so we can
  contemplate solving (3) epoch by epoch in the time direction.
  
  To simplify the problem, we will assume that the temporal basis
  element {tau[j,u]} is at least C0 continuous and supported on two
  consecutive time steps, centered on epoch {j}. That is, its support
  is the interval from {tj-tStep} through {tj+tStep}, where {tj} is
  the time of epoch {j} and {tStep} is the time step. We also assume
  that each element {pi[s,v]} is supported on a single time step, from
  epoch {s-1} to epoch {s}.
  
  With this assumption, each equation of (3) involves coefficients
  from only two consecutive epochs, {c[i,j-1,u]} and {c[i,j,u]}.
  Moreover, since the space-time elements (approximators and gauges)
  in each epoch are time-translates of those in the previous epoch,
  the coefficients of {M} repeat epoch after epoch. 
  
  More precisely, {M} has a periodic block-band structure
    {
        N[1]  N[0]  0     0     0     0     ... 
        
        ...   ...   ...   ...   ...   ...   ...
           
        ...   N[1]  N[0]  0     0     0     ...
        
        ...   0     N[1]  N[0]  0     0     ...
        
        ...   0     0     N[1]  N[0]  0     ...
  
        ...   ...   ...   ...   ...   ...   ...
 
    }
    
  where each block {N[k]}, for {k=0..1}, is the { m·q × m·q } matrix
  of the scalar products of the differentiated basis elements in epoch
  {j-k} against the gauge functions in epoch {j}, i.e.
  
    { N[k][rv,iu] == <D(phi[i,j-k,u]) | theta[r,j,v]> }

  (Note that the block {N[k]} does not depend on {j}.) Therefore, at
  each stage of this computation, we advance by one time step.
  Therefore, we can solve the system above one epoch at a time. At
  each stage {j} we solve the equation

    { N[0] a[0] == d - N[1] a[1] }  (6)
  
  where each {a[k]}, called a /frame/, is the subvector of {m·q}
  coefficients from {c} relative to epoch {j-k}, i.e. {a[k][iu] =
  c[i,j-k,u]}, and {d} is the subvector of {b} that corresponds to the
  gauge functions of epoch {j}, that is {d[rv] = b[r,j,v]}.
  
  To start process (6), the client must specify the first frame {a[0]}
  associated with the initial epoch {j = 0}. Then, at each stage we
  update the current state by doing {a[1] = a[0]}, {j = j+1}, and
  we compute the current frame {a[0]} for epoch {j} from the known
  previous frame {a[1]} of epoch {j-1}.

  If the time basis elements are continuous of order {z}, then each
  frame {a[0]} associated with an epoch {j} defines the value of the
  solution, and all time derivatives to order {z}, at the epoch's time
  {tj}. Thus, given only the frame {a[0]}, we can extract from it {z}
  space-only spherical functions that shows the solution at time {tj}
  and its time derivatives up to order {z}.
*/

#define _GNU_SOURCE

#include <SPTriang.h>
#include <SPMatrix.h>
#include <SPVector.h>
#include <SPTriang.h>
#include <SPFunction.h>
#include <SPSpline.h>
#include <SPIntegral.h>
#include <SPTimeSpaceProcFunction.h>
#include <SPApprox.h>
#include <SPTimeSpaceFuncMap.h> 
#include <SPOptions.h>
#include <SPPlot.h> 
#include <SPPlotOptions.h> 
#include <SP1DSpline.h>
#include <SPRange.h>
#include <SPH3.h>
#include <SPTimeSpaceSolveG1C0Aux.h>

#include <udg_pulse.h>
#include <rn.h>
#include <affirm.h>
#include <nat.h>
#include <bool.h>
#include <filefmt.h>
#include <fget.h>
#include <nget.h>

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <values.h>
#include <limits.h>

#define PROG_USAGE \
    PROG_NAME " \\\n" \
    "  -rhsName NAME \\\n" \
    "  -cDiffuse NUM -cDrift NUM -cDecay NUM \\\n" \
    "  -iniName NAME -solName NAME \\\n" \
    "  -timeIni NUM -timeStep NUM -nSteps NUM\\\n" \
    "  -spcBasis NAME [ -matName NAME ] \\\n" \
    SPSys_LinOptionsHelp " \\\n" \
    SPSys_GenOptionsHelp " \\\n" \
    "  [-minVal NUM] [ -smpOrder SPACEORD TIMEORD ] \\\n" \
    "  [ -verbose ] \\\n" \
    "  -outName NAME \\\n" \
    "  [ -plotEvery NUM ] [ -writeEvery NUM ] [ -printEvery NUM ] \\\n" \
    SPPlotOptions_FunctionHelp " \n"

typedef struct Options
  { char *outName;          /* Prefix for output filenames. */
    /* Parameters of the problem: */
    char *rhsName;          /* Name of right-hand-side function {TMap}. */
    double cDiffuse;        /* Diffusion coeff {F} of equation. */
    double cDrift;          /* Angular speed {R} of drift term. */
    double cDecay;          /* Decay coeff {L} of equation. */
    char *iniName;          /* Name of initial state coeff file. */
    char *solName;          /* Name of true solution. */
    /* The spatial aproximation basis: */
    char *spcBasis;         /* Prefix for spatial basis filenames. */
    char *matName;          /* Prefix for spatial matrix filenames (default {spcBasis}). */
    /* The temporal aproximation basis: */
    double timeIni;         /* Starting time. */
    double timeStep;        /* Ending time. */
    double nSteps;          /* Number of time steps. */
    /* Parameters for PDE integration: */
    bool_t verbose;         /* TRUE prints right-hand-side vectors, etc.. */
    int smpOrderSpace;      /* Sampling order for Gaussian quadrature in space. */
    int smpOrderTime;       /* Number of samples for Gaussian quadrature in time. */
    double minVal;          /* Threshold for matrix factor cleanup. */ 
    SPSys_LinOptions_t lso; /* Parameters for linear system solving */
    SPSys_GenOptions_t gso; /* Parameters for solving the non-linear system. */
    /* Plotting/writing options: */
    int printEvery;         /* Print soln and error every this many frames (0 = none). */
    int plotEvery;          /* Plot soln and error every this many frames (0 = none). */
    int writeEvery;         /* Write one every this many frames. */
    SPPlotOptions_t plt;    /* Plot format and style options. */
  } Options;

Options *GetOptions(int argn, char **argc);

SPTimeSpaceFunction *GetTrueSolution(char *solName);
  /* Reads the true solution from file "{solName}.tfn". */

/* ODE INTEGRATION */

/* In the following procedures, frames {a[0..w-1]} are interpreted
  relative to the space-time basis {stbas}. The arguments {N[0..w-1]} are
  the matrices of system (6), as explained under {GetMatrices}. If
  {sol} is not NULL, the space component defined by frame {a[0]} is
  compared with {sol(timeIni)}. */

void LoadInitialFrame(char *frmName, SPVector a, bool_t verbose);
  /* Reads, from file "{frmName}.cof", a frame of the initial state
    for some epoch {j}. Stores the frame into the vector {a} (that is,
    sets {a[iu] = c[i,j,u]} for all {i,u}). */

char *FrameName(char *outName, int epoch);
  /* Returns a name for the frame number {epoch}, "{outName}-{epoch}".
    However the name is "{outName}-ini" if {epoch} is -1, and 
    "{outName}-fin" if {epoch} is {MAXINT}. */

void OutputFrame
  ( SPPlot_Stream *fps,
    Options *o, 
    int j,
    double tj,
    double tStep,
    SPVector a, 
    STBasis *stbas, 
    Triangulation *tri,
    SPTimeSpaceFunction *sol,
    bool_t writeIt,
    bool_t printIt,
    bool_t plotIt
  );
  /* Outputs the frame {a} as appropriate. Assumes that {a} is
    centered at epoch {j}, i.e. its support is the interval {tj -
    tStep} to {tj + tStep}.
    
    If {writeIt} is true, the frame {a} is written to disk as
    "{o.outName}-{j}.cof". If {printIt} is true, a snapshot of the
    function at epoch {j} (time {tj}) is evaluated at selected points
    of the sphere, compared to the true solution {sol}, and the
    results are printed to {stderr}. If {plotIt} is true, that
    snapshot is plotted to the figure stream {fps}.  */
  
void PrintMaxErrorValues
  ( int j,
    double tj,
    SPFunction *app,
    SPTimeSpaceFunction *sol,
    Triangulation *tri, 
    double *gMax, 
    double *sMax, 
    double *eMax, 
    double *eAvg
  );
  /* Computes and prints the maximum absolute values {gMax,sMax} of
    the static approximation {app}, and of the true solution {sol} at
    time {tj}, respectively; and the maximum and root-mean-square error
    {eMax,eAvg} between the approximation and {sol}. */

SPPlot_Stream *GetPlotStream(char *outName, SPPlotOptions_t *po);
  /* Opens an output Postscript stream with name prefix {outName}
    configured according to the options {po}. */

int main(int argn, char **argc)
  { Options *o = GetOptions(argn, argc);
    SPPlotOptions_t *po = &(o->plt);
    double tStep = o->timeStep;
    SPIntegral_SetDefaultSamplingOrder(o->smpOrderSpace);

    int smpCheckSpace = 40; /* Space sampling order when checking error. */
    int smpCheckTime  = 10; /* Time samples per step when checking error. */

    /* Get spatial basis: */
    Basis sBas = SPApprox_ReadBasis(txtcat(o->spcBasis, ".bas"));
    int sDim = sBas.ne;    /* Dimension {n} of space basis. */
    Triangulation *tri = SPSpline_BasisTriangulation(sBas);

    /* Get parameters of temporal basis (Hermite-style 2-step): */
    udg_pulse_family_t tFam = udg_pulse_family(udg_PK_H, 0, 1);
    int tDim = tFam.nmp;  /* Number {q} of temporal elements per epoch. */
    affirm(tDim == 1, "wrong num of pulses per epoch");
    affirm(udg_pulse_max_supp_count(&tFam) == 2, "wrong time pulse width");
    
    /* Compose the spatio-temporal basis: */
    STBasis stbas = (STBasis) { sBas, tFam };
    int frameDim = tDim * sDim; /* Num coeffs in each frame. */
    
    SPTimeSpaceFuncMap RHS = SPTimeSpaceFuncMap_FromName(o->rhsName);

    SPTimeSpaceFunction *sol = GetTrueSolution(o->solName);
    
    SPVector a[2]; /* {a[k]} is the coeff vector of frame {j - k}. */
    int k;
    for (k = 0; k < 2; k++) { a[k] = double_vec_new(frameDim); }
  
    /* Obtain the system's matrices: */
    SPMatrix N[2];  /* The matrices {N[0..2-1]}. */ 
    SPTimeSpaceG1C0_GetMatrices
      ( o->matName, 
        o->cDiffuse, o->cDrift, o->cDecay, 
        tStep, &stbas, N, 
        o->verbose
      );

    /* Obtain the factors of {N[0]}, if necessary: */
    SPMatrix N0L, N0R; /* Left and right factors of {N[0]}. */ 
    SPVector N0D;
    SPSys_ComputeNeededFactors(N[0], o->lso.mth, o->minVal, &N0L, &N0D, &N0R);
    
    FILE *errWr = open_write(txtcat(o->outName, ".erp"), TRUE);
    fprintf(errWr, "# %6s %22s   %22s %22s   %22s %22s\n",  
      "epoch", "time", 
      "max(app-sol)", "rms(app-sol)", 
      "max(resid)",   "rms(resid)"
    );

    /* Fix perspective parameters, if needed: */
    SPPlot_FixView(&(po->obs), 0.0, &(po->upp), NULL, 0.0, NULL, TRUE);

    /* Open the figure stream: */
    SPPlot_Stream *fps = GetPlotStream(o->outName, po);
    
    int j = 0; /* Index of current epoch. */

    double errMax = 0; /* Maximum error {app(t)-sol(t)} seen. */
    double errSum = 0; /* Sum of mean square error over all steps. */
    
    double resMax = 0; /* Maximum residue {lhs(t)-rhs(t)} seen. */
    double resSum = 0; /* Sum of mean square residue over all steps. */
    
    while(TRUE)
      { 
        /* Main iteration: compute (or read) the frame {a[0]} for epoch {j},
          that is, {a[0][iu] = c[i,j,u]}. */
        
        /* Compute the epoch's time {tj}: */
        double tj = o->timeIni + j*tStep;
        fprintf(stderr, "### epoch t[%d] = %8.2f #################\n", j, tj);
        
        if (j == 0)
          { /* Load frame {j} of the initial state: */
            LoadInitialFrame(o->iniName, a[0], FALSE);
          }
        else
          { /* At this point we know the previous frame {a[1]}, that is, the
              coefficients {a[1][iu] = c[i,j-1,u]}. Compute the next frame
              {a[0]} from {a[1]} by solving the Galerking constraints in the 
              interval between {tj-tStep} and {tj}: */
            SPTimeSpaceG1C0_ComputeNextFrame
              ( tj, tStep, a, N, N0L, N0D, N0R,
                RHS, &stbas, tri, 
                o->smpOrderTime,
                &(o->lso), &(o->gso), o->verbose
              );
            /* Check error in {[tj-tStep __ tj]}: */
            double errMaxStep = 0, errAvgStep = 0;
            double resMaxStep = 0, resAvgStep = 0;
            SPTimeSpaceG1C0_CheckSolution
              ( a, tj, tStep, &stbas, tri,
                o->cDiffuse, o->cDrift, o->cDecay, RHS, sol,
                smpCheckTime, smpCheckSpace,
                &errMaxStep, &errAvgStep,
                &resMaxStep, &resAvgStep
              );
            fprintf(errWr, "  %6d %22.14e   %22.14e %22.14e   %22.14e %22.14e\n",  
              j, tj, errMaxStep, errAvgStep, resMaxStep, resAvgStep
            );
            /* Update global error and square sum: */
            if (errMaxStep > errMax) { errMax = errMaxStep; }
            if (resMaxStep > resMax) { resMax = resMaxStep; }
            errSum += errAvgStep*errAvgStep;
            resSum += resAvgStep*resAvgStep;
          }
        if (o->verbose)
          { int i;
            fprintf(stderr, "  a[0] = (\n");
            for (i = 0; i < a[0].ne; i++) 
              { fprintf(stderr, "    %22.14e\n", a[0].e[i]); }
            fprintf(stderr, " )\n");
          }

        /* Output the new frame {a[0]} if appropriate: */
        bool_t wrt = ((j == 0) || (j == o->nSteps) || (j % o->writeEvery == 0));
        bool_t prn = ((j == 0) || (j == o->nSteps) || (j % o->printEvery == 0));
        bool_t plt = ((j == 0) || (j == o->nSteps) || (j % o->plotEvery == 0));
        OutputFrame(fps, o, j, tj, tStep, a[0], &stbas, tri, sol, wrt, prn, plt);

        /* Are we done? */
        if (j >= o->nSteps) { break; }
        
        /* Prepare for next iteration: */
        { SPVector atmp = a[1]; a[1] = a[0]; a[0] = atmp; }
        j++;
      }
      
    fprintf(stderr, "### end of integration ###################\n"); 
    /* Force output of final frame, with name "{o.outName}-fin": */
    double timeFin = o->timeIni + o->nSteps*tStep;
    OutputFrame(fps, o, INT_MAX, timeFin, tStep, a[0], &stbas, tri, sol, TRUE, TRUE, TRUE);
    
    /* Report max and average error: */
    double errAvg = sqrt(errSum/o->nSteps);
    fprintf(stderr, "errMax = %+9.1e errAvg = %+9.1e\n", errMax, errAvg);

    /* Report max and average residue: */
    double resAvg = sqrt(resSum/o->nSteps);
    fprintf(stderr, "resMax = %+9.1e resAvg = %+9.1e\n", resMax, resAvg);
    
    fclose(errWr);
    pswr_close_stream(fps);
    return 0;
  }
  
void LoadInitialFrame(char *frmName, SPVector a, bool_t verbose)
  { 
    char *fileName = jsprintf("%s.cof", frmName);
    FILE *rd = open_read(fileName, TRUE);
    SPVector ar = SPVector_Read(rd);
    affirm(ar.ne == a.ne, "size mismatch");
    int i;
    for (i = 0; i < a.ne; i++) { a.e[i] = ar.e[i]; }
    free(ar.e);
    free(fileName);
  }
  
SPPlot_Stream *GetPlotStream(char *outName, SPPlotOptions_t *po)
  {
    return SPPlot_NewStream
      ( po->eps, outName, po->paperSize, po->figSize, po->projLonLat, po->caption.ne );
  }

void OutputFrame
  ( SPPlot_Stream *fps,
    Options *o, 
    int j,
    double tj,
    double tStep,
    SPVector a, 
    STBasis *stbas, 
    Triangulation *tri,
    SPTimeSpaceFunction *sol,
    bool_t writeIt,
    bool_t printIt,
    bool_t plotIt
  )
  { char *frmName = FrameName(o->outName, j);
    SPPlotOptions_t *po = &(o->plt);
    
    /* Write frame coefficients to disk: */
    if (writeIt)
      { char *fileName = NULL;
        char *fileName = jsprintf("%s.cof", frmName);
        FILE *wr = open_write(fileName, TRUE);
        SPVector_Write(wr, a);
        fclose(wr);
        free(fileName);
      }

    if (printIt || plotIt) 
      {
        /* Generate snapshot of computed sol at the frame's central epoch {tj}: */
        SPFunction *g = NULL; 

        fprintf(stderr, "building snapshot for t = %.4f...\n", tj); 
        g = SPTimeSpaceG1C0_TimeSliceSingle(a, tj, tStep, stbas, tri, 0, FALSE);

        /* Compare with true solution, get max values and errors: */
        double gMax, sMax, eMax, eAvg;
        fprintf(stderr, "summary for epoch %d\n", j); 
        PrintMaxErrorValues(j, tj, g, sol, tri, &gMax, &sMax, &eMax, &eAvg);
        /* fprintf(errWr, "  %6d %6d %16.12f\n",  j, iter, eAvg); */
        if (plotIt)
          { 
            /* Choose plot range for funcs: */
            double gsMax = (gMax > 2.0*sMax ? gMax : sMax);
            double fMax = po->autoRange*gsMax + (1 - po->autoRange)*po->fRange;
          
            /* Get plot mesh size in scene units: */
            double relMeshSize = po->meshSize/(po->figSize/2);

            auto double solt(S2Point *s); /* Evaluates {sol} at time {tj}. */
            double solt(S2Point *s)
              { return sol->m->eval(sol, s, tj); }

            if (po->radius < 1.0) { fprintf(stderr, "** \"-radius\" ignored\n"); }
            pswr_new_canvas(fps, NULL);
            char *frmTag = fmt_int(j, 6);
            SPApprox_PlotFunctionAndError
              ( fps, frmTag, g, solt, 
                tri, relMeshSize, (tri != NULL), (! po->eps),
                fMax, eMax,
                &(po->obs), &(po->upp),
                &(po->caption),
                j, tj
              );
            free(frmTag);
          }
      }

    free(frmName);
  }
  
char *FrameName(char *outName, int epoch)
  { char *name = NULL;
  if (epoch < 0) 
    { char *name = jsprintf("%s-ini", outName); }
  else if (epoch == INT_MAX)
    { char *name = jsprintf("%s-fin", outName); }
  else
    { char *name = jsprintf("%s-%06d", outName, epoch); }    
    return name;
  }

SPTimeSpaceFunction *GetTrueSolution(char *solName)
  { SPTimeSpaceFunction *f = (SPTimeSpaceFunction *)SPTimeSpaceProcFunction_FromName(solName);
    if (f == NULL) 
      { fprintf(stderr, "Unknown solution = \"%s\"\n", solName);
        affirm(FALSE, "aborted");
      }
    return f;
  }

Options *GetOptions(int argn, char **argc)
  { Options *o = notnull(malloc(sizeof(Options)), "no mem");
    SPOptions_Parser_t *pp = SPOptions_NewParser(stderr, argn, argc);

    SPOptions_SetUsage(pp, PROG_USAGE);

    SPOptions_GetKeyword(pp, "-rhsName");                               
    o->rhsName = SPOptions_GetNext(pp);  

    SPOptions_GetKeyword(pp, "-cDiffuse");                               
    o->cDiffuse = SPOptions_GetNextDouble(pp, -DBL_MAX, +DBL_MAX);  

    SPOptions_GetKeyword(pp, "-cDrift");                               
    o->cDrift = SPOptions_GetNextDouble(pp, -DBL_MAX, +DBL_MAX);  

    SPOptions_GetKeyword(pp, "-cDecay");                               
    o->cDecay = SPOptions_GetNextDouble(pp, -DBL_MAX, +DBL_MAX);  

    SPOptions_GetKeyword(pp, "-iniName");                               
    o->iniName = SPOptions_GetNext(pp);  

    SPOptions_GetKeyword(pp, "-solName");                               
    o->solName = SPOptions_GetNext(pp);  

    SPOptions_GetKeyword(pp, "-timeIni");
    o->timeIni = SPOptions_GetNextDouble(pp, -DBL_MAX, DBL_MAX);

    SPOptions_GetKeyword(pp, "-timeStep");
    o->timeStep = SPOptions_GetNextDouble(pp, 1.0e-10, 1.0e+10);

    SPOptions_GetKeyword(pp, "-nSteps");
    o->nSteps = SPOptions_GetNextInt(pp, 1, 1000000);

    SPOptions_GetKeyword(pp, "-spcBasis");                               
    o->spcBasis = SPOptions_GetNext(pp);  
       
    if (SPOptions_TestKeyword(pp, "-matName"))
      { o->matName = SPOptions_GetNext(pp); }
    else
      { o->matName = o->spcBasis; }

    if (SPOptions_TestKeyword(pp, "-minVal"))
      { o->minVal = SPOptions_GetNextDouble(pp, 0.0, DBL_MAX); }
    else 
      { o->minVal = 0.0; }
    
    /* Parameters for the linear equation system solver: */
    o->lso = SPSys_LinOptionsParse(pp, FALSE);
    if (o->lso.mth == SPSys_LM_NONE)
      { SPOptions_Error(pp, "must specify a linear solution method"); }
    else if (o->lso.mth == SPSys_LM_Cholesky) 
      { SPOptions_Error(pp, "Cholesky not allowed"); }

    /* Parameters for the global (non-linear) equation system solver: */
    o->gso = SPSys_GenOptionsParse(pp, FALSE);
      
    SPOptions_GetKeyword(pp, "-smpOrder");
    o->smpOrderSpace = SPOptions_GetNextInt(pp, 1, INT_MAX);
    o->smpOrderTime = SPOptions_GetNextInt(pp, 1, INT_MAX);

    SPOptions_GetKeyword(pp, "-outName");                               
    o->outName = SPOptions_GetNext(pp);  
    
    o->verbose = SPOptions_TestKeyword(pp, "-verbose");
                                                 
    /* Plotting options: */

    if (SPOptions_TestKeyword(pp, "-plotEvery"))
      { o->plotEvery = SPOptions_GetNextInt(pp, 0, INT_MAX); }
    else
      { o->plotEvery = (o->nSteps + 19)/20; }

    if (SPOptions_TestKeyword(pp, "-writeEvery"))
      { o->writeEvery = SPOptions_GetNextInt(pp, 0, INT_MAX); }
    else
      { o->writeEvery = (o->nSteps + 199)/200; }

    if (SPOptions_TestKeyword(pp, "-printEvery"))
      { o->printEvery = SPOptions_GetNextInt(pp, 0, INT_MAX); }
    else
      { o->printEvery = (o->nSteps + 199)/200; }
    
    o->plt = SPPlotOptions_FunctionParse(pp);
    
    SPOptions_Finish(pp);
    return o;
  }

void PrintMaxErrorValues
  ( int j,
    double tj,
    SPFunction *app,
    SPTimeSpaceFunction *sol,
    Triangulation *tri, 
    double *gMax, 
    double *sMax, 
    double *eMax, 
    double *eAvg
  )
  {
    auto double solValue(S2Point *p);
    auto double appValue(S2Point *p);
    auto double error(S2Point *p);
    auto double errorSqr(S2Point *p);
    
    double solValue(S2Point *p) { return sol->m->eval(sol, p, tj); }

    double appValue(S2Point *p) { return app->m->eval(app, p); }

    double error(S2Point *p) 
      { return app->m->eval(app, p) - sol->m->eval(sol, p, tj); }

    double errorSqr(S2Point *p) 
      { double d = app->m->eval(app, p) - sol->m->eval(sol, p, tj);
        return d*d;
      }

    (*gMax) = SPRange_OnSphere(appValue, 120);
    fprintf(stderr, "max(fabs(app)) = %9.1e\n", *gMax);
    (*sMax) = SPRange_OnSphere(solValue, 120);
    fprintf(stderr, "max(fabs(sol)) = %9.1e\n", *sMax);
    (*eMax) = SPRange_OnSphere(error, 120);
    fprintf(stderr, "max(fabs(app-sol)) = %9.1e\n", *eMax);
    (*eAvg) = sqrt(SPIntegral_OnSphere(errorSqr, tri)); 
    fprintf(stderr, "rms(app-sol) = %9.1e\n", *eAvg);
  }

