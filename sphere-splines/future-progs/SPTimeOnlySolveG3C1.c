/* Solves a time-only differential equation with g=3, c=1 basis. */
/* Last edited on 2023-02-12 07:55:45 by stolfi */

#define PROG_NAME "SPTimeOnlySolveG3C1"

/*
  This program solves a time-domain ordinary differential equation
  
    { D(f)(t) == TMap(f(t), t) }              (1)
    
  where {f} is the function to be determined, which depends on time
  {t}; {D} is a linear operator that depends on {f} and its
  derivatives; and {TMap} is a given function.
  
  An example is the damped harmonic oscillator equation, where the 
  left-hand side is:
  
    { M·TDD(f)(t) + R·TD(f)(t) + K·f(t) } (2)
  
  where 
  
    {TD} denotes single derivative with respect to time;
    {TDD} denotes double derivative with respect to time;
    
    {M} is the mass;
    {R} is the friction coefficient;
    {K} is the spring coefficient;
  
  Here the right-hand side {TMap(f(t), t)} is the
  external driving force.
  
  PETROV-GALERKIN APPROACH
  
  We look for an approximation {h} to the solution {f}, consisting of
  a linear combination of /finite time elements/,
  
    { h(p) = SUM { c[j,u] · tau[j,u](t) } } 
   
  where each {c[j,u]} is a real coefficient (to be determined), and
  each finite element {tau[j,u]} is a known function of bounded support
  in time. Each element is associated with an /epoch/ (special time)
  identified by an index {j} in {0..n-1}, but there may be several
  elements for the same epoch, distinguished by an extra index {u} in
  {0..q-1}.
  
  After Petrov and Galerkin, we look for the coefficients {c[j,u]}
  that anihilate scalar products
  
    { <D(h) - H | pi[s,v]> }
    
  where {H} is the function {t -> TMap(h(t), t)}, and {pi}
  is another family of finite space-time elements, the /gauge
  functions/, which are `dense enough' in a sense that will become 
  clear later on. Like the basis elements {tau}, each gauge function
  {pi[s,v]} is associaeted with some epoch {s} in {0..n-1},
  and there are {q} functions per epoch, identified by 
  the index {v} in {0..q-1}.
  
  The Galerkin criterion reduces to solving the  non-linear system
  
    { M c == b }        (3)
    
  where {M} is a fixed matrix
  
    { M[ju,sv] == < D(tau[j,u]) | pi[s,v] > }
        
  and {b} is a time-varying vector

    { b[sv] == <H | pi[s,v]> }.  (4)
           
  Here we use the convention that {ju = j*q + u} and {sv = s*q + v}.
  Note that this system is not trivial since the right-hand side {b}
  depends on the approximate solution {h}, and therefore on the
  coefficient vector {c}, usually in a non-linear way. So (3) must be
  solved iteratively.

  PROGRESSIVE SOLUTION
  
  System (3) is usually too big to solve at once. Fortunately the
  boundary conditions are of the initial-value type so we can
  contemplate solving (3) epoch by epoch in the time direction.
  
  To simplify the problem, we will assume that the temporal basis
  element {tau[j,u]} is at least C1 continuous and supported on two
  consecutive time steps, centered on epoch {j}. That is, its support
  is the interval from {tj - tStep} through {tj + tStep}, where {tj} is
  the time of epoch {j} and {tStep} is the time step. We also assume
  that each element {pi[s,v]} is supported on a single time step, from
  epoch {s-1} to epoch {s}.
  
  With this assumption, each equation of (3) involves coefficients
  from only two consecutive epochs, {c[j-1,u]} and {c[j,u]}.
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
    
  where each block {N[k]}, for {k=0..1}, is the {q × q} matrix of scalar
  products of the differentiated basis elements in epoch {j-k} against
  the gauge functions in epoch {j}, i.e.
  
    { N[k][v][u] == <D(tau[j-k,u]) | pi[j,v]> }

  (Note that the block {N[k]} does not depend on {j}.) At each stage
  of this computation, we advance by one time step. Therefore, we can
  solve the system above one epoch at a time. At each stage {j} we
  solve the equation

    { N[0] a[0] == d - N[1] a[1] }  (6)
  
  where each {a[k]}, called a /frame/, is the vector of {q}
  coefficients in epoch {j-k}, i.e. {a[k][u] = c[j-k,u]}, and {d} is
  the subvector of {b} that corresponds to the gauge functions of
  epoch {j}, that is {d[v] = b[j,v]}.
  
  To start process (6), the client must specify the first frame {a[0]}
  associated with the initial epoch {j = 0}. Then, at each stage we
  update the current state by doing {a[1] = a[0]}, {j = j+1}, and
  we compute the current frame {a[0]} for epoch {j} from the known
  previous frame {a[1]} of epoch {j-1}.

  If the time basis elements are continuous to order {z}, then each
  frame {a[0]} associated with an epoch {j} defines the value of the
  solution, and all time derivatives to order {z}, at the epoch's time
  {tj}. Thus, given only the frame {a[0]}, we can evaluate the
  approximate solution at time {tj}, and its derivatives up to order {z}.
*/

#define _GNU_SOURCE

#include <SPVector.h>
#include <SPTimeOnlyProcFunction.h>
#include <SPTimeOnlyFuncMap.h> 
#include <SPOptions.h>
#include <SP1DSpline.h>
#include <SPTimeOnlySolveG3C1Aux.h>

#include <udg_pulse.h>
#include <rn.h>
#include <r2x2.h>
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
    "  -cMass NUM -cFriction NUM -cSpring NUM \\\n" \
    "  -iniName NAME -solName NAME \\\n" \
    "  -time NUM NUM -nSteps NUM\\\n" \
    "  [ -maxIter NUM ] [ -absTol NUM ] [ -relTol NUM ]\\\n" \
    "  [ -smpOrder NUM ] \\\n" \
    "  -outName NAME \\\n" \
    "  [ -writeEvery NUM ] [ -printEvery NUM] \\\n" \
    "  [ -plotSteps NUM] \\\n" \
    "  [ -verbose ] \n"

typedef struct Options
  { char *outName;          /* Prefix for output filenames. */
    /* Parameters of the problem: */
    char *rhsName;          /* Name of right-hand-side function {TMap}. */
    double cMass;           /* Mass coefficient {M} of equation. */
    double cFriction;       /* Friction coefficient {R} of equation. */
    double cSpring;         /* Spring constant {K} of equation. */
    char *iniName;          /* Name of initial state coeff file. */
    char *solName;          /* Name of true solution. */
    /* The temporal approximation basis: */
    double timeIni;         /* Starting time. */
    double timeStep;        /* Time step. */
    double nSteps;          /* Number of time steps. */
    /* Parameters for ODE integration: */
    bool_t verbose;         /* TRUE prints right-hand-side vectors, etc.. */
    int smpOrder;           /* Triangle sampling order for spherical integrals. */
    SPSys_GenOptions_t gso; /* Parameters for solving the non-linear system. */
    /* Output options: */
    int printEvery;         /* Print one every this many frames (0 = none). */
    int writeEvery;         /* Write one every this many frames. */
    int plotSteps;          /* Number of plot sub-steps per time step (0 = no plot). */
  } Options;

Options *GetOptions(int argn, char **argc);

SPTimeOnlyFunction *GetTrueSolution(char *solName);
  /* Reads the true solution from file "{solName}.tfn". */

/* ODE INTEGRATION */

/* In the following procedures, frames {a[0..w-1]} are interpreted
  relative to the time basis {tbas}. The arguments {N[0..w-1]} are
  the matrices of system (6), as explained under {GetMatrices}. If
  {sol} is not NULL, the space component defined by frame {a[0]} is
  compared with {sol(timeIni)}. */

void LoadInitialFrame(char *iniName, int j, SPVector a, bool_t verbose);
  /* Reads frame {j} of the initial state into the vector {a} (that
    is, sets {a[u] = c[j,u]} for all {i,u}), from file
    "{iniName}-{j}.tst". Ignores the frame time {tj} 
    specified in the file. */

char *FrameName(char *outName, int epoch);
  /* Returns a name for the frame number {epoch}, "{outName}-{epoch}".
    However the name is "{outName}-ini" if {epoch} is -1, and 
    "{outName}-fin" if {epoch} is {MAXINT}. */

void OutputFrame
  ( Options *o, 
    int j,
    double tj,
    double tStep,
    SPVector a, 
    TBasis *tbas, 
    SPTimeOnlyFunction *sol,
    bool_t writeIt,
    bool_t printIt
  );
  /* Outputs the frame {a} as appropriate. Assumes that {a} is
    centered at epoch {j}, i.e. its support is the interval {tj -
    tStep} to {tj + tStep}.
    
    If {writeIt} is true, the frame {a} is written to disk as
    "{o.outName}-{j}.tst". If {printIt} is true, the 
    function value at epoch {j} (time {tj}) is evaluated,
    compared to the true solution {sol}, and the
    results are printed to {stderr}.  */
  
void PrintMaxErrorValues
  ( int j,
    double tj,
    double gVal,
    SPTimeOnlyFunction *sol,
    double *gAbs, 
    double *sAbs, 
    double *eVal
  );
  /* Computes and prints the absolute values {gMax,sMax} of
    the approximation {app}, and of the true solution {sol(tj)} at
    time {tj}, respectively; and the signed error
    {eVal} between {app} and {sol(tj)}. */

void ReadFrame(FILE *rd, SPVector a, double *tj);
  /* Reads a frame from {rd}, stores it in {a}, and stores the
    corresponding epoch time in {*tj}. */

void WriteFrame(FILE *wr, SPVector a, double tj);
  /* Writes a frame {a} to {wr}. Assumes that the corresponding
    epoch time is {tj}. */

int main(int argn, char **argc)
  { Options *o = GetOptions(argn, argc);
    
    int smpCheck = 50; /* Time samples per step when checking error. */

    /* Get parameters of temporal basis (Hermite-style 2-step): */
    double tStep = o->timeStep;
    udg_pulse_family_t tFam = udg_pulse_family(udg_PK_H, 1, 3);
    int tDim = tFam.nmp;  /* Number {q} of temporal elements per epoch. */
    affirm(tDim == 2, "wrong num of pulses per epoch");
    affirm(udg_pulse_max_supp_count(&tFam) == 2, "wrong time pulse width");
    
    /* Compose the temporal basis: */
    TBasis tbas = (TBasis) { tFam };
    int frameDim = tDim; /* Num coeffs in each frame. */
    
    SPTimeOnlyFuncMap RHS = SPTimeOnlyFuncMap_FromName(o->rhsName);

    SPTimeOnlyFunction *sol = GetTrueSolution(o->solName);
    
    SPVector a[2]; /* {a[k]} is the coeff vector of frame {j - k}. */
    int k;
    for (k = 0; k < 2; k++) { a[k] = double_vec_new(frameDim); }
  
    /* Obtain the system's matrices: */
    r2x2_t N[2];  /* The matrices {N[0..2-1]}. */ 
    SPTimeOnlyG3C1_GetMatrices
      ( o->cMass, o->cFriction, o->cSpring, 
        tStep, &tbas, N, 
        o->verbose
      );

    FILE *pltWr = NULL;
    if (o->plotSteps > 0)
      {
        pltWr = open_write(txtcat(o->outName, ".plt"), TRUE);
        fprintf(pltWr, "# %16s",  "t");
        fprintf(pltWr, "  %16s %16s %16s",  "app(t)", "TD(app(t))", "TDD(app(t))");
        fprintf(pltWr, "  %16s %16s %16s",  "LHS(app(t))", "RHS(app(t),t)", "LHS-RHS");
        fprintf(pltWr, "  %16s %16s",  "sol(t)", "err(t)");
        fprintf(pltWr, "\n");
      }
    int j = 0; /* Index of current epoch. */
    
    double errMax = 0; /* Maximum error {app(t)-sol(t)} seen. */
    double resMax = 0; /* Maximum residue {lhs(t)-rhs(t)} seen. */
    while(TRUE)
      { 
        /* Main iteration: compute (or read) the frame {a[0]} for epoch {j},
          that is, {a[0][u] = c[j,u]}. */
        
        /* Compute the epoch's time {tj}: */
        double tj = o->timeIni + j*tStep;
        fprintf(stderr, "### epoch t[%d] = %8.2f #################\n", j, tj);
        
        if (j == 0)
          { /* Load frame {j} of the initial state: */
            LoadInitialFrame(o->iniName, j, a[0], o->verbose);
          }
        else
          { /* At this point we know the previous frame {a[1]}, that is, the
              coefficients {a[1][u] = c[j-1,u]}. Compute the next frame
              {a[0]} from {a[1]} by solving the Galerking constraints in the 
              interval between {tj-tStep} and {tj}: */
            SPTimeOnlyG3C1_ComputeNextFrame
              ( tj, tStep, a, N, 
                RHS, &tbas, 
                o->smpOrder, 
                &(o->gso), o->verbose
              );
            if (o->plotSteps > 0)
              { /* Plot solution in {[tj-tStep __ tj]}: */
                SPTimeOnlyG3C1_PlotSolution
                  ( pltWr, 
                    a, tj, tStep, &tbas, 
                    o->cMass, o->cFriction, o->cSpring, RHS, sol,
                    o->plotSteps
                  );
              }
            /* Check error in {[tj-tStep __ tj]}: */
            SPTimeOnlyG3C1_CheckSolution
              ( a, tj, tStep, &tbas, 
                o->cMass, o->cFriction, o->cSpring, RHS, sol,
                smpCheck, &errMax, &resMax
              );
          }
        if (o->verbose)
          { int i;
            fprintf(stderr, "  a[0] = (");
            for (i = 0; i < a[0].ne; i++) 
              { fprintf(stderr, " %16g", a[0].e[i]); }
            fprintf(stderr, " )\n");
          }
        /* Output the new frame {a[0]} if appropriate: */
        bool_t wrt = ((j == 0) || (j == o->nSteps) || (j % o->writeEvery == 0));
        bool_t prn = ((j == 0) || (j == o->nSteps) || (j % o->printEvery == 0));
        OutputFrame(o, j, tj, tStep, a[0], &tbas, sol, wrt, prn);

        /* Are we done? */
        if (j >= o->nSteps) { break; }
        
        /* Prepare for next iteration: */
        { SPVector atmp = a[1]; a[1] = a[0]; a[0] = atmp; }
        j++;
      }
    fprintf(stderr, "### end of integration ###################\n"); 
    /* Force output of final frame, with name "{o.outName}-fin": */
    double timeFin = o->timeIni + o->nSteps*tStep;
    OutputFrame(o, INT_MAX, timeFin, tStep, a[0], &tbas, sol, TRUE, TRUE);
    
    fprintf(stderr, "errMax = %16g resMax = %16g\n", errMax, resMax);
    
    if (o->plotSteps > 0) { fclose(pltWr); }
        
    return 0;
  }

void LoadInitialFrame(char *iniName, int j, SPVector a, bool_t verbose)
  { 
    char *frmName = FrameName(iniName, j);
    char *fileName = jsprintf("%s.tst", frmName);
    FILE *rd = open_read(fileName, TRUE);
    double t0;
    ReadFrame(rd, a, &t0);
    free(fileName); free(frmName);
  }
  
void OutputFrame
  ( Options *o, 
    int j,
    double tj,
    double tStep,
    SPVector a, 
    TBasis *tbas, 
    SPTimeOnlyFunction *sol,
    bool_t writeIt,
    bool_t printIt
  )
  { char *frmName = FrameName(o->outName, j);
    
    /* Write frame coefficients to disk: */
    if (writeIt)
      { char *fileName = NULL;
        char *fileName = jsprintf("%s.tst", frmName);
        FILE *wr = open_write(fileName, TRUE);
        WriteFrame(wr, a, tj);
        fclose(wr);
        free(fileName);
      }

    if (printIt) 
      {
        /* Compute value of approximate solution at time {tj}: */
        double aVal = a.e[0]; 
        /* Compute value of true solution at time {tj}: */
        double sVal = sol->m->eval(sol, tj);
        /* Compare solutions: */
        double eVal = aVal - sVal;
        
        fprintf(stderr, "iter %6d", j);
        fprintf(stderr, " app = %16.12f", aVal);
        fprintf(stderr, " sol = %16.12f", sVal);
        fprintf(stderr, " app-sol = %16.12f\n", eVal);
      }
    free(frmName);
  }
  
#define SPTimeOnlyFunctionFrame_FileFormat "2005-08-18"

void ReadFrame(FILE *rd, SPVector a, double *time)
  { int N, i;
    filefmt_read_header(rd, "SPTimeOnlyFunctionFrame", SPTimeOnlyFunctionFrame_FileFormat);
    fget_skip_formatting_chars(rd);
    (*time) = nget_double(rd, "time"); fget_eol(rd);
    N = nget_int32(rd, "coeffs"); fget_eol(rd);
    affirm(N == a.ne, "wrong frame size");
    for (i = 0; i < N; i++)
      { a.e[i] = fget_double(rd);
        fget_eol(rd);
      }
    filefmt_read_footer(rd, "SPTimeOnlyFunctionFrame");
  }

void WriteFrame(FILE *wr, SPVector a, double tj)
  { int i;
    filefmt_write_header(wr, "SPTimeOnlyFunctionFrame", SPTimeOnlyFunctionFrame_FileFormat);
    fprintf(wr, "time = %.6f\n", tj);
    fprintf(wr, "coeffs = %d\n", a.ne);
    for (i = 0; i < a.ne; i++)
      { fprintf(wr, "%22.16e\n", a.e[i]); }
    filefmt_write_footer(wr, "SPTimeOnlyFunctionFrame");
    fflush(wr);
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

SPTimeOnlyFunction *GetTrueSolution(char *solName)
  { SPTimeOnlyFunction *f = (SPTimeOnlyFunction *)SPTimeOnlyProcFunction_FromName(solName);
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

    SPOptions_GetKeyword(pp, "-cMass");                               
    o->cMass = SPOptions_GetNextDouble(pp, -DBL_MAX, +DBL_MAX);  

    SPOptions_GetKeyword(pp, "-cFriction");                               
    o->cFriction = SPOptions_GetNextDouble(pp, -DBL_MAX, +DBL_MAX);  

    SPOptions_GetKeyword(pp, "-cSpring");                               
    o->cSpring = SPOptions_GetNextDouble(pp, -DBL_MAX, +DBL_MAX);  

    SPOptions_GetKeyword(pp, "-iniName");                               
    o->iniName = SPOptions_GetNext(pp);  

    SPOptions_GetKeyword(pp, "-solName");                               
    o->solName = SPOptions_GetNext(pp);  

    SPOptions_GetKeyword(pp, "-timeIni");
    o->timeIni = SPOptions_GetNextDouble(pp, -DBL_MAX, DBL_MAX);

    SPOptions_GetKeyword(pp, "-timeStep");
    o->timeStep = SPOptions_GetNextDouble(pp, -DBL_MAX, DBL_MAX);

    SPOptions_GetKeyword(pp, "-nSteps");
    o->nSteps = SPOptions_GetNextInt(pp, 1, 1000000);

    SPOptions_GetKeyword(pp, "-outName");                               
    o->outName = SPOptions_GetNext(pp);  
    
    if (SPOptions_TestKeyword(pp, "-maxIter"))
      { o->gso.maxIter = SPOptions_GetNextInt(pp, 0, INT_MAX); }
    else
      { o->gso.maxIter = 50; }

    if (SPOptions_TestKeyword(pp, "-relTol"))
      { o->gso.relTol = SPOptions_GetNextDouble(pp, 0.0, DBL_MAX); }
    else
      { o->gso.relTol = 0.0; }

    if (SPOptions_TestKeyword(pp, "-absTol"))
      { o->gso.absTol = SPOptions_GetNextDouble(pp, 0.0, DBL_MAX); }
    else
      { o->gso.absTol = 0.0; }

    SPOptions_GetKeyword(pp, "-smpOrder");
    o->smpOrder = SPOptions_GetNextInt(pp, 1, INT_MAX);

    o->verbose = SPOptions_TestKeyword(pp, "-verbose");
                                                 
   /* Output options: */

    if (SPOptions_TestKeyword(pp, "-writeEvery"))
      { o->writeEvery = SPOptions_GetNextInt(pp, 0, INT_MAX); }
    else
      { o->writeEvery = (o->nSteps + 199)/200; }

    if (SPOptions_TestKeyword(pp, "-printEvery"))
      { o->printEvery = SPOptions_GetNextInt(pp, 0, INT_MAX); }
    else
      { o->printEvery = (o->nSteps + 199)/200; }
    
    if (SPOptions_TestKeyword(pp, "-plotSteps"))
      { o->plotSteps = SPOptions_GetNextInt(pp, 0, INT_MAX); }
    else
      { o->plotSteps = 0; }
    
    SPOptions_Finish(pp);
    return o;
  }

