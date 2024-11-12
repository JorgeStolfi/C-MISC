/* Solves a space-time differential equation on S^2. */
/* Last edited on 2023-10-15 03:25:44 by stolfi */

#define PROG_NAME "SPTimeSpaceSolveGen"

/* ??? FUTURE WORK */

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
  
    { TD(f)(p) - F·SLap(f)(p) + R·V(s)¤SGrd(f)(p) - L·f(p) } (2)
  
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
  
  We look for an approximation {h} to the solution {f}, consisting of
  a linear combination of /finite space-time elements/,
  
    { h(p) = SUM { c[iju] · phi[iju](p) } } 
   
  where each {c[iju]} is a real coefficient (to be determined), and
  each finite element {phi[iju]} is a known function of bounded support
  in time and space. (The element index is written {iju} to emphasize
  that the element grid stretches both in space {i} and time {j},
  possibly with several elements {u} for each time step.)
  
  After Galerkin, we look for the coefficients {c[iju]} that anihilate
  scalar products
  
    { <D(h) - H | theta[rsv]> }
    
  where {H} is the function {p -> TMap(h(p), p)}, and {theta}
  is another family of finite space-time elements, the /gauge
  functions/, which are `dense enough' in a sense that will become 
  clear later on.
  
  The Galerkin criterion reduces to solving the  non-linear system
  
    { M c == b }        (3)
    
  where {M} is a fixed matrix
  
    { M[iju,rsv] == < D(phi[iju]) | theta[rsv] > }
        
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
  
  System (3) is usually too big to solve at once. Fortunately the
  boundary conditions are of the initial-value type so we can
  contemplate solving (3) epoch by epoch in the time direction.
  
  At each stage of this computation, we advance by one grid step in
  the time direction. Namely, after stage {j} of the computation we
  assume that all coefficients {c[i,j',u]} with {j' < j} are known. We
  then compute the next layer of unknown coefficients {c[r,j,v]},
  using only the equations of {M} that involve these coefficients and
  do not involve any coefficients {c[i,j',u]} with {j' > j}.
  
  To simplify the indexing, we make an assumption that the support of
  each temporal basis element {tau[j,u]} consists of some number
  {wTau} of consecutive time steps, ending at epoch {j}. That is, its
  support is the interval from {tj-wTau*tstep} through {tj}, where
  {tj} is the time of epoch {j}. Similarly, we assume that each
  element {pi[s,v]} is supported on {wPi} consecutive time steps,
  ending at epoch {s}.

  With this assumption, each equation of (3) involves coefficients
  from only {w = wTau+wPi-1} consecutive epochs
  {c[i,j,u],c[i,j-1,u],... c[i,j-(w-1),u]}. Moreover, since the
  space-time elements (approximators and gauges) in each epoch are
  time-translates of those in the previous epoch, the coefficients of
  {M} repeat epoch after epoch. More precisely, {M} has a periodic
  block-band structure
  
    {
        ... 
           
        ...   0     N[w-1]  N[w-2]  ...     ...  N[0]   0      0     0   ...
        
        ...   0     0       N[w-1] N[w-2]   ...  N[1]   N[0]   0     0   ...
        
        ...   0     0       0      N[w-1]   ...  N[2]   N[1]   N[0]  0   ...
  
        ... 
 
    }
    
  where each block {N[k]}, for {k=0..w-1}, is the { m·q × m·q } matrix
  of the scalar products of the differentiated basis elements in epoch
  {j-k} against the gauge functions in epoch {j}, i.e.
  
    { N[k][iu,rv] == <D(phi[i,j-k,u]) | theta[r,j,v]> }
  
  (Note that the block {N[k]} does not depend on {j}.) Therefore, at
  each stage {j} we solve the equation

    { N[0] a[0] == d - N[1] a[1] - N[2] a[2] ... - N[w-1] a[w-1] }  (6)
  
  where each {a[k]}, called a /frame/, is the subvector of {c}
  consisting of the coefficients in epoch {j-k}, i.e. {a[k][iu] =
  c[i,j-k,u]}, and {d} is the subvector of {c} that corresponds to the
  gauge functions of epoch {j}, that is {d[rv] = c[r,j,v]}.
  
  A set of {m} consecutive frames, where {m >= wTau}, defines
  the solution over {m-(wTau-1)} consecutive time intervals. If
  {t} is the final time of the latest frame in the set, the interval of
  definition begins at time {t - tstep*(wTau-1)}.
  
  A set of {w-1} consecutive frames is a /state/ of the simulation.
  Thus, at each stage we compute the next frame {a[0]} from the
  current state (frames {a[1..w-1]}, and update the current state by
  doing {a[1..w-1] = a[0..w-2]}.
  
  To start process (6), the client must specify the first {w-1} frames
  {a[1..w-1]} --- i.e., the initial state.
  
  For a fixed approximation order, one can usually trade the number of
  frames in the state, {w-1}, against the number of temporal elements
  per epoch {q}, while maintaining the approximation order.
*/

#define _GNU_SOURCE

#include <SPTriang.h>
#include <SPMatrix.h>
#include <SPVector.h>
#include <SPSys.h>
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
#include <SPTimeOnlySolveG3C1Aux.h>

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

typedef struct Options
  { char *outName;       /* Prefix for output filenames. */
    /* Parameters of the problem: */
    char *rhsName;          /* Name of right-hand-side function {TMap}. */
    double cDiffuse;        /* Diffusion coeff {F} of equation. */
    double cDrift;          /* Angular speed {R} of drift term. */
    double cDecay;          /* Decay coeff {L} of equation. */
    char *iniName;          /* Name of initial state coeff file. */
    char *solName;          /* Name of true solution. */
    /* The spatial aproximation basis: */
    char *spcName;          /* Prefix for spatial basis filenames. */
    char *matName;          /* Prefix for spatial matrix filenames (default {spcName}). */
    /* The temporal aproximation basis: */
    double timeIni;         /* Starting time. */
    double timeFin;         /* Ending time. */
    double nSteps;          /* Number of time steps. */
    double timeKind;        /* Kind of temporal basis elements. */
    double timeDegree;      /* Degree of temporal splines. */
    double timeCont;        /* Continuity class of temporal splines. */
    /* Parameters for PDE integration: */
    bool_t verbose;         /* TRUE prints right-hand-side vectors, etc.. */
    int smpOrderSpace;      /* Sampling order for Gaussian quadrature in space. */
    int smpOrderTime;       /* Number of samples for Gaussian quadrature in time. */
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

void GetPDEMatrices
  ( char *matName, 
    double cDiffuse, 
    double cDrift,
    double cDecay, 
    double tStep,
    STBasis *stbas,
    SPMatrix N[]
  );
  /* Obtains the matrices {N[0],N[1]} of the linear system (6) for the
    diffusion-decay equation (1,2). The arguments {cDiffuse},
    {cDrift}, and {cDecay} are the constants {F}, {R} and {L} of
    equation (2), respectively. The matrices are allocated as needed.
    
    When computing scalar products of sphere-time elements, the
    space-domain factor is obtained from the scalar product matrices
    {<sigma[i]|sigma[r]>}, {<V¤SGrd(sigma[i])|sigma[r]>}, and
    {<SLap(sigma[i])|sigma[r]>}, which are read from files
    "{matName}-ev.mat", "{matName}-vg.mat", "{matName}-sl.mat",
    respectively. */
  
/* ODE INTEGRATION */

/* In the following procedures, frames {a[0..w-1]} are interpreted
  relative to the space-time basis {stbas}. The arguments {N[0..w-1]} are
  the matrices of system (6), as explained under {GetPDEMatrices}. If
  {sol} is not NULL, the space component defined by frame {a[0]} is
  compared with {sol(timeIni)}. */

double LoadInitialFrame(char *iniName, int j, SPVector a);
  /* Reads frame {j} of the initial state into the vector {a} (that
    is, sets {a[iu] = c[i,j,u]} for all {i,u}), from file
    "{iniName}-{j}.tst". Returns the extinction time {tj} of the 
    frame, specified in the file. */

char *FrameName(char *outName, int epoch);
  /* Returns a name for the frame number {epoch}, "{outName}-{epoch}".
    However the name is "{outName}-ini" if {epoch} is -1, and 
    "{outName}-fin" if {epoch} is {MAXINT}. */

void OutputState
  ( SPPlot_Stream *fps,
    Options *o, 
    int j,
    double tj,
    double tstep,
    int w,
    SPVector a[], 
    STBasis *stbas, 
    Triangulation *tri,
    SPTimeSpaceFunction *sol,
    bool_t writeIt,
    bool_t printIt,
    bool_t plotIt
  );
  /* Outputs the state {a[0..w-2]} ({w-1} frames) as appropriate.
    Assumes that frame {a[k]} refers to epoch {j-k}, i.e. its support
    ends at time {tj - k*tstep}. If {writeIt} is true, each frame {a[k]}
    is written to disk as "{o.outName}-{j-k}.tst". If {printIt} is true, a
    snapshot of the function at epoch {j-wTau+1} is evaluated
    at selected points of the sphere, compared to the true solution
    {sol}, and the results are printed to {stderr}. If {plotIt} is
    true, that snapshot is plotted to the stream {fps}.  */
  
void PrintMaxErrorValues
  ( int epoch,
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

void ReadFrame(FILE *rd, SPVector a, double *tj);
  /* Reads a frame from {rd}, stores it in {a}, and stores the
    corresponding epoch time in {*tj}. */

void WriteFrame(FILE *wr, SPVector a, double tj);
  /* Writes a frame {a} to {wr}. Assumes that the corresponding
    epoch time is {tj}. */

SPPlot_Stream *GetPlotStream(char *outName, SPPlotOptions_t *po);
  /* Opens an output Postscript stream with name prefix {outName}
    configured according to the options {po}. */

int main(int argn, char **argc)
  { Options *o = GetOptions(argn, argc);
    SPPlotOptions_t *po = &(o->plt);
    SPIntegral_SetDefaultSamplingOrder(o->smpOrderSpace);

    /* Get spatial basis: */
    Basis sBas = SPApprox_ReadBasis(txtcat(o->spcName, ".sBas"));
    int sDim = sBas.ne;    /* Dimension {n} of space basis. */
    Triangulation *tri = SPSpline_BasisTriangulation(sBas);

    /* Get parameters of temporal basis: */
    udg_pulse_family_t tFam = udg_pulse_family(o->timeKind, o->timeCont, o->timeDegree);
    int tDim = tFam.nmp;  /* Number {q} of temporal elements per epoch. */
    
    /* Compute the time step: */
    double tStep = (o->timeFin - o->timeIni)/o->nSteps;
    
    /* Compose the spatio-temporal basis: */
    STBasis stbas = (STBasis) { sBas, tFam };
    int frameDim = tDim * sDim; /* Num coeffs in each frame. */
    
    /* Number of frames that appear in each equation: */
    int w = 2*udg_pulse_max_supp_count(&tFam) - 1;
    
    SPTimeSpaceFuncMap RHS = SPTimeSpaceFuncMap_FromName(o->rhsName);

    SPTimeSpaceFunction *sol = GetTrueSolution(o->solName);
    
    SPVector a[w]; /* {a[k]} is the coeff vector of frame {j - k}. */
    int k;
    for (k = 0; k < w; k++) { a[k] = double_vec_new(frameDim); }
  
    /* Obtain the system's matrices: */
    SPMatrix N[w];  /* The matrices {N[0..w-1]}. */ 
    GetPDEMatrices(o->matName, o->cDiffuse, o->cDrift, o->cDecay, tStep, &stbas, N);

    /* Obtain the factors of {N[0]}, if necessary: */
    SPMatrix N0L, N0R; 
    SPSys_ComputeNeededFactors(&(N[0]), &(o->lso), &N0L, &N0R);
    
    /* FILE *errWr = open_write(txtcat(o->outName, ".erp"), TRUE); */
    /* fprintf(errWr, "# %6s %6s %16s\n",  "epoch", "iter", "rmsError"); */

    /* Fix perspective parameters, if needed: */
    SPPlot_FixView(&(po->obs), 0.0, &(po->upp), NULL, 0.0, NULL, TRUE);

    /* Open the figure stream: */
    SPPlot_Stream *fps = GetPlotStream(o->outName, po);
    
    int j = 1; /* Index of current epoch. */
    while(TRUE)
      { /* At this point we know the frame {a[k]}, that is, the
          coefficients {a[k][iu] = c[i,j-k,u]}, for {k=1..min(w,j)-1}.
          In this interation we get frame {a[0}, that is, {a[0][iu] =
          c[i,j,u]}. Not that at the first iteration ({j=1}) we
          have no frames. */
          
        double alpha = ((double)j)/((double)o->nSteps);
        double tj = (1.0 - alpha)*o->timeIni + alpha*o->timeFin;
        fprintf(stderr, "### epoch t[%d] = %8.2f #################\n", j, tj);
        if (j < w - 1)
          { /* Load frame {j} of the initial state: */
            LoadInitialFrame(o->iniName, j, a[0]);
          }
        else
          { /* Compute next frame from current one: */
            ComputeNextFrame
              ( tj, a, N, N0L, RHS, &stbas, tri, 
                o->lso, o->gso, o->verbose
              );
          }
        if (j >= wTau-1)
          { int jm = j-(wTau-1);
            /* We have enough frames to compute the function at epoch {jm}. */
            /* Output if appropriate, with name "{o.outName}-{j}": */
            bool_t wrt = ((jm < w) || (jm == o->nSteps) || ((jm - w) % o->writeEvery == 0));
            bool_t prn = ((jm < w) || (jm == o->nSteps) || ((jm - w) % o->printEvery == 0));
            bool_t plt = ((jm < w) || (jm == o->nSteps) || ((jm - w) % o->plotEvery == 0));
            OutputState(fps, o, j, tj, tstep, a[0], w, &stbas, tri, sol, wrt, prn, plt);
          }
        
        /* Are we done? */
        if (j >= o->nSteps) { break; }
        
        /* Prepare for next iteration: */
        SPVector atmp = a[w-1];
        for (k = w-1; k > 0; k--) { a[k] = a[k-1]; }
        a[0] = atmp;
        j++;
      }
    fprintf(stderr, "### end of integration ###################\n"); 
    /* Force output of final frame, with name "{o.outName}-fin": */
    OutputState(fps, o, INT_MAX, o->timeFin, a[0], &stbas, tri, sol, TRUE, TRUE, TRUE);
    
    /* fclose(errWr); */
    pswr_close_stream(fps);
    return 0;
  }
  
void LoadInitialFrame(char *iniName, int j, SPVector a)
  { 
    char *frmTag = fmt_int(j, 6);
    char *fileName = NULL;
    asprintf(&fileName, "%s-%s.tst", iniName, frmTag);
    FILE *rd = open_read(fileName, TRUE);
    double t0;
    ReadFrame(rd, a, &t0);
    free(frmTag); free(fileName);
  }
  
void GetPDEMatrices
  ( char *matName, 
    double cDiffuse, 
    double cDrift,
    double cDecay, 
    STBasis *stbas,
    SPMatrix N[]
  )
  { 
    /* 
      We must build { N[k][iu,rv] == <D(phi[i,j-k,u]) | theta[r,j,v]> } 
      that is 

        { < TD(phi) - F*SLap(phi') + v¤SGrd(phi') - L*f | theta[r,j,v]> }

      where {F = cDiffuse}, {V = cDrift}, {L = cDecay}, and {phi'} is
      short for {phi[i,j-k,u] = sigma[i](e)·tau[j-k,u](t)}. Therefore
      we can write
      
        { N[k] == NTD[k] + F*NSL[k] + V·NSG[k] - L*NF[k] }
        
      where 
        
        { NTD[k][iu,rv]
            := < TD(phi[i,j-k,u])(t) | theta[r,j,v] >
            == < TD(sigma[i]·tau[j-k,u]) | rho[r]·pi[j,v] >
            == < sigma[i]·TD(tau[j-k,u]) | rho[r]·pi[j,v] >
            == < sigma[i] | rho[r] > · < TD(tau[j-k,u]) | pi[j,v] > 
            == SE[i,r] · TDW[k,u,v]
            
          NSL[k][iu,rv] 
            := < SLap(phi[i,j-k,u]) | theta[r,j,v] > 
            == < SLap(sigma[i]·tau[j-k,u]) | rho[r]·pi[j,v] >
            == < SLap(sigma[i])·tau[j-k,u] | rho[r]·pi[j,v] >
            == < SLap(sigma[i]) | rho[r] > · < tau[j-k,u] | pi[j,v] > 
            == SL[i,r] · TW[k,u,v]
            
          NSG[k][iu,rv] 
            := < v¤SGrd(phi[i,j-k,u]) | theta[r,j,v] > 
            == < v¤SGrd(sigma[i]·tau[j-k,u]) | rho[r]·pi[j,v] >
            == < v¤SGrd(sigma[i])·tau[j-k,u] | rho[r]·pi[j,v] >
            == < (v/V)¤SGrd(sigma[i]) | rho[r] > · < tau[j-k,u] | pi[j,v] > 
            == SV[i,r] · TW[k,u,v]
            
          NF[k][iu,rv] 
            := < phi[i,j-k,u] | theta[r,j,v] > 
            == < sigma[i]·tau[j-k,u] | rho[r]·pi[j,v] >
            == < sigma[i] | rho[r] > · < tau[j-k,u] | pi[j,v] > 
            == SE[i,r] · TW[k][u,v]
        }
        
     and
     
        { SE[i,r] := < sigma[i] | rho[r] >
          SL[i,r] := < SLap(sigma[i]) | rho[r] >
          SV[i,r] := < V¤SGrd(sigma[i]) | rho[r] >
          
          TW[k][u,v] := < tau[-k,u] | pi[0,v] >
          TDW[k][u,v] := < TD(tau[-k,u]) | pi[0,v] >
        }
     
     If the velocity field {v} is known beforehand, the matrices
     {SE,SL,SV} can be precomputed. The matrices {TW[k]} and {TDW[k]}
     are small and can be quickly computed.
    */

    Basis sBas = stbas->sBas;        /* Spatial basis. */
    int sDim = sBas.ne;            /* Dimension of space basis. */
    int tDeg = stbas->tFam.g;        /* Degree of temporal elements. */
    int tCont = stbas->tFam.c;       /* Continuity of temporal elements. */
    int tDim = stbas->tFam.nmp;      /* Number of temporal elements per epoch. */
    int colsN = sDim*tDim, rowsN = sDim*tDim;
    N[0] = SPMatrix_Null(colsN, rowsN); 
    N[1] = SPMatrix_Null(colsN, rowsN);

    SPMatrix SE = SPApprox_ReadMatrix(txtcat(matName, "-ev.mat"));
    SPMatrix SL = SPApprox_ReadMatrix(txtcat(matName, "-sl.mat"));
    SPMatrix SV = SPApprox_ReadMatrix(txtcat(matName, "-vg.mat"));
    
    /* Computing the matrices: */
    int u, v, k;
    for (u = 0; u < tDim; u++)
      { for (v = 0; v < tDim; v++)
          { for (k = 0; k < 2; k++)
              { /* Compute the time-domain scalar prods {TW[k,u,v],TDW[k,u,v]}: */
                double TWkuv = TimeBasisProduct(tDeg, tCont, k,u,v);
                double TDWkuv = TimeBasisDiffProduct(tDeg, tCont, k,u,v);
                /* Compute the submatrix {N[k][iu,rv]} for current {u,v}: */
                double SEcoef = TDWkuv - cDecay*TWkuv;
                double SLcoef = cDiffuse*TWkuv;
                double SVcoef = cDrift*TWkuv;
                SPMatrix Wkuv = SPMatrix_Mix(SEcoef, SE, SLcoef, SL);
                SPMatrix Zkuv = SPMatrix_Mix(1.0, Wkuv, SVcoef, SV);
                /* Insert {Zkuv} in {N[k]}: */
                SPMatrix Mk = SPMatrix_MixBlock(1.0, N[k], 1.0, Zkuv, u*sDim, v*sDim);
                free(N[k].ents.e); 
                N[k] = Mk;
                free(Wkuv.ents.e); 
                free(Zkuv.ents.e);
              }
          }
      }
  }
  
SPPlot_Stream *GetPlotStream(char *outName, SPPlotOptions_t *po)
  {
    int nCap = (po->eps ? 0 : 1); /* Number of caption lines. */
    return SPPlot_NewStream
      ( po->eps, outName, po->paperSize, po->figSize, po->projLonLat, nCap );
  }

void OutputState
  ( SPPlot_Stream *fps,
    Options *o, 
    int j,
    double tj,
    SPVector a,
    STBasis *stbas, 
    Triangulation *tri,
    SPTimeSpaceFunction *sol,
    bool_t writeIt,
    bool_t printIt,
    bool_t plotIt
  )
  { 
    char *frmTag = fmt_int(j, 6);
    SPPlotOptions_t *po = &(o->plt);
    
    /* Write frame coefficients to disk: */
    if (writeIt)
      { char *fileName = NULL;
        asprintf(&fileName, "%s-%s.tst", o->outName, frmTag);
        FILE *wr = open_write(fileName, TRUE);
        WriteFrame(wr, a, tj);
        fclose(wr);
        free(fileName);
      }

    if (printIt || plotIt) 
      {
        /* Generate snapshot of computed sol at the frame's central epoch {tj}: */
        SPFunction *g = NULL; 

        fprintf(stderr, "building snapshot for t = %.4f...\n", tj); 
        g = BuildStaticFunction(a, stbas);

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
            SPApprox_PlotFunctionAndError
              ( fps, frmTag, g, solt, 
                tri, relMeshSize, (tri != NULL), TRUE,
                fMax, eMax,
                &(po->obs), &(po->upp)
              );
          }
      }

    free(frmTag);
  }
  
#define SPTimeSpaceFunctionFrame_FileFormat "2003-08-18"

void ReadFrame(FILE *rd, SPVector a, double *time)
  { int N, i;
    filefmt_read_header(rd, "SPTimeSpaceFunctionFrame", SPTimeSpaceFunctionFrame_FileFormat);
    fget_skip_formatting_chars(rd);
    (*time) = nget_double(rd, "time"); fget_eol(rd);
    N = nget_int32(rd, "coeffs"); fget_eol(rd);
    affirm(N == a.ne, "wrong frame size");
    for (i = 0; i < N; i++)
      { a.e[i] = fget_double(rd);
        fget_eol(rd);
      }
    filefmt_read_footer(rd, "SPTimeSpaceFunctionFrame");
  }

void WriteFrame(FILE *wr, SPVector a, double tj)
  { int i;
    filefmt_write_header(wr, "SPTimeSpaceFunctionFrame", SPTimeSpaceFunctionFrame_FileFormat);
    fprintf(wr, "time = %.6f\n", tj);
    fprintf(wr, "coeffs = %d\n", a.ne);
    for (i = 0; i < a.ne; i++)
      { fprintf(wr, "%22.16e\n", a.e[i]); }
    filefmt_write_footer(wr, "SPTimeSpaceFunctionFrame");
    fflush(wr);
  }

char *FrameName(char *outName, int epoch)
  { char *name = NULL;
    asprintf(&name, "%s-%06d", outName, epoch);
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

    SPOptions_SetUsage(pp, 
      PROG_NAME " \\\n"
      "  -rhsName NAME \\\n"
      "  -cDiffuse NUM -cDrift NUM -cDecay NUM \\\n"
      "  -iniName NAME -solName NAME \\\n"
      "  -time NUM NUM -nSteps NUM\\\n"
      "  -timeBasis KIND CONT DEGREE\\\n"
      "  -spcName NAME [ -matName NAME ] \\\n"
      SPSys_LinOptionsHelp " \\\n"
      SPSys_GenOptionsHelp " \\\n"
      "  [ -smpOrder NUM NUM ] \\\n"
      "  -outName NAME \\\n"
      "  [ -plotAll ] [ -plotFinal ] [ -verbose ] \\\n"
      SPPlotOptions_FunctionHelp " \n"
    );

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

    SPOptions_GetKeyword(pp, "-time");
    o->timeIni = SPOptions_GetNextDouble(pp, -DBL_MAX, DBL_MAX);
    o->timeFin = SPOptions_GetNextDouble(pp, o->timeIni, DBL_MAX);

    SPOptions_GetKeyword(pp, "-nSteps");
    o->nSteps = SPOptions_GetNextInt(pp, 1, 1000000);

    SPOptions_GetKeyword(pp, "-timeBasis");
    o->timeKind = SPOptions_GetNextInt(pp, udg_pulse_kind_MIN, udg_pulse_kind_MAX);
    o->timeCont = SPOptions_GetNextInt(pp, 0, 2);
    o->timeDegree = SPOptions_GetNextInt(pp, 0, 7);

    SPOptions_GetKeyword(pp, "-spcName");                               
    o->spcName = SPOptions_GetNext(pp);  
       
    if (SPOptions_TestKeyword(pp, "-matName"))
      { o->matName = SPOptions_GetNext(pp); }
    else
      { o->matName = o->spcName; }

    SPOptions_GetKeyword(pp, "-outName");                               
    o->outName = SPOptions_GetNext(pp);  
    
    o->verbose = SPOptions_TestKeyword(pp, "-verbose");
                 
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

    /* Plotting options: */

    if (SPOptions_TestKeyword(pp, "-plotEvery"))
      { o->plotEvery = SPOptions_GetNextInt(pp, 0, INT_MAX); }
    else
      { o->plotEvery = (o->nSteps + 19)/20; }

    if (SPOptions_TestKeyword(pp, "-writeEvery"))
      { o->writeEvery = SPOptions_GetNextInt(pp, 0, INT_MAX); }
    else
      { o->writeEvery = (o->nSteps + 199)/200; }

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
    fprintf(stderr, "max(fabs(app)) = %16.12f\n", *gMax);
    (*sMax) = SPRange_OnSphere(solValue, 120);
    fprintf(stderr, "max(fabs(sol)) = %16.12f\n", *sMax);
    (*eMax) = SPRange_OnSphere(error, 120);
    fprintf(stderr, "max(fabs(app-sol)) = %16.12f\n", *eMax);
    (*eAvg) = sqrt(fabs(SPIntegral_OnSphere(errorSqr, tri))); 
    fprintf(stderr, "rms(app-sol) = %16.12f\n", *eAvg);
  }

SPFunction *BuildStaticFunction(SPVector a, STBasis *stbas)
  {
    /* Get parameters of temporal basis: */
    affirm(stbas->tFam.g = 2*stbas->tFam.c + 1, "unimplemented degree/cont combination");
    int tDim = stbas->tFam.c + 1;  /* Number {q} of temporal elements per epoch. */
    
    /* Compute coeff vector {c} of spatial basis: */
    SPFunction_Basis sbas = stbas->sBas;
    double c[sbas.ne];
    int i;
    for (i = 0; i < sbas.ne; i++)
      { /* Evaluate and add all time basis elements at central time: */
        double s = 0.0;
        int u;
        for (u = 0; u < tDim; u++)
          { int iu = sbas.ne*u + i;
            s += a.e[iu] * TimeBasisEvalCenter(stbas->tFam.c, u); 
          }
        c[i] = s;
      }
      
    /* Combine the basis elements: */
    return SPFunction_LinComb(c, sbas);
  }
