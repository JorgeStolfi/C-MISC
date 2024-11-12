#define PROG_NAME "tsp_enum"
#define PROG_DESC "compute cost statistics of hamiltoninan tours or cycle covers"
#define PROG_VERS "1.1"

/* Copyright © 2004 by the State University of Campinas (UNICAMP). */
/* See the copyright, authorship, and warranty notice at end of file. */
/* Last edited on 2023-03-31 04:13:56 by stolfi */

#define PROG_HELP \
  "  " PROG_NAME " \\\n" \
  "  [-seed NUM] \\\n" \
  "  -nv NUM [-tour BOOL] \\\n" \
  "  -distr DISTRNAME DIM \\\n" \
  "  [-npmax NUM] [-nsmax NUM] \\\n" \
  "  [-sdim NUM] \\\n" \
  "  -outName OUTNAME"

#define PROG_INFO \
  "NAME\n" \
  "  " PROG_NAME " - " PROG_DESC "\n" \
  "\n" \
  "SYNOPSIS\n" \
  PROG_HELP "\n" \
  "\n" \
  "DESCRIPTION\n" \
  "  The program enumerates hamiltonian cycles or cycle covers in a weighted graph. See the comments below for futher information.\n" \
  "\n" \
  "OPTIONS\n" \
  "  -seed {NUM}\n" \
  "  -nv {NUM}\n" \
  "  -tour {BOOL}\n" \
  "  -distr {DISTRNAME {DIM}\n" \
  "  -npmax {NUM}\n" \
  "  -nsmax {NUM}\n" \
  "  -sdim {NUM}\n" \
  "  -outName {OUTNAME}\n" \
  "\n" \
  "FILES\n" \
  "  ???\n" \
  "\n" \
  "SEE ALSO\n" \
  "  A doctor, if symptoms persist.\n" \
  "\n" \
  "BUGS\n" \
  "  Ants, flies, cockroaches, etc..\n" \
  "\n" \
  "AUTHOR\n" \
  "  Jorge Stolfi, DCC-IMECC-UNICAMP <stolfi@dcc.unicamp.br>\n" \
  "\n" \
  "MODIFICATION HISTORY\n" \
  "  Created 2004-dec-22 by Jorge Stolfi, IC-UNICAMP."

/* Must define _GNU_SOURCE in order to get {asprintf} */
#define _GNU_SOURCE
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <values.h>
#include <math.h>

#include <bool.h>
#include <sign.h>
#include <affirm.h>
#include <intsort.h>
#include <nget.h>
#include <fget.h>
#include <jsstring.h>
#include <jsrandom.h>
#include <jsfile.h>

#include <tsp_lib.h>
#include <tsp_fit.h>
#include <tsp_cost.h>

#define NV_MAX (256)
/* Maximum number of vertices. */

#define NP_MAX (INT32_MAX/2 + 1)
/* Maximum number of perms that we can generate. */

#define SHOW_PERM (TRUE)
/* Define this as TRUE to print the tour as well as its cost. */

#define DEBUG TRUE
/* Define this as TRUE to print various diagnostics. */

typedef struct options_t 
  { int32_t seed;        /* Seed for random number gnerators. */
    int32_t nv;          /* Number of vertices. */
    bool_t tour;     /* TRUE enumerates tours, FALSE enumerates spins. */
    int32_t npmax;       /* Max number of perms to generate. */
    int32_t nsmax;       /* Max size of sub-sample for asymptotic analysis. */
    Distr_t distr;   /* Type of distribution for edge costs. */
    int32_t d;           /* Parameter for edge cost distribution. */
    int32_t sdim;        /* Value of {Sdim} for model fitting; 0 to be determined. */
    char *outName;   /* Prefix for output file names. */
  } options_t;

void enum_perms(int32_t nv, bool_t tour, int32_t npmax, int32_t *np, vtx_t **p);
/* Enumerates all or some permutations {v[0..nv-1]} of the integers {0..nv-1}.
  
  If {tour} is true, interprets each perm as a /tour/: a Hamiltonian
  circuit of {K_nv}, namely the connected 1-1-regular digraph with arcs
  {v[i]->v[(i+1)%nv]} for all {i}. In that case, considers only those
  tours with {v[0] = 0}, thus avoiding tours that differ only by a a
  cyclic shift; and with {v[1] < v[nv-1]}, thus avoiding tours that
  differ only by reversal.
  
  If {tour} is false, considers each perm as a /spin/ or /derangement/: a cycle cover,
  that is, the (possibly disconnected) 1-1-regular digraph with arcs
  {i->v[i]} with {v[i] != i}, for all {i}.
  
  If {npmax} is nonzero and less than the total number of valid perms,
  generates a random sample of {npmax} perms (possibly with
  repetitions), instead of all of them.
  
  In either case, the number of permutations actually generated is
  returned in {*np}, and the permutations are returned in the array
  {**p}, which is allocated by the procedure. Element {i} of
  permutation {k} is stored in {(*p)[k*nv+i]}, for {k} in {0..(*np)-1}
  and {i} in {0..nv-1}. */
    
void process_global_model(options_t *o, int32_t np, int32_t *ix, vtx_t *p, double *Z);
/* Fits the global sphere-slice model to the costs {Z[0..np-1]}, 
  then writes the actual and model costs, and the permutations {p},
  to a file called "{o->outName}-glob.dat".  The entries are
  listed in the order specified by the indices {ix[0..np-1]}. */

void process_asymptotic_model(options_t *o, int32_t ns, sign_t dir, int32_t np, int32_t *ix, vtx_t *p, double *Z);
/* Fits the asymptotic power model to the {ns} lowest costs (if
  {dir=+1}) or largest costs (if {dir=-1}, then writes the actual and
  model costs of those elements to a file called
  "{o->outName}-asym-{xx}.dat" where {xx} is "lo" or "hi", depending
  on {dir}. Assumes that {ix[0..np-1]} are the indices of the elements
  {0..np-1}, in increasing order of actual cost. */

void write_plot_data(FILE *wr, int32_t nv, int32_t np, int32_t *ix, double *Z, double *TM, double *ZM, vtx_t *p);
/* Writes to {wr} the data for a set of {np} perms.

  For each index {k} in {0..np-1}, writes a line with its actual rank
  {j}, index {k}, actual relative rank {T[k] = (j+0.5)/np}, actual
  cost {Z[k]}, model-predicted relative rank {TM[k]}, model-predicted
  costs {ZM[k]}, and the perm itself {p[k*nv..(k+1)*nv-1]}.

  Assumes that {k=ix[j]} is the index into {Z,TM,ZM} of the perm with
  rank {j} in increasing order of the costs, that is, {Z[ix[j]]} is
  increasing with {j}. The items are written in order of increasing
  {j}.
  
  If {ix} is NULL, assumes {ix[j]=j} for all {j}; if {TM} is NULL,
  assumes {TM[k] = T[k]}; if {ZM} is NULL, assumes {ZM[k]=Z[k]}. */

void write_perm(FILE *wr, int32_t nv, vtx_t *v);
/* Writes the perm {v[0..nv-1]} to file {wr}. */
  
int32_t compute_num_tours(int32_t nv, bool_t sym);
/* Computes the number of tours on {nv} vertices, counting only tours
  that start at 0. If {sym=TRUE} and {nv >= 2}, considers only
  tours {v[0..nv-1]} where {v[1] <= v[nv-1]}. If the result would
  exceed {NP_MAX}, returns {NP_MAX+1}. */

int32_t compute_num_spins(int32_t nv);
/* Computes the number of derangements on {nv} vertices, that is,
  permutations of {0..nv-1} that have no fixed point. If the result
  would exceed {NP_MAX}, returns {NP_MAX+1}. */

int32_t compute_num_quasi_spins(int32_t nv);
/* Number of perms of {0..nv-1} which have at least one cycle and no
  fixed points, except possibly {v[0]=0}. If the result would exceed
  {NP_MAX}, returns {NP_MAX+1}. */

int32_t main (int32_t argc, char **argv);
options_t *get_options(int32_t argc, char **argv);
double parse_double(int32_t *argn, int32_t argc, char **argv, double lo, double hi);
void data_error(int32_t line, char *msg);
void arg_error(char *msg, int32_t argn, char *arg);

FILE *open_plot_file(char *tag, options_t *o, char *title);

int32_t main (int32_t argc, char **argv)
  { 
    options_t *o = get_options(argc, argv);
    
    double c[o->nv*o->nv];

    srandom(o->seed);
    gen_arc_cost_matrix(o->nv, o->distr, o->d, c);

    /* Enumerate tours: */
    vtx_t *p;
    int32_t np;
    enum_perms(o->nv, o->tour, o->npmax, &np, &p);

    /* Compute their costs: */
    double Z[np];
    compute_perm_costs(o->nv, o->tour, c, np, p, Z);
    
    /* Index-sort of costs: */
    int32_t ix[np];
    sort_costs(np, ix, Z);

    /* Fit global sphere-slice model: */
    process_global_model(o, np, ix, p, Z);
    
    /* Fit asymptotic power model: */
    int32_t ns = 2*(int32_t)sqrt(np);
    if (ns > np/2) { ns = np/2; }
    if (ns > o->nsmax) { ns = o->nsmax; }
    process_asymptotic_model(o, ns, +1, np, ix, p, Z);
    process_asymptotic_model(o, ns, -1, np, ix, p, Z);
    
    return 0;
  }

void sort_costs(int32_t np, int32_t *ix, double *Z)
  { 
    auto int32_t cmp(int32_t x, int32_t y);
    /* Compares indices {x,y} according to cost {Z}. */

    int32_t cmp(int32_t x, int32_t y)
      { if (Z[x] < Z[y]) 
          { return -1; }
        else if (Z[x] > Z[y])
          { return +1; }
        else 
          { return 0; }
      }

    int32_t j; 
    for(j = 0; j < np; j++) { ix[j] = j; }
    isrt_heapsort(ix, np, cmp, +1);
  }

void process_global_model(options_t *o, int32_t np, int32_t *ix, vtx_t *p, double *Z)
  { 
    /* Fit global model: */
    int32_t Sdim; double Zmid, Zrad;
    int32_t useSdim = o->sdim;
    fit_global_model(np, ix, Z, o->nv, o->tour, useSdim, &Sdim, &Zmid, &Zrad);
    /* Compute ranks {TM} predicted by global model: */
    double TM[np];
    compute_global_model_ranks(np, Z, TM, Sdim, Zmid, Zrad);
    /* Open plot file: */
    FILE *wr = open_plot_file("glob", o, "Global model");
    /* Write model parameters: */
    fprintf(wr, "Sdim = %4d\n", Sdim);
    fprintf(wr, "Zmid = %24.16e\n", Zmid);
    fprintf(wr, "Zrad = %24.16e\n", Zrad);
    /* Write plot data: */
    write_plot_data(wr, o->nv, np, ix, Z, TM, NULL, p);
    fclose(wr);
  }

void process_asymptotic_model(options_t *o, int32_t ns, sign_t dir, int32_t np, int32_t *ix, vtx_t *p, double *Z)
  { 
    char *hilo = (dir > 0 ? "lo" : "hi");
    char *title = NULL;
    asprintf(&title, "asymptotic model for %d out of %d samples at %s end", ns, np, hilo);
    
    fprintf(stderr, "%s: Fitting %s\n", __FUNCTION__, title); 
    int32_t nv = o->nv;
    demand(ns <= np, "bad ns");
    /* Extract subset to analyze: */
    double ZS[ns];
    vtx_t ps[ns*nv];
    int32_t j;
    for (j = 0; j < ns; j++)
      { int32_t k = (dir > 0 ? ix[j] : ix[np - 1 - j]);
        ZS[j] = (dir > 0 ? Z[k] : -Z[k]);
        vtx_t *v = &(p[k*nv]);
        vtx_t *vs = &(ps[j*nv]);
        int32_t i;
        for (i = 0; i < nv; i++) { vs[i] = v[i]; }
      }
    /* Fit asymptotic model: */
    double Expt, Zref, Coef;
    double useExpt = (o->sdim <= 0 ? 0.0 : ((double)o->sdim + 1)/2);
    fit_asymp_model(ns, ZS, np, nv, useExpt, &Expt, &Zref, &Coef);
    /* Compute costs according to asymptotic model: */
    double ZT[ns];
    compute_asymp_model_costs(ns, ZT, np, Expt, Zref, Coef);
    /* Write plot file: */
    char *tag = txtcat("asym-", hilo);
    FILE *wr = open_plot_file(tag, o, title);
    free(tag); free(title);
    /* Write model parameters: */
    fprintf(wr, "Expt = %6.1f\n", Expt);
    fprintf(wr, "Zref = %24.16e\n", Zref);
    fprintf(wr, "Coef = %24.16e\n", Coef);
    /* Write plot data: */
    write_plot_data(wr, o->nv, ns, NULL, ZS, NULL, ZT, ps);
    fclose(wr);
  }

int32_t compute_num_tours(int32_t nv, bool_t sym)
  {
    if (nv == 0)
      { return 0; }
    else if (nv == 1)
      { return 1; }
    else if (nv == 2)
      { return 1; }
    else if (nv == 3)
      { return (sym ? 1 : 2); }
    else
      { int32_t nt1 = compute_num_tours(nv-1, sym);
        if(nt1 > NP_MAX/(nv-1)) { return NP_MAX+1; }
        return nt1*(nv-1); 
      }
  }

int32_t compute_num_spins(int32_t nv)
  {
    if (nv == 0)
      { return 1; }
    else if (nv == 1)
      { return 0; }
    else 
      { int32_t nq1 = compute_num_quasi_spins(nv-1);
        if (nq1 > NP_MAX/(nv-1)) { return NP_MAX+1; }
        return nq1*(nv-1); 
      }
  }

int32_t compute_num_quasi_spins(int32_t nv)
  /* Number of perms of {0..nv-1} which have at least one 
    cycle and no loops except at vertex {0}. */
  {
    if (nv == 0)
      { return 0; }
    else if (nv == 1)
      { return 1; }
    else
      { /* Count quasi-spins that HAVE a fixed point at 0: */
        int32_t ns1 = compute_num_spins(nv-1);
        /* Count quasi-spins that DO NOT HAVE a fixed point at 0: */
        int32_t nq1 = compute_num_quasi_spins(nv-1);
        if (nq1 > (NP_MAX - ns1)/(nv-1)) { return NP_MAX+1; }
        return ns1 + nq1*(nv-1); 
      }
  }

void enum_perms(int32_t nv, bool_t tour, int32_t npmax, int32_t *np, vtx_t **p)
  {
    demand(nv >= 0, "invalid number of vertices");
    demand(nv <= NV_MAX, "too many vertices");
    
    /* Eventually this will be a parameter: */
    bool_t sym = TRUE; /* If {TRUE}, ignores tours that differ by reversal. */

    /* Compute number {np} of perms to generate. */
    if (npmax != 0) 
      { demand(npmax <= NP_MAX, "asking for too many perms");
        (*np) = npmax; 
      }
    else 
      { if (tour)
          { (*np) = compute_num_tours(nv, sym); }
        else
          { (*np) = compute_num_spins(nv); }
        demand((*np) <= NP_MAX, "there are too many perms");
      }
    
    /* Allocate perm table: */
    (*p) = (vtx_t *)notnull(malloc((*np)*nv*sizeof(vtx_t)), "no mem");
    
    vtx_t v[nv];
    
    auto void do_enum(int32_t k);
      /* If {npmax == 0}, generates all valid perms of {v} that start with
        {v[0..k-1]}, else generates a random one of those perms. Either
        way, stores the enumeradted perm(s) in the {p} array starting at
        {p[ngen*nv]}, and increments {ngen}. Upon exit, {v} is
        restored to its original contents. */

    int32_t ngen = 0; /* Counts output perms. */

    void do_enum(int32_t k)
      { 
        if ((!tour) && (k > 0) && (v[k-1] == k-1))
          { /* Skip spins with loops: */
            return;
          }
        else if (k == nv-1) 
          { /* Perm is complete. */
            /* If costs are symmetric, ignore "decreasing" tours: */
            if (tour && sym && (k > 0) && (v[k] < v[1])) { return; }
            /* Store perm in {p} array: */
            assert (ngen < (*np));
            { int32_t j; 
              vtx_t *pp = &((*p)[ngen*nv]); 
              for (j = 0; j < nv; j++) { pp[j] = v[j]; }
            }
            ngen ++; 
          }
        else if (npmax == 0)
          { /* Try all remaining choices for {v[k]}, recurse: */
            int32_t i;
            for (i = k; i < nv; i++)
              { vtx_t t = v[k]; v[k] = v[i]; v[i] = t;
                do_enum(k+1);
                t = v[k]; v[k] = v[i]; v[i] = t;
              }
          }
        else
          { /* Pick a random choice for {v[k]}, recurse: */
            int32_t i = k + int32_abrandom(0, nv - k-1);
            vtx_t t = v[k]; v[k] = v[i]; v[i] = t;
            do_enum(k+1);
            t = v[k]; v[k] = v[i]; v[i] = t;
          }
      }

    /* Initial perm (identity): */
    for (int32_t i = 0; i < nv; i++) { v[i] = (vtx_t)i; }

    /* Enumerate desired perms: */
    while (ngen < (*np))
      { /* generate all perms or or one perm: */
        if (tour) 
          { /* Enumerate all tours with {v[0] = 0}: */
            do_enum(1);
          }
        else
          { /* Enumerate all derangements: */
            do_enum(0);
          }
      }

    fprintf(stderr, "generated %d valid perms\n", ngen);
    assert(ngen == (*np));
  }

void write_plot_data(FILE *wr, int32_t nv, int32_t np, int32_t *ix, double *Z, double *TM, double *ZM, vtx_t *p)
  { int32_t j;
    for (j = 0; j < np; j++)
      { /* Index of element with actual rank j: */
        int32_t k = (ix == NULL ? j : ix[j]);
        /* Actual relative rank: */
        double Tk = ((double)j + 0.5)/((double) np);
        /* Actual cost: */
        double Zk = Z[k];
        /* Model-predicted relative rank: */
        double TMk = (TM == NULL ? Tk : TM[k]);
        /* Model-predicted cost: */
        double ZMk = (ZM == NULL ? Zk : ZM[k]);
        /* Write data: */
        fprintf(wr, "%10d", j);
        fprintf(wr, " %9d", k);
        fprintf(wr, " %18.10e", Tk);
        fprintf(wr, " %18.10e", Zk);
        fprintf(wr, " %18.10e", TMk);
        fprintf(wr, " %18.10e", ZMk);
        if (SHOW_PERM) 
          { /* Write permutation: */
            vtx_t *v = &(p[k*nv]); 
            write_perm(wr, nv, v);
          }
        /* Finish line: */
        fprintf(wr, "\n");
      }
    fflush(wr);
  }

void write_perm(FILE *wr, int32_t nv, vtx_t *v)
  {
    int32_t i;
    for (i = 0; i < nv; i++) { fprintf(wr, " %d", v[i]); } 
  }

FILE *open_plot_file(char *tag, options_t *o, char *title)
  {
    char *fname = NULL;
    affirm(asprintf(&fname, "%s-%s.dat", o->outName, tag) > 0, "no mem");
    FILE *wr = open_write(fname, TRUE);
    free(fname);
    
    fprintf(wr, "# %s\n", title);
    write_params(wr, o->seed, o->nv, o->tour, o->npmax, o->distr, o->d);
    return wr;
  }

options_t *get_options(int32_t argc, char **argv)
  {
    options_t *o = notnull(malloc(sizeof(options_t)), "no mem"); 
    
    int32_t argn = 1;
    
    /* Preliminary defaults: */
    o->seed = 4615;
    o->nv = -1;
    o->tour = TRUE;
    o->npmax = 0;
    o->nsmax = 200;
    o->distr = DISTR_GEN_USPHE;
    o->d = -1;
    o->sdim = 0;
    o->outName = NULL;
    
    /* Parse command line options: */

    while ((argn < argc) && (argv[argn][0] == '-') && (argv[argn][1] != '\0'))
      { if (strcmp(argv[argn], "-seed") == 0)
          { argn++;
            o->seed = (int32_t)parse_double(&argn, argc, argv, 1.0, ((double)INT32_MAX));
          } 
        else if (strcmp(argv[argn], "-nv") == 0)
          { argn++;
            o->nv = (int32_t)parse_double(&argn, argc, argv, 1.0, ((double)NV_MAX));
          }
        else if (strcmp(argv[argn], "-tour") == 0)
          { argn++;
            o->tour = ((int32_t)parse_double(&argn, argc, argv, 0.0, 1.0) != 0);
          }
        else if (strcmp(argv[argn], "-npmax") == 0)
          { argn++;
            o->npmax = (int32_t)parse_double(&argn, argc, argv, 0.0, ((double)NP_MAX));
          }
        else if (strcmp(argv[argn], "-nsmax") == 0)
          { argn++;
            o->nsmax = (int32_t)parse_double(&argn, argc, argv, 0.0, ((double)NP_MAX));
          }
        else if (strcmp(argv[argn], "-distr") == 0)
          { argn++;
            o->distr = distr_from_name(argv[argn]); argn++;
            o->d = (int32_t)parse_double(&argn, argc, argv, 0.0, ((double)INT32_MAX));
          }
        else if (strcmp(argv[argn], "-sdim") == 0)
          { argn++;
            o->sdim = (int32_t)parse_double(&argn, argc, argv, 0.0, ((double)NV_MAX));
          }
        else if (strcmp(argv[argn], "-outName") == 0)
          { argn++;
            o->outName = argv[argn]; argn++;
          }
        else 
          { arg_error("invalid option", argn, argv[argn]); exit(1); } 
      }

    if (argn != argc) { arg_error("excess arguments", argn, argv[argn]); exit(1); } 
    
    /* Complete defaults: */
    if (o->nv <= 0) 
      { arg_error("must define \"-nv\"", -1, NULL); }
    if (o->outName == NULL) 
      { arg_error("must define \"-outName\"", -1, NULL); }
    return o;
  }

double parse_double(int32_t *argn, int32_t argc, char **argv, double lo, double hi)
  { char *rest;
    double x;
    if (*argn >= argc) 
      { arg_error("argument value is missing", (*argn)-1, argv[(*argn)-1]); }
    x = strtod(argv[*argn], &rest);
    if (*rest != '\0') 
      { arg_error("invalid number", *argn, argv[*argn]); } 
    if ((x < lo) || (x > hi)) 
      { arg_error("out of range", *argn, argv[*argn]); }
    ++(*argn);
    return x;
  }
  
void data_error(int32_t line, char *msg)
  {
    fprintf(stderr, "%s:%d: **%s\n", "-", line, msg);
    exit(1);
  }

void arg_error(char *msg, int32_t argn, char *arg)
  {
    fprintf(stderr, "%s: **", PROG_NAME);
    if ((argn >= 0) && (arg != NULL))
      { fprintf(stderr, " argv[%d] = %s:", argn, arg); }
    fprintf(stderr, " %s\n", msg);
    fprintf(stderr, "usage: %s \\\n%s", PROG_NAME, PROG_HELP);
    exit(1);
  }     
   
