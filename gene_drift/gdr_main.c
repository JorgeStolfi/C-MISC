#define PROG_NAME "gdr_main"
#define PROG_VERS "2.0"

/* Analyzes surviving asexual lineage ratios female/male. */

/* Last edited on 2023-11-25 16:36:07 by stolfi */

#define gdr_main_C_COPYRIGHT \
  "Duh?"

#define PROG_HELP \
  " " PROG_NAME " \\\n" \
  "    -yStart {yStart} \\\n" \
  "    -yStop {yStop} \\\n" \
  "    -iniSize {iniSize} \\\n" \
  "    -finSize {finSize} \\\n" \
  "    -nRuns {nRuns} \\\n" \
  "    -outPrefix {outPrefix} \\\n" \
  "    -tag {tag} \\\n" \
  "    -demoParms 0 {demoParms[0]} \n" \
  "    -demoParms 1 {demoParms[1]} \n" \
  "  where each {demoParms[s]} is \\\n" \
  "    " gdr_demo_options_HELP ""

#define PROG_INFO \
  "SYNOPSIS\n" \
  "  " PROG_HELP "\n" \
  "\n" \
  "DESCRIPTION\n" \
  "  Simulates the history of a asexually reproduced traits in" \
  " a bisexual population over many years, where each individual" \
  " has a variable number of children with a specific probability" \
  " distribution depending on age and sex.  The simulation is" \
  " repeated many times to obtain reliable statistics.\n" \
  "\n" \
  "  The main parameters are the initial population size, the" \
  " demographic parameters of each sex, the overall" \
  " population growth rate, and the" \
  " number of repetitions of the experiment.\n" \
  "\n" \
  "  The simulation spans a total of {ny} years, numbered" \
  " from {yStart} (most remote in the past) to {yStop} (most recent).  Individuals" \
  " are identified by unique /indices/ starting from 0, in chrono order of birth.\n" \
  "\n" \
  "  The variable of interest in each simulation is the number {nsl(y)} of" \
  " lineages whose founder is alive in year {y} and has at least one descendant in the last" \
  " year {yStop}.\n" \
  "\n" \
  "  The descent is assumed to be asexual (one parent for each" \
  " individual), so this program could be used to simulate heredity" \
  " of the Y chromosome or of mitochondial DNA in humans.\n" \
  "\n" \
  "   Currently the program does not simulate mutations.  It is implicitly" \
  " assumed that each individual is genetically distinct, and that the year" \
  " of least common ancestry of two present-day individuals can be" \
  " determined by comparison of the relevant genomic sequences.\n" \
  "\n" \
  "   The program also does not simulate geographic spread and diffusion, assuming" \
  " that it is not significant for asexually reproduced traits.\n" \
  "\n" \
  "OUTPUTS\n" \
  "  The output files are:\n" \
  "\n" \
  "    \"{outPrefix}-{tag}-{s}-nlins.txt\" for each sex {s} (0 or 1)\n" \
  "      " gdr_lineage_counts_median_range_file_INFO "\n" \
  "\n" \
  "    \"{outPrefix}-{tag}-rsamp.txt\"\n" \
  "      " gdr_lineage_sex_ratio_samples_file_INFO "\n" \
  "\n" \
  "    \"{outPrefix}-{tag}-ratio.txt\"\n" \
  "      " gdr_lineage_sex_ratio_median_range_file_INFO "\n" \
  "\n" \
  "    \"{outPrefix}-{tag}-{s}-probs.tex\"\n" \
  "      " gdr_demo_child_distr_tex_table_file_INFO "\n" \
  "\n" \
  "    \"{outPrefix}-{tag}-parms.tex\"\n" \
  "      " gdr_main_tex_parms_file_INFO "\n" \
  "\n" \
  "OPTIONS\n" \
  "  -yStart {yStart}\n" \
  "  -yStop {yStop}\n" \
  "    These mandatory arguments specify the initial and final year of the simulation. The values may be negative.\n" \
  "\n" \
  "  -iniSize {iniSize}\n" \
  "  -finSize {finSize}\n" \
  "    These mandatory arguments specify the maximum cohort" \
  " size (number of individuals born) on years {yStart} and" \
  " {yStop}, respectively. The max cohort size at intermediate" \
  " years will be interpolated exponentially between" \
  " these two values.\n" \
  "\n" \
  "  -nRuns {nRuns}\n" \
  "    This mandatory argument specifies the number of simulation" \
  " runs to be executed for each sex {s}.\n" \
  "\n" \
  "  -outPrefix {outPrefix}\n" \
  "    This mandatory argument specifies a common prefix for all output file names.\n" \
  "\n" \
  "  -tag {tag}\n" \
  "    This mandatory argument specifies a tag that can be used to" \
  " distinguish separate runs of the program.  It must be a single" \
  " upercase letter. It is also used in the names of macros" \
  " written to the TeX parameters file \"{outPrefix}-{tag}-parms.tex\".\n" \
  "\n" \
  "  -demoParms {s} {demoParms[s]} \n" \
  "    This mandatory argument specifies the key parameters" \
  " of the demographic parameters of each sex. It must appear" \
  " exactly twice, with {s} 0 (male) and 1 (female).  The" \
  " syntax of the parameters is\n" \
  "\n" \
  gdr_demo_options_INFO "\n" \
  "\n" \
  "Duh9?"
    
#define gdr_sample_plots_MAX 20
  /* Max number of sample plot of f:m ratio. */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#include <bool.h>
#include <cmp.h>
#include <jsmath.h>
#include <jsfile.h>
#include <jsrandom.h>
#include <argparser.h>
#include <in.h>
#include <rn.h>
#include <vec.h> 

#include <gdr_demo.h>
#include <gdr_sim.h>
#include <gdr_lineage.h>

typedef gdr_demo_parms_t gdr_chops_t; /* For short. */

typedef struct gdr_main_options_t
  { int32_t iniSize;            /* Number of individuals in initial population. */
    int32_t finSize;            /* Desired number of individuals in final population. */
    int32_t yStart;             /* First year covered by simulation. */
    int32_t yStop;              /* Last year covered by simulation. */
    gdr_chops_t *demoParms[2]; /* Demographic parameters from command line. */
    int32_t nRuns;              /* How many times to repeat the simulation. */
    char *outPrefix;            /* Prefix for output files. */
    char *tag;                  /* Tag identifying the run of the program. */
  } gdr_main_options_t;
  /* The command line options.  */
   
/* DATA STRUCTURES AND NOTATIONS

  See also the comments in {gdr_demo.h}, {gdr_sim.h}, {gdr_lineage.h}.

  SEXES 
  
  The sex of an individual is represented by an integer either 0
  (sometimes called 'male') and 1 ('female'). 
  
  YEARS
  
  Externally the years are numbered {o.yStart..o.yStop}, but internally
  they are numbered {0..ny-1}t. So internal year {y} means external year
  {yStart+y} in files. Unless said otherwise, in the following comments
  the internal numbering is used. */
  
void gdr_main_write_tex_parms_file
  ( gdr_main_options_t *o,
    double_vec_t cProb[],
    double_vec_t cFreq[],
    int32_t nrSample,  /* Number of sample ratio plots. */
    double rsl_sv[]
  );
  /* Writes the TeX parameter file "{outPrefix}-{tag}-parms.tex" as 
    described by {gdr_main_tex_parms_file_INFO}. */

#define gdr_main_tex_parms_file_INFO \
  "A file suitable for input in TeX papers. Defines several TeX" \
  " macros \"\\parm{TAG}{NAME}\" where" \
  " {TAG} is the one-letter tag specified in the command line and {NAME} is" \
  " \"nRuns\", \"nGens\", \"iniSize\", or \"finSize\"."

void gdr_main_save_counts
  ( int32_t ny, 
    int32_t vec[], 
    int32_t r, 
    int32_t nr, 
    int32_t sv[]
  );
  /* Saves counts {vec[0..ny-1]} into column {r} of table
    {sv[0..ny-1,0..nr-1]}. */

void gdr_main_dump_saved_counts(char *name, int32_t s, int32_t ny, int32_t nr, int32_t sv[]);
  /* The argument {sv} must be an array of {ny} rows and {nr} columns, stored by rows.
    Prints its contents on {stderr}, with title "{name}_sv[{s}]". */

#define gdr_cohort_size_MAX 2000
#define gdr_num_years_MAX 20000
#define gdr_nRuns_MAX 1000000
  /* Limits to (hopefully) avoid overflow and out-of-memory crash. */

void gdr_main_check_limits(int32_t ny, int32_t nr, int32_t nbt[]);
  /* Checks whether {ny,nr,nbt[0..ny-1]} are not too big. */

gdr_main_options_t *gdr_parse_options(int32_t argc, char **argv);
  /* Parses the command line arguments. */
  
int32_t main(int32_t argc, char **argv);

/* IMPLEMENTATIONS */

int32_t main(int32_t argc, char **argv)
  {
    bool_t debug = TRUE;
    if (debug) { fprintf(stderr, "> %s\n", __FUNCTION__); }
    
    gdr_main_options_t *o = gdr_parse_options(argc, argv);
    srand(417); srandom(417);
    
    int32_t ny = o->yStop - o->yStart + 1;  /* Number of years in history. */
    int32_t nr = o->nRuns;                  /* Number of iterations of experiment. */

    /* Compute the desired size {nbt[0..ny-1]} of each year: */
    int32_t *nbt = gdr_sim_compute_target_year_cohort_sizes(ny, o->iniSize, o->finSize);
    
    /* Sanity limits check: */
    gdr_main_check_limits(ny, nr, nbt);
    
    /* Population and lineage counts for a single run: */
    int32_t pop[ny];     /* {pop[y]} is the population for year {y}, for some sex {s} and run {r}. */ 
    int32_t nsl[ny];     /* {nsl[y]} is the surv lineage count for year {y}, for some sex {s} and run {r}. */ 

    /* Population and lineage counts tables for all runs and both sexes: */
    int32_t *nsl_sv[2];     /* {nsl_sv[s][r,y]} stores {nsl[y]} for each sex {s}, run {r}, year {y}. */ 
    int32_t *pop_sv[2];     /* {pop_sv[s][r,y]} stores {pop[y]} for each sex {s}, run {r}, year {y}. */ 
    
    double_vec_t cProb[2];  /* Nominal child probabilities for each sex {s}. */
    double_vec_t cFreq[2];  /* Observed child count frequencies for each sex {s}. */
    for (int32_t s = 0; s < 2; s++) 
      { 
        fprintf(stderr, "  ### sex %d ##################################################\n", s);
        gdr_demo_parms_t *dmp = o->demoParms[s];
        cProb[s] = gdr_demo_compute_child_count_distr(dmp);
        gdr_sim_state_t *st = gdr_sim_state_new
          (dmp->fMin, dmp->fMax, &(cProb[s]), dmp->cPrec, ny, nbt);

        gdr_demo_show_distr("target child count", s, st->cProb, NULL);
    
        int64_vec_t cCount = int64_vec_new(0);
        int32_t cMax_obs = -1; /* Max actual child count seen. */
        
        pop_sv[s] = in_alloc(ny*nr); 
        nsl_sv[s] = in_alloc(ny*nr); 
        
        /* Safe year interval to count children: */
        int32_t ybrMin = 2*st->fMax;            /* Try to avoid the initial transient. */  
        int32_t ybrMax = ny - 1 - st->fMax - 1; /* Avoid final period with unrestricted births. */

        for (int32_t r = 0; r < nr; r++) 
          { /* Run the simulation once for {st.ny} years: */
            if (debug) { fprintf(stderr, "=== run %d ===\n", r); }
            gdr_sim_run(st);
            gdr_sim_accumulate_child_count_histogram(st, ybrMin, ybrMax, &(cCount), &(cMax_obs));
            gdr_sim_compute_populations(st, pop);
            gdr_main_save_counts(ny, pop, r, nr, pop_sv[s]);

            gdr_lineage_compute_surviving_counts(st, nsl);
            gdr_main_save_counts(ny, nsl, r, nr, nsl_sv[s]);
          }
        int64_vec_trim(&(cCount), cMax_obs+1);
        cFreq[s] = gdr_demo_freqs_from_counts(cCount.ne, cCount.e);
        
        if (debug && (ny < 1000) && (nr < 10))
          { gdr_main_dump_saved_counts("pop", s, ny, nr, pop_sv[s]);
            gdr_main_dump_saved_counts("nsl", s, ny, nr, nsl_sv[s]);
          }
       
        gdr_lineage_write_counts_median_range
          ( o->outPrefix, o->tag, s, ny, o->yStart, nr, nsl_sv[s], pop_sv[s] );

        gdr_demo_show_distr("child count", s, st->cProb, &(cFreq[s]));
        gdr_demo_write_child_distr_tex_table
          ( o->outPrefix, o->tag, s, 
            st->cProb, st->cPrec, 
            &(cFreq[s]), st->cPrec+2
          );

        free(cCount.e);
        gdr_sim_state_free(st);
        
        fprintf(stderr, "  ############################################################\n");
      }
    double *rsl_sv = rn_alloc(ny*nr); 
    gdr_lineage_compute_sex_ratio_table(ny, nr, nsl_sv[0], nsl_sv[1], rsl_sv);
    int32_t nrSample = (nr < gdr_sample_plots_MAX ? nr : gdr_sample_plots_MAX); /* How many sample plots to write. */
    gdr_lineage_write_sex_ratio_samples
      ( o->outPrefix, o->tag, ny, o->yStart, nr, nrSample, rsl_sv );
    gdr_lineage_write_sex_ratio_median_range
      ( o->outPrefix, o->tag, ny, o->yStart, nr, rsl_sv );
    gdr_main_write_tex_parms_file
      ( o, cProb, cFreq, nrSample, rsl_sv );

    for (int32_t s = 0; s < 2; s++) 
      { free(nsl_sv[s]); 
        free(cProb[s].e);
        free(cFreq[s].e);
      }
    free(rsl_sv);
    if (debug) { fprintf(stderr, "< %s\n", __FUNCTION__); }
    return 0;
  }

void gdr_main_save_counts
  ( int32_t ny, 
    int32_t vec[], 
    int32_t r, 
    int32_t nr, 
    int32_t sv[]
  )
  {
    demand((r >= 0) && (r < nr), "invalid column index {r}");
    for (int32_t y = 0; y < ny; y++) 
      { sv[y*nr + r] = vec[y]; }
  }

void gdr_main_dump_saved_counts(char *name, int32_t s, int32_t ny, int32_t nr, int32_t sv[])
  {
    fprintf(stderr, "  --- %s_sv[%d] ---\n", name, s);
    for (int32_t y = 0; y < ny; y++)
      { fprintf(stderr, "     %4d ", y);
        for (int32_t r = 0; r < nr; r++)
          { fprintf(stderr, " %8d", sv[y*nr + r]); }
        fprintf(stderr, "\n");
      }
    fprintf(stderr, "  --------------\n");
  }

void gdr_main_check_limits(int32_t ny, int32_t nr, int32_t nbt[]) 
  {
    bool_t debug = FALSE;
    if (debug) { fprintf(stderr, "  > %s\n", __FUNCTION__); }
    
    int32_t nbtMax = INT32_MAX/ny/100;
    if (debug) { fprintf(stderr, "    ny = %d nr = %d nbtMax = %d\n", ny, nr, nbtMax); }
    demand(ny <= gdr_lineage_num_years_MAX, "too many years");
    demand(nr <= gdr_lineage_num_runs_MAX, "too many runs");
    demand(nr <= INT32_MAX/ny/2, "too many table entries");
    for (int32_t y = 0; y < ny; y++) 
      { if (debug) { fprintf(stderr, "    nbt[%d] = %d\n", y, nbt[y]); }
        demand(nbt[y] <= nbtMax, "too many births");
      }

    if (debug) { fprintf(stderr, "  < %s\n", __FUNCTION__); }
  }

/* OUTPUT */

void gdr_main_write_tex_parms_file
  ( gdr_main_options_t *o,
    double_vec_t cProb[],
    double_vec_t cFreq[],
    int32_t nrSample,  /* Number of sample ratio plots. */
    double rsl_sv[]
  )
  {
    int32_t ny = o->yStop - o->yStart + 1;
    int32_t nr = o->nRuns;
    char *tag = o->tag;
    
    /* Memorable child count for males: */
    int32_t cBigM;      /* A large number of sons for a man. */
    int32_t cBigNumM;   /* Number of men in 1000 with {chBigM} or more sons. */
    gdr_demo_get_interesting_features(&(cProb[0]), &cBigM, &cBigNumM);
    
    /* Interesting ratio plot features: */
    int32_t gcPlat; /* Year (reverse chrono) when the ratio flattens out. */
    double rslPlat;   /* Sub-maximum of the ratio. */
    gdr_lineage_find_notable_ratio_parameters(ny, nr, rsl_sv, &gcPlat, &rslPlat);

    /* Open the file: */
    char *fname = NULL;
    asprintf(&fname, "%s-%s-parms.tex", o->outPrefix, tag);
    FILE *wr = open_write(fname, TRUE);
    free(fname);
    
    auto void wrdef(char *name, double dval, int32_t ival);
      /* Writes the TeX definition of "\parms{tag}{name}" as
        the given {dval}, or {ival} if {dval=NAN}. */
      
    /* From the command line: */
    wrdef("cAlphaM", o->demoParms[0]->cAlpha, 0);
    wrdef("cAlphaF", o->demoParms[1]->cAlpha, 0);
    wrdef("cPrec",   NAN, o->demoParms[0]->cPrec);
    assert(o->demoParms[1]->cPrec == o->demoParms[0]->cPrec);
    wrdef("nRuns",   NAN, o->nRuns);
    wrdef("yStart",  NAN, o->yStart), 
    wrdef("yStop",   NAN, o->yStop), 
    wrdef("iniSize", NAN, o->iniSize);
    wrdef("finSize", NAN, o->finSize);
    
    /* Internal parameters: */
    wrdef("ny",      NAN, ny), 
    wrdef("nrSample", NAN, nrSample);
    
    /* Max child counts: */
    wrdef("cProbMaxM",  NAN, cProb[0].ne - 1);
    wrdef("cProbMaxF",  NAN, cProb[1].ne - 1);
    wrdef("cFreqMaxM",  NAN, cFreq[0].ne - 1);
    wrdef("cFreqMaxF",  NAN, cFreq[1].ne - 1);

    /* Notable child counts: */
    wrdef("chBigM",    NAN, cBigM);
    wrdef("chBigNumM", NAN, cBigNumM);
    
    /* Notable features of ratio plot: */
    wrdef("gcPlat",   NAN, gcPlat);
    wrdef("rslPlat",  rslPlat, -1);

    fclose(wr);
    return;
    
    void wrdef(char *name, double dval, int32_t ival)
      { 
        fprintf(wr, "\\newcommand{\\parms%s%s}", tag, name);
        if (isnan(dval))
          { fprintf(wr, "{%d}", ival); }
        else
          { fprintf(wr, "{%.14g}", dval); }
        fprintf(wr, "\n");
      }
  }

gdr_main_options_t *gdr_parse_options(int32_t argc, char **argv)
  { 
    argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, PROG_HELP);
    argparser_set_info(pp, PROG_INFO);
    argparser_process_help_info_options(pp);
    
    gdr_main_options_t *o = notnull(malloc(sizeof(gdr_main_options_t)), "no mem");

    argparser_get_keyword(pp, "-iniSize");
    o->iniSize = (int32_t)argparser_get_next_int(pp, 1, gdr_cohort_size_MAX);

    argparser_get_keyword(pp, "-finSize");
    o->finSize = (int32_t)argparser_get_next_int(pp, 1, gdr_cohort_size_MAX);

    argparser_get_keyword(pp, "-yStart");
    o->yStart = (int32_t)argparser_get_next_int(pp, -100000, +3000);

    argparser_get_keyword(pp, "-yStop");
    int32_t yStopMax = o->yStart + gdr_num_years_MAX - 1;
    o->yStop = (int32_t)argparser_get_next_int(pp, o->yStart+2, yStopMax);

    /* The keyword "-demoParms" must appear exactly twice, once for each sex: */
    for (int32_t s = 0; s < 2; s++) { o->demoParms[s] = NULL; }
    for (int32_t ks = 0; ks < 2; ks++)
      { argparser_get_keyword(pp, "-demoParms");
        int32_t s = (int32_t)argparser_get_next_int(pp, 0, 1);
        if (o->demoParms[s] != NULL) 
          { argparser_error(pp, "multiple \"-demoParms\" for the same sex"); }
        o->demoParms[s] = gdr_demo_options_parse(pp);
      }
    if (o->demoParms[0]->cPrec != o->demoParms[1]->cPrec)
      { argparser_error(pp, "the {cPrec} parameter must be the same for both sexes"); }

    argparser_get_keyword(pp, "-nRuns");
    o->nRuns = (int32_t)argparser_get_next_int(pp, 1, gdr_nRuns_MAX);

    argparser_get_keyword(pp, "-outPrefix");
    o->outPrefix = argparser_get_next(pp);

    argparser_get_keyword(pp, "-tag");
    o->tag = argparser_get_next(pp);
    if (strlen(o->tag) != 1)
      { argparser_error(pp, "\"-tag\" must be just one letter"); }
    if ((strcmp(o->tag, "A") < 0) || (strcmp(o->tag, "Z") > 0)) 
      { argparser_error(pp, "\"-tag\" must be in [A-Z]"); }

    argparser_skip_parsed(pp);
    argparser_finish(pp);
    return o;
  }
