#define PROG_NAME "fit_btc_bubble_model"
#define PROG_DESC "Fits a multi-bubble model to historic bitcoin prices"
#define PROG_VERS "1.0"

#define fit_btc_bubble_model_C_COPYRIGHT \
  "Copyright © 2015 by the State University of Campinas (UNICAMP)"

/* Last edited on 2015-05-02 21:46:24 by stolfilocal */

#define PROG_HELP \
  "  " PROG_NAME " \\\n" \
  "    -inParms {IN_PARMS_FILE} \\\n" \
  "    [ -inRangeParms {IN_LO_PARMS_FILE} {IN_HI_PARMS_FILE} ] \\\n" \
  "    [ -maxLSQIters {MAX_LSQ_ITER} ] \\\n" \
  "    [ -maxNLIters {MAX_NL_ITER} ] \\\n" \
  "    [ -evalRange {EVAL_DT_INI} {EVAL_DT_FIN} ] \\\n" \
  "    -hrad {HRAD} \\\n" \
  "    [ -truncate {TRUNC_LEVEL} ] \\\n" \
  "    -outPrefix {OUT_PREFIX} \\\n" \
  "    " argparser_help_info_HELP " \\\n" \
  "    < {IN_PRICES_FILE} \\\n" \
  "    > {OUT_PRICES_FILE}"

#define PROG_INFO \
  "NAME\n" \
  "  " PROG_NAME " - " PROG_DESC "\n" \
  "\n" \
  "SYNOPSIS\n" \
  PROG_HELP "\n" \
  "\n" \
  "DESCRIPTION\n" \
  "  The program reads from {IN_PRICES_FILE} a series of daily bitcoin prices, and from {IN_PARMS_FILE} some parameters of a multi-bubble model for those prices.  The program adjusts the parameters of the model so that the prices predicted by the model are a best fit for the input prices.  It then writes to standard output a file with the input and model prices, as well as the contribution of each bubble, and also a file \"{OUT_PREFIX}_parms.txt\" with the adjusted values of the parameters.\n" \
  "\n" \
  "INPUT PRICE FILE\n" \
  "  The file {IN_PRICES_FILE} must begin with a line \"days = {ND}\", followed by {ND} lines, one for each day in the modeled time interval.  In each line, leading whitespace is ignored, and the field separator is one or more consecutive whitespace characters.  Comments and blank lines are not allowed.  The first field of each line must be an ISO date \"{YYYY}-{MM}-{DD}\".  The dates are assumed to be in increasing order and consecutive, with no gaps.  The second column must be the mean price for that day.  If the price is zero, the program assumes that there is no price data for that day.\n" \
  "\n" \
  "PARAMETERS FILE\n" \
  btc_bubble_parms_read_INFO "\n" \
  "\n" \
  "  The {COEF} fields in the input file {IN_PARMS_FILE} are ignored.  In the output file, they are set to the coefficients in the fitted linear compbination.\n" \
  "\n" \
  "MAIN OUTPUT FILE\n" \
  "  The main output file has one line for each day in the series, in the format\n" \
  "\n" \
  "    \"{DATE} {TIME} | {P_IN} | {P_MOD} | {RES} | {RAT} | {COMP[0]} ... | {COMP[nb-1]}\"\n" \
  "\n" \
  "  The time field is always \"00:00:00\".  The fields {P_IN} and {P_MOD} are the input and model prices.  The field {RES} is the residual price {P_IN - P_MOD}, while {RAT} is the ratio {P_IN/P_MOD} \n" \
  "\n" \
  "OPTIMIZATION\n" \
  "  The basic fitting procedere adjusts only the coefficients of the linear combiantion of all bubbles, by minimization of the sum of squared absolute differences between the model and given proces.  This basic step uses a robust least squares method with iterative de-emphasis of outliers.\n" \
  "\n" \
  "  Optionally, the program can adjust also the other parameters -- end of rally date, width of plateau, rally rate, decay rate --  of one or more bubbles.  The discrete parameters (dates and durations) are adjusted by exhaustive enumeration. (Beware of combinatorial explosion!) The continuous parameters are adjusted by the divided-edge simplex method of non-linear optimization.  The criterion here is minimization of the sum of squares of the log of the ratio of model and gien prices.\n" \
  "\n" \
  "  In this case, the program read two additional bubble parameter files {IN_LO_PARMS_FILE} and {IN_HI_PARMS_FILE}, that define the lower and upper bounds of those parameters for the optimization.  Every parameter value {P_LO} in {IN_LO_PARMS_FILE} must be less than or equal to the correspondng parameter value {P_HI} inf {IN_HI_PARMS_FILE}. The corresponding value {P} in {IN_PARMS_FILE} must be between {P_LO} and {P_HI}.  If {P_LO < P_HI}, then the parameter {P} will be adjusted by the non-linear fitting procedure.\n" \
  "\n" \
  "  Note that the start of decay {DT_INI_DN} is not considered an independent parameter for this purpose, but being derived from the end-of-rally date {DT_FIN_UP} and width of the plateau {WD_PLAT} (difference between {DT_FIN_UP} and {DT_INI_DN}).  Thus if {DT_INI_DN} is different in the files {IN_LO_PARMS_FILE} and {IN_HI_PARMS_FILE}, but the width of the platau is the same, only the date {DT_FIN_UP} will be optimized.\n" \
  "\n" \
  "OPTIONS\n" \
  "  -hrad {HRAD}\n" \
  "    This mandatory argument specifies the half-width of the smoothing Hann window to be applied to the bubble functions before fitting.  The window will span {2*HRAD + 1} samples.  The input prices are assumed to have been smoothed with the same {HRAD} parameter.  If {HRAD} is zero, there will be no smoothing.\n" \
  "\n" \
  "  -inParms {IN_PARMS_FILE}\n" \
  "    This mandatory argument specifies the name of the file that contains the nominal values for the bubble model parameters.\n" \
  "\n" \
  "  -maxLSQIters {MAX_LSQ_ITER}\n" \
  "    This optional argument specifies the maximum number of iterations of the robust least-squares coefficient fitting procedure.  If it is omitted or zero, a simple least squares fit will be used.\n" \
  "\n" \
  "  -maxNLIters {MAX_NL_ITER}\n" \
  "    This optional argument specifies the maximum number of iterations of the non-linear optimization procedure that adjusts the peak dates and rates of the bubbles.  If it is omitted or zero, the non-linear optimization will be suppressed.\n" \
  "\n" \
  "  -inRangeParms {IN_LO_PARMS_FILE} {IN_HI_PARMS_FILE}\n" \
  "    This argument specifies the names of two files that contain the lower bounds and upper bounds or the model parameters, for their non-linear optimization.  The files are read only if {MAX_NL_ITER} is positive; if {MAX_NL_ITER} is zero, this argument is optional, and is ignored if present.\n" \
  "\n" \
  "  -evalRange {EVAL_DT_INI} {EVAL_DT_FIN}\n" \
  "    This optional argument specifies a range of dates to be considered when evaluating the quality of fit for the non-linear fitting procedure. If omitted, the full input date range is considered for that purpose.  The least-squares fitting of the bubble amplitudes always considers the full data range.\n" \
  "\n" \
  "  -truncate {TRUNC_LEVELS}\n" \
  "    This optional argument specifies a relevance threshold for bubbles in output file (a non-negative fraction).  A bubble component will be set to zero in the output file if its value is less than {TRUNC_LEVEL} times the total model price.  If this parameter is not specified, the programs assumes {TRUNC_LEVEL = 0}, i.e. the entire bubble will be written out.\n" \
  "\n" \
  "  -outPrefix {OUT_PREFIX}\n" \
  "    This mandatory parameter specifies the prefix for all output files except {OUT_PRICES_FILE}.\n" \
  "\n" \
  "DOCUMENTATION OPTIONS\n" \
  argparser_help_info_HELP_INFO "\n" \
  "\n" \
  "SEE ALSO\n" \
  "  gawk(1).\n" \
  "\n" \
  "AUTHOR\n" \
  "  Created 2015-03-29 by Jorge Stolfi, IC-UNICAMP.\n" \
  "\n" \
  "MODIFICATION HISTORY\n" \
  "  2015-04-18 Added non-linear optimization.\n" \
  "  2015-04-19 Split out library {libbtc}.\n" \
  "  2015-04-21 Added \"-truncate\" option.\n" \
  "  2015-04-22 Added \"-evalRange\" option.\n" \
  "  2015-04-29 Added {COEF} field in bubble params file.\n" \
  "  2015-04-29 Making width of plateau instead of start of decay an indep parameter.\n" \
  "\n" \
  "WARRANTY\n" \
  argparser_help_info_NO_WARRANTY "\n" \
  "\n" \
  "RIGHTS\n" \
  "  " fit_btc_bubble_model_C_COPYRIGHT ".\n" \
  "\n" \
  argparser_help_info_STANDARD_RIGHTS

#define _GNU_SOURCE
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <bool.h>
#include <affirm.h>
#include <nget.h>
#include <fget.h>
#include <jsfile.h>
#include <lsq_robust.h>
#include <argparser.h>

#include <btc_bubble_t.h>
#include <btc_price_series_read.h>
#include <btc_date_lookup.h>
#include <btc_bubble_parms_read.h>
#include <btc_bubble_nl_opt_check_bubble_parms_in_range.h>
#include <btc_bubble_compute_basis.h>
#include <btc_bubble_nl_opt_adjust_parameters.h>
#include <btc_bubble_parms_write.h>

/* COMMAND-LINE OPTIONS */

typedef struct options_t  /* Program arguments parsed from the command line. */
  { char* inParms;      /* Name of input parameters file. */
    int maxLSQIters;    /* max number of iterations in the robust least-squares fitting procedure. */
    int maxNLIters;     /* Max number of iterations in the non-linear parameter fitting procedure. */
    char* inLoParms;    /* Name of input file with lower bounds of the parameters. */
    char* inHiParms;    /* Name of input file with upper bounds of parameters. */
    int hrad;           /* Half-width of smoothing window. */
    double truncate;    /* Truncation level for bubbles in output file. */
    char* eval_dt_ini;  /* Low ISO date of period to be consiered for NL fit quality evaluation (or NULL). */
    char* eval_dt_fin;  /* High ISO date of period to be consiered for NL fit quality evaluation (or NULL). */
    char* outPrefix;    /* Prefix for names of output files. */
  } options_t;

#define fbbm_MAX_BUBBLES 100
  /* max number of bubbles expected in parameter file. */

/* INTERNAL PROTOTYPES */

options_t *parse_options(int argc, char **argv);
  /* Parses the command line arguments and packs them as an {options_t}. */

int main(int argc,char** argv);
    
void fbbm_write_all_prices
  ( FILE* wr, 
    int nd, 
    char* dt[], 
    double ap[], 
    int nb, 
    btc_bubble_t bp[], 
    double bval[],
    double truncate
  );
  /* Writes to {wr} the results of the fit. Namely, the dates
    {dt[0..nd-1]}, the input prices {ap[0..nd-1]}, the model prices,
    residual and ratio, then the bubble functions {bval[0..nd*nb-1]}
    scaled by the respective coefficients {bp[0..nb-1].coef}. Each
    scaled bubble value is replaced by zero if it is less than
    {truncate} times the model price. */

/* IMPLEMENTATIONS */

int main(int argc, char** argv)
  {
    options_t *o = parse_options(argc, argv);
    
    /* Read the price series: */
    int nd = 0; /* Number of days in series. */
    char** dt = NULL;    /* Dates in analyzed interval, in ISO format, indexed {0..nd-1}. */
    double* ap = NULL; /* Dail mean prices, indexed {0..nd-1}. */
    btc_price_series_read(stdin, &nd, &dt, &ap);
    
    /* Read the initial bubble parameters: */
    int nb = 0; /* Number of bubbles in model. */
    btc_bubble_t* bp = NULL; /* Bubble parametes, indexed {0..nb-1}. */
    btc_bubble_parms_read(o->inParms, nd, dt, &nb, &bp);
    
    /* Get the lower and upper bubble parameters for non-linear fitting: */
    btc_bubble_t* bp_lo = NULL; /* Lower bounds for the bubble parametes, indexed {0..nb-1}. */
    btc_bubble_t* bp_hi = NULL; /* Upper bounds for the bubble parametes, indexed {0..nb-1}. */
    if (o->maxNLIters > 0)
      { int nb_lo = 0;
        btc_bubble_parms_read(o->inLoParms, nd, dt, &nb_lo, &bp_lo);
        demand(nb_lo == nb, "Wrong number of bubbles in lower bounds file");
        int nb_hi = 0;
        btc_bubble_parms_read(o->inHiParms, nd, dt, &nb_hi, &bp_hi);
        demand(nb_hi == nb, "Wrong number of bubbles in upper bounds file");
        btc_bubble_nl_opt_check_bubble_parms_in_range(nb, bp_lo, bp, bp_hi);
      }
      
    /* Get the range to consider in NL optimization: */
    int eval_id_ini = (o->eval_dt_ini == NULL? 0 : btc_date_lookup(nd, dt, o->eval_dt_ini));
    int eval_id_fin = (o->eval_dt_fin == NULL? nd-1 : btc_date_lookup(nd, dt, o->eval_dt_fin));
    
    /* Compute the bubble function basis: */
    double* bval = notnull(malloc(nd*nb*sizeof(double)), "no mem"); /* Smoothed values of all bubbles. */
    
    /* Compute weights for least squares so as to simulate logscale weight: */
    double* wt = notnull(malloc(nd*sizeof(double)), "no mem"); /* Weight of each point. */
    int id;
    for (id = 0; id < nd; id++) { wt[id] = (ap[id] == 0 ? 0.0 : 1.0/ap[id]); }

    /* Adjust the parameters: */
    btc_bubble_nl_opt_adjust_parameters
      ( nd, dt, ap, wt, nb, bp_lo, bp, bp_hi, o->hrad, 
        o->maxLSQIters, o->maxNLIters, 
        eval_id_ini, eval_id_fin,
        o->outPrefix, bval
      );

    /* Write adjusted parameters: */
    btc_bubble_parms_write(o->outPrefix, "-adj", nd, dt, nb, bp);
    
    /* Write model and component bubbles: */
    fbbm_write_all_prices(stdout, nd, dt, ap, nb, bp, bval, o->truncate);
    
    return 0;
  }

void fbbm_write_all_prices
  ( FILE *wr, 
    int nd, 
    char* dt[], 
    double ap[], 
    int nb, 
    btc_bubble_t bp[], 
    double bval[],
    double truncate
  )
  { 
    int id, jb;
    for (id = 0; id < nd; id++)
      { double api = ap[id];
        double* bvali = &(bval[nb*id]);
        /* Recompute the model price: */
        double mpi = 0.0;
        for (jb = 0; jb < nb; jb++)
          { mpi += bp[jb].coef * bvali[jb]; }
        /* Compute the residual and ratio: */
        double resi = (api == 0 ? 0.0 : api - mpi);
        double rati = (api == 0 ? 1.0 : api/mpi);
        /* Print everything: */
        fprintf(wr, "%s 00:00:00 | %18.5f | %18.5f", dt[id], api, mpi);
        fprintf(wr, " | %18.5f | %12.6f", resi, rati);
        for (jb = 0; jb < nb; jb++)
          { double bvj = bp[jb].coef * bvali[jb];
            if (bvj < truncate*mpi) { bvj = 0.0; }
            fprintf(wr, " | %18.5f", bvj);
          }
        fprintf(wr, "\n");
      }
    fflush(wr);     
  }

options_t *parse_options(int argc, char **argv)
  {
    /* Initialize argument parser: */
    argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, PROG_NAME " version " PROG_VERS ", usage:\n" PROG_HELP);
    argparser_set_info(pp, PROG_INFO);
    argparser_process_help_info_options(pp);
    
    /* Allocate the command line argument record: */
    options_t *o = (options_t *)malloc(sizeof(options_t)); 
    
    /* Parse keyword parameters: */
    
    argparser_get_keyword(pp, "-inParms");
    o->inParms = argparser_get_next(pp);
    
    if (argparser_keyword_present(pp, "-maxLSQIters"))
      { o->maxLSQIters = (int)argparser_get_next_int(pp, 1, 100); }
    else
      { o->maxLSQIters = 0; }
    
    if (argparser_keyword_present(pp, "-maxNLIters"))
      { o->maxNLIters = (int)argparser_get_next_int(pp, 1, 100); }
    else
      { o->maxNLIters = 0; }
    
    if (argparser_keyword_present(pp, "-inRangeParms"))
      {  o->inLoParms = argparser_get_next(pp);
         o->inHiParms = argparser_get_next(pp);
      }
    else
      { if (o->maxNLIters > 0)
          { argparser_error(pp, "parameter \"-inRangeParms\" required for non-linear optimization"); }
        o->inLoParms = NULL; 
        o->inHiParms = NULL; 
      }
    
    argparser_get_keyword(pp, "-hrad");
    o->hrad = (int)argparser_get_next_int(pp, 0, 127);
    
    if (argparser_keyword_present(pp, "-truncate"))
      { o->truncate = argparser_get_next_double(pp, 0.00, 1.00); }
    else
      { o->truncate = 0.00; }
    
    if (argparser_keyword_present(pp, "-evalRange"))
      { o->eval_dt_ini = argparser_get_next(pp);
        o->eval_dt_fin = argparser_get_next(pp);
      }
    else
      { o->eval_dt_ini = NULL;
        o->eval_dt_fin = NULL;
      }
    
    argparser_get_keyword(pp, "-outPrefix");
    o->outPrefix = argparser_get_next(pp);
    
    /* Check for spurious arguments: */
    argparser_finish(pp);
    
    return o;
  }
