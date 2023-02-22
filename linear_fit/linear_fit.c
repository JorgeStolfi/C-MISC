#define PROG_NAME "linear_fit"
#define PROG_DESC "Performs multivariate linear correlation anaysis on a set of data."
#define PROG_VERS "2013-12-15"

#define linear_fit_C_COPYRIGHT \
  "Copyright © 2013 by the State University of Campinas (UNICAMP)"
/* Last edited on 2023-02-12 11:15:50 by stolfi */
    
#define PROG_HELP \
  "  " PROG_NAME " \\\n" \
  "    -terms {NX} \\\n" \
  "    [ -unitTerm {UNIT} ] \\\n" \
  "    [ -weighted {WTFLAG} ] \\\n" \
  "    [ -varNames {NAME}.. ] \\\n" \
  "    [ -writeFormula {FORMFILE} ] \\\n" \
  "    [ -verbose {VERBOSE} ] \\\n" \
  "    [ -format {FORMAT} ] \\\n" \
  "    " argparser_help_info_HELP " \\\n" \
  "    < {INFILE} \\\n" \
  "    > {OUTFILE}"

#define PROG_INFO \
  "NAME\n" \
  "  " PROG_NAME " - " PROG_DESC "\n" \
  "\n" \
  "SYNOPSIS\n" \
  PROG_HELP "\n" \
  "\n" \
  "DESCRIPTION\n" \
  "  The program reads from {stdin} a sequence of tuples\n" \
  "\n" \
  "     {ID[i]} {Z[i]} {X[i,0]} ... {X[i,NX-1]}\n" \
  "\n" \
  "  one per line, where {ID[i]} is an arbitrary string without blanks (a tuple ID) and the other" \
  " elements are {double}s.\n" \
  "\n" \
  "  It then fits a linear function of the variables {X[i,k]} to {Z[i]}, by least squares, namely\n" \
  "\n" \
  "    {Y[i] = C[0]*X[0] + C[1]*X[1] + ... + C[NX-1]*X[NX-1]}\n" \
  "\n" \
  "  where {C[0],C[1],...,C[NX-1]} are real coefficients.\n" \
  "\n" \
  "  Finally, it writes" \
  " to {stdout} the triplets\n" \
  "\n" \
  "    {ID[i]} {Z[i]} {Y[i]}\n" \
  "\n" \
  "  with the original values {Z[i}} and the corresponding" \
  " fitted values {Y[i]}, one triplet per line.   The formula is also printed to {stderr}.\n" \
  "\n" \
  "  In the input file, blank lines and comments starting with \"#\" to the end" \
  " of the line ignored.\n" \
  "\n" \
  "  Optionally, each input line may include a non-negative weight {W[i]} after the" \
  " target value {Z[i]}:\n" \
  "\n" \
  "     {ID[i]} {Z[i]} {W[i]} {X[i,0]} ... {X[i,NX-1]}\n" \
  "\n" \
  " If there are no explicit weights, the effect is the same" \
  " as specifying the same weight for all records.  Only the relative magnitudes" \
  " of the weights are relevant; multiplying all weights by the same positive" \
  " factor will not change the result.\n" \
  "\n" \
  "  The coefficients {C[J]} are chosen so as to minimize the weighted sum of the" \
  " squared errors, namely {S = SUM{i : W[i]*(Z[i]-Y[i])^2}}. The points with positive" \
  " weight must be sufficiently varied, otherwise the minimum is not unique and" \
  " the coefficients found may be meaningless.\n" \
  "\n" \
  "  Optionally (see the \"-unutTerm\" argument below), the program will fit an affine function, that is, a linear" \
  " function of the {X[i]} plus a constant term {C[NX]}:\n" \
  "\n" \
  "     {Y[i] = C[0]*X[0] + C[1]*X[1] + ... + C[NX-1]*X[NX-1]  +  C[NX]}\n" \
  "\n" \
  "OPTIONS\n" \
  "  -terms {NX}\n" \
  "    This mandatory argument is the number {NX} of independent variables {X[i,0],...X[i,NX-1]} on" \
  " each input tuple.\n" \
  "\n" \
  "  -termNames {NAME}\n" \
  "    This optional keyword is followed by {NX} strings that will be used as variable names in printouts.\n" \
  "\n" \
  "  -unitTerm {UNIT}\n" \
  "    This optional argument specifies whether an extra term {X[i,NX]} with" \
  " value 1 and name \"1\" should be implicitly appended ({UNIT = 1 or \"T\") or" \
  " not ({UNIT = 0 or \"F\") to the given terms {X[0..NX-1]}.   If omitted, it" \
  " defaults to \"-unitTerm F\".\n" \
  "\n" \
  "    With this option, the program fits an affine function of the {X[i]}, instead" \
  " of a purely linear one.  In theory, the same effect could" \
  " be achieved by adding a column of \"1\" values at the end of each line.  In" \
  " practice, besides being more conveninent, the \"-unitTerm\" is handled in a way" \
  " that generaly gives a more accurate result.  The two approaches are mutually" \
  " exclusive, however: \"-unitTerm\" should not be specified if one of" \
  " the independent variables {X[i,j]} is constant for all {j}. \n" \
  "\n" \
  "  -weighted {WTFLAG}\n" \
  "    This optional argument specifies whether the weights {W[i]} are present" \
  " in the input file ({WTFLAG = 1} or \"T\"), or are absent and should be assumed" \
  " to be all 1.0 ({WTFLAG = 0} or \"F\").  In the first case, each" \
  " weight {W[i]} (which must be non-negative) should be specified between" \
  " the target value {Z[i]} and the first independent" \
  " variable {X[i,0]}.  If this parameter is omitted, the program assumes \"-weighted F\".\n" \
  "\n" \
  "  -writeFormula {FORMFILE}\n" \
  "    If this optional parameter is present, the fitted formula is written to" \
  " file called {FORMFILE}.  The first line will have the number of" \
  " variables {NX}, then there will be {NX} lines each with an index {k} (from" \
  " 0 to {NX-1}) and the coefficient {C[k]}.  Finally, two lines with mean and" \
  " the standard deviation of the residual {Z[i]-Y[i]}.\n" \
  "\n" \
  "  -format {FORMAT}\n" \
  "    This optional parameter specifies the format to use to print the values" \
  " of the given data {Z[i]} and the fitted data {Y[i]} in the {stdout} output" \
  " file.  It should be a '%' format specification as would be accepted" \
  " by {printf}, suitable for formatting {double} values, like 'f', 'g', or 'e'.  If not" \
  " specified, the program assumes \"-format '%24.16e'\".\n" \
  "\n" \
  "  -verbose {VERBOSE}\n" \
  "    This optional argument asks for progress and disgnostic ouptut to be printed" \
  " to {stderr} ({VERBOSE = 1} or \"T\"), or suppressed ({VERBOSE = 0} or \"F\"). If this" \
  " parameter is omitted, the program assumes \"-verbose F\".\n" \
  "\n" \
  "DOCUMENTATION OPTIONS\n" \
  argparser_help_info_HELP_INFO "\n" \
  "\n" \
  "SEE ALSO\n" \
  "  gawk(1) gnuplot(1).\n" \
  "\n" \
  "AUTHOR\n" \
  "  Created 2013-12-15 by Jorge Stolfi, IC-UNICAMP.\n" \
  "\n" \
  "MODIFICATION HISTORY\n" \
  "  2014-03-17 by J.Stolfi: added \"-weighted\" option\n" \
  "  2013-12-20 by J.Stolfi: added \"-writeFormula\" option\n" \
  "  2013-12-20 by J.Stolfi: compute mean and deviation of residual.\n" \
  "  2019-12-30 by J.Stolfi: changed type of ID from int to string.\n" \
  "  2023-01-07 by J.Stolfi: added term names and printout sorted by abs coeff.\n" \
  "  2023-01-07 by J.Stolfi: added \"-verbose\" option.\n" \
  "  2023-02-07 by J.Stolfi: added \"-unitTerm\" and \"-format\" options.\n" \
  "\n" \
  "WARRANTY\n" \
  argparser_help_info_NO_WARRANTY "\n" \
  "\n" \
  "RIGHTS\n" \
  "  " linear_fit_C_COPYRIGHT ".\n" \
  "\n" \
  argparser_help_info_STANDARD_RIGHTS

#define stringify(x) strngf(x)
#define strngf(x) #x

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#include <affirm.h>
#include <argparser.h>
#include <jsfile.h>
#include <vec.h>
#include <bool.h>
#include <fget.h>
#include <rn.h>
#include <rmxn.h>
#include <gauss_elim.h>

typedef struct lif_options_t
  { int32_t terms;       /* Number of known variables in formula. */
    char **termNames;    /* Names of independent variables, or printout. */
    bool_t unitTerm;     /* If true, adds an extra term "1" internally. */
    char *writeFormula;  /* File where to write the formula, or {NULL}. */
    bool_t weighted;     /* TRUE iff there is an explicit weight after the target value {Z[i]}. */
    char *format;        /* The format for output {Z[i]} and {Y[i]} values. */
    bool_t verbose;      /* Print progress report and disgnostic output. */
  } lif_options_t;
  /* Arguments from command line. */
     
void lif_read_data(FILE *rd, bool_t weighted, int32_t *NZp, char ***IDp, double **Zp, double **Wp, int32_t NX, double **Xp);
  /* Reads the tuples {ID[i],Z[i],W[i],X[i,0],..,X[i,NX-1]} from {stdin}, in the format expected by the program
    into a newly allocated vectors {Z[0..NZ-1],W[0..NZ-1]} and an array {X[0..NZ*NX-1]}, where {NZ} 
    is the number of input tuples.  The array {X} conceptually has {NZ} rows and {NX} columns.
    It is linearized by rows; element {X[i*NX+k]} is {X[i,k]} for {k} in {0..NX-1}.  
    Returns the number of samples {NZ} in {*NZp} and the arrays {ID,Z,W,X} 
    {*IDp,*Zp,*Wp,*Xp}.
    
    If {weighted} is true, each weight {W[i]} is read from the input line, just after 
    the target value {Z[i]}.  If {weighted} is false, {W[i]} is set to 1 for all {i}.  */
   
void lif_accum_system(double Zi, double Wi, int32_t NX, double Xi[], double *A, double *B);
  /* Updates the least squares system {*A,*B} with one more observation,
    consisting of the given value {Zi = Z[i]} with weight {Wi = W[i]} and the corresponding
    values {Xi[0..NX-1]} of the independent variables. The matrix {A}
    must have {NX} rows and columns, and the vector {B} must have {NX}
    elements. */
 
void lif_build_model
  ( int32_t NZ, 
    double Z[], 
    double W[], 
    int32_t NX, 
    double X[], 
    bool_t unitTerm, 
    int32_t NC, 
    double C[], 
    bool_t verbose
  );
  /* Takes a list {Z[0..NZ-1]} of data samples with weights {W[0..NZ-1]}, and a linearized {NZ} by
    {NX} array {X[0..NZ*NX-1]} of the corresponding values of the
    independent variables. Builds a linear model that gives an
    approximation {Y[i]} to each {Z[i]} by a combination of the
    independent variables {X[i,0..NX-1]}.
    
    The model is defined by a column vector {C[0..NC-1]} of {NC}
    coefficients that are determined by the procedure. 
    
    If {unitTerm} is false, {NC} should be {NX}, and the model is {Y =
    X*C} where {Y[0..NZ-1]} is the vector of predicted values.
    
    If {unitTerm} is true, {NC} must be {NX+1}, and the model is {Y =
    X*C[0..NX-1] + C[NX]}. */
 
void lif_apply_model(int32_t NZ, int32_t NX, double X[], bool_t unitTerm, int32_t NC, double C[], double Y[]);
  /* Takes a series {Z[0..NZ-1]} of data samples, and the coefficients
    {C[0..NC-1]} of a linear predictor. Computes the approximation values
    {Y[0..NZ-1]}. See {lif_build_model} for the meaning of {NC,C,Y}. */
  
void lif_residual_stats(int32_t NZ, double Z[], double W[], double Y[], double *avgP, double *devP);
  /* Computes the mean {*avgP} and deviation {*devP} of the residual {Z[i]-Y[i]}.
    weighted by {W[i]}.  No correction is made for the number of fitted parameters. */
  
void lif_write_model(char *fname, int32_t NX, char *tName[], bool_t unitTerm, int32_t NC, double C[], double avg, double dev);
  /* Writes the fitted model {C[0..NC-1]} and the residual stats {ave,dev} to file "{fname}" in a machine-readable format.
    The term names {tname[0..NX-1]} are printed as #-comments. */
  
void lif_print_model(FILE *wr, int32_t NX, char *tName[], bool_t unitTerm, int32_t NC, double C[], double avg, double dev);
  /* Prints the fitted model {C[0..NC-1]} and the residual stats {ave,dev} to {wr}, in a human-readable format.
    The coefficients are sorted by decreasing absolute value. */

void lif_write_data(FILE *wr, int32_t NZ, char *ID[], double Z[], double Y[], char *fmt);
  /* Writes to {wr} the pairs {ID[i],Z[i],Y[i]} for {I} in {0..NZ-1}. */

lif_options_t *lif_parse_options(int32_t argc, char **argv);

int32_t main (int32_t argc, char **argv)
  {
    lif_options_t *o = lif_parse_options(argc, argv);
    bool_t verbose = o->verbose;
    int32_t NX = o->terms;
    bool_t unitTerm = o->unitTerm;
    int32_t NZ = -1; /* Number of input samples. */
    char **ID = NULL; /* Input data point IDs {ID[0..NZ-1]}. */
    double *Z = NULL; /* Input data values are {Z[0..NZ-1]}. */
    double *W = NULL; /* Weights of data records are {W[0..NZ-1]}. */
    double *X = NULL; /* Input independent variable values are {X[0..NZ*NX-1]}. */
    lif_read_data(stdin, o->weighted, &NZ, &ID, &Z, &W, NX, &X);
    if (verbose) { fprintf(stderr, "read %d data samples {Z[i],X[i,j]}\n", NZ); }

    /* Build the model: */
    int32_t NC = NX + (unitTerm ? 1 : 0);
    double C[NX]; /* Fitted model coefficients. */
    lif_build_model(NZ, Z, W, NX, X, unitTerm, NC, C, verbose);

    /* Apply the model: */
    if (verbose) { fprintf(stderr, "computing fitted values {Y[i]}\n"); }
    double *Y = notnull(malloc(NZ*sizeof(double)), "no mem");
    lif_apply_model(NZ, NX, X, unitTerm, NC, C, Y);
    
    /* Apply the model: */
    if (verbose) { fprintf(stderr, "computing residual stats\n"); }
    double avg, dev; /* Statistics of residual. */
    lif_residual_stats(NZ, Z, W, Y, &avg, &dev);
    
    if (verbose) 
      { lif_print_model(stderr, NX, o->termNames, unitTerm, NC, C, avg, dev); }
    
    if (o->writeFormula != NULL) 
      { if (verbose) { fprintf(stderr, "writing formula to file\n"); }
        lif_write_model(o->writeFormula, NX, o->termNames, unitTerm, NC, C, avg, dev);
      }
      
    lif_write_data(stdout, NZ, ID, Z, Y, o->format);
    
    return 0;
  }

void lif_read_data(FILE *rd, bool_t weighted, int32_t *NZp, char ***IDp, double **Zp, double **Wp, int32_t NX, double **Xp)
  {
    int32_t NZ0 = 2048; /* Initial allocation; may be expanded. */
    string_vec_t ID = string_vec_new(NZ0);
    double_vec_t Z = double_vec_new(NZ0);
    double_vec_t W = double_vec_new(NZ0);
    double_vec_t X = double_vec_new(NX*NZ0);

    /* Read the input samples: */
    int32_t NZ = 0; /* Input values are {Z.e[0..NZ-1]}. */
    while (TRUE)
      { bool_t ok = fget_test_comment_or_eol(rd, '#');
        if (ok) { continue; }
        if (fget_test_eof(rd)) { break; } 
        string_vec_expand(&ID,NZ);
        double_vec_expand(&Z,NZ);
        double_vec_expand(&W,NZ);
        double_vec_expand(&X,(NZ+1)*NX-1);
        ID.e[NZ] = fget_string(rd);
        Z.e[NZ] = fget_double(rd);
        if (weighted)
          { W.e[NZ] = fget_double(rd);
            demand(W.e[NZ] >= 0, "weights must be non-negative");
          }
        else
          { W.e[NZ] = 1.0; }
        double *Xi = &(X.e[NZ*NX]); /* Row of {X} for this input line. */
        for (int32_t k = 0; k < NX; k++) { Xi[k] = fget_double(rd); }
        NZ++;
        fget_comment_or_eol(rd, '#');
      }
    string_vec_trim(&ID,NZ);
    double_vec_trim(&Z,NZ);
    double_vec_trim(&W,NZ);
    double_vec_trim(&X,NZ*NX);
      
    (*NZp) = NZ;
    (*IDp) = ID.e;
    (*Zp) = Z.e;
    (*Wp) = W.e;
    (*Xp) = X.e;
   }
  
void lif_build_model
  ( int32_t NZ, 
    double Z[], 
    double W[], 
    int32_t NX, 
    double X[], 
    bool_t unitTerm, 
    int32_t NC, 
    double C[], 
    bool_t verbose
  )
  {
    double Xave[NX];     /* If {unitTerm}, average value of each {X[i]}. */
    double Zave = NAN;   /* If {unitTerm}, average value of {Z}. */
    if (unitTerm)
      { demand(NC == NX + 1, "{NC} should be {NX+1}");
        if (verbose) { fprintf(stderr, "removing variable averages...\n"); }
        /* Compute the weighted sums of {X[k],Z} in {Xave[k],Zave}: */
        for (int32_t k = 0; k < NX; k++) { Xave[k] = 0; }
        Zave = 0;
        double sumW = 0;
        for (int32_t i = 0; i < NZ; i++)
          { double Wi = W[i];
            demand(Wi >= 0, "invalid weight {W[i]}");
            double *Xi = &(X[i*NX]); 
            for (int32_t k = 0; k < NX; k++) 
              { Xave[k] += Wi*Xi[k]; }
            Zave += Wi*Z[i];
            sumW += Wi;
          }
        /* Convert weighted sums to weighted averages: */
        if (sumW > 0)
          { for (int32_t k = 0; k < NX; k++) { Xave[k] /= sumW; }
            Zave /= sumW;
          }
        /* Subtract averages from all variables: */
        for (int32_t i = 0; i < NZ; i++)
          { double *Xi = &(X[i*NX]); 
            for (int32_t k = 0; k < NX; k++) { Xi[k] -= Xave[k]; }
            Z[i] -= Zave;
          }
      }
    else
      { demand(NC == NX, "{NC} should be {NX}"); }

    /* Build the least squares system: */
    if (verbose) { fprintf(stderr, "building linear system matrices...\n"); }
    double *A = rmxn_alloc(NX, NX);
    double B[NX];
    rn_zero(NX, B);
    rmxn_zero(NX, NX, A);
    
    /* Accumulate the matrix and vector of the least squares system: */
    for (int32_t i = 0; i < NZ; i++)
      { double *Xi = &(X[i*NX]); 
        lif_accum_system(Z[i], W[i], NX, Xi, A, B);
      }
    if (verbose) { rmxn_gen_print2(stderr, NX,  NX, A,  1, B,  "%+18.9f", "  ","\n  ","\n", "[ "," "," ]", "  "); }
    
    if (verbose) { fprintf(stderr, "solving system...\n"); }
    double tiny = 1.0e-8;
    int32_t rank = gsel_solve(NX, NX, A, 1, B, C, tiny);
    if (rank < NX) { fprintf(stderr, "!! warning: system rank = %d\n", rank); }
    if (verbose) { rmxn_gen_print(stderr, NX, 1, C, "%+18.9f", "  ","\n  ","\n", "[ "," "," ]"); }
    
    if (unitTerm)
      { if (verbose) { fprintf(stderr, "adding variable averages to indep term...\n"); }
        double Cunit = Zave;
        for (int32_t k = 0; k < NX; k++) { Cunit -= C[k]*Xave[k]; }
        C[NC-1] = Cunit;
      }
    
    free(A);
  }
  
void lif_apply_model(int32_t NZ, int32_t NX, double X[], bool_t unitTerm, int32_t NC, double C[], double Y[])
  { 
    assert(NC == NX + (unitTerm ? 1 : 0));
    for (int32_t i = 0; i < NZ; i++)
      { double *Xi = &(X[i*NX]); 
        double sum = 0;
        for (int32_t k = 0; k < NX; k++) { sum += C[k]*Xi[k]; }
        if (unitTerm) { sum += C[NX]; }
        Y[i] = sum;
      }
  }

void lif_residual_stats(int32_t NZ, double Z[], double W[], double Y[], double *avgP, double *devP)
  {
    double sumWD = 0;
    double sumW = 0;
    for (int32_t i = 0; i < NZ; i++) { sumWD += W[i]*(Z[i] - Y[i]); sumW += W[i]; }
    double avg = sumWD/sumW;
    demand(sumW > 0, "total weight is zero");
    double sumWD2 = 0;
    for (int32_t i = 0; i < NZ; i++) { double Di = (Z[i] - Y[i]) - avg;  sumWD2 += W[i]*Di*Di; }
    double dev = sqrt(sumWD2/sumW);
    (*avgP) = avg;
    (*devP) = dev;
  }

void lif_write_model(char *fname, int32_t NX, char *tName[], bool_t unitTerm, int32_t NC, double C[], double avg, double dev)
  {
    assert(NC == NX + (unitTerm ? 1 : 0));
    FILE *wr = open_write(fname, TRUE);
    fprintf(wr, "%d\n", NX);
    for (int32_t k = 0; k < NX; k++) 
      { fprintf(wr, "%3d %+24.16e", k, C[k]); 
        if (tName != NULL)
          { fprintf(wr, " # %s", tName[k]); }
        fprintf(wr, "\n");
      }
    if (unitTerm) { fprintf(wr, "%3d %+24.16e # 1\n", NX, C[NX]); }
    fprintf(wr, "    %+24.16e\n", avg);
    fprintf(wr, "    %+24.16e\n", dev);
    fclose(wr);
  }
  
void lif_print_model(FILE *wr, int32_t NX, char *tName[], bool_t unitTerm, int32_t NC, double C[], double avg, double dev)
  {
    assert(NC == NX + (unitTerm ? 1 : 0));

    /* Index sort of {C[0..NX-1]} by absolute value of coefficient (leave {C[NX]} at end): */
    int32_t ix[NC];
    for (int32_t k = 0; k < NC; k++) { ix[k] = k; }
    for (int32_t k = 0; k < NX; k++) 
      { /* Largest coeffs are {C[ix[i]]} for {i} in {0..k-1}. */
        /* Find largest unsorted coeff: */
        int32_t imax = k;
        for (int32_t i = k+1; i < NX; i++)
          { if (fabs(C[ix[i]]) > fabs(C[ix[imax]])) { imax = i; } }
        if (imax != k)
          { int32_t t = ix[k]; ix[k] = ix[imax]; ix[imax] = t; }
      }
  
    fprintf(wr, "fitted model:\n");
    for (int32_t k = 0; k < NC; k++) 
      { int32_t i = ix[k];
        fprintf(wr, "  %+14.9f", C[i]); 
        if (i < NX)
          { fprintf(wr, " * X[%2d]", i);
            if (tName != NULL) { fprintf(wr, " # %s", tName[i]); }
          }
        else
          { fprintf(wr, "%8s", "");
            fprintf(wr, " # 1");
          }
        fprintf(wr, "\n"); 
      }
    fprintf(wr, "  %+14.9f         # avg residual\n", avg);
    fprintf(wr, "  %+14.9f * RND() # dev residual\n", dev);
  }
    
void lif_accum_system(double Zi, double Wi, int32_t NX, double Xi[], double *A, double *B)
  { 
    for (int32_t k = 0; k < NX; k++)
      { B[k] += Wi*Zi*Xi[k];
        for (int32_t j = 0; j < NX; j++) { A[k*NX + j] += Wi*Xi[k]*Xi[j]; }
      }
  }
  
void lif_write_data(FILE *wr, int32_t NZ, char *ID[], double Z[], double Y[], char *fmt)
  {
    for (int32_t i = 0; i < NZ; i++)
      { fprintf(wr, "%s ", ID[i]);
        fprintf(wr, fmt, Z[i]);
        fprintf(wr, " ");
        fprintf(wr, fmt, Y[i]);
        fprintf(wr, "\n");
      }
    fflush(wr);
  }

#define lif_terms_MAX 50
  /* Maximum number of terms (group averages) in model. */

lif_options_t *lif_parse_options(int32_t argc, char **argv)
  {
    /* Initialize argument parser: */
    argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, PROG_NAME " version " PROG_VERS ", usage:\n" PROG_HELP);
    argparser_set_info(pp, PROG_INFO);
    argparser_process_help_info_options(pp);

    /* Allocate the command line argument record: */
    lif_options_t *o = notnull(malloc(sizeof(lif_options_t)), "no mem");

    /* Parse keyword parameters: */

    argparser_get_keyword(pp, "-terms");
    o->terms = (int32_t)argparser_get_next_int(pp, 1, lif_terms_MAX);
    
    if (argparser_keyword_present(pp, "-termNames"))
      { int32_t NX = o->terms;
        o->termNames = notnull(malloc(NX*sizeof(char*)), "no mem");
        for (int32_t k = 0; k < NX; k++) 
          { char *name = argparser_get_next_non_keyword(pp);
            o->termNames[k] = name;
          }
      }
    else
      { o->termNames = NULL; }

    if (argparser_keyword_present(pp, "-unitTerm"))
      { o->unitTerm = argparser_get_next_bool(pp); }
    else
      { o->unitTerm = FALSE; };

    if (argparser_keyword_present(pp, "-weighted"))
      { o->weighted = argparser_get_next_bool(pp); }
    else
      { o->weighted = FALSE; };
    
    if (argparser_keyword_present(pp, "-writeFormula"))
      { o->writeFormula = argparser_get_next_non_keyword(pp); }
    else
      { o->writeFormula = NULL; }

    if (argparser_keyword_present(pp, "-format"))
      { o->format = argparser_get_next_non_keyword(pp);
        /* Minimal sanity check: */
        if (strchr(o->format, '%') == NULL) { argparser_error(pp, "no '%' in format"); }
      }
    else
      { o->format = "%24.16e"; };

    if (argparser_keyword_present(pp, "-verbose"))
      { o->verbose = argparser_get_next_bool(pp); }
    else
      { o->verbose = FALSE; };
    
    /* Check for spurious arguments: */
    argparser_finish(pp);

    return o;
  }
