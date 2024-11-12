#define PROG_NAME "tseries_aggregate"
#define PROG_DESC "Performs principal component analysis on an EEG dataset."
#define PROG_VERS "2013-12-16"

#define tseries_aggregate_C_COPYRIGHT \
  "Copyright © 2013 by the State University of Campinas (UNICAMP)"
/* Last edited on 2023-10-14 21:08:17 by stolfi */
    
#define PROG_HELP \
  "  " PROG_NAME " \\\n" \
  "    -terms {NX} \\\n" \
  "    [ -group { none | binary | fibonacci } ] \\\n" \
  "    [ -constant {X0} ] \\\n" \
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
  "  The program reads from {stdin} a sequence of {NZ} pairs {ID[i],Z[i]}, one per" \
  " line, for {i} in {0..NZ-1}; where {ID[i]} is an arbitrary integer identifier and {Z[i]} is" \
  " a floating-point number.\n" \
  "\n" \
  "  The program then writes to {stdout} a set of tuples {ID[i],Z[i],X[i,1],X[i,2],...,X[1,NX]}" \
  " where each value {X[i,k]} is an average of some values {Z[j]} preceding {Z[i]}.  This file" \
  " is useful, for example, as input to {linear_fit} to find a linear predictor for the time series {Z}.\n" \
  "\n" \
  "  In the simplest case, the variables {X[i,1],X[i,2],...,X[i,NX]} may" \
  " be the {NX} samples {Z[i-1],Z[i-2],...,Z[i-NH]} that precede {Z[i]}, where {NH=NX}.  Optionally," \
  " the program may compute {X[i,1],X[i,2],...X[i,NX]} from a" \
  " larger set of {NH} previous input samples {Z[i-1],Z[i-2],...,Z[i-NH]}, by grouping" \
  " them into {NX} blocks of various sizes, and averaging each block.\n" \
  "\n" \
  "  The program may also insert a constant term {X[i,0]} before the other historical" \
  " variables {X[i,1],X[i,2],...,X[1,NX]}in the output file.\n" \
  "\n" \
  "  In the input file, blank lines and comments starting with \"#\" to the end" \
  " of the line ignored. There is an output line for each input line, except for" \
  " the first {NH} input lines for which there are not enough previous values to compute the group averages.\n" \
  "\n" \
  "OPTIONS\n" \
  "  -terms {NX}\n" \
  "    This mandatory argument is the number {NX} of output variables" \
  " (input samples or block averages) {X[i,1],X[i,2],...X[i,NX]} to produce for each input sample {Z[i]}.\n" \
  "\n" \
  "  -group { none | binary | fibonacci }\n" \
  "    This optional argument specifies whether and how the input samples should" \
  " be condensed into block averages to yield the variables {X[i,k]}.  Specifically:\n" \
  "\n" \
  "      * If the argument is \"none\" (the default), then each group is a single" \
  " input sample, namely {X[i,k]} will be {Z[i-k]} for {k=1,2,...,NX}.\n" \
  "\n" \
  "      * If the argument is \"binary\", then {NH=2^NX-1}, and the previous input" \
  " values are grouped into blocks whose sizes increase by a factor of 2, namely" \
  " 1,2,4,8,...,{2^{NX-1}}.  That is, {X[i,1]=Z[i-1]}, {X[i,2]=(Z[i-2]+Z[i-3])/2}," \
  " {X[i,3]=(Z[i-4]+Z[i-5]+Z[i-6]+Z[i-7])/4}, and so on.\n" \
  "\n" \
  "      * If the argument is \"fibonacci\", then the blocks sizes are the first" \
  " {NX} nonzero Fibonacci numbers, namely 1, 1, 2, 3, 5, 8, 13, .... Thus," \
  " {X[i,1]=Z[i-1]}, {X[i,2]=Z[i-2]}, {X[i,3]=(Z[i-3]+Z[i-4])/2}, {X[i,4]=(Z[i-5]+Z[i-6]+Z[i-7])/3}, and" \
  " so on.\n" \
  "\n" \
  "  -constant {X0}\n" \
  "    This optional argument requests that the constant term {X[i,0]=X0} be inserted in each output" \
  " line before all averages {X[i,1],...X[i,NX]}.\n" \
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
  "  2013-12-16 by J.Stolfi: split off from {linear_fit}.\n" \
  "\n" \
  "WARRANTY\n" \
  argparser_help_info_NO_WARRANTY "\n" \
  "\n" \
  "RIGHTS\n" \
  "  " tseries_aggregate_C_COPYRIGHT ".\n" \
  "\n" \
  argparser_help_info_STANDARD_RIGHTS

#define stringify(x) strngf(x)
#define strngf(x) #x

/* !!! Should take the output format as parameter !!! */
     
#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <values.h>

#include <affirm.h>
#include <argparser.h>
#include <vec.h>
#include <bool.h>
#include <rn.h>
#include <fget.h>

typedef enum
  { tsag_group_NONE,
    tsag_group_BINARY,
    tsag_group_FIBONACCI
  } tsag_group_t;
  /* A grouping criterion. */

typedef struct tsag_options_t
  { int terms;           /* Number of known variables in formula. */
    tsag_group_t group;  /* Group averaging criterion. */
    double constant;     /* Constant term to insert on each output line, or {NAN} if none. */
  } tsag_options_t;
  /* Arguments from command line. */
     
void tsag_define_groups(int NX, tsag_group_t group, int gsz[], int ixg[]);
  /* Sets {gsz[k-1]} to the number of input samples in averaging group {k},
    for {k} in {1..NX}.  Also sets {ixg[0..NX]} such that the input samples 
    in group {k} when predicting {Z[i]} are {Z[i-r]} where {r} is in {ixg[k-1]+1..ixg[k]}. */
   
void tsag_read_data(FILE *rd, int *NXp, int **IDp, double **Zp);
  /* Reads the pairs {ID[i],Z[i]} from {stdin}, in the format expected by the program,
    to EOF, and saves them into newly allocated vectors {ID[0..NZ-1]}, {Z[0..NZ-1]}. 
    Returns the number of samples {NZ} in {*NXp} and the vector {ID,Z} in {*IDP,*Zp}. */
  
void tsag_write_data(FILE *wr, int NZ, int ID[], double Z[], double X0, int NX, int gsz[], int ixg[]);
  /* Writes to {wr} the tuples {ID[i],Z[i],X[i,1],X[i,2],...,X[1,NX]} for {i} in {NH..NZ-1},
    where {NH=ixg[NX]}.  If {X0} is not {NAN}, writes it before {X[i,1]}. */

void tsag_update_group_sums(int NX, double S[], int ixg[], double Z[], int i);
  /* Updates the group sums {S[0..NX-1]} after processing element {Z[i]}. 
    Assumes that group {k} consists of elements {Z[i-r]} where
    {r} is in {ixg[k-1]+1..ixg[k]}, and its sum is stored in {S[k-1]}, for {k} in {1..NX}.
    Updates the sums assuming that all groups shift forward by one position. */

tsag_options_t *tsag_parse_options(int argc, char **argv);

int main (int argc, char **argv)
  {
    bool_t verbose = TRUE;
    
    tsag_options_t *o = tsag_parse_options(argc, argv);
    int NX = o->terms;
    int gsz[NX];    /* The size of group {k} is {gsz[k-1]}, for {k} in {1..NX}. */
    int ixg[NX+1];  /* The input samples in group {k} are {Z[i-r]} where {r} is in {ixg[k-1]+1..ixg[k]}. */
    tsag_define_groups(NX, o->group, gsz, ixg);
    int NH = ixg[NX]; /* Number of previous samples needed to compute each {X} tuple. */
    
    int NZ = -1; /* Number of input samples. */
    int *ID = NULL;
    double *Z = NULL; /* Input values are {Z[0..NZ-1]}. */
    tsag_read_data(stdin, &NZ, &ID, &Z);
    if (verbose) { fprintf(stderr, "read %d data samples\n", NZ); }

    if (verbose) { fprintf(stderr, "writing %d tuples\n", NZ - NH); }
    double X0 = o->constant; /* Constant term, or {NAN}. */
    tsag_write_data(stdout, NZ, ID, Z, X0, NX, gsz, ixg);
    
    return 0;
  }

void tsag_define_groups(int NX, tsag_group_t group, int gsz[], int ixg[])
  {
    assert(NX >= 0);
    ixg[0] = 0;
    int k;
    for (k = 0; k < NX; k++)
      { 
        switch(group)
          {
            case tsag_group_NONE:
              gsz[k] = 1;
              break;
            case tsag_group_BINARY:
              gsz[k] = (k < 1 ? 1 : 2*gsz[k-1]);
              break;
            case tsag_group_FIBONACCI:
              gsz[k] = (k < 2 ? 1 : gsz[k-2] + gsz[k-1]);
              break;
            default:
              assert(FALSE);
          }
        ixg[k+1] = ixg[k] + gsz[k];
      }
  }

void tsag_read_data(FILE *rd, int *NXp, int **IDp, double **Zp)
  {
    int NZ0 = 2048; /* Initial allocation; may be expanded. */
    int_vec_t ID = int_vec_new(NZ0);
    double_vec_t Z = double_vec_new(NZ0);
    int NZ = 0; /* Input values are {Z.e[0..NZ-1]}. */
    
    /* Read the input samples: */
    while (TRUE)
      { bool_t ok = fget_test_comment_or_eol(rd, '#', NULL);
        if (ok) { continue; }
        if (fget_test_eof(rd) { break; }
        int_vec_expand(&ID,NZ);
        double_vec_expand(&Z,NZ);
        ID.e[NZ] = fget_int32(rd);
        Z.e[NZ] = fget_double(rd);
        fget_comment_or_eol(rd, '#', NULL);
        NZ++;
      }
    int_vec_trim(&ID,NZ);
    double_vec_trim(&Z,NZ);
      
    (*NXp) = NZ;
    (*IDp) = ID.e;
    (*Zp) = Z.e;
   }
  
void tsag_write_data(FILE *wr, int NZ, int ID[], double Z[], double X0, int NX, int gsz[], int ixg[])
  {
    int NH = ixg[NX]; /* Number of {Z} values used in all groups. */
    double S[NX]; /* The sum of group {k} is {S[k-1]}, for {k} in {1..NX}. */
    rn_zero(NX, S);
    int i;
    for (i = 0; i < NZ; i++)
      { if (i >= NH)
          { fprintf(wr, "%d", ID[i]);
            fprintf(wr, "  %.8g ", Z[i]);
            if (! isnan(X0)) { fprintf(wr, "  %.8g ", X0); }
            int k;
            for (k = 1; k <= NX; k++)
              { double Zik = S[k-1]/gsz[k-1];
                fprintf(wr, " %.8g", Zik); 
              }
            fprintf(wr, "\n"); 
          }
        tsag_update_group_sums(NX, S, ixg, Z, i);
      }
    fflush(wr);
  }

void tsag_update_group_sums(int NX, double S[], int ixg[], double Z[], int i)
  {
    int naddmax = 16; /* Max samples to add directly. */
    int k;
    for (k = 1; k <= NX; k++)
      { int nxt = i-ixg[k-1];  /* {Z[nxt]} is the oldest sample earlier than group {k}. */
        int fin = i-ixg[k];    /* {Z[fin]} is the oldest sample currently in group {k}. */
        assert(nxt <= i);
        assert(fin < nxt);
        if (nxt < 0)
          { /* Group is entirely before first sample: */
            S[k-1] = 0;
          }
        else if (nxt - fin == 1)
          { /* Single sample {Z[nxt-1]} in group is just replaced by {Z[nxt]}: */
            if (nxt >= 0) { S[k-1] = Z[nxt]; }
          }
        else if (nxt - fin <= naddmax)
          { /* Just add {Z[fin+1..nxt]}: */
            int j;
            double sum = 0;
            for (j = nxt; (j > fin) && (j > 0); j--) { sum += Z[j]; }
            S[k-1] = sum;
          }
        else
          { /* Too many samples to add; {Z[fin]} goes out, {Z[nxt]} goes in: */
            if (fin >= 0) { S[k-1] -= Z[fin]; }
            S[k-1] += Z[nxt];
          }
      }
  }

#define tsag_terms_MAX 30
  /* Maximum number of terms (group averages) in model. */

tsag_options_t *tsag_parse_options(int argc, char **argv)
  {
    /* Initialize argument parser: */
    argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, PROG_NAME " version " PROG_VERS ", usage:\n" PROG_HELP);
    argparser_set_info(pp, PROG_INFO);
    argparser_process_help_info_options(pp);

    /* Allocate the command line argument record: */
    tsag_options_t *o = notnull(malloc(sizeof(tsag_options_t)), "no mem");

    /* Parse keyword parameters: */

    argparser_get_keyword(pp, "-terms");
    o->terms = (int)argparser_get_next_int(pp, 0, tsag_terms_MAX);

    if (argparser_keyword_present(pp, "-group"))
      { if (argparser_keyword_present_next(pp, "none"))
          { o->group = tsag_group_NONE; }
        else if (argparser_keyword_present_next(pp, "binary"))
          { o->group = tsag_group_BINARY; }
        else if (argparser_keyword_present_next(pp, "fibonacci"))
          { o->group = tsag_group_FIBONACCI; }
        else 
          { argparser_error(pp, "invalid grouping strategy"); }
      }
    else
      { o->group = tsag_group_NONE; }

    if (argparser_keyword_present(pp, "-constant"))
      { o->constant = argparser_get_next_double(pp, -DBL_MAX, +DBL_MAX); }
    else
      { o->constant = NAN; }
      
    /* Check for spurious arguments: */
    argparser_finish(pp);

    return o;
  }
