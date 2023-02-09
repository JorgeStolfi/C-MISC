#define PROG_NAME "fit_piecewise"
#define PROG_DESC "Piecewise polynomial fit and prediction to a time series."
#define PROG_VERS "2013-12-15"

#define fit_piecewise_C_COPYRIGHT \
  "Copyright © 2013 by the State University of Campinas (UNICAMP)"
/* Last edited on 2015-10-18 02:56:10 by stolfilocal */
    
/* !!! Provide optional format strings for {X}, {W}, {Z}. !!! */
/* !!! Use {(X[i]-XC)/XR} instead of just {(X[i]-XC)/XR} as the poly argument. !!! */
/* !!! Use Bernstein or orthogonal polynomials insted of {x^u}. */

#define PROG_HELP \
  "  " PROG_NAME " \\\n" \
  "    -maxDegree {MAXDEG} \\\n" \
  "    -penalty break {PTBREAK} skip {PTSKIP} \\\n" \
  "    [ -weighted {WTFLAG} ] \\\n" \
  "    [ -writeFunction {FORMFILE} ] \\\n" \
  "    [ -resample {NR} {XINI} {XSTEP} {SAMPLEFILE} ] \\\n" \
  "    [ -predict {XDELTA} {PREDFILE} ] \\\n" \
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
  "  The program reads from {stdin} a sequence of {NP >= 2} tuples {ID[i],X[i],Z[i],W[i]}, one per" \
  " line, where {ID[i]} is an arbitrary integer (the tuple's ID), {X[i]} is an argument value, {Z[i]} is" \
  " a target valuem and {W[i]} is a non-negative weight, for {i} in {0..NP-1}.  The program" \
  " then fits a piecewise polynomial function {F} of maximum degree {MAXDEG} to those values, by" \
  " least squares and dynamic programming.  Finally, it writes to {stdout} the" \
  " tuples {ID[i],X[i],Z[i],W[i],Y[i]}, where {Y[i]} is the fitted function's" \
  " value corresponding to {X[i]}, one tuple per line.\n" \
  "\n" \
  "  In the input file, blank lines and comments starting with \"#\" to the end" \
  " of the line ignored.  The lines must be sorted so that the {X[i]} are strictly" \
  " increasing.  If the spacing between the arguments {X[i]} is too" \
  " irregular, the fitted polynomials may be quite strange.\n" \
  "\n" \
  "   The formulas for each segment of the fitted function are printed" \
  " to {stderr}, and, optionally to a separate file \"{FORMFILE}\".  The program" \
  " can be asked also to produce a file \"{SAMPLEFILE}\" with the fitted function {F} re-sampled" \
  " at regularly spaced arguments.  Finally, the program can be asked to try to" \
  " predict the future value of the fitted function {F(X[i]+XDELTA)} from" \
  " the data tuples {0..i} only, and write these predictions to a file \"{PREDFILE}\".\n" \
  "\n" \
  "\n" \
  "  The weight {W[i]} may be omitted in the input lines, in which case the effect is the same" \
  " as specifying the same weight for all records.  Only the relative magnitudes" \
  " of the weights are relevant; multiplying all weights by the same positive" \
  " factor will not change the result.\n" \
  "\n" \
  "FITTED FUNCTION\n" \
  "  The breaks between two consecutive pieces always occurs at /potential breaks/" \
  " which are midpoint between two successive argument values, {B[i]=(X[i-1]+X[i])/2}.  Additional" \
  " potential breaks {B[0]} and {B[NP]} are are assumed before {X[0]} and after {X[NP-1]}," \
  " namely {B[0] = ((3*X[0]-X[1])/2} and {B[NP] = (3*X[NP-1]-X[NP-2])/2}.  Note that each" \
  " input argument {X[i]} lies between {B[i]} and {B[i}1]}. The /actual breaks/ of the fitted" \
  " function {F()} are a sublist {B[r[k]]} of these potential breaks, for {k} in {0..NS} where {NS} is" \
  " the number of segments of {F()}.  These" \
  " actual breaks always include {B[r[0]] = B[0]} and {B[r[NS]] = B[NP]}.\n" \
  "\n" \
  "  For each {k} in {0..NS-1}, segment {k} of the fitted function is a polynomial {FS[k]()} defined" \
  " over the interval {D[k]} between two potential breaks {B[r[k]]} and {B[r[k+1]]}, with {r[k+1] > r[k]}.  The" \
  " degree of the polynomial fitted to that segment is {g[k] = min{MAXDEG,r[k+1]-r[k]-2}}; if {g[k] < 0}, then" \
  " there is no polynomial and any points {X[i]} in the domain {D[k]} are simply treated as outliers, independently" \
  " of their values {Z[i]}.\n" \
  "\n" \
  "SEGMENT FITTING\n" \
  "  Once the initial and final breaks {B[r[k]]} and {B[r[k+1]]} are chosen, segment {FS[k]()} of" \
  " the function is determined by the implied degree {g[k]} and the points {(X[i],Z[i])} and the" \
  " weight {W[i]} with {r[k] <= i < r[k+1]}.\n" \
  "\n" \
  "  If {g[k]} is non-negative, the segment {FS[k]()} is determined by fitting a polynomial of" \
  " and the chosen degree {g[k]} to those points, bu the weighted least squares method.  The fitting" \
  " process returns for each point a probability {PROUT[i]} that the value {Z[i]} is an /outlier/, a value" \
  " that cannot be explained as the trend value {FS[k](X[i])} plus some Gaussian noise.  The effective" \
  " weight of each point {i} used in the least squares fitting is {WEF[i] = W[i]*sqrt(1-PROUT[i])}, so" \
  " points that look like outliers have less influence on the fit, even if {W[i]} is maximum.\n" \
  "\n" \
  "  In particular, if {g[k]<0}, then {FS[k]()} is not defined (its value is always {NAN}), and all points in that" \
  " range (usually only one) are assumed to be outliers, with {PROUT[i] = 1} and thus {WEF[i] = 0}.\n" \
  "\n" \
 "PENALTY AND BREAK CHOOSING\n" \
  "   The /fitting penalty/ of a segment is defined as the weighted sum of the squared" \
  " deviations from the fitted trend, namely" \
  " {QFIT[k] = SUM{ WEF[i]*(Z[i] - FS[k](X[i]))^2 : i \\in r[k]..r[k+1]-1}. In addition," \
  " that segment has an assigned /outlier penalty/ defined as" \
  " {QOUT[k] = PTSKIP*SUM{ W[i]*PROUT[i] : i \\in r[k]..r[k+1]-1}.  In particular," \
  " if {g[k]} is negative, {QFIT[k]} is zero and {QOUT[k]} is {PTSKIP} times the sum of the weights {W[i]}.\n" \
  "\n" \
  "   Each potential break {B[i]} also has a /break penalty/ {QBRK[i]}.  For {i}" \
  " in {1..NP-1}, {QBRK[i]} is {PTBREAK*WHAF[i]} where {WHAF[i]} is the harmonic" \
  " mean of the weights of the data points adjacent to the break, namely {W[i-1]} and" \
  " {W[i-1]}.  The break penalties {QBRK[0]} and {QBRK[NP]} for {BR[0]} and" \
  " {BR[NP]} are zero.\n" \
  "\n" \
  "   The /total penalty/ for the function {F()} is defined" \
  " as {Q(F) = SUM{ QFIT[k] + QSKP[k] + QBRK[r[k]] : k \\in 0..NS-1 }}. The program chooses" \
  " the breaks so as to minimize this pernalty.\n" \
  "\n" \
  "OPTIONS\n" \
  "  -maxDegree {MAXDEG}\n" \
  "    This mandatory argument specifies the maximum degree {MAXDEG} of each segment of the fitted function.\n" \
  "\n" \
  "  -penalty break {PTBREAK} skip {PTSKIP}\n" \
  "    This mandatory argument specifies the {PTBREAK} and {PTSKIP} parameters for" \
  " computing the total penalty of a candidate function {F()}.\n" \
  "\n" \
  "  -weighted {WTFLAG}\n" \
  "    This optional argument specifies whether the weights {W[i]} are present" \
  " in the input file ({WTFLAG = 1} or \"T\"), or are absent and should be assumed" \
  " to be all 1.0 ({WTFLAG = 0} or \"F\").  In the first case, each" \
  " weight {W[i]} (which must be non-negative) should be specified between" \
  " the target value {Z[i]} and the first independent" \
  " variable {X[i,0]}.  If this parameter is omitted, the program assumes \"weighted F\".\n" \
  "\n" \
  "  -writeFunction {FORMFILE}\n" \
  "    If this optional parameter is present, the fitted function is written to" \
  " file called \"{FORMFILE}\", in addition to {stderr}.  The first line will have the number of" \
  " segments {NS}, then there will be {NS} lines each with an index {k} from" \
  " 0 to {NS-1}.  Each line will have the format\n" \
  "      \"{r[k]} {B[r[k]]} {r[k+1]} {B[r[k+1]]} {XC[k]} {QFIT[k]} {QOUT[k]} {g[k]} {C[k,0]}.. {C[k,g[k]]}\"\n" \
  "    where {XC[k] = (B[r[k]] + B[r[k+1]])/2}, and ({C[k,t]} is the coefficient of {(X-XC[k])^t} in {FS[k]()}.   If" \
  " this argument is omitted, the formula is written only to {stderr}.\n" \
  "\n" \
  "  -resample {NR} {XINI} {XSTEP} {SAMPLEFILE}\n" \
  "    If this optional argument is present, the program will evaluate the fitted " \
  " function {F()} at the sample points {XR[j] = XINI+ j*XSTEP} for {j} in {0..NR-1}, and" \
  " write the results to file \"{SAMPLEFILE}\".  Each line will have the" \
  " format \"{j} {XR[j]} {F(XR[j])}\".  Whenever the argument {XR[j]} is outside" \
  " the domain {[B[0]..B[NP]]}, the value of {F(XR[j])} will be {NAN}.  If this" \
  " argument is omitted, no resampling is done.\n" \
  "\n" \
  "  -predict {XDELTA} {PREDFILE}\n" \
  "    If this optional argument is present, the program will try to predict the" \
  " function value at argument {XF[i] = X[i] + XDELTA} for each input argument {X[i]}, using" \
  " only the partial best fit {F'[i]()} that considers only data tuples {0..i}.  The result" \
  " is written to file \"{PREDFILE}\"; each line will have the" \
  " format \"{i} {X[i]} {XF[i]} {F'[i](XF[i]})}.  Note that {XF[i]} may" \
  " be greater than the next break {B[i+1]}, in particular" \
  " greater than {B[NP]}.\n" \
  "\n" \
  "DOCUMENTATION OPTIONS\n" \
  argparser_help_info_HELP_INFO "\n" \
  "\n" \
  "SEE ALSO\n" \
  "  gawk(1) gnuplot(1).\n" \
  "\n" \
  "AUTHOR\n" \
  "  Created 2014-03-17 by Jorge Stolfi, IC-UNICAMP.\n" \
  "\n" \
  "MODIFICATION HISTORY\n" \
  "  2014-03-17 by J.Stolfi: created\n" \
  "\n" \
  "WARRANTY\n" \
  argparser_help_info_NO_WARRANTY "\n" \
  "\n" \
  "RIGHTS\n" \
  "  " fit_piecewise_C_COPYRIGHT ".\n" \
  "\n" \
  argparser_help_info_STANDARD_RIGHTS

#define stringify(x) strngf(x)
#define strngf(x) #x

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#include <affirm.h>
#include <argparser.h>
#include <jsfile.h>
#include <jsmath.h>
#include <vec.h>
#include <bool.h>
#include <fget.h>
#include <rn.h>
#include <in.h>
#include <rmxn.h>

typedef struct fpc_penalty_options_t
  { double ptBreak;  /* Penalty for a break in the trend.  */
    double ptSkip;   /* Penalty factor for treating a point as an outlier.  */
  } fpc_penalty_options_t;
  /* Penalty parameters. */

typedef struct fpc_predict_options_t
  { double xDelta;    /* Argument interval from prediction made to prediction valid for.  */
    char *fileName;   /* File where to write the predictions. */
  } fpc_predict_options_t;
  /* Prediction parameters. */

typedef struct fpc_resample_options_t
  { int  NR;          /* Number of sampling points.  */
    double XIni;      /* First sample.  */
    double XStep;     /* Sampling step.  */
    char *fileName;   /* File where to write the predicted values. */
  } fpc_resample_options_t;
  /* Resampling parameters. */

typedef struct fpc_options_t
  { int terms;             /* Number of known variables in formula. */
    int maxDegree;            /* Max polynomial degree.*/
    bool_t weighted;       /* TRUE iff there is an explicit weight after the target value {Z[i]}. */
    char *writeFunction;   /* File where to write the formula, or {NULL}. */
    fpc_penalty_options_t *penalty; 
    fpc_resample_options_t *resample; 
    fpc_predict_options_t *predict; 
  } fpc_options_t;
  /* Arguments from command line. */

fpc_options_t *fpc_parse_options(int argc, char **argv);
fpc_penalty_options_t *fpc_parse_penalty_options(argparser_t *pp);
fpc_resample_options_t *fpc_parse_resample_options(argparser_t *pp);
fpc_predict_options_t *fpc_parse_predict_options(argparser_t *pp);
  /* Procedures to parse the command line options. */

void fpc_read_data(FILE *rd, bool_t weighted, int *NPp, int **IDp, double **Xp, double **Zp, double **Wp);
  /* Reads the tuples {ID[i],X[i],Z[i],W[i]} from {stdin}, in the format expected by the program
    into a newly allocated vectors {X[0..NP-1],Z[0..NP-1],W[0..NP-1]}, where {NP} 
    is the number of input tuples.  Returns the number of tuples {NP} in {*NPp} and the arrays {ID,X,Z,W} 
    {*IDp,*Xp,*Zp,*Wp}.
    
    If {weighted} is true, each weight {W[i]} is read from the input line, just after 
    the target value {Z[i]}.  If {weighted} is false, {W[i]} is set to 1 for all {i}.  */
       
typedef struct fpc_segment_t
  {
    int degree;            /* Degree of fitted segment, or {-1} if undefined. */
    double QFIT;      /* Fitting penalty. */
    double QOUT;      /* Outlier penalty. */
    double *coeff;    /* Coefficient vector (with {degree+1} elements). */
  } fpc_segment_t;
  /* Describes a potential segment of the fitted function, spanning the gap between breaks {i}
    and {j} (supplied externally). */

void fpc_choose_potential_breaks(int NP, double X[], double B[]);
  /* Computes the potential breaks {B[0..NP]} from the original arguments {X[0..NP-1]}, such that 
    {B[i] < X[i] < B[i+1]} for all valid {i}. */
   
void fpc_fit_function
  ( int NP, 
    double X[], 
    double Z[], 
    double W[], 
    int gMax, 
    fpc_penalty_options_t *pto, 
    double B[], 
    fpc_segment_t *A[], 
    int *NSP, 
    int **rP
  );
  /* Fits a piecewise polynomial function to the data {X[0..NP-1],Z[0..NP-1]) with data
    point weights {W[0..NP-1]}, where each segment has degree at most {gMax}, with the
    penatly options {pto} and potential breaks {B[0..NP]}.

    As part of its work, the preocedure partially fills the triangular tableau
    {A[0..NA-1]} with potential segments of the function. Namely, entry in line {i} and
    column {j} of the tableau, with {0<=i<j<=NP}, if not null, will be the best fitting
    polynomial for the data between breaks {B[i]} and {B[j]}
    (the triples {X[t],Z[t],W[t]} with {t} in {i..j-1}).
    
    The fitted function {F} is defined by the number of segments {NS} and the list
    {r[0..NS]} of the indices of the actual breaks, both chosen by the procedure. Namely,
    for each {k} in {0..NS-1]}, segment {k} of {F} will extend from {B[r[k]]} to
    {B[r[k+1]]}, and {0 = r[0] <= r[k] < r[k+1] <= r[NS] = NP}. The segment will be
    described by the entry of the tableau in row {r[k]} and column {r[k+1]}, which will be
    non-null.  The procedure allocates the vector {r} and returns {NS,r} in {*NSP,*rP}. */
    
fpc_segment_t *fpc_fit_segment(int NP, double X[], double Z[], double W[], int degree, double B[], int i, int j, double ptSkip);
  /* Finds the best-fitting segment to the data points {X[i..j-1],Z[i..j-1]} with weights {W[i..j-1]}, assumed
    to lie between the breaks {B[i]} and {B[j]}.   Uses {ptSkip*W[u]*PROUT[u]} as the outlier penalty. */

double fpc_break_penalty_weight(int NP, double W[], int j);
  /* Returns the relative penalty for using break number {j}, assumed to lie between data points {j-1} 
    and {j}. The result is the harmonic mean of the weights of those two points; it is zero if
    either point does not exist (that is, if {j} is {0} or {NP}). */

void fpc_write_function(char *fName, int NP, double B[], fpc_segment_t *A[], int NS, int r[]);
  /* Writes to a file called {fName} a machine-readable description of the fitted function, in the format 
    documented above. */
 
void fpc_print_function(FILE *wr, int NP, double B[], fpc_segment_t *A[], int NS, int r[], double avg, double dev);
  /* Writes to {wr} a human-readable description of the fitted function. */
   
double fpc_eval_function(double XArg, int NP, double B[], fpc_segment_t *A[], int NS, int r[], int *kHintP);
  /* Evaluates the function {F} defined by the breaks {B[0..NP]}, the tableau {A}, and the
    index vector {r[0..NS]} at the argument value {XArg}. Will return {NAN} if {XArg <
    B[0]} or {XArg > B[NP]}, or if {XArg} is in the domain of a segment that is undefined.
    
    If {XArg} coincides with a break {B[i]}, the result is the average of the values given
    by the two segments on either side; except that if one one of them is {NAN}, the
    result is the other one.
    
    The parameter {kHinP} is passed on to the function {fpc_locate_segment} (quod videt)
    to locate the proper segment(s) {FS[k]} of {F}. */

int fpc_locate_segment(double XArg, int NP, double B[], int NS, int r[], int *kHintP);
  /* Expects that {B[0..NP]} is a list of strictly increasing breaks, and {r[0..NS]} is 
    a list of strictly increasing indices with {r[0]==0} and {r[NS]==NP}. 
    Returns an index {k} such {B[r[k]] <= XArg <= B[r[k+1]]}.
    
    Returns {-1} if there is no such {k}, that is, if {Xarg < B[0]} or {XArg > B[NP]}. The
    result is ambiguous if {Xarg} coincides with some break {B[r[k]]} and this {k} is
    neither 0 nor {NS-1}.
    
    The procedure does a binary search of {XArg} among the breaks {B[0..NP]} to find the
    index {k}. The argument {kHint}, if not null is assumed to be the address of a
    variable whose value is a hint to the value of {k}. If {*kHint} is in {0..NS-1}, the
    procedure will try {k=*kHint} and the two adjacent values, before doing the binary
    search. If the procedure succeeds, it updates {*kHintP} with {k}, otherwise leaves it
    unchanged. */

double fpc_eval_segment_by_index(double XArg, int k, int NP, double B[], fpc_segment_t *A[], int NS, int r[]);
  /* Evaluates at argument {XArg} the segment {FS[k]} of the function {F} defined by
   {B[0..NP]}, the tableau {A}, and the break index list {r[0..NS]}. 
   
   Returns {NAN} {k} is not in {0..NS-1}, or if the selected segment is undefined at
   {XArg}. If neither condition is true, the segment's underlying polynomial is evaluated
   normally, even if {XArg} is outside the domain of the segment, namely the closed
   interval {[ B[r[k]] _ B[r[k+1]] ]}, */

double fpc_eval_segment(double x, fpc_segment_t *FS);
  /* Evaluates the segment {FS}, whose nominal domain is {D = [Bi _ Bj]}, at the 
    domai-normalized argument {x}. The original domain of the segment corresponds
    to the normalized domain {[-1 _ +1]}, but the underlying function is
    usually evaluated even when {x} lies outside it. */ 

double fpc_normalize_arg(double XArg, double Bi, double Bj);
  /* Applies to {XArg} an affine (shift and scale) map that takes the interval {[Bi _ Bj]}
    to the interval {[-1 _ +1]}. */

void fpc_write_segment(FILE *wr, int NP, int i, int j, double B[], fpc_segment_t *A[]);
  /* Writes to {wr} a description of one segment of the fitted funtion {F}, spanning
     breaks {B[i]} to {B[j]}.  The segment's description is taken from the tableau {A}. */

void fpc_print_segment(FILE *wr, int NP, int k, int i, int j, double B[], fpc_segment_t *A[]);
  /* Writes to {wr} a human-readable description of one segment of the fitted funtion {F}, spanning
     breaks {B[i]} to {B[j]}.  The segment's description is taken from the tableau {A}. */

void fpc_resample_function(fpc_resample_options_t *rso, int NP, double B[], fpc_segment_t *A[], int NS, int r[]);
  /* Writes to a file called {rso.fileName} the values of the fitted function at {rso.NR} equally spaced 
    points {rso.XIni + k*rso.XStep} for {k} in {0..rso.NR-1}. */
    
void fpc_make_predictions
  ( int NP, 
    int ID[],
    double X[], 
    double Z[], 
    double W[], 
    int gMax,
    fpc_penalty_options_t *pto, 
    double B[], 
    fpc_segment_t *A[],
    fpc_predict_options_t *pro 
  );
  /* For each {i} in {0..NP-1}, writes to a file called {pro.fileName} the predicted 
    values values of {F'[i](X[i] + pro.xDelta)}, where {F'[i]()} is the function 
    fitted to the data points {X[t],Z[t],W[t]} with {t} in {0..i} only.
    
    The last segment of each partial function {F'[i]()} is fitted with slightly different rules
    than the others.  In particular, its max degree is {min(gMax,m-1)}, where {m}
    is the number of data points spanned by that last segment, instead of 
    {min(gMax,m-2)} as in the other segments. */



void fpc_eval_function_original(int NP, double X[], double B[], fpc_segment_t *A[], int NS, int r[], double Y[]);
  /* For each {i} in {0..NP-1}, computes the value {Y[i] = F(X[i])} of the fitted function described 
    by {B[0..NP],A[0..NA-1],NS,r[0..NS]}.
    
    The last segment of each partial function {F'[i]()} is fitted with slightly different rules
    than the others.  In particular, its max degree is {min(gMax,m-1)}, where {m}
    is the number of data points spanned by that last segment, instead of 
    {min(gMax,m-2)} as in the other segments. */

void fpc_write_data(FILE *wr, int NP, int ID[], double X[], double Z[], double W[], double Y[]);
  /* Writes to {wr} the original data {ID[0..NP-1],X[0..NP-1],Z[0..NP-1],W[0..NP-1]}, the 
    values {Y[0..NP-1]} assumed to be {F(X[0..NP-1])}, and the residuals {Z[i] - Y[i]}. */

void fpc_residual_stats(int NP, double Z[], double W[], double Y[], double *avgP, double *devP);
  /* Returns in {*avgP} and {*devP} the weighted mean and the weighted standard deviation
    of the residuals {Z[i] - Y[i]}, computed with weights {W[i]}. No correction is
    made for the number of fitted parameters. */

fpc_segment_t *fpc_segment_alloc(int degree);
  /* Allocates a descriptor for a segment of the specified {degree}, including the
    {.coeffs} vector. The penalties {.QFIT,.QOUT} are set to {NAN}. The degree must be
    {-1} (for the everywhere-{NAN} function) or non-negative. */

void fpc_segment_free(fpc_segment_t *FS);
  /* Reclaims the space used by {FS}, including the coefficient vector {FS->coeff}. */

fpc_segment_t **fpc_tableau_alloc(int NP);
  /* Allocates a tableau of segments for {NP} data points, that is, {NP+1} potential breaks.
    The tableau is stored by rows as a vector with {choose(NP+1,2)} entries,
    which are all set to {NULL} initially. */
    
void fpc_tableau_free(int NP, fpc_segment_t ** A);
  /* Reclaims the space used by {A}, including all its segment descriptors and the respective 
    coefficient vectors {FS->coeff}. */

fpc_segment_t **fpc_get_segment_address(int NP, fpc_segment_t *A[], int i, int j);
  /* Returns the address of the element in row {i} and column {j} of the triangular tableau {A}.
    Requires {0 <= i < j <= NP}. */

void fpc_set_segment(int NP, fpc_segment_t *A[], int i, int j, fpc_segment_t *FS);
  /* Stores {FS} into element in row {i} and column {j} of the triangular tableau {A}.
    Requires {0 <= i < j <= NP}. */

fpc_segment_t *fpc_get_segment(int NP, fpc_segment_t *A[], int i, int j);
  /* Returns the pointer stored in row {i} and column {j} of the triangular tableau {A}.
    Requires {0 <= i < j <= NP}. */
  
void fpc_fit_polynomial(int NP, int i, int j, double X[], double Z[], double W[], double B[], int degree, double C[]);
/* Fits a polynomial of the specified {degree} to points {X[u],Z[u]}, weighted by {W[u]}, for {u} in {i..j-1}.
  Expects {0 <= i < j <= NP}, and {B[i] < X[u] < B[j]} for all such {u}.
  
  The coefficients are returned in {C[0..degree]}. The polynomial is
  {FS(XANY) = SUM{ C[t]*(XANY-XC)^t : t \in 0..degree }}, where {XC = (B[i]+B[j])/2}. */

void fpc_eval_polynomial_basis(double x, int degree, int NH, double H[]);
  /* Fills {C[0..NH-1]} with the values of the polynomial basis of the given {degree}
    at the given domain-normalized argument {x}. Currently the basis is {x^t} for {t} in {0..degree},
    so it requires {NH = degree+1}.  */

void fpc_accum_system(double Zu, double Wu, int NH, double Hu[], double *M, double *B);
  /* Updates the least squares system {*M,*B} with one more observation,
    consisting of the given value {Zu = Z[u]} with weight {Wu = W[u]} and the corresponding
    values {Hu[0..NH-1]} of the basis functions at the same argument {X[u]}, for some {u} in {0..NP-1}.
    The matrix {M} must have {NH} rows and columns, and the vector {B} must have {NH}
    elements. */

int main (int argc, char **argv)
  {
    bool_t verbose = TRUE;
    
    fpc_options_t *o = fpc_parse_options(argc, argv);
    int NP = -1; /* Number of input samples. */
    int *ID = NULL; /* Input id numbers {ID[0..NP-1]}. */
    double *X = NULL; /* Input argument values are {X[0..NP-1]}. */
    double *Z = NULL; /* Input data values are {Z[0..NP-1]}. */
    double *W = NULL; /* Weights of data records are {W[0..NP-1]}. */
    fpc_read_data(stdin, o->weighted, &NP, &ID, &X, &Z, &W);
    if (verbose) { fprintf(stderr, "read %d data samples\n", NP); }

    /* Build the model: */
    double *B = rn_alloc(NP+1); /* Potential breaks. */
    fpc_choose_potential_breaks(NP, X, B);
    
    fpc_segment_t **A = fpc_tableau_alloc(NP); /* Potential segments. */
    int NS = -1;   /* Number of segments in fitted formula. */
    int *r = NULL; /* Indices of selected breaks. */
    fpc_fit_function(NP, X, Z, W, o->maxDegree, o->penalty, B, A, &NS, &r);

    if (o->writeFunction != NULL)
      { /* Write formula to disk: */
        fpc_write_function(o->writeFunction, NP, B, A, NS, r);
      }
      
    if (o->resample != NULL)
      { /* Resample the formula: */
        fpc_resample_function(o->resample, NP, B, A, NS, r);
      }
      
    if (o->predict != NULL)
      { /* Make incremental predictions: */
        fpc_make_predictions(NP, ID, X, Z, W, o->maxDegree, o->penalty, B, A, o->predict);
      }
      
    /* Evaluate the function at the original arguments and write to {stdout}: */
    double *Y = notnull(malloc(NP*sizeof(double)), "no mem");
    fpc_eval_function_original(NP, X, B, A, NS, r, Y);
    fpc_write_data(stdout, NP, ID, X, Z, W, Y);

    double avg, dev; /* Statistics of residual. */
    fpc_residual_stats(NP, Z, W, Y, &avg, &dev);
    
    if (verbose) { fpc_print_function(stderr, NP, B, A, NS, r, avg, dev);  }
    
    free(B);
    fpc_tableau_free(NP, A);
    free(ID);
    free(X);
    free(Z);
    free(W);

    return 0;
  }

void fpc_read_data(FILE *rd, bool_t weighted, int *NPp, int **IDp, double **Xp, double **Zp, double **Wp)
  {
    int NP0 = 2048; /* Initial allocation; may be expanded. */
    int_vec_t ID = int_vec_new(NP0);
    double_vec_t X = double_vec_new(NP0);
    double_vec_t Z = double_vec_new(NP0);
    double_vec_t W = double_vec_new(NP0);

    /* Read the input samples: */
    int NP = 0; /* Input values are {Z.e[0..NP-1]}. */
    while (TRUE)
      { fget_skip_spaces(rd);
        int c = fgetc(rd);
        if (c == EOF) 
          { break; }
        else if (c == '#')
          { do { c = fgetc(rd); } while ((c != '\n') && (c != EOF)); 
            if (c == EOF) { break; }
            continue;
          }
        else if (c == '\n')
          { continue; }
        else
          { ungetc(c, rd); 
            int_vec_expand(&ID,NP);
            double_vec_expand(&X,NP);
            double_vec_expand(&Z,NP);
            double_vec_expand(&W,NP);
            ID.e[NP] = fget_int32(rd);
            X.e[NP] = fget_double(rd);
            Z.e[NP] = fget_double(rd);
            if (weighted)
              { W.e[NP] = fget_double(rd);
                demand(W.e[NP] >= 0, "weights must be non-negative");
              }
            else
              { W.e[NP] = 1.0; }
            NP++;
            fget_comment_or_eol(rd, '#');
          }
      }
    int_vec_trim(&ID,NP);
    double_vec_trim(&X,NP);
    double_vec_trim(&Z,NP);
    double_vec_trim(&W,NP);
      
    (*NPp) = NP;
    (*IDp) = ID.e;
    (*Xp) = X.e;
    (*Zp) = Z.e;
    (*Wp) = W.e;
   }
  
double fpc_eval_function(double XArg, int NP, double B[], fpc_segment_t *A[], int NS, int r[], int *kHintP)
  {
    /* Find the index {k1} of a segment whose domain contains {XArg}, and the break indices {i1,j1}: */
    int k1 = fpc_locate_segment(XArg, NP, B, NS, r, kHintP);
    if ((k1 < 0) || (k1 >= NS)) { return NAN; }
    double val1 = fpc_eval_segment_by_index(XArg, k1, NP, B, A, NS, r);
    /* If {Xarg} is at a segment boundary, get the other segment {FS2}: */
    int k2 = -1;
    if ((XArg == B[r[k1]]) && (k1 >= 1))
      { k2 = k1 - 1; }
    else if ((XArg == B[r[k1+1]]) && (k1 <= NP-2))
      { k2 = k1 + 1; }
    if (k2 == -1)
      { return val1; }
    else
      { double val2 = fpc_eval_segment_by_index(XArg, k2, NP, B, A, NS, r);
        /* Return the average of the non-{NAN} values: */
        if (isnan(val1))
          { return val2; }
        else if (isnan(val2))
          { return val1; }
        else
          { return 0.5*val1 + 0.5*val2; }
      }
  } 

int fpc_locate_segment(double XArg, int NP, double B[], int NS, int r[], int *kHintP)
  {
    /* Overall domain: */
    if ((XArg < B[0]) || (XArg > B[NP])) { return -1; }
    if (NS == 0) { return -1; }
    /* Initialize the range {k0..k1-1} of the segment index {k}: */
    int k0 = 0;
    int k1 = NS;
    int kHint = (kHintP == NULL ? -1 : (*kHintP));
    if ((kHint >= 0) && (kHint < NS))
      { /* Consider the hint for {k1}: */
        if ((kHint >= 1) && (XArg <= B[r[kHint]]))
          { k1 = kHint; }
        else if ((kHint < NS) && (XArg <= B[r[kHint+1]]))
          { k1 = kHint+1; }
        /* Consider the hint for {k0}: */
        if ((kHint < NS) && (XArg >= B[r[kHint]]))
          { k0 = kHint; }
        else if ((kHint > 0) && (XArg >= B[r[kHint-1]]))
          { k0 = kHint+1; }
      }
    assert((0 <= k0) && (k0 < k1) && (k1 <= NS));
    assert((XArg >= B[r[k0]]) && (XArg <= B[r[k1]]));
    if (k1 - k0 > 1) 
      { /* Consider the first and last intervals: */
        if (XArg <= B[r[k0+1]]) 
          { k1 = k0+1; }
        else if (XArg >= B[r[k1-1]])
          { k0 = k1-1; }
        else
          { /* Binary search: */
            while (k1 - k0 >= 2)
              { int km = (k0 + k1)/2;
                if (XArg <= B[r[km]]) { k1 = km; } else { k0 = km; }
              }
          }
      }
    assert(k1 - k0 == 1);
    assert((XArg >= B[r[k0]]) && (XArg <= B[r[k1]]));
    if (kHintP != NULL) { (*kHintP) = k0; }
    return k0;
  }

double fpc_eval_segment_by_index(double XArg, int k, int NP, double B[], fpc_segment_t *A[], int NS, int r[])
  { if ((k < 0) || (k >= NS)) { return NAN; }
    /* Get the breaks {i,j} spanned by segment {k}: */
    int i = r[k];
    int j = r[k+1];
    assert(B[i] < B[j]); /* Strictly smaller. */
    /* Get the segment {FS} spanning breaks {i,j}: */
    fpc_segment_t *FS = fpc_get_segment(NP, A, i, j);
    /* Evaluate the {FS} at {XArg}, save value in {val}: */
    double x = fpc_normalize_arg(XArg, B[i], B[j]);
    double val = fpc_eval_segment(x, FS);
    return val;
  }
  
double fpc_normalize_arg(double XArg, double Bi, double Bj)
  { 
    double R = 0.5*Bj - 0.5*Bi;
    double XC = 0.5*Bi + 0.5*Bj;
    assert((Bi < XC) && (XC < Bj));
    assert(R > 0);
    return (XArg - XC)/R;
  }

double fpc_eval_segment(double x, fpc_segment_t *FS)
  { 
    if (FS == NULL){ return NAN; }
    int degree = FS->degree; /* Degree of polynomial. */
    assert(degree >= -1);
    if (degree < 0)
      { return NAN; }
    else
      { int NH = degree + 1; /* Number of coefficients. */
        double H[NH];
        fpc_eval_polynomial_basis(x, degree, NH, H);
        double sum = 0;
        int t;
        for (t = 0; t < NH; t++) { sum += FS->coeff[t]*H[t]; }
        return sum;
      }
  }
    
void fpc_eval_polynomial_basis(double x, int degree, int NH, double H[])
  { 
    assert(NH == degree+1);
    double h = 1.0;
    int t;
    for (t = 0; t <= degree; t++) { H[t] = h; h = h * x; } 
  }
  
void fpc_eval_function_original(int NP, double X[], double B[], fpc_segment_t *A[], int NS, int r[], double Y[])
  { 
    int i;
    int kHint = -1;
    for (i = 0; i < NP; i++) { Y[i] = fpc_eval_function(X[i], NP, B, A, NS, r, &kHint); }
  }

void fpc_choose_potential_breaks(int NP, double X[], double B[])
  {
    assert(NP >= 2);
    B[0] = X[0] - 0.5*(X[1] - X[0]);
    int i;
    for (i = 1; i < NP; i++) { B[i] = 0.5*(X[i-1] + X[i]); }
    B[NP] = X[NP-1] + 0.5*(X[NP-1] - X[NP-2]);
  }

fpc_segment_t *fpc_segment_alloc(int degree)
  { 
    fpc_segment_t *FS = notnull(malloc(sizeof(fpc_segment_t)), "no mem");
    assert(degree >= -1);
    int NH = degree + 1;
    (*FS) = (fpc_segment_t) 
      { .degree = degree,
        .QFIT = NAN,
        .QOUT = NAN,
        .coeff = (NH == 0 ? NULL : rn_alloc(NH))
      };
    return FS;
  }

void fpc_segment_free(fpc_segment_t *FS)
  { if (FS != NULL)
      { free(FS->coeff); free(FS); }
  }

fpc_segment_t **fpc_tableau_alloc(int NP)
  { int NA = NP*(NP+1)/2;  /* Number of potential segments {choose(NP+1,2)}. */
    fpc_segment_t **A = notnull(malloc(NA*sizeof(fpc_segment_t *)), "no mem");
    int ij;
    for (ij = 0; ij < NA; ij++) { A[ij] = NULL; }
    return A;
  }
    
void fpc_tableau_free(int NP, fpc_segment_t ** A)
  { if (A != NULL)
      { int NA = NP*(NP+1)/2;  /* Number of potential segments {choose(NP+1,2)}. */
        int ij;
        for (ij = 0; ij < NA; ij++) { if (A[ij] != NULL) { fpc_segment_free(A[ij]); } }
        free(A);
      }
  }

void fpc_fit_function
  ( int NP, 
    double X[], 
    double Z[], 
    double W[], 
    int gMax, 
    fpc_penalty_options_t *pto, 
    double B[], 
    fpc_segment_t *A[], 
    int *NSP, 
    int **rP
  )
  {
    /* Optimal segment chains from {B[0]} to each {B[i]}: */
    double *QBEST = rn_alloc(NP+1); /* {QBEST[i]} is penalty of the opt chain that ends at {B[i]}, incl. break pen.. */
    int *iPrev = in_alloc(NP+1);    /* {iPrev[i]} is the index of the previous break that chain, or {-1}. */
    int *iQueue = in_alloc(NP+1);   /* Queue of best previous breaks. */ 
    
    /* Start with the empty chain from {B[0]} to {B[0]}: */
    QBEST[0] = 0;   /* No penalty for the break at {B[0]}. */
    iPrev[0] = -1;  /* No previous break. */
    iQueue[0] = 0;  /* The only break so far is {B[0]}. */
    
    /* Now compute the best chains to all other {B[j]}: */
    int j;
    for (j = 1; j <= NP; j++)
      { QBEST[j] = +INF; /* No chains considered yet. */
        iPrev[j] = -1;   /* Just in case. */
        /* Try all values {i} of the previous break point: */
        /* As this point, {iQueue[0..j-1]} are the indices {0..j-1} sorted by decreasing {QBEST}. */
        int s;
        for (s = j - 1; s >= 0; s--)
          { /* Get the next best candidate {i}: */
            int i = iQueue[s];
            /* Consider a segment from {B[i]} to {B[j]}: */
            if (QBEST[i] < QBEST[j])
              { int NPi = j - i; /* Points spanned by segment. */
                int degree = ((i == 0) || (j == NP) ? NPi - 1 : NPi - 2);
                if (degree >  gMax) { degree = gMax; }
                fpc_segment_t *FS = fpc_fit_segment(NP, X, Z, W, degree, B, i, j, pto->ptSkip);
                fpc_set_segment(NP, A, i, j, FS);
                double QTOTi = QBEST[i] + FS->QFIT + FS->QOUT;
                if (QTOTi < QBEST[j]) { QBEST[j] = QTOTi; iPrev[j] = i; }
              }
          }
        /* Add the break penalty at {B[j]}, if any: */
        if (j < NP)
          { double QBREAKj = pto->ptBreak * fpc_break_penalty_weight(NP, W, j);
            QBEST[j] += QBREAKj;
          }
        /* Insert {j} in {iQueue}: */
        s = j;
        while ((s > 0) && (QBEST[iQueue[s-1]] < QBEST[j])) { iQueue[s] = iQueue[s-1]; s--; }
        iQueue[s] = j;
      }
    /* Now gather the optimal chain from {B[0]} to {B[NP]}: */
    int NS = 0;
    int i = NP;
    while (i > 0) { NS++; i = iPrev[i]; }
    assert((NP == 0) || (NS > 0));
    int *r = in_alloc(NS+1);
    r[NS] = NP;
    int k = NS;
    while (k > 0) { k--; r[k] = iPrev[r[k+1]]; }
    /* Cleanup: */
    free(QBEST);
    free(iPrev);
    free(iQueue);
    /* Return it: */
    (*NSP) = NS;
    (*rP) = r;
  }

void fpc_residual_stats(int NP, double Z[], double W[], double Y[], double *avgP, double *devP)
  {
    int i;
    double sumWD = 0;
    double sumW = 0;
    for (i = 0; i < NP; i++) { sumWD += W[i]*(Z[i] - Y[i]); sumW += W[i]; }
    double avg = sumWD/sumW;
    demand(sumW > 0, "total weight is zero");
    double sumWD2 = 0;
    for (i = 0; i < NP; i++) { double Di = (Z[i] - Y[i]) - avg;  sumWD2 += W[i]*Di*Di; }
    double dev = sqrt(sumWD2/sumW);
    (*avgP) = avg;
    (*devP) = dev;
  }

void fpc_resample_function(fpc_resample_options_t *rso, int NP, double B[], fpc_segment_t *A[], int NS, int r[])
  {
    FILE *wr = open_write(rso->fileName, TRUE);
    int kHint = -1; /* Hint for segment search. */
    int NR = rso->NR;
    int j;
    for (j = 0; j < NR; j++)
      { double Xj = rso->XIni + j*rso->XStep;
        double Zj = fpc_eval_function(Xj, NP, B, A, NS, r, &kHint);
        fprintf(wr, "%6d %g %g\n", j, Xj, Zj);
      }
    fclose(wr);
  }

void fpc_write_data(FILE *wr, int NP, int ID[], double X[], double Z[], double W[], double Y[])
  {
    int i;
    for (i = 0; i < NP; i++)
      { fprintf(wr, "%d %g %g %g  %g %g\n", ID[i], X[i], Z[i], W[i], Y[i], Z[i] - Y[i]); }
    fflush(wr);
  }

void fpc_write_function(char *fName, int NP, double B[], fpc_segment_t *A[], int NS, int r[])
  {
    FILE *wr = open_write(fName, TRUE);
    fprintf(wr, "%d\n", NS);
    int k;
    for (k = 0; k < NS; k++) { fpc_write_segment(wr, NP, r[k], r[k+1], B, A); }
    fclose(wr);
  }
  
void fpc_write_segment(FILE *wr, int NP, int i, int j, double B[], fpc_segment_t *A[])
  {
    assert((0 <= i) && (i < j) && (j <= NP));
    /* Get the domain interval's center {XC}: */
    double XC = (B[i] + B[j])/2;
    assert((B[i] < XC) && (XC < B[j]));
    fprintf(wr, "%5d %20g %5d %20g %20g", i, B[i], j, B[j], XC); 
    /* Get segment data: */
    fpc_segment_t *FS = fpc_get_segment(NP, A, i, j);
    assert(FS != NULL);
    fprintf(wr, "  %20g %20g %3d", FS->QFIT, FS->QOUT, FS->degree);
    int u;
    for (u = 0; u <= FS->degree; u++) { fprintf(wr, " %20g", FS->coeff[u]); }
    fprintf(wr, "\n");
  }
      
void fpc_print_function(FILE *wr, int NP, double B[], fpc_segment_t *A[], int NS, int r[], double avg, double dev)
  {
    fprintf(wr, "fitted function:\n");
    fprintf(wr, "number of segments = %d\n", NS);
    int k;
    for (k = 0; k < NS; k++) { fpc_print_segment(wr, NP, k, r[k], r[k+1], B, A); }
    fprintf(wr, "weighted mean residual = %+14.9f\n", avg);
    fprintf(wr, "weighted root mean squared residaual = %+14.9f\n", dev);
    fflush(wr);
  }
      
void fpc_print_segment(FILE *wr, int NP, int k, int i, int j, double B[], fpc_segment_t *A[])
  {
    assert((0 <= i) && (i < j) && (j <= NP));
    /* Get the domain interval's center {XC}: */
    double XC = (B[i] + B[j])/2;
    assert((B[i] < XC) && (XC < B[j]));
    fprintf(wr, "segment %d spans from break B[%d] = %g to break B[%d] = %g\n", k, i, B[i], j, B[j]);
    /* Get segment data: */
    fpc_segment_t *FS = fpc_get_segment(NP, A, i, j);
    assert(FS != NULL);
    fprintf(wr, "  center XC = %g", XC);
    fprintf(wr, "  fit penalty QFIT = %g  outlier penalty = %g  degree = %3d\n", FS->QFIT, FS->QOUT, FS->degree);
    int u;
    for (u = 0; u <= FS->degree; u++) 
      { fprintf(wr, "  %+14.9f", FS->coeff[u]);
        if (u > 0) { fprintf(wr, "*(X-XC)"); }
        if (u > 1) { fprintf(wr, "^%d", u); }
        fprintf(wr, "\n");
      }
  }
 
fpc_segment_t **fpc_get_segment_address(int NP, fpc_segment_t *A[], int i, int j)
  {
    assert((0 <= i) && (i < j) && (j <= NP));
    int k = i*(2*NP - 1 - i)/2 + j;
    return &(A[k]);
  }

void fpc_set_segment(int NP, fpc_segment_t *A[], int i, int j, fpc_segment_t *FS)
  {
    fpc_segment_t **Aij = fpc_get_segment_address(NP, A, i, j);
    (*Aij) = FS;
  }

fpc_segment_t *fpc_get_segment(int NP, fpc_segment_t *A[], int i, int j)
  { 
    fpc_segment_t **Aij = fpc_get_segment_address(NP, A, i, j);
    return (*Aij);
  }

double fpc_break_penalty_weight(int NP, double W[], int j)
  {
    if ((j == 0) || (j == NP)) { return 0.0; }
    double W0 = W[j-1];
    double W1 = W[j];
    if ((W0 == 0) || (W1 == 0)) { return 0.0; }
    return 2.0*W0*W1/(W0 + W1);
  }

void fpc_fit_polynomial(int NP, int i, int j, double X[], double Z[], double W[], double B[], int degree, double C[])
  {
    bool_t verbose = TRUE;
    
    int NH = degree+1;
    if (NH == 0) { return; }
    
    /* Get the domain interval's center {XC}: */
    assert((0 <= i) && (i < j) && (j <= NP));
    assert(j - i >= NH);

    /* The least squares system: */
    double *M = rmxn_alloc(NH, NH);  /* Matrix of linear system. */
    double *R = rn_alloc(NH);        /* Right-hand side of system. */
    rn_zero(NH, R);
    rmxn_zero(NH, NH, M);
    
    /* Accumulate the matrix and vector of the least squares system: */
    double *Hu = rn_alloc(NH);        /* Basis values. */
    int u;
    for (u = i; i < j; u++)
      { double xu = fpc_normalize_arg(X[i], B[i], B[j]);
        fpc_eval_polynomial_basis(xu, degree, NH, Hu);
        fpc_accum_system(Z[u], W[u], NH, Hu, M, R);
      }
      
    if (verbose)
      { /* Print the Least Squares system: */
        int r, s;
        for (r = 0; r < NH; r++)
          { for (s = 0; s < NH; s++) 
              { fprintf(stderr, " %12.9f", M[r*NH + s]); }
            fprintf(stderr, "  ");
            fprintf(stderr, " %12.9f", R[r]);
            fprintf(stderr, "\n");
          }
      }
    
    /* Solve the least squares system: */
    double *L = rmxn_alloc(NH,NH);
    (void)rmxn_inv(NH, M, L);
    rmxn_map_col(NH, NH, L, R, C);
    
    free(L);
    free(Hu);
    free(M);
    free(R);
  }

void fpc_accum_system(double Zu, double Wu, int NH, double Hu[], double *M, double *R)
  { 
    int j, k;
    for (k = 0; k < NH; k++)
      { R[k] += Wu*Zu*Hu[k];
        for (j = 0; j < NH; j++) { M[k*NH + j] += Wu*Hu[k]*Hu[j]; }
      }
  }
  
#define fpc_maxDegree_MAX 6
  /* Maximum "-maxDegree" argument allowed on command line. */

#define fpc_penalty_break_MAX 1.0e+100
#define fpc_penalty_skip_MAX 1.0e+100
  /* Maximum penalties allowed on command line. */
  
#define fpc_resample_NR_MAX ((1<<30)-1)
  /* Maximum number of resampling points. */

#define fpc_X_TINY 1.0e-100
#define fpc_X_HUGE 1.0e+100
  /* Min and max arguments allowed on command line. */

fpc_options_t *fpc_parse_options(int argc, char **argv)
  {
    /* Initialize argument parser: */
    argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, PROG_NAME " version " PROG_VERS ", usage:\n" PROG_HELP);
    argparser_set_info(pp, PROG_INFO);
    argparser_process_help_info_options(pp);

    /* Allocate the command line argument record: */
    fpc_options_t *o = notnull(malloc(sizeof(fpc_options_t)), "no mem");

    /* Parse keyword parameters: */

    argparser_get_keyword(pp, "-maxDegree");
    o->maxDegree = (int)argparser_get_next_int(pp, 1, fpc_maxDegree_MAX);
    
    argparser_get_keyword(pp, "-penalty");
    o->penalty = fpc_parse_penalty_options(pp);
    
    o->weighted = argparser_keyword_present(pp, "-weighted");
    
    if (argparser_keyword_present(pp, "-writeFunction"))
      { o->writeFunction = argparser_get_next_non_keyword(pp); }
    else
      { o->writeFunction = NULL; }

    if (argparser_keyword_present(pp, "-resample"))
      { o->resample = fpc_parse_resample_options(pp); }
    else
      { o->resample = NULL; }

    if (argparser_keyword_present(pp, "-predict"))
      { o->predict = fpc_parse_predict_options(pp); }
    else
      { o->predict = NULL; }

    /* Check for spurious arguments: */
    argparser_finish(pp);

    return o;
  }

fpc_penalty_options_t *fpc_parse_penalty_options(argparser_t *pp)
  {
    /* Allocate the command line argument record: */
    fpc_penalty_options_t *o = notnull(malloc(sizeof(fpc_penalty_options_t)), "no mem");
    o->ptBreak = NAN;
    o->ptSkip = NAN;

    /* Parse keyword parameters after "-penalty": */
    while (TRUE)
      { if (argparser_keyword_present_next(pp, "break"))
          { o->ptBreak = argparser_get_next_double(pp, 0.0, fpc_penalty_break_MAX); }
        else if (argparser_keyword_present_next(pp, "skip"))
          { o->ptSkip = argparser_get_next_double(pp, 0.0, fpc_penalty_skip_MAX); }
        else
          { break; }
      }
    if (isnan(o->ptBreak)) { argparser_error(pp, "missing \"break\" penalty parameter"); }
    if (isnan(o->ptSkip)) { argparser_error(pp, "missing \"skip\" penalty parameter"); }
    return o;
  }

fpc_resample_options_t *fpc_parse_resample_options(argparser_t *pp)
  {
    /* Allocate the command line argument record: */
    fpc_resample_options_t *o = notnull(malloc(sizeof(fpc_resample_options_t)), "no mem");
    o->NR = (int)argparser_get_next_int(pp, 1, fpc_resample_NR_MAX);
    o->XIni = argparser_get_next_double(pp, -fpc_X_HUGE, +fpc_X_HUGE);
    o->XStep = argparser_get_next_double(pp, 1.0e-13*fabs(o->XIni), fpc_X_HUGE);
    o->fileName = argparser_get_next_non_keyword(pp);
    return o;
  }

fpc_predict_options_t *fpc_parse_predict_options(argparser_t *pp)
  {
    /* Allocate the command line argument record: */
    fpc_predict_options_t *o = notnull(malloc(sizeof(fpc_predict_options_t)), "no mem");
    o->xDelta = argparser_get_next_double(pp, fpc_X_TINY, fpc_X_HUGE);
    o->fileName = argparser_get_next_non_keyword(pp);
    return o;
  }

