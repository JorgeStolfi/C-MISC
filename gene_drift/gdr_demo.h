/* Probability distribution tools for asexual lineage studies. */
#ifndef gdr_demo_H
#define gdr_demo_H
/* Last edited on 2023-06-01 07:28:10 by stolfi */

#define gdr_demo_H_COPYRIGHT \
  "Duh?"

#define _GNU_SOURCE
#include <stdint.h>

#include <bool.h>
#include <vec.h> 
#include <argparser.h> 

#define gdr_demo_child_count_MAX 200
  /* Max child count of an individual. */

#define gdr_demo_age_MAX 200
  /* Max age of an individual. */
 
typedef struct gdr_demo_parms_t {
    int32_t cMax;      /* Max children count. */
    double cAlpha;     /* Decay ratio for child count distribution. */
    int32_t cPrec;     /* Decimal fraction digits to round probs to. */
    int32_t fMin;      /* Minimum age to have children. */
    int32_t fMax;      /* Maximum age to have children. */
  } gdr_demo_parms_t;
  /* A set of /demographic parameters/ that define the key probabilities 
    for simulation of a specific sex {s}.  See {gdr_demo_parms_throw_children}.

    All parameters consider only children with the same sex as the parent.  */   

gdr_demo_parms_t *gdr_demo_parms_new(void);
  /* Allocates a new {gdr_demo_parms_t} record.  Integer fields are set to zero.
    Distribution tables are set to {NULL}. */
 
void gdr_demo_parms_free(gdr_demo_parms_t *dmp);
  /* Frees the storage used by {*dmp}, including internal tables. */
  
double_vec_t gdr_demo_compute_child_count_distr(gdr_demo_parms_t *dmp);
  /* Given the demographic parameters {dmp} for some sex {s},
    returns a vector {cProb} such that {cProb.e[c]}, for {c} in
    {0..dmp->cMax}, is the probability that an individual of sex {s} will
    have {c} children of the same sex.
    
    Each probability will be rounded to the nearest integer multiple of
    multiples of {unit = 10^{-dmp.cPrec}}, but set to {unit} if it would
    round to zero. */

void gdr_demo_show_parms(char *title, int32_t s, gdr_demo_parms_t *dmp);
  /* Writes to {stderr} the given demographic parameters {dmp},
    with the given title.  If {s} is not negative, also shows {s} on the title. */

void gdr_demo_throw_children
  ( double_vec_t *cProb,
    int32_t fMin,
    int32_t fMax,
    int32_t *cNum_P, 
    int32_t cAge[]
  );
  /* Generates a bunch of children for some individual {i} of some sex
    {s} according to the demographic parameters {cProb,fMin,fMax}.
    Considers only children with the same sex {s}. 
    
    The chosen number {cNum} of children will be in {0..cMax} where
    {cMax = cProb.ne-1}. It is returned in {*cNum_P}. The ages of {i} at
    which the children will be born, all in the range {fMin .. fMax},
    are returned in {cAge[0..cNum-1]}. The vector {cAge} must have at
    least {cMax} elements.
    
    The procedure first chooses {cNum} according to the distribution
    {cProb}. Then it picks {cNum} ages in the range {fMin .. fMax}, with
    uniform distribution -- all distinct, if possible. */

void gdr_demo_show_distr(char *title, int32_t s, double_vec_t *iProb, double_vec_t *iFreq);
  /* Writes to {stderr} the distribution {iProb}, with the given title. 
    If {iFreq} is not {NULL}, assumes it is the observed frequencies
    while {iProb} are the theoretical probabilities, and shows 
    both side by side.  If {s} is not negative, also shows {s} on the title. */

void gdr_demo_exponential_distr(int32_t iMin, int32_t iMax, double alpha, double_vec_t *iProb);
  /* Sets {iProb.e[i]} to zero for {i} in {0..iMin-1}, to {alpha^{i-iMin}} for {i}
    in {iMin..iMax}. Truncates {iProb} to {iMax+1} elements.  */

void gdr_demo_distr_normalize_sum(double val, double_vec_t *iProb);
  /* Scales the distribution {iProb} so that its sum is {val}. */

void gdr_demo_distr_round_off(double unit, double_vec_t *iProb);
  /* Rounds off all entries of {iProb} to integer multiples of {unit}.
    Trims off any trailing elements that are zero after rounding. */

void gdr_demo_distr_check(double_vec_t *cProb, int32_t cMax, double avg, int32_t prec);
  /* Checks whether the distribution {cProb} has exactly {cMax+1} entries, all in {[0_1]},
    and its sum and mean are (approximately) 1.0 and {avg}, assuming that they were rounded
    to {prec} decimal fraction digits.  The check of the mean is skipped if {avg} is {NAN}. */
   
void gdr_demo_get_interesting_features(double_vec_t *cProb, int32_t *cBig_P, int32_t *cBigNum_P);
  /* Picks a largish number of children {cBig} less than {cMax = cProb.ne-1}
    and computes the probability of an individual having {cBig} or more children of the
    same sex, namely the sum of {cProb.e[c]} for {c} in {cBig..cMax}. This
    probability is then expressed as an integer number of individuals per 10000,
    rounded down. Returns the results in {*cBig_P} and {*cBigNumP}, respectively. */
    
void gdr_demo_increment_count(int64_vec_t *rCount, int32_t *rMax_P, int32_t r, int32_t amt);
  /* Assumes that {rCount.e[r']} is a set of counts of events identified
    by integers {r'} in {0..rMax}, where {rMax} is the value of
    {*rMax_P} on entry. Adds {amt} to the count {rCount.e[r]}. updating
    {*rMax_P} and expanding {rCount} as needed. */

double_vec_t gdr_demo_freqs_from_counts(int32_t nr, int64_t rCount[]);
  /*  Assumes that {rCount[r]} is a count of events identified by 
    integer {r}, for {r} in {0..nr-1}.  Returns a vector 
    {rFreq} with the corresponding fractional frequencies.
    Namely {rFreq.e[r]} is {rCount[r]} divided by the total count.
    
    Trailing zero elements of {rCount} are ignored and omitted from 
    the frequency vector {rFreq}. */

double_vec_t gdr_demo_reverse_cumul_distr (double_vec_t *iProb);
  /* Assumes {iProb} contains the probabilities of a discrete
    distribution. Returns a vector {iCumu} with same size as {iProb}
    with the reverse cumulative distribution. Namely, 
    {iCumu.e[c] = SUM{iProb.e[k] : k >= c}}. */
   
double gdr_demo_find_quantile(int32_t n, double x[], double frac);
  /* Assumes that the list {x[0..n-1]} is sorted in increasing order,
    interpolated linearly (affinely) as a function {L(z)} for {z} in {[0
    _ n]}. Returns the value of that function for argument {z=frac*n}
    where {frac} is a fraction in {[0_1]}.
    
    Specifically, if {z <= 0.5} returns {x[0]}, else if {z >= n-0.5}
    returns {x[n-1]}, else returns {interp(x[k], x[k+1], z-k-0.5)}
    {k=floor(z-0.5)} and {interp(a,b,f)} is the affine interpolation at
    abscissa {f} between points {(0,a)} and {(1,b)}. */

void gdr_demo_sort_doubles(int32_t n, double x[]);
  /* Sorts the {x[0..n-1]} in increasing order. */

void gdr_demo_write_child_distr_tex_table
  ( char *outPrefix, 
    char *tag,
    int32_t s,
    double_vec_t *cProb,
    int32_t cProbPrec,
    double_vec_t *cFreq,
    int32_t cFreqPrec
  );
  /* Writes to "{outPrefix}-{tag}-{s}-probs.tex" a TeX file with the
    child count distribution {cProb} for sex {s}, and the corresponding actual
    statistics {cFreq} seen in the simulations, as described in
    {gdr_demo_child_distr_tex_table_file_INFO}.
    
    The values in {cProb} and {cFreq} are written with {cProbPrec} and
    {cFreqPrec} decimal fraction digits, respectively. */
 
#define gdr_demo_child_distr_tex_table_file_INFO \
  "A file suitable for input in TeX papers, with a table of child bearing" \
  " probabilities per age, for the sex {s}. The file contains a complete macro" \
  " call \"\\begin{tabular}...\end{tabular}\". The table body" \
  " has one line for each age {a} (in years), with the format\n" \
  "\n" \
  "        {a} & {CCD[s,c]} & {CCF[s,c]} & {CCD_CUM[s,c]} & {CCF_CUM[s,c]} \\\\ \n" \
  "\n" \
  "      where {CCD[s,c]} is the nominal probability that an individual of sex {s}" \
  " will have {c} children of the same sex; {CCF[s,c]} is the actual frequency" \
  " of each child count {c} observed in simulations; and {CCD_CUM} and {CCF_CUM} are" \
  " the complemented cumulative distributions, namely the  the probability" \
  " and frequency of {c} or more children.  The frequencies are" \
  " expected to differ from the nominal probabilities due to" \
  " sampling error and population size adjustment, and consider" \
  " only the individuals who were born and died in the simulated year range."

gdr_demo_parms_t *gdr_demo_options_parse(argparser_t *pp);
  /* Parses the command line arguments that specify the 
    demographic parameters, as specified in 
    {gdr_demo_options_HELP} and {gdr_demo_options_INFO}. */
    
#define gdr_demo_options_HELP \
  " {cMax} {cAlpha} {cPrec} {fMin} {fMax}"

#define gdr_demo_options_INFO \
  "      " gdr_demo_options_HELP "\n" \
  "\n" \
  "    where {cMax} is the max number of children of" \
  " the same sex, {cAlpha} is the decay ratio for the" \
  " child count distribution, {cPrec} is the rounding" \
  " precision for the same, and {fMin} and {fMax} is" \
  " the min and max age for having children."
 
#endif

