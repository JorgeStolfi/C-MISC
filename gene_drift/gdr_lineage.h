/* Tools for surviving lineages counts and and ratios. */
#ifndef gdr_lineage_H
#define gdr_lineage_H

/* Last edited on 2023-06-01 13:43:55 by stolfi */

#define gdr_lineage_H_COPYRIGHT \
  "Duh?"

#define _GNU_SOURCE
#include <stdint.h>

#include <bool.h>
#include <vec.h> 

#include <gdr_demo.h>
#include <gdr_sim.h>

#define gdr_lineage_cohort_size_MAX 90000
#define gdr_lineage_num_years_MAX 10000
#define gdr_lineage_num_runs_MAX 100000

/* 
  The following concepts appply to a simulation state {st} (a
  {cdr_sim_state_t} record):
  
  LINEAGE
  
  The /lineage/ {L[i]} of an individual {i} is the set of all
  individuals of the same sex {s} as {i}, including {i} itself, which
  descend from {i} entirely through children of sex {s} (i.e. strictly
  patrilinear or matrilinear descent).  The individual {i} is the /founder/ of 
  the lineage.
  
  Note that two distinct individuals born on the same year will have
  disjoing lineages.
  
  SURVIVING LINEAGE
  
  A lineage is said to be /surviving/ if it includes at least one
  individual which is actually alive on the last year {st.ny-1} of the
  simulation. 
  
  This concept is defined only after the simulation has progressed to
  that final year (that is, when {st.yCur == st.ny-1}), so that every
  individual which could be alive on year {st.ny-1} has been
  generated and expanded.
  
  SURVIVING LINEAGE COUNT
  
  A surviving lineage is said to /count for/ some year {y} if its founder {i} 
  was born on or before {y}, and (1) it is still alive at year {ny-1},
  or (2) is not alive at year {ny-1} but has at least one child born 
  ustrictly after {y} whose lineage is surviving.
  
  The /surviving lineage count/ {nsl[y]} for a year {y} is the number of
  surviving lineages that are countable on {y}. */
     
void gdr_lineage_compute_surviving_counts(gdr_sim_state_t *st, int32_t nsl[]);
  /* The parameter {nsl} should be an array of {st.ny} integers. At
    the end of the simulation, for each year {y} in {0..st.ny-1}, stores
    into {nsl[y]} the surviving lineage count for year {y}. */    

void gdr_lineage_compute_sex_ratio_table
  ( int32_t ny, 
    int32_t nr, 
    int32_t nsl_sv0[],
    int32_t nsl_sv1[],
    double rsl_sv[]
  );
  /* The parameters {nsl_sv0}, {nsl_sv1}, and {rsl_sv} must be arrays of
    {ny} rows and {nr} columns, linearized by rows.
    
    The tables {nsl_sv0} and {nsl_sv1} must have the number of surviving
    lineages per year and run for sexes 0 and 1, respectively, as
    computed by {gdr_lineage_compute_counts_table}.
    
    The procedure sets {rsl_sv[y,r]} to the female-to-male surviving
    lineage ratios {nsl_sv1[y,r]/nsl_sv0[y,r]}. */
    
void gdr_lineage_find_notable_ratio_parameters
  ( int32_t ny, 
    int32_t nr, 
    double rsl_sv[],
    int32_t *gcPlat_P,
    double *rslPlat_P
  );
  /* Assumes that {rsl_sv[0..ny-1,0..nr-1]} is the table of surviving
    lineage ratios per year and run, as created by
    {gdr_lineage_compute_sex_ratio_table}.
    
    For each year {y} in {0..ny-1}, extracts from row {y} of that
    table the median {rslMed[y]}.
    
    Then finds the approximate value {rslPlat} of the median when it is
    close to reaching its max abs value and the year {gcPlat}
    (counting back from 0 = most recent) where that value is reached. */

/* OUTPUT */

void gdr_lineage_write_counts_median_range
  ( char *outPrefix,
    char *tag,
    int32_t s,
    int32_t ny, 
    int32_t yStart,
    int32_t nr,
    int32_t nsl_sv[],
    int32_t pop_sv[]
  );
  /* Assumes that {pop[y]} is the number of individuals actually
    alive on each year {y} in {0..ny-1}. Assumes that{nsl_sv} is an
    arrray of {nr} columns and {ny} rows, linearized by rows; where
    {nsl_sv[y,r]} is the number of surviving lineages in year {y} of run
    {r}.
    
    Writes to a file "{outPrefix}-{tag}-{s}-nlins.txt" the statistical
    summary (median and range) of {nsl_sv} for each year {y}, as described in
    {gdr_lineage_counts_median_range_file_INFO}. Also writes the data to
    {stderr}. */

#define gdr_lineage_counts_median_range_file_INFO \
  "Has one line for each simulated year {y} in {yStart..yStop}, with a statistical" \
  " summary of the number of surviving lineages {nsl[y]} over" \
  " all repetitions of the experiment for sex {s}.  Namely, each line" \
  " has the following fields:\n" \
  "\n" \
  "        {y} {plo95(y)} {pmed(y)} {phi95(y)} {flo95(y)} {fmed(y)} {fhi95(y)}\n" \
  "\n" \
  "      Here {pmed(y)} is the median of the population sizes" \
  " (counts of individuals alive) on" \
  " year {y}, {pop_sv[y,0..nr-1]}, over all runs; and {plo95}" \
  " and {phi95} are the low and high ends of the 95% range (that is, the" \
  " 2.5% and 97.5% percentiles) of those populations.\n" \
  "\n" \
  "      An individual is considered /alive/ on every" \
  " year from its birth year up to but not" \
  " including its death year.\n" \
  "\n" \
  "      Similarly, the" \
  " quantities {fmed(y),flo95(y),fhi95(y)}} are the median and 95% range" \
  " of the percentages {pct_sv[y,r]=100*nsl_sv[y,r]/pop_sv[y,r]} over all" \
  " runs {r} in {0..nr-1}."
  
void gdr_lineage_write_sex_ratio_samples
  ( char *outPrefix,
    char *tag,
    int32_t ny, 
    int32_t yStart, 
    int32_t nr, 
    int32_t nrSample,
    double rsl_sv[]
  );
  /* The parameter {rsl_sv[r,y]} should have the of female-to-male
    surviving lineage ratios per year {y} in {0..ny-1} and run {r}
    in {0..nr-1}, as computed by {gdr_lineage_compute_sex_ratio_table}.
    The procedure writes the ratios for {nrSample} different runs
    to file "{outPrefix}-{tag}-rsamp.txt",
    as described in {gdr_lineage_sex_ratio_samples_file_INFO}. */

#define gdr_lineage_sex_ratio_samples_file_INFO \
  "Contains several plots of {rsl(y)}, the ratio of the female surviving" \
  " lineage count to the male one, corresponding to several separate runs" \
  " of the experiment.  The plots are separated by blank lines. Each plot" \
  " has one line for each year {y} in the range {yStart..yStop}, with" \
  " the format \"{y} {kr} {rslp(y)}\" where {kr} is the plot" \
  " index counting from 0, and {rslp(y)} the ratio {rsl(y)}, slightly" \
  " perturbed to avoid overlapping plots." 

void gdr_lineage_write_sex_ratio_median_range
  ( char *outPrefix,
    char *tag,
    int32_t ny, 
    int32_t yStart, 
    int32_t nr, 
    double rsl_sv[]
  );
  /* The parameter {rsl_sv[r,y]} should have the of female-to-male
    surviving lineage ratios per year {y} in {0..ny-1} and run {r}
    in {0..nr-1}, as computed by {gdr_lineage_compute_sex_ratio_table}.
    The procedure writes the median and range of those ratios for each
    year {y} to the file "{outPrefix}-{tag}-ratio.txt", as described in
    {gdr_lineage_sex_ratio_median_range_file_INFO}. */


#define gdr_lineage_sex_ratio_median_range_file_INFO \
  "Has one line for each simulated year {y} in {yStart..yStop}, with a statistical" \
  " summary of the ratio {rsl(y)} of the female surviving lineage count to" \
  " the male one, over" \
  " all repetitions of the experiment.  Namely, each line" \
  " has the following fields:\n" \
  "\n" \
  "        {y} {rsl_lo95[y]} {rsl_med[y]} {rsl_hi95[y]}\n" \
  "\n" \
  "      where {rsl_med(y(} is the median of {rsl(y)} over" \
  " all {nr} runs, and {rsl_lo95(y)} and {rsl_hi95(y)} are" \
  " the low and high ends of the 95% range" \
  " (that is, the 2.5% and 97.5% percentiles of those values)." 

#endif

