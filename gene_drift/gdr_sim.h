/* Basic demography parameters and tools for asexual lineage studies. */
#ifndef gdr_sim_H
#define gdr_sim_H
/* Last edited on 2023-06-01 10:33:17 by stolfi */

#define gdr_sim_H_COPYRIGHT \
  "Duh?"

#define _GNU_SOURCE
#include <stdint.h>

#include <bool.h>
#include <vec.h> 
#include <argparser.h> 

#include <gdr_demo.h>
  
/* 
  CONCEPTS AND NOTATION
  
    YEARS

    Each year covered by the simulation is represented as an integer
    {y}, possibly negative, in in chronological order.

    INDIVIDUALS

    Each individual in a simulation run is identified with an integer,
    starting from 0. All individuals are assumed to be of the same sex {s}.
    
    The key events on the simulated life of an individual are its birth,
    its death, and the generation of zero or more child individuals.
    Each of those events happens on a determinate year.
    
    The goal of the simulation is to generate a description of the
    population and of those events for each year {y} in some range {0..ny-1}. 
  */

/* SIMULATION STATE */

typedef struct gdr_sim_state_t {
    /* Fixed simulation parameters: */
    int32_t fMin;         /* Min child bearing age. */
    int32_t fMax;         /* Max child bearing age. */
    double_vec_t *cProb;  /* Target child count distribution. */
    int32_t cPrec;        /* Decimal fraction digits of {cProb}. */
    int32_t ny;           /* Years of interest are {0..ny-1}. */
    int32_t *nbt;         /* Target number of births per year {0..ny-1}. */
    /* Current state counters: */
    int32_t yCur;         /* Year for which the state is up to date. */
    int32_t ni;           /* Number of individuals created so far. */
    /* Individual properties, indexed by individual {0..ni-1}: */
    int32_vec_t ybr;      /* Year of birth. */
    int32_vec_t ydt;      /* Year of death, or {-1}. */
    int32_vec_t nch;      /* Number of children, or {-1}. */
    int32_vec_t lch;      /* Index of last child, or {-1}. */
    int32_vec_t par;      /* Parent individual, or {-1}. */
    int32_vec_t psy;      /* Another individual with same birth year, or {-1}. */
    /* Cohort tables (``yearbook''), indexed by year {0..ni-1}: */
    int32_t *lib;         /* Last individual born on each year. */
    int32_t *nbr;         /* Number of individuals born on each year. */
  } gdr_sim_state_t;
  /* A record that describes the state during a simulation run.

    The fields {fMin,fMax,cProb,cPrec,ny,nbt} are general demographic
    parameters that govern the simulation. They are defined when the
    state is created, and are not changed afterwards.
    
    At the beginning of each major iteration of the simulation loop, the
    following /state invariants/ hold:
    
    CURRENT YEAR
    
    The /current year/ {yCur} is in the range {0..ny-1}. The population
    history for years {0 .. yCur} has been computed and will not change.
    
    NUMBER OF INDIVIDUALS
    
    The field {ni} is the number of individuals created so far in 
    the simulation. They are numbered {0 .. ni-1}. They include all those who 
    will ever be years {0 .. yCur}, but possibly others.
    
    YEARS OF BIRTH AND DEATH

    Each individual {i} in {0 .. ni-1} is assigned a /year of birth/ {ybr[i]}
    and a /year of death/ {ydt[i]}, with {ybr[i] <= ydt[i]}.
    Both may be negative or greater than {yCur}, even greater
    than {ny-1}.

    Individuals with {ybr[i] <= yCur} are said to have been /expanded/,
    while those with {ybr[i] > yCur} are only /scheduled/ for year
    {ybr[i]}. Every individual
    begins its existence as scheduled, and may later be expanded
    when {yCur} reaches ist birth year.
    
    Once individual {i} is expanded, the year of death {ydt[i]} is
    defined, and all children it is supposed to ever have will have been
    created and added to the state. These children will normally be only
    scheduled, since their years of birth will be greater than {yCur}.
    
    (However, during initialization or some children may be born on or
    before year {yCur}, in which case they will be recursively expanded
    too.)
    
    If individual {i} is still scheduled, none of its children will have
    been created yet, not even just scheduled, and its year of death
    is still undefined; so {ydt[i]} is {-1} by convention. 
    
    NUMBER OF CHILDREN AND LAST CHILD
    
    When an individual {i} is expanded, it is assigned a /children
    count/ {nhc[i]}, which is the number of children created (0 or
    more), whether they are expanded or just scheduled. 
    
    The number of children of a scheduled individual {i} is undefined.
    In that case {nch[i]} is temporarily set to {-1}.  This condition
    marks {i} as still unexpanded.
    
    Expansion also also assigns a /last child index/ {lch[i]}, the
    largest of the indices of its children. If {i} is still only
    scheduled ({nch[i]=-1}), or got no children in the expansion
    ({nch[i]=0}), then {lch[i]} is set to {-1}.
    
    The parameters {nch[i]} and {lch[i]} may increase later by twinning
    during cohort size adjustment, see below.

    ACTUAL LIFE

    An expanded individual is considered to be /alive/ from the year or birth
    up to but /not/ including the year of death.
    
    That is, individual {i} is alive on some year {y} if {ybr[i] <= y}
    and {ydt[i] > y}.
    
    The actual life is undefined before the individual has been expanded.
    
    PARENT POINTERS

    For every individual {i}, expanded or scheduled, the index of the
    parent is saved in {par[i]}, if known.  This will be the case for all
    individuals with {ybr[i]} positive or zero. 
    
    If {ybr[i]} is negative, the parent may not have been created 
    in the simulation; in which case {par[i]} is conventionally set to {-1}.

    COHORTS AND YEARBOOK
    
    The /cohort/ of a year {y} is the set of all individuals created
    so far (expanded or scheduuled) with birth year {y}.
    
    The cohort of each year {y} in {0..ny-1} is organized in the 
    state as a linked list. The element {st.lib[y]} is the index 
    of the last individual added to the cohort, or {-1} if no 
    individual was born on {y}.  

    For any individual {i} on that cohort, {psy.e[i]} the individual
    preceding {i} on the list, or {-1} if {i} is was the first to be
    created with birth year {y}.  
     
    The number of individuals in the cohort (the length of the list) is
    {st.nbr[y]}, the /birth count/ for year {y}.
    
    Note that {st.nbr[y]} and {st.lib[y]} do not extist if the year {y}
    is negative or greater than {st.ny-1}. Individuals born in such
    years are not tracked by those tables.
    
    The tables {st.lib}, {st.nbr}, and {st.psy} are said to be 
    the /yearbook/. */
 
gdr_sim_state_t *gdr_sim_state_new
  ( int32_t fMin,
    int32_t fMax,
    double_vec_t *cProb,
    int32_t cPrec,
    int32_t ny,
    int32_t *nbt
  );
  /* Allocates a new {gdr_sim_state_t} record {st}.
  
    The parameters {fMin,fMax,cProb,cPrec,ny,nbt} are saved in {st} as given. 
    Note that the {cProb} vector descriptor and the {nbt} table are shared, not copied.
    The {cPrec} parameter is assumed to be the number of decimal fraction digits
    that the {cProb} probabilities were rounded to.
   
    The number of individuals {st.ni} is set to zero. The current year
    field {st.yCur} is left undefined.  The vectors
    {ybr,ydt,nch,lch,par,psy} are created but their cointents and sizes are left
    undefined.    
    
    The cohort tables {st.lib} and {st.nbr} are allocated for {ny} elements. 
    Elements {st.lib[0..ny-1]} are set to {-1} (empty list) and {st.nbr[0..ny-1]}
    are set to zero. */
 
void gdr_sim_state_free(gdr_sim_state_t *st);
  /* Frees the storage used by {*st}, including internal tables,
    except the target child count distribution {st.cProb}
    and the target cohort size table {st.nbt}. */

int32_t *gdr_sim_compute_target_year_cohort_sizes(int32_t ny, int32_t iniSize, int32_t finSize);
  /* Creates a table {nbt[0..ny-1]} with the target number of
    individuals to be born on each year, suitable for
    {gdr_sim_state_new}. Currently generates exponential growth from
    {iniSize} to {finSize}, rounded to the nearest integer. */

/* SIMULATION FUNCTIONS */

void gdr_sim_run(gdr_sim_state_t *st);
  /* Simulates the evolution of a population over {ny} consecutive years
    {0..ny-1}, according to the demographic parameters in {st} and the
    target cohort sizes {st.nbt[0..ny-1]}.
    
    The procedure consider only individuals of a single sex {s}, thus
    it assumes that reproduction is asexual.
    
    When the procedure returns, the the main loop invariant above will
    hold, with the current year {st.yCur} at {ny-1}. */

void gdr_sim_check_final_state(gdr_sim_state_t *st);
  /* Checks if the state {st}, which should have {st.yCur=ny-1},
    satisfies the main loop invariant and is consistent with the 
    demography parameters in {st} and the target year birth counts
    {st.nbt[0..ny-1]}. */    

/* STATISTICS GATHERING */
  
void gdr_sim_accumulate_child_count_histogram
  ( gdr_sim_state_t *st, 
    int32_t ybrMin, 
    int32_t ybrMax, 
    int64_vec_t *cCount,
    int32_t *cMax_obs_P
  );
  /* For every expanded individual {i} in the state {st} which has year
    of birth in the range {ybrMin..ybrMax}, increments {cCount.e[c]},
    where {c} is its actual child count {st->nch[i]} (including twins
    created by cohort adjustment). The range {ybrMin..ybrMax} must be a subset
    of {0..st.ny-1}
    
    On entry, assumes that the maximum child count seen so far is {cm =
    *cMax_obs_P}, and that {cCount.e[c]} is defined for {c} in {0..cm}.
    If no individuals have been counted yet, {*cMax_obs_P} must be {-1}.
    The procedure will expand the vector {cCount} and update
    {*cMax_obs_P} as necessary. */

int64_vec_t gdr_sim_compute_child_count_histogram(gdr_sim_state_t *st, int32_t ybrMin, int32_t ybrMax);
  /* Returns a vector {cCount} such that {cCount.e[c]} is the number of
    expanded individuals in the state {st} who were born and died in the
    years {ybrMin..ybrMax} and have exactly {c} children. The range
    {ybrMin..ybrMax} must be a subset of {0..st.ny-1}. */

void gdr_sim_compare_child_count_histograms(char *title, int64_vec_t *cCount_a, int64_vec_t *cCount_b);
  /* Prints the child count histograms {cCount_a,cCount_b} side by side to {stderr}. */

void gdr_sim_compute_populations(gdr_sim_state_t *st, int32_t pop[]);
  /* The parameter {pop} should be an array of {st.ny} integers. At
    the end of the simulation, for each year {y} in {0..ny-1}, stores
    into {pop[y]} the number of individuals actually alive on year {y}.
    That is, which are expanded and have {st.ybr[i] <= y < st.ydt[i]}. */    

void gdr_sim_write_cohort_sizes
  ( char *outPrefix, 
    char *tag, 
    int32_t yStart,
    int32_t ny, 
    int32_t nbr[]
  );
  /* Writes to file "{outPrefix}-{tag}-births.txt" the cohort sizes {nbr[0..ny-1]},
    one per line, in the format "{yStart+y} {nbr[y]}".
    
    ??? Add effective alive individuals ???
    */

#endif

