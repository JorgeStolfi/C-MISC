/* See {grd_sim.h} */
/* Last edited on 2023-06-03 15:27:19 by stolfi */

#define gdr_sim_C_COPYRIGHT \
  "Duh?"

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

/* INTERNAL PROTOTYPES */

void gdr_sim_trim_state(gdr_sim_state_t *st);
  /* To be called at the end of a simulation run. Trims the vectors
    {ybr,ydt,nch,lch,par,psy} of {st} to exactly {st->ni} elements. */

void gdr_sim_create_initial_population(gdr_sim_state_t *st); 
  /* Sets the number of individuals {st.ni} to zero,
    sets the current year {st.yCur} to zero, and creates an
    initial population.
    
    The initial population will have exactly {st->nbt[0]} individuals
    born on year 0, all of them expanded and with children.
    
    The procedure may also create additional expanded individuals that
    are born in negative years and are still alive on year 0, as well 
    as individuals scheduled on positive years because they are children 
    of expanded individuals. */

void gdr_sim_one_year(gdr_sim_state_t *st);
  /* Simulates one year of the population's evolution.  
    
    The procedure assumes that the main loop invariant is true on entry
    and that {st.yCur} is in {0..ny-2}.  It simulates the evolution from
    year {st.yCur} to {st.yCur+1}.
    
    First it adjust the state so that the cohort of year {st.yCur+1} has
    the target size, {st->nbt[st.yCur+1]}. Then it expands all
    individuals in that cohort, while ensuring that the target cohort
    sizes of the following years are not exceeded. Then it increments
    {st.yCur}.  This should restore the main loop invariant. */

int32_t gdr_sim_add_individual(gdr_sim_state_t *st, int32_t ybri, int32_t pari);
  /* 
    Assumes that {st} satisfies the main loop invariant, except that 
    there may be scheduled (unexpanded) individuals born on or before {st->yCur}.
  
    Adds a new scheduled individual to the state {st}, with birth date {ybri}
    and parent pointer {pari}, and returns its index {i}. 
    
    Specifically, sets {st.ybr.e[i] = ybri}, {st.par.e[i] = pari},
    {st.ydt.e[i] = -1}, {st.nch.e[i] = -1}, and {st.lch.e[i] = -1}. The
    value {pari} should be {-1} if the parent was not created, otherwise
    it must be an expanded individual in {0..st.ni-1}; in this case, the
    child count {st.nch.e[pari]} is incremented, and {st->lch[pari]} is
    set to {i}.
    
    The year {ybri} may be negative, but if {ybri <= st.yCur} the
    individual {i} will have to be expanded in order to restore the main
    loop invariant.
    
    The returned index {i} will be the value of {st.ni} on entry; on
    exit, {st.ni} will have been incremented by one.
    
    If {ybri} is in the range {0 .. st.ny-1}, procedure also adds {i} to
    the yearbook by setting {st.lib[ybri]} and {st.psy.e[i]}, and
    increments {st.nbr[ybri]}. The procedure fails if the latter would
    exceed {st->nbt[i]}. */

void gdr_sim_expand_individual 
  ( gdr_sim_state_t *st,
    int32_t i,
    int32_t cNum,
    int32_t cAge[]
  );
  /* Assumes that the state {st} satifies the loop invariant except that
    there may be one or more individuals {j} (including {i}) that are
    only scheduled, not expanded, even though they have
    {st.ybr[j]<=st.yCur}.
    
    This procedure expands the individual {i}, by scheduling up to
    {cNum} children {j} to be born on years {st.ybr[j] = st.ybr[i] +
    cAge[k]} for {k} in {0..cNum-1}.
    
    If any of those children end up with {st.ybr[j] <= st.yCur}, the
    procedure recursively expands them, by choosing parameters
    {cNum,cAge} for each of them according to the demographic
    parameters in {st}.
    
    The procedure also sets the year of death {st.ydt[i]}, currently
    to {st.ybr[i] + st.fMax + 1}.
    
    The procedure also adds to the yearbook every new individual created
    this way, scheduled or expanded, whose birth year is in {0..st.ny-1}.
    
    In any case, the creation of a new individual {j} with birth year
    {ybrj} in {0..st.ny-1} is inhibited if it would cause {st.nbr[ybrj]}
    to exceed {st.nbt[ybrj]}. No such restriction applies if {ybrj} is
    outside {0..st.ny-1}. */

void gdr_sim_adjust_cohort_size(gdr_sim_state_t *st, int32_t y);
  /* If the year {y} is in the range {st.yCur+1 .. st.ny-1},
    adjusts the state so that {st.nbr[y]} (the number of individuals
    born on year {y}) is exactly the target number {st.nbt[y]}. 
    
    The procedure works by by repeatedy chosing, at random, one
    individual {i} born on year {y}, adding a scheduled "twin"
    {j} with same parent and birth year as {j}, until
    the target {st.nbt[y]} is met.
    
    Any individuals added with birth year {st->yCur} or less are expanded.
    Children of {j} are not added if they would cause their target coort sizes 
    to be exceeded. 
    
    The procedure fails if {y} is outside the range {st->yCur+1 ..st->ny-1}, or if {st.nbr[y]}
    is zero (since there is no one to twin) or already exceeds
    {st.nbt[y]}.
    
    If {y} is {0..st.ny-1}, the procedure has no effect. */

void gdr_sim_check_main_loop_invariant(gdr_sim_state_t *st, int32_t yCur);
  /* Checks whether the state {st} satisfies the main loop invariant 
    with {st->yCur == yCur}. */
    
/* EXPORTED FUNCTIONS */
 
gdr_sim_state_t *gdr_sim_state_new
  ( int32_t fMin,
    int32_t fMax,
    double_vec_t *cProb,
    int32_t cPrec,
    int32_t ny,
    int32_t *nbt
  )
  { 
    demand(ny >= 3*fMax, "too few years, too long life");
    
    gdr_sim_state_t *st = (gdr_sim_state_t *)notnull(malloc(sizeof(gdr_sim_state_t)), "no mem");

    st->fMin = fMin;
    st->fMax = fMax;
    st->cProb = cProb;
    st->cPrec = cPrec;
    st->ny = ny;   
    st->nbt = nbt;

    st->ni = 0;  /* Number of individuals created so far. */

    /* Will be indexed {0..ni-1}: */
    st->ybr = int32_vec_new(1000); /* Year of birth. */
    st->ydt = int32_vec_new(1000); /* Year of death. */
    st->nch = int32_vec_new(1000); /* Number of children. */
    st->lch = int32_vec_new(1000); /* Index of last child. */
    st->par = int32_vec_new(1000); /* Parent individual. */
    st->psy = int32_vec_new(1000); /* Previous indiv in cohort. */
    
    /* Will be indexed {0..ny-1}: */
    st->lib = (int32_t *)notnull(malloc(ny*sizeof(int32_t)), "no mem"); /* First indiv in cohort. */
    st->nbr = (int32_t *)notnull(malloc(ny*sizeof(int32_t)), "no mem"); /* Size of cohort. */
    for (int32_t y = 0; y < ny; y++)
      { st->lib[y] = -1; 
        st->nbr[y] = 0;
      }
    
    st->yCur = INT32_MIN; /* Last expanded year. */
    
    return st;
  }
 
void gdr_sim_trim_state(gdr_sim_state_t *st)
  { int32_vec_trim(&(st->ybr), st->ni);
    int32_vec_trim(&(st->ydt), st->ni);
    int32_vec_trim(&(st->nch), st->ni);
    int32_vec_trim(&(st->lch), st->ni);
    int32_vec_trim(&(st->par), st->ni);
    int32_vec_trim(&(st->psy), st->ni);
  }

void gdr_sim_state_free(gdr_sim_state_t *st)
  {
    free(st->ybr.e);
    free(st->ydt.e);
    free(st->nch.e);
    free(st->lch.e);
    free(st->par.e);
    free(st->psy.e);

    free(st->lib);
    free(st->nbr);
    
    free(st);
  }

int32_t *gdr_sim_compute_target_year_cohort_sizes(int32_t ny, int32_t iniSize, int32_t finSize)
  { 
    demand((iniSize > 0) && (finSize > 0), "year sizes must be positive");
    int32_t *nbt = (int32_t *)notnull(malloc(ny*sizeof(int32_t)), "no mem");
    double rat = ((double)finSize)/((double)iniSize);
    for (int32_t y = 0; y < ny; y++)
      { nbt[y] = (int32_t)floor(0.5 + iniSize*exp(y*log(rat)/ny)); }
    return nbt;
  }

void gdr_sim_run(gdr_sim_state_t *st)
  { bool_t debug = FALSE;
    if (debug) { fprintf(stderr, "  > %s\n", __FUNCTION__); }

    int32_t ny = st->ny;
    assert(ny >= 0);

    /* Initial population: */
    gdr_sim_create_initial_population(st);
    assert(st->yCur == 0);
    assert(st->nbr[0] == st->nbt[0]);
    
    /* Subsequent years: */
    for (int32_t y = 1; y < ny; y++)
      { /* Simulate evolution from year {y-1} to year {y}: */
        gdr_sim_one_year(st);
      }
    if (debug) { fprintf(stderr, "    created %d individuals\n", st->ni); }
    gdr_sim_trim_state(st);

    if (debug) { fprintf(stderr, "  < %s\n", __FUNCTION__); }
  }
       
void gdr_sim_create_initial_population(gdr_sim_state_t *st)
  { bool_t debug = FALSE;
    if (debug) { fprintf(stderr, "   > %s\n", __FUNCTION__); }

    int32_t cMax = st->cProb->ne - 1;
    
    /* Clear all individual tables: */
    st->ni = 0;
    for (int32_t y = 0; y < st->ny; y++) { st->nbr[y] = 0; st->lib[y] = -1; }
    
    /* Make year 0 the current one: */
    st->yCur = 0;
    
    if (debug) { fprintf(stderr, "     trying to get %d births on year 0\n", st->nbt[0]); }
    
    while (st->nbr[0] < st->nbt[0])
      { if (debug) { fprintf(stderr, "       still only %d births\n", st->nbr[0]); }
        /* Generate an individual, chosing birth date {ybri} and children {cNum,cAge[0..cNum-1]}: */
        int32_t cNum;
        int32_t cAge[cMax+1];
        gdr_demo_throw_children(st->cProb, st->fMin, st->fMax, &cNum, cAge);
        if (debug) { fprintf(stderr, "       generated guy with cNum = %d\n", cNum); }
        if (cNum == 0)
          { /* Better not waste slots with it. */
            if (debug) { fprintf(stderr, "         ignored\n"); }
          }
        else 
          { /* Find max {aMax} of {cAge[0..cNum-1]}: */
            if (debug) { fprintf(stderr, "         cAge ="); }
            int32_t aMax = -1;
            for (int32_t ka = 0; ka < cNum; ka++) 
              { if (debug) { fprintf(stderr, " %d", cAge[ka]); }
                if (cAge[ka] > aMax) { aMax = cAge[ka]; }
              }
            int32_t ybri = 0 - int32_abrandom(0, aMax); 
            if (debug) { fprintf(stderr, " ybr = %d\n", ybri); }
            /* Will be born at year 0 or will have a child at year 0 or later. */
            int32_t i = gdr_sim_add_individual(st, ybri, -1);
            assert((i >= 0) && (i < st->ni));
            assert(st->ybr.e[i] == ybri);
            assert(st->ydt.e[i] == -1);
            assert(st->par.e[i] == -1);
            assert(st->nch.e[i] == -1);
            assert(st->lch.e[i] == -1);
            /* Recursively expand individual {i}: */
            gdr_sim_expand_individual(st, i, cNum, cAge);
            /* Must have been expanded: */
            assert(st->nch.e[i] >= 0); 
            assert(st->ydt.e[i] >= ybri); 
            assert((st->nch.e[i] == 0) == (st->lch.e[i] == -1));
          }
      }

    if (debug) { fprintf(stderr, "   < %s\n", __FUNCTION__); }
  }

void gdr_sim_one_year(gdr_sim_state_t *st)
  { bool_t debug = FALSE;
    if (debug) { fprintf(stderr, "    > %s\n", __FUNCTION__); }

    if (debug) { fprintf(stderr, "      current year %d, %d individuals created so far\n", st->yCur, st->ni); }
    
    /* Make sure that the cohort of year {yCur+1} has the target size: */
    int32_t y1 = st->yCur + 1;
    assert((y1 >= 0) && (y1 < st->ny));
    int32_t nhGoal = st->nbt[y1]; /* Ideal cohort size. */
    if (debug) { fprintf(stderr, "      trying to adjust cohort of year %d to %d births\n", y1, nhGoal); }
    gdr_sim_adjust_cohort_size(st, y1); 
    
    /* Advance the current year: */
    st->yCur = y1;
    
    /* Expand the individuals scheduled to be born on year {yCur}: */
    /* It is safe to scan the list of cohort {yCur} because it will not change: */
    int32_t cMax = st->cProb->ne -1;
    int32_t i = st->lib[y1];
    while (i != -1)
      { assert((i >= 0) && (i < st->ni)); /* Individual must exist. */
        assert(st->ybr.e[i] == y1);   /* It must be born on year {y1}. */
        assert(st->ydt.e[i] == -1);   /* It must be scheduled not expanded. */
        assert(st->nch.e[i] == -1);   /* It must be scheduled not expanded. */
        assert(st->lch.e[i] == -1);   /* It must be scheduled not expanded. */
        /* Choose its children: */
        int32_t cNum;
        int32_t cAge[cMax+1];
        gdr_demo_throw_children(st->cProb, st->fMin, st->fMax, &cNum, cAge);
        gdr_sim_expand_individual(st, i, cNum, cAge);
        /* It must be expanded now: */
        assert(st->ydt.e[i] >= st->ybr.e[i]);   /* It must be expanded now. */
        assert(st->nch.e[i] >= 0);   /* It must be expanded now. */
        assert((st->nch.e[i] == 0) == (st->lch.e[i] == -1));
        /* Get next individual in cohort: */
        i = st->psy.e[i]; 
      }

    if (debug) { fprintf(stderr, "     < %s\n", __FUNCTION__); }
  }
            
int32_t gdr_sim_add_individual(gdr_sim_state_t *st, int32_t ybri, int32_t pari)
  {
    int32_t i = st->ni;
    
    int32_vec_expand(&(st->ybr), i); st->ybr.e[i] = ybri;
    int32_vec_expand(&(st->ydt), i); st->ydt.e[i] = -1;
    int32_vec_expand(&(st->par), i); st->par.e[i] = pari;
    int32_vec_expand(&(st->nch), i); st->nch.e[i] = -1; /* Mark as not expanded yet. */
    int32_vec_expand(&(st->lch), i); st->lch.e[i] = -1; /* No children yet. */

    if (pari != -1)
      { assert((pari >= 0) && (pari < i)); 
        /* Parent must be expanded. */
        assert(st->nch.e[pari] >= 0); 
        assert((st->nch.e[pari] == 0) == (st->lch.e[pari] == -1));
        /* Assign paternity: */
        st->nch.e[pari]++;   /* Indiv {pari} got one more child. */
        st->lch.e[pari] = i; /* It is now the last child of {pari}. */
      }

    st->ni++;

    /* Add to yearbook: */
    if ((ybri >= 0) && (ybri < st->ny))
      { int32_vec_expand(&(st->psy), i); 
        st->psy.e[i] = st->lib[ybri];
        st->lib[ybri] = i;
        st->nbr[ybri]++;
        assert(st->nbr[ybri] <= st->nbr[ybri]); /* Caller beware. */
      }

    return i;
  }

void gdr_sim_expand_individual 
  ( gdr_sim_state_t *st,
    int32_t i,
    int32_t cNum,
    int32_t cAge[]
  )  
  { bool_t debug = FALSE;
    if (debug) { fprintf(stderr, "      > %s\n", __FUNCTION__); }

    if (debug) { fprintf(stderr, "        expanding %d", i); }
    demand((i >= 0) && (i < st->ni), "invalid individual {i}");
    
    if (debug) 
      { fprintf(stderr, " ybr = %d ydt = %d par = %d", st->ybr.e[i], st->ydt.e[i], st->par.e[i]); 
        fprintf(stderr, " nch = %d lch = %d\n", st->nch.e[i], st->lch.e[i]);
      }
    demand(st->ybr.e[i] <= st->yCur, "can't expand {i}, too early");
    demand(st->ydt.e[i] == -1, "is already expanded?");
    demand(st->nch.e[i] == -1, "is already expanded?");
    demand(st->lch.e[i] == -1, "inconsistent {lch}");

    /* Define year of death: */
    st->ydt.e[i] = st->ybr.e[i] + st->fMax + 1; /* Fixed lifespan for now. */

    /* Mark {i} as expanded: */
    st->nch.e[i] = 0;
    
    /* Schedule the children: */
    if (debug) { fprintf(stderr, "        trying to add %d children\n", cNum); }
    int32_t j0 = st->ni;  /* Index of first child of {i}, if any. */
    int32_t j1 = j0 - 1;  /* Index of last child of {i}, or {-1} if none. */
    for (int32_t kc = 0; kc < cNum; kc++)
      { int32_t ybrj = st->ybr.e[i] + cAge[kc]; /* Year when child should be born. */
        bool_t ok = (ybrj < 0) || (ybrj >= st->ny) || (st->nbr[ybrj] < st->nbt[ybrj]);
        if (ok)
          { /* Can add child to the cohort {ybrj}: */
            int32_t j = gdr_sim_add_individual(st, ybrj, i);
            if (debug) { fprintf(stderr, "          added child %d with ybr = %d\n", j, ybrj); }
            assert(j == j1 + 1); /* The children of {i} must end in consecutive positions. */
            j1 = j;              /* Remember last child actually created. */
            assert((j >= 0) && (j < st->ni));     /* The child {j} must exist. */
            assert(st->nch.e[i] == (j1 - j0 + 1));  /* We created that many children. */
            assert(st->lch.e[i] == j1);             /* And {j1} was the last one. */
            
            assert((st->ybr.e[j] == ybrj) && (st->par.e[j] == i)); /* Insertion must have succeeded. */
            assert(st->nch.e[j] == -1); /* Child {j} must still be unexpanded. */
            assert(st->ydt.e[j] == -1); /* Child {j} must still be unexpanded. */
            assert(st->lch.e[j] == -1); /* Child {j} must still be unexpanded. */
          }
        else
          { 
            if (debug) { fprintf(stderr, "          suppressed a child with ybr = %d due to cohort limit\n", ybrj); }
          }
      }
    if (debug) { fprintf(stderr, "        added %d children %d .. %d\n", j1-j0+1, j0, j1); }
    int32_t nch_act = j1 - j0 + 1;
    assert(st->nch.e[i] == nch_act);
    assert(st->lch.e[i] == (nch_act == 0 ? -1 : j1));

    /* Expand recursively any children born on or before year {st->yCur}: */
    int32_t cMax = st->cProb->ne - 1;
    for (int32_t j = j0; j <= j1; j++)
      { if (st->ybr.e[j] <= st->yCur)
          { /* Must expand {j} too: */
            assert(st->ydt.e[j] == -1); /* Child {j} must still be unexpanded. */
            assert(st->nch.e[j] == -1); /* Child {j} must still be unexpanded. */
            assert(st->lch.e[j] == -1); /* Child {j} must still be unexpanded. */
            int32_t cNum_j;
            int32_t cAge_j[cMax+1];
            gdr_demo_throw_children(st->cProb, st->fMin, st->fMax, &cNum_j, cAge_j);
            gdr_sim_expand_individual(st, j, cNum_j, cAge_j);
            /* Child {j} must be expanded now: */
            assert(st->ydt.e[j] >= st->ybr.e[j]); 
            assert(st->nch.e[j] >= 0); 
            assert((st->nch.e[j] == 0) == (st->lch.e[j] == -1));
          }
      }

    if (debug) { fprintf(stderr, "      < %s\n", __FUNCTION__); }
  }
 
void gdr_sim_adjust_cohort_size(gdr_sim_state_t *st, int32_t y)
  {
    bool_t debug = FALSE;
    if (debug) { fprintf(stderr, "      > %s\n", __FUNCTION__); }

    demand((y > st->yCur) && (y < st->ny), "cannot adjust cohort of this year");

    int32_t nh = st->nbr[y]; /* Size of cohort. */
    int32_t nhGoal = st->nbt[y]; /* Ideal size. */
    
    if (nh == nhGoal)
      { /* Nothing to do: */ return; }
    else if (nh == 0)
      { fprintf(stderr, "        cohort of year %d is empty, cannot increase it\n", y); }
    else if (nh >st->nbt[y])
      { fprintf(stderr, "        cohort of year %d has %d births, cannot shrink to %d\n", y, nh, nhGoal); }
    else
      { if (debug) { fprintf(stderr, "        cohort of year %d has %d births, twinning to %d\n", y, nh, nhGoal); }
        /* Collect the cohort individuals {h[0..nh-1]}: */
        int32_t h[nh]; /* Indices of individuals in cohort {y}. */
        int32_t i = st->lib[y]; /* Next individual in cohort. */
        for (int32_t kh = 0; kh < nh; kh++)
          { assert((i >= 0) & (i < st->ni)); /* Individual must exist. */
            assert(st->ybr.e[i] == y); /* Must be of the right cohort. */
            h[kh] = i;
            i = st->psy.e[i];
          }
        assert(i == -1); /* List length must be exactly {nh}. */
        
        /* Create twins until target is reached: */
        while (st->nbr[y] <nhGoal)
          { /* Pick a random individual {i} to twin: */
            int32_t kh = int32_abrandom(0, nh-1);
            int32_t i = h[kh];
            /* Add his twin {j}: */
            int32_t pari = st->par.e[i];
            gdr_sim_add_individual(st, y, pari);
          }
      }
    if (debug) { fprintf(stderr, "      < %s\n", __FUNCTION__); }
  }

void gdr_sim_accumulate_child_count_histogram
  ( gdr_sim_state_t *st, 
    int32_t ybrMin, 
    int32_t ybrMax, 
    int64_vec_t *cCount, 
    int32_t *cMax_obs_P
  )
  {
    bool_t debug = FALSE;
    if (debug) { fprintf(stderr, "      > %s\n", __FUNCTION__); }
    
    demand((ybrMin >= 0) && (ybrMax < st->ny), "invalid year range");
    
    int32_t ni = st->ni;
      
    /* Increment {cCount.e[c]} with num of individuals with {c} children: */
    if (debug) { fprintf(stderr, "        looking at years %d .. %d\n", ybrMin, ybrMax); }
    int32_t niCt = 0;
    int32_t cm = (*cMax_obs_P);
    if (debug) { fprintf(stderr, "        cm = %d cCount.ne = %d\n", cm, cCount->ne); }
    assert(cm <= (int32_t)cCount->ne-1);
    for (int32_t i = 0; i < ni; i++)
      { int32_t nchi = st->nch.e[i];
        if (nchi != -1)
          { /* Individual is expanded: */
            int32_t ybri = st->ybr.e[i];
            int32_t ydti = st->ydt.e[i];
            assert(nchi >= 0);
            assert(ydti >= ybri);
            if (debug && (nchi != 0)) { fprintf(stderr, "          i = %d ybr = %d ydt = %d nch = %d\n", i, ybri, ydti, nchi); }
            if ((ybri >= ybrMin) && (ybri <= ybrMax))
              { /* Count individual: */
                /* Note that {nchi} may be greater than current {cCount.ne-1}. */
                gdr_demo_increment_count(cCount, &cm, nchi, 1);
                niCt++;
              }
          }
      }
    if (debug) { fprintf(stderr, "       counted %d individuals, cMax = %d\n", niCt, cm); }

    (*cMax_obs_P) = cm;

    if (debug) { fprintf(stderr, "      < %s\n", __FUNCTION__); }
  }
  
int64_vec_t gdr_sim_compute_child_count_histogram(gdr_sim_state_t *st, int32_t ybrMin, int32_t ybrMax)
  { 
    int64_vec_t cCount = int64_vec_new(0);
    int32_t cMax_obs = -1;
    gdr_sim_accumulate_child_count_histogram(st, ybrMin, ybrMax, &(cCount), &(cMax_obs));
    assert(cMax_obs < (int32_t)cCount.ne);
    /* Remove any zero trail (paranoia): */
    while ((cMax_obs >= 0) && (cCount.e[cMax_obs] == 0)) { cMax_obs--; }
    int64_vec_trim(&(cCount), cMax_obs+1);
    return cCount;
  }

void gdr_sim_compare_child_count_histograms(char *title, int64_vec_t *cCount_a, int64_vec_t *cCount_b)
  { int32_t c_max = (int32_t)imax(cCount_a->ne, cCount_b->ne) - 1;
    fprintf(stderr, "        comparing histograms of child counts %d..%d %s:\n", 0, c_max, title);
    for (int32_t c = 0; c <= c_max; c++)
      { int64_t cta = (c < (int32_t)cCount_a->ne ? cCount_a->e[c] : 0);
        int64_t ctb = (c < (int32_t)cCount_b->ne ? cCount_b->e[c] : 0);
        fprintf(stderr, "          %3d %6ld %6ld\n", c, cta, ctb);
      }
  }

void gdr_sim_compute_populations(gdr_sim_state_t *st, int32_t pop[])
  {
    int32_t ny = st->ny;
    int32_t ni = st->ni;
    demand(st->yCur == ny-1, "state is not final"); 
    
    /* Clear population counts: */
    for (int32_t y = 0; y < ny; y++) { pop[y] = 0; }
    
    /* Scan individuals: */
    for (int32_t i = 0; i < ni; i++)
      { int32_t nchi = st->nch.e[i];
        if (nchi != -1)
          { /* Individual is expanded: */
            int32_t ybri = st->ybr.e[i]; /* Year of birth. */
            int32_t ydti = st->ydt.e[i]; /* Year of death. */
            assert(ydti != -1);
            assert(ydti >= ybri);
            if (nchi > 0)
              { /* Consistency of parenting dates: */
                int32_t j = st->lch.e[i]; /* Last child of {i}. */
                assert(j != -1); /* Since it is expanded and has children. */
                assert((j >= 0) && (j < ni));
                assert(st->ybr.e[j] <= st->ydt.e[i]); /* Can't have children after death. */
              }
            /* Count it on years it is actually alive: */
            int32_t y0 = (ybri < 0 ? 0 : ybri);
            int32_t y1 = (ydti >= ny ? ny-1 : ydti);
            for (int32_t y = y0; y <= y1; y++) { pop[y]++; }
          }
      }
  }
  
void gdr_sim_write_cohort_sizes
  ( char *outPrefix, 
    char *tag, 
    int32_t ny, 
    int32_t yStart, 
    int32_t nbr[]
  )
  {
    char *fname = NULL;
    asprintf(&fname, "%s-%s-births.txt", outPrefix, tag);
    FILE *wr = open_write(fname, TRUE);
    free(fname);
    
    for (int32_t y = 0; y < ny; y++)
      { fprintf(wr, "%+6d %10d\n", yStart+y, nbr[y]); }
    fclose(wr);
  }

void gdr_sim_check_main_loop_invariant(gdr_sim_state_t *st, int32_t yCur)
  {
    bool_t debug = FALSE;
    if (debug) { fprintf(stderr, "  > %s\n", __FUNCTION__); }

    demand(st != NULL, "state is {NULL}");
    demand(st->cProb != NULL, "{cProb} is {NULL}");
    demand(st->fMin > 0, "invalid {fMin}");
    demand(st->fMax >= st->fMin, "invalid {fMax}");
    
    int32_t ny = st->ny;
    demand(ny >= 2, "bad {ny}");
    
    demand(st->nbt != NULL, "{nbt} is {NULL}");
    for (int32_t y = 0; y < ny; y++) 
      { demand(st->nbt[y] >= 1, "invalid {nbt[y]}"); }
    
    int32_t cMax = st->cProb->ne - 1;
    gdr_demo_distr_check(st->cProb, cMax, 1.0, st->cPrec);
    
    demand((0 <= st->yCur) && (st->yCur < st->ny), "{yCur} out of range");
    demand(st->yCur == yCur, "wrong {yCur}");

    int32_t ni = st->ni;
    demand(ni > 0, "invalid {ni}");
    demand(st->ybr.ne >= ni, "invalid {ybr.ne}");
    demand(st->ydt.ne >= ni, "invalid {ydt.ne}");
    demand(st->nch.ne >= ni, "invalid {nch.ne}");
    demand(st->lch.ne >= ni, "invalid {lch.ne}");
    demand(st->par.ne >= ni, "invalid {par.ne}");
    demand(st->psy.ne >= ni, "invalid {psy.ne}");
    int32_t *nch_act = in_alloc(ni); /* Actual child counts. */
    int32_t *lch_act = in_alloc(ni); /* Actual last child index. */
    for (int32_t i = 0; i < ni; i++) { nch_act[i] = 0; lch_act[i] = -1; }
    int32_t ni_main = 0; /* Count of individuals born in {0..ny-1}. */
    for (int32_t i = 0; i < ni; i++)
      { int32_t ybri = st->ybr.e[i];
        int32_t ydti = st->ydt.e[i];
        int32_t nchi = st->nch.e[i];
        if (debug) { fprintf(stderr, "    checking i = %d ybr = %d ydt = %d nch = %d\n", i, ybri, ydti, nchi); }
        if (ybri <= yCur)
          { /* Must be expanded: */
            demand(nchi >= 0, "indiv {i} should have {nch} non-negative");
            demand(ybri <= ydti, "death before birth");
            demand((nchi == 0) == (st->lch.e[i] == -1), "inconsistent {lch[i]}");
          }
        else
          { /* Must be scheduled only: */
            demand(ydti == -1, "indiv {i} should have {ydt=-1}");
            demand(nchi == -1, "indiv {i} should have {nch=-1}");
            demand(st->lch.e[i] == -1, "indiv {i} should have {lch=-1}");
          }
        if ((ybri >= 0) && (ybri < ny)) { ni_main++; }
        int32_t pari = st->par.e[i];
        if (pari == -1)
          { /* Individual {i} must have been created at initialization: */
            demand(ybri <= 0, "individual should have a parent");
          }
        else
          { /* Check parent, increment its children: */
            demand((pari >= 0) && (pari < i), "invalid {par[i]}");
            int32_t ybrp = st->ybr.e[pari]; /* Birth year of parent. */
            int32_t ydtp = st->ydt.e[pari]; /* Death year of parent. */
            demand(ybri <= ydtp, "parent was dead");
            demand(ybri >= ybrp + st->fMin, "parent too young");
            demand(ybri <= ybrp + st->fMax, "parent too old");
            nch_act[pari]++;
            lch_act[pari] = i;
          }
      }

    /* Check children counts: */
    for (int32_t i = 0; i < ni; i++)
      { int32_t nchi = st->nch.e[i];
        if (nchi == -1)
          { /* Unexpanded, should have no children: */
            demand(nch_act[i] == 0, "unexpanded has children");
            demand(lch_act[i] == -1, "unexpanded has last child");
          }
        else
          { /* Expanded: */
            demand(nchi == nch_act[i], "inconsistent {nch,par}");
            demand(st->lch.e[i] == lch_act[i], "inconsistent {lch,par}");
          }
      }
      
    /* Check yearbook: */
    int32_t ni_book = 0; /* Individuals in yearbook. */
    for (int32_t y = 0; y < ny; y++)
      { int32_t nbry = st->nbr[y] ;
        demand((nbry >= 0) && (nbry <= ni_main), "invalid {nbr[y]}");
        int32_t nbr_obs = 0; /* Count of indivs in list of year {y}. */
        int32_t i = st->lib[y];
        while (i != -1)
          { demand((i >= 0) && (i < ni), "invalid indiv in yearbook list");
            demand(st->ybr.e[i] == y, "birth year inconsistency in yearbook");
            nbr_obs++;
            i = st->psy.e[i];
          }
        demand(nbry == nbr_obs, "incorrect {nbr[y]}");
        ni_book += nbry;
      }
    demand(ni_book == ni_main, "missing indivs in yearbook");

    free(nch_act);
    free(lch_act);
    
    if (debug) { fprintf(stderr, "  < %s\n", __FUNCTION__); }
  }

void gdr_sim_check_final_state(gdr_sim_state_t *st)
  { 
    gdr_sim_check_main_loop_invariant(st, st->ny-1);
  }
