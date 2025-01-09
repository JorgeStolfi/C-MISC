/* See {gdr_lineage.h} */

/* Last edited on 2023-06-03 15:28:04 by stolfi */

#define gdr_lineage_C_COPYRIGHT "Duh"
    
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
#include <in.h>
#include <rn.h>
#include <vec.h> 

#include <gdr_demo.h>
#include <gdr_sim.h>
#include <gdr_lineage.h>

/* INTERNAL PROTOTYPES */

void gdr_lineage_count_surviving_descendants(gdr_sim_state_t *st, int32_t sud[], int32_t lsc[]);
  /* The vector {sud} must have {st->ni} elements. For each individual
    {i} in {0..st.ni-1}, stores in {sud[i]} the number of descendants of
    {i} which are actually alive on year {st.ny-1}.
    
    Also stores into {lsc[i]} the index of the child of {i} with surviving descendants
    that has largest birth year. If there is no such child, sets {lsc[i]} to {-1}.
    
    The state {st} must be final (with {st.yCur==st.ny-1}). */

void gdr_write_quantiles(FILE *wr, int32_t nr, double vals[], int32_t prec);
      /* Writes to {wr} the quantiles {frac[0..nq-1]} of the numbers {vals[0..nr-1]}. */

/* IMPLEMENTATIONS */

void gdr_lineage_compute_surviving_counts(gdr_sim_state_t *st, int32_t nsl[])
  {
    bool_t debug = FALSE;
    if (debug) { fprintf(stderr, "  > %s\n", __FUNCTION__); }

    gdr_sim_check_final_state(st);

    int32_t ny = st->ny;
    int32_t ni = st->ni;
    if (debug) { fprintf(stderr, "    ny = %d ni = %d\n", ny, ni); }

    /* Counts of surviving descendants per individual: */
    if (debug) { fprintf(stderr, "    counting descendants\n"); }
    int32_t *sud = in_alloc(ni); /* Number of descendants of each individual alive on {ny-1}. */
    int32_t *lsc = in_alloc(ni); /* Last-born child of each individual that has surviving descendants, or {-1}. */
    gdr_lineage_count_surviving_descendants(st, sud, lsc);

    /* Clear {nsl}: */
    if (debug) { fprintf(stderr, "    initializing {nsl}\n"); }
    for (uint32_t y = 0;  y < ny; y++) { nsl[y] = 0; }
    
    /* Count individuals with surviving lineages for the years that they count for: */
    if (debug) { fprintf(stderr, "    scanning individuals\n"); }
    for (uint32_t i = 0;  i < ni; i++) 
      { int32_t nchi = st->nch.e[i]; /* Number of children. */
        int32_t ybri = st->ybr.e[i]; /* Year of birth. */
        int32_t ydti = st->ydt.e[i]; /* Year of birth. */
        if (nchi == -1)
          { /* Still unexpanded, must be past the end: */
            assert(ybri >= ny); 
            assert(ydti == -1); 
          }
        else 
          { /* Individual is expanded: */ 
            assert(ybri <= st->yCur); 
            assert(ydti >= ybri);
            int32_t ycti; /* Last year for which {i} counts: */
            if ((ydti <= 0) || (ybri >= ny))
              { /* Individual life is out of range, does not count */
                ycti = ybri - 1; /* So that {ybri..ycti} is empty. */
              }
            else if (ydti > ny-1)
              { /* The individual itself is surviving.  Its whole life counts: */
                ycti = ydti - 1;
              }
            else if (sud[i] == 0)
              { /* Individual life overlaps {0..ny-1} but his lineage is not surviving: */
                ycti = ybri - 1; /* So that {ybri..ycti} is empty. */
              }
            else 
              { /* Individual has a surviving lineage: */
                /* Since he is not itself surviving, it must have children: */
                assert(nchi > 0); 
                /* Find age at which last child with surviving lineage is born: */
                int32_t lsci = lsc[i];
                assert(lsci != -1);
                assert((lsci > i) && (lsci < ni));
                /* The lineage of {i} counts from birth of {i} to just before that of {lsci}: */
                ycti = st->ybr.e[lsci] - 1;
              }

            /* Count {i} for the years {ybri..ycti} inter {0..ny-1}: */
            int32_t y0 = (ybri < 0 ? 0 : ybri);
            int32_t y1 = (ycti >= ni ? ni-1 : ycti);
            for (int32_t y = y0; (y <= y1); y++) { nsl[y]++; }
          }
      }

    free(sud);
    free(lsc);
    
    if (debug) { fprintf(stderr, "  < %s\n", __FUNCTION__); }
  }

void gdr_lineage_count_surviving_descendants(gdr_sim_state_t *st, int32_t sud[], int32_t lsc[])
  { 
    bool_t debug = FALSE;
    if (debug) { fprintf(stderr, "    > %s\n", __FUNCTION__); }
    
    int32_t ni = st->ni;
    int32_t ny = st->ny;
    demand(st->yCur == ny - 1, "state is not final");

    /* Set {sud[i]=1} for individuals actually alive on last year, others zero: */
    if (debug) { fprintf(stderr, "      initializing\n"); }
    for (uint32_t i = 0;  i < ni; i++) 
      { bool_t alive;
        int32_t ybri = st->ybr.e[i];  /* Year of birth. */
        int32_t ydti = st->ydt.e[i];  /* Year of death. */
        int32_t nchi = st->nch.e[i];  /* Number of children. */
        if (debug) { fprintf(stderr, "      i = %d ybr = %d ydt = %d nch = %d\n", i, ybri, ydti, nchi); }
        if (nchi != -1)
          { /* Individual is expanded. */
            assert((ydti != -1) && (ydti >= ybri));
            /* Individual {i} is actually alive up to {ydt[i]-1}: */
            alive = ((ybri <= ny-1) && (ydti > ny-1));
          }
        else
          { /* Individual {i} is not expanded, must be born after {ny-1}: */
            assert(ybri >= ny);
            assert(ydti == -1);
            /* Individual {i} was not actually alive on {ny-1} : */
            alive = FALSE;
          }
        sud[i] = (alive ? 1 : 0);
        lsc[i] = -1;
      }
      
    /* Propagate {sud} backwards through {par}, set {lsc}: */
    if (debug) { fprintf(stderr, "      propagating\n"); }
    for (int32_t i = ni-1; i >= 0; i--)
      { int32_t p = st->par.e[i];
        if (p != -1)
          { assert((p >= 0) && (p < i));
            assert(st->ybr.e[p] < st->ybr.e[i]);
            sud[p] += sud[i];
            int32_t lscp = lsc[p];
            if (lscp != -1)
              { assert ((lscp > p) && (lscp < ni));
                if (st->ybr.e[i] > st->ybr.e[lscp]) 
                  { lsc[p] = i; }
              }
            else
              { lsc[p] = i; }
          }
      }
      
    if (debug) { fprintf(stderr, "    < %s\n", __FUNCTION__); }
  }
  
void gdr_lineage_compute_sex_ratio_table
  ( int32_t ny, 
    int32_t nr, 
    int32_t nsl_sv0[],
    int32_t nsl_sv1[],
    double rsl_sv[]
  )
  {
    for (uint32_t y = 0;  y < ny; y++)
      { for (uint32_t r = 0;  r < nr; r++)
          { int32_t k = y*nr + r;
            double nsl0k = (double)(nsl_sv0[k]); assert(nsl0k > 0.0);
            double nsl1k = (double)(nsl_sv1[k]); assert(nsl1k > 0.0);
             rsl_sv[k] = nsl1k/nsl0k;
           }
       }
  }

void gdr_lineage_find_notable_ratio_parameters
  ( int32_t ny, 
    int32_t nr, 
    double rsl_sv[],
    int32_t *gcPlat_P,
    double *rslPlat_P
  )
  { 
    /* Compute the median {rslMed[0..ny-1]} for all years, find true max: */
    double rslMed[ny];    /* Medians. */
    double rslMax = 0.0;  /* Median with max abs value. */
    for (uint32_t y = 0;  y < ny; y++)
      { /* Get sex lineage fractions of year {y}, sorted: */
        double rslg[nr]; /* Copy of row {y} of table. */
        for (uint32_t r = 0;  r < nr; r++) { rslg[r] = rsl_sv[nr*y + r]; }
        gdr_demo_sort_doubles(nr, rslg);
        double med = gdr_demo_find_quantile(nr, rslg, 0.50);
        rslMed[y] = med;
        if (fabs(med) > rslMax) { rslMax = med; }
      }
    assert(rslMax != 0.0);
    
    /* Find the first and last reverse gens of the plateau: */
    double rslPlat = 0.90 * rslMax; /* Plateau value.*/
    if (rslPlat >= 10)
      { rslPlat = floor(rslPlat + 0.5); }
    else
      { rslPlat = floor(rslPlat/0.1 + 0.5)*0.1; }
    
    int32_t gcPlat = -1;
    /* Scan years in reverse chron order (most recent first): */
    for (int32_t y = ny-1; y >= 0; y--)
      { int32_t gc = ny - 1 - y; /* Complemented year (0 = last). */
        if (fabs(rslMed[y]) >= rslPlat) 
          { if (gcPlat == -1) { gcPlat = gc; } }
      }
    (*gcPlat_P) = gcPlat;
    (*rslPlat_P) = rslPlat;
  }

/* OUTPUT */

#define gdr_quantiles_NUM  3
#define gdr_quantiles_FRACTIONS  0.02500,  0.50000,  0.97500
  /* The quantiles of values that are to be written out to disk by
    {gdr_lineage_write_counts_median_range} and
    {gdr_lineage_write_sex_ratio_median_range}. */

void gdr_lineage_write_counts_median_range
  ( char *outPrefix,
    char *tag,
    int32_t s,
    int32_t ny, 
    int32_t yStart,
    int32_t nr,
    int32_t nsl_sv[],
    int32_t pop_sv[]
  )
  { char *fname = NULL;
    char *fname = jsprintf("%s-%s-%d-nlins.txt", outPrefix, tag, s);
    FILE *wr = open_write(fname, TRUE);
    free(fname);

    /* Run quantile fractions and their column headers: */ 
    for (uint32_t y = 0;  y < ny; y++)
      { 
        /* Write external year: */
        fprintf(wr, "%+6d", yStart + y);
        
        /* Write quantiles of population: */
        double pops[nr]; /* Copy of row {y} of {pop_sv}. */
        for (uint32_t r = 0;  r < nr; r++) 
          { pops[r] = (double)pop_sv[nr*y + r]; }
        gdr_write_quantiles(wr, nr, pops, 1);
        
        /* Get lineage counts of year {y}, sorted: */
        double psls[nr]; /* Row {y} of {nsl_sv} as percentage of {pop_sv}. */
        for (uint32_t r = 0;  r < nr; r++) 
          { psls[r] = 100*(double)nsl_sv[nr*y + r]/(double)pop_sv[nr*y + r]; }
        gdr_write_quantiles(wr, nr, psls, 6);

       fprintf(wr, "\n");
      }
    fclose(wr);
  }
  
void gdr_lineage_write_sex_ratio_samples
  ( char *outPrefix,
    char *tag,
    int32_t ny, 
    int32_t yStart, 
    int32_t nr, 
    int32_t nrSample,
    double rsl_sv[]
  )
  {
    char *fname = jsprintf("%s-%s-rsamp.txt", outPrefix, tag);
    FILE *wr = open_write(fname, TRUE);
    free(fname);
  
    /* Choose the indices {which[0..nrSample]}:*/
    int32_t which[nr]; /* Indices of runs to write. */
    for (uint32_t r = 0;  r < nr; r++) { which[r] = r; }
    for (uint32_t q = 0;  q < nrSample; q++) 
      { int32_t j = int32_abrandom(q, nr-1);
        if (j != q) { int32_t t = which[q]; which[q] = which[j]; which[j] = t; }
      }

    for (uint32_t q = 0;  q < nrSample; q++)
      { int32_t r = which[q];
        assert((r >= 0) && (r <nr));
        for (uint32_t y = 0;  y < ny; y++)
          { double rsl_gr = rsl_sv[y*nr + r];
            /* Perturbe values slightly: */
            rsl_gr += dabrandom(-0.01, +0.01);
            fprintf(wr, "%+6d %3d %10.3f\n", yStart + y, q, rsl_gr);
          }
        fprintf(wr, "\n"); /* Blank line between plots. */
      }
    fclose(wr);
  }

void gdr_lineage_write_sex_ratio_median_range
  ( char *outPrefix,
    char *tag,
    int32_t ny, 
    int32_t yStart, 
    int32_t nr, 
    double rsl_sv[]
  )
  { char *fname = NULL;
    char *fname = jsprintf("%s-%s-ratio.txt", outPrefix, tag);
    FILE *wr = open_write(fname, TRUE);
    free(fname);

    for (uint32_t y = 0;  y < ny; y++)
      { 
        /* Write external year: */
        fprintf(wr, "%+6d ", yStart + y);

        /* Get sex lineage fractions of year {y}, sorted: */
        double rsls[nr]; /* Copy of row {y} of table. */
        for (uint32_t r = 0;  r < nr; r++) { rsls[r] = rsl_sv[nr*y + r]; }
        gdr_write_quantiles(wr, nr, rsls, 3);
        fprintf(wr, "\n");
      }
    fclose(wr);
  }

void gdr_write_quantiles(FILE *wr, int32_t nr, double vals[], int32_t prec)
  {  
    int32_t nq = gdr_quantiles_NUM;
    double frac[gdr_quantiles_NUM] = { gdr_quantiles_FRACTIONS };

    gdr_demo_sort_doubles(nr, vals);
    fprintf(wr, " ");
    for (uint32_t q = 0;  q < nq; q++) 
      { double nsl_qt = gdr_demo_find_quantile(nr, vals, frac[q]);
        fprintf(wr, " %10.*f", prec, nsl_qt);
      }
  }
