/* See {tmaze_gen.h} */
/* Last edited on 2024-12-21 11:31:43 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <assert.h>

#include <bool.h>
#include <affirm.h>
#include <frgb.h>
#include <dgraph.h>
#include <jsrandom.h>

#include <tmaze.h>

#include <tmaze_ba.h>
#include <tmaze_ba_graph.h>

#include <tmaze_bb.h>
#include <tmaze_bb_graph.h>

#include <tmaze_gen.h>

char *tmaze_gen_family_name_lc(tmaze_family_t fam)
  {
    if (fam == tmaze_family_BRASILIA)
      { return "brasilia"; }
    else if (fam == tmaze_family_BLIPBLOP)
      { return "blipblop"; }
    else
      { demand(FALSE, "unknown family"); }
  }

void tmaze_gen_random(tmaze_t *M, tmaze_family_t fam, int seed)
  {
    if (fam == tmaze_family_BRASILIA)
      { tmaze_ba_random(M, seed); }
    else if (fam == tmaze_family_BLIPBLOP)
      { tmaze_bb_random(M, seed); }
    else
      { demand(FALSE, "unknown family"); }
  }

tmaze_t tmaze_gen_from_string(char *string, tmaze_family_t fam, bool_t torus)
  { 
    if (fam == tmaze_family_BRASILIA)
      { return tmaze_ba_from_string(string, torus); }
    else if (fam == tmaze_family_BLIPBLOP)
      { return tmaze_bb_from_string(string, torus);}
    else
      { demand(FALSE, "unknown family"); }
  }

tmaze_pattern_t tmaze_gen_pattern_from_string(char *string, tmaze_family_t fam, bool_t torus, bool_t sub, int val)
  { 
    if (fam == tmaze_family_BRASILIA)
      { return tmaze_ba_pattern_from_string(string, torus, sub, val); }
    else if (fam == tmaze_family_BLIPBLOP)
      { return tmaze_bb_pattern_from_string(string, torus, sub, val); }
    else
      { demand(FALSE, "unknown family"); }
  }
  
void tmaze_gen_print(FILE *wr, tmaze_t *M, tmaze_family_t fam)
  {
    if (fam == tmaze_family_BRASILIA)
      { tmaze_ba_print(wr, M); }
    else if (fam == tmaze_family_BLIPBLOP)
      { tmaze_bb_print(wr, M); }
    else
      { demand(FALSE, "unknown family"); }
  }
  
void tmaze_gen_plot_maze
  ( char *dir,
    char *name,
    tmaze_t *M, 
    tmaze_family_t fam,
    bool_t cells,  /* TRUE paints the background of each tile. */
    bool_t grids,  /* TRUE paints the cell grid. */
    double tsize,  /* Cell size (mm) */
    bool_t curved  /* If TRUE uses curved roads, FALSE uses polygonal ones. */
  )
  {
    if (fam == tmaze_family_BRASILIA)
      { if (curved) 
          { fprintf(stderr, "\"curved\" style not available for Brasilia mazes\n"); }
        tmaze_ba_plot_maze(dir, name, M, cells, grids, tsize);
      }
    else if (fam == tmaze_family_BLIPBLOP)
      { tmaze_bb_plot_maze(dir, name, M, cells, grids, tsize, curved); }
    else
      { demand(FALSE, "unknown family"); }
  }

dgraph_t tmaze_gen_graph_make(tmaze_t *M, tmaze_family_t fam)
  {
    if (fam == tmaze_family_BRASILIA)
      { return tmaze_ba_graph_make(M); }
    else if (fam == tmaze_family_BLIPBLOP)
      { return tmaze_bb_graph_make(M); }
    else
      { demand(FALSE, "unknown family"); }
  }

tmaze_size_t tmaze_gen_tot_size(tmaze_t *M, tmaze_family_t fam)
  {
    if (fam == tmaze_family_BRASILIA)
      { return tmaze_ba_tot_size(M); }
    else if (fam == tmaze_family_BLIPBLOP)
      { return tmaze_bb_tot_size(M); }
    else
      { demand(FALSE, "unknown family"); }
  }

tmaze_size_t tmaze_gen_max_comp_size(tmaze_t *M, tmaze_family_t fam)
  {
    if (fam == tmaze_family_BRASILIA)
      { return tmaze_ba_max_comp_size(M); }
    else if (fam == tmaze_family_BLIPBLOP)
      { return tmaze_bb_max_comp_size(M); }
    else
      { demand(FALSE, "unknown family"); }
  }

void tmaze_gen_graph_count_components_by_size
  ( dgraph_t *G, 
    tmaze_t *M, 
    tmaze_family_t fam, 
    tmaze_size_t *msP,
    tmaze_comp_count_t **ctP
  )
  {
    if (fam == tmaze_family_BRASILIA)
      { tmaze_ba_graph_count_components_by_size(G, M, msP, ctP); }
    else if (fam == tmaze_family_BLIPBLOP)
      { tmaze_bb_graph_count_components_by_size(G, M, msP, ctP); }
    else
      { demand(FALSE, "unknown family"); }
  }

double tmaze_gen_predicted_comp_count(tmaze_t *M, tmaze_family_t fam, tmaze_size_t size)
  {
    if (fam == tmaze_family_BRASILIA)
      { return tmaze_ba_graph_predicted_comp_count(M, size); }
    else if (fam == tmaze_family_BLIPBLOP)
      { return tmaze_bb_graph_predicted_comp_count(M, size); }
    else
      { demand(FALSE, "unknown family"); }
  }

void tmaze_gen_comp_stats_gather
  ( tmaze_t *M, 
    tmaze_family_t fam, 
    int trials, 
    int seed,
    tmaze_size_t *msP, 
    double **ct_avg_smallP,
    double **ct_dev_smallP, 
    double **ct_avg_largeP,
    double **ct_dev_largeP
  )
  {
    /* We need at least two trials to estimate the deviation: */
    demand(trials >= 2, "too few trials");
    
    /* Get themax component size: */
    tmaze_size_t max_size = tmaze_gen_max_comp_size(M, fam);
    tmaze_size_t tot_size = tmaze_gen_tot_size(M, fam);

    /* Allocate tables for sums of counts and sums of counts squared: */
    double *ct_sum1[2];  /* {ct_sum1[c]} is sum of {ct[s]}  for comps of size {s} and class {c}. */
    double *ct_sum2[2];  /* {ct_sum2[c]} is sum of {ct[s]^2} for comps of size {s} and class {c}. */
    tmaze_size_t size;
    int class; /* 0 = small, 1 = large */
    for (class = 0; class <= 1; class++)
      { ct_sum1[class] = notnull(malloc((max_size+1)*sizeof(double)), "no mem");
        ct_sum2[class] = notnull(malloc((max_size+1)*sizeof(double)), "no mem");
        for (size = 0; size <= max_size; size++)
          { ct_sum1[class][size] = ct_sum2[class][size] = 0; }
      }
    
    int it;
    for (it = 0; it < trials; it++)
      {
        /* Fill the maze with random tiles: */
        tmaze_gen_random(M, fam, seed);
        
        /* Get the underlying graph {G}: */
        dgraph_t G = tmaze_gen_graph_make(M, fam);
        assert(G.rows == G.cols);
        
        /* Get the connected components as a function of size: */
        tmaze_size_t ms_loc;  /* Max possible size of a component in this instance. */
        tmaze_comp_count_t *ct_loc; /* Count cmponents by size. */
        tmaze_gen_graph_count_components_by_size (&G, M, fam, &ms_loc, &ct_loc);
        assert(ms_loc <= max_size);
        
        /* Find the maximum and total component size in this trial: */
        tmaze_size_t max_size_loc = 0;
        tmaze_size_t tot_size_loc = 0;
        for (size = 0; size <= ms_loc; size++)
          { tmaze_comp_count_t ct = ct_loc[size];
            if (ct != 0) max_size_loc = size;
            tot_size_loc += size*ct;
          }
        /* Check additivity of component sizes: */
        assert(tot_size_loc == tot_size);
        
        /* Accumulate the instance counts onto the total counts: */
        for (size = 0; size <= ms_loc; size++)
          { /* Find out the {class} (small or large) of the comps with this size: */
            class = (size == max_size_loc ? 1 : 0);
            tmaze_comp_count_t ct = ct_loc[size];
            ct_sum1[class][size] += ct;
            ct_sum2[class][size] += ct*ct;
          }
        
        /* Free instance-specific storage: */
        free(ct_loc);
        dgraph_trim(&G, 0);

        /* Prepare for the next trial: */
        seed += 17;
      }

    /* Compute avg and dev for each class and size (reuse {ct_sum1,ct_sum2}): */
    for (class = 0; class <= 1; class++)
      { for (size = 0; size <= max_size; size++)
          { double avg = ct_sum1[class][size]/trials;
            double var = fmax(0, (ct_sum2[class][size] - trials*avg*avg)/(trials - 1));
            double dev = sqrt(var);
            ct_sum1[class][size] = avg;
            ct_sum2[class][size] = dev;
          }
      }
    
    /* Return results:*/
    (*msP) = max_size;
    if (ct_avg_smallP != NULL) { (*ct_avg_smallP) = ct_sum1[0]; } else { free(ct_sum1[0]); } 
    if (ct_dev_smallP != NULL) { (*ct_dev_smallP) = ct_sum2[0]; } else { free(ct_sum2[0]); } 
    if (ct_avg_largeP != NULL) { (*ct_avg_largeP) = ct_sum1[1]; } else { free(ct_sum1[1]); } 
    if (ct_dev_largeP != NULL) { (*ct_dev_largeP) = ct_sum2[1]; } else { free(ct_sum2[1]); } 
  }

void tmaze_gen_print_comp_size_counts
  ( FILE *wr, 
    tmaze_t *M, 
    tmaze_family_t fam,
    tmaze_size_t ms,
    tmaze_comp_count_t ct[]
  )
  {
    /* Get the maze's size and the max component size: */
    int tot_size = tmaze_gen_tot_size(M, fam);
    int max_size = tmaze_gen_max_comp_size(M, fam);

    /* Print the component counts: */
    fprintf(wr, "# maze dimensions = %d × %d (%d tiles)\n", M->nx, M->ny, M->nx*M->ny);
    fprintf(wr, "# maze topology = %s\n", (M->torus ? "torus" : "open"));
    fprintf(wr, "# total size of maze = %d\n", tot_size);
    fprintf(wr, "# max possible component size = %d\n", max_size);
    fprintf(wr, "# \n");
    fprintf(wr, "# %9s %11s %11s %11s\n", "SIZE", "CT", "TSZ", "ECT");
    fprintf(wr, "# \n");
    tmaze_size_t max_size_obs = 0;
    tmaze_size_t tot_size_obs = 0;
    tmaze_size_t size;
    for (size = 0; size <= ms; size++)
      { /* Get data for this line: */
        tmaze_comp_count_t oct = ct[size];  /* Count of components with this size. */
        tmaze_size_t tsz = oct*size;        /* Total size of those components. */
        /* Get the theoretical prediction for {act}: */
        double ect = tmaze_gen_predicted_comp_count(M, fam, size);

        /* keep track of total and max component size: */
        if (oct != 0) { max_size_obs = size; }
        tot_size_obs += tsz;

        /* Print line if significant: */
        if ((oct != 0) || ((! isnan(ect)) && (ect != 0)))
          { fprintf(wr, "  %9d %11d %11d", size, oct, tsz);
            fprintf(wr, " %11.1f", ect);
            fprintf(wr, "\n"); 
          }
      }
    fprintf(wr, "# max component size seen = %d\n", max_size_obs);
    fprintf(wr, "# total component size = %d\n", tot_size_obs);
    fprintf(wr, "# fraction of maze covered = %8.6f\n", ((double)tot_size_obs)/((double)tot_size));
  }

void tmaze_gen_print_comp_size_distr
  ( FILE *wr, 
    tmaze_t *M,
    tmaze_family_t fam,
    tmaze_size_t ms,
    double ct_avg[],
    double ct_dev[]
  )
  {
    /* Get the maze's size and the max component size: */
    int tot_size = tmaze_gen_tot_size(M, fam);
    int max_size = tmaze_gen_max_comp_size(M, fam);

    /* Print some maze attributes: */
    fprintf(wr, "# maze dimensions = %d × %d (%d tiles)\n", M->nx, M->ny, M->nx*M->ny);
    fprintf(wr, "# maze topology = %s\n", (M->torus ? "torus" : "open"));
    fprintf(wr, "# total size of maze = %d\n", tot_size);
    fprintf(wr, "# max possible component size = %d\n", max_size);
    fprintf(wr, "# \n");
    fprintf(wr, "# %9s %11s %11s %11s %11s %11s %11s\n", "SIZE", "ACT", "DCT", "TSZ", "ACT_CUM", "TSZ_CUM", "ECT");
    fprintf(wr, "# \n");
 
    /* Print the table entries: */
    double act_cum = 0;  /* Cumulative count of components. */
    double tsz_cum = 0;  /* Cumulative count of vertices in those components. */
    tmaze_size_t max_size_obs = 0; /* Max component size with nonzero count. */
    tmaze_size_t size;
    for (size = 0; size <= ms; size++)
      { 
        /* Get data for this line: */
        double act = ct_avg[size];  /* Avg of count of components with this size. */
        double dct = ct_dev[size];  /* Dev of count of components with this size. */
        double tsz = act*size;      /* Total size of those components. */
        
        /* get the theoretical prediction {ect} for {act}: */
        double ect = tmaze_gen_predicted_comp_count(M, fam, size);  
        
        /* Update cumulative counts: */
        act_cum += act;
        tsz_cum += tsz;
        
        /* Keep track of max component size: */
        if (act != 0) { max_size_obs = size; }
          
        /* Print line if significant: */
        if ((act != 0) || ((! isnan(ect)) && (ect != 0)))
          { fprintf(wr, "  %9d %11.1f %11.1f %11.1f", size, act, dct, tsz);
            fprintf(wr, " %11.1f %11.1f", act_cum, tsz_cum);
            fprintf(wr, " %11.1f", ect);
            fprintf(wr, "\n"); 
          }
      }
    fprintf(wr, "# max component size seen = %d\n", max_size_obs);
    fprintf(wr, "# average covered fraction of maze = %8.6f\n", tsz_cum/((double)tot_size));
    
    fflush(wr);
  }    
