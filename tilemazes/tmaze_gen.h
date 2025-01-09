#ifndef tmaze_gen_H
#define tmaze_gen_H

/* Tools for mazes of a generic given family. */
/* Last edited on 2024-12-21 11:31:47 by stolfi */

/* !!!{Change {tmaze_t} to include a {family} field, then merge this file into {tmaze.h}.} */

#include <stdint.h>
#include <limits.h>
#include <stdio.h>

#include <bool.h>
#include <frgb.h>
#include <dgraph.h>
#include <jsrandom.h>

#include <tmaze.h>

typedef enum
  { tmaze_family_BRASILIA = 0, /* The Brasilia Airport maze family. */
    tmaze_family_BLIPBLOP = 1  /* The Blip Blop maze family. */
  } tmaze_family_t;
  /* The known maze families. */
    
char *tmaze_gen_family_name_lc(tmaze_family_t fam);
  /* A string with the name of maze family {fam}, in lowercase. */

void tmaze_gen_random(tmaze_t *M, tmaze_family_t fam, int seed);    
  /* Intializes the random generator with {srandom(seed)} and then
    fills {M} with randomly chosen blip or blop tiles. */

tmaze_t tmaze_gen_from_string(char *string, tmaze_family_t fam, bool_t torus);
  /* Creates a tile maze {M} of family {fam} with given {M.torus} attribute
    whose size and tiles are specified by {string}, in a family-dependent
    encoding (one character per cell). The cells are filled row-by-row,
    from TOP TO BOTTOM and left to right. The character ';' is a 
    row terminator. */

tmaze_pattern_t tmaze_gen_pattern_from_string(char *string, tmaze_family_t fam, bool_t torus, bool_t sub, int val);
  /* Creates a tile maze pattern {M} of family {fam} with given
    attributes {M.torus,M.sub,M.val} whose size and tilesets are
    specified by {string}, in a family-dependent encoding (one
    character per cell). The cells are filled row-by-row, from TOP TO
    BOTTOM and left to right. The character ';' is a row terminator;
    '*' is any tile of the family, 'x' is no tiles, and ';' is a row
    terminator. */

/* PRINTING */

void tmaze_gen_print(FILE *wr, tmaze_t *M, tmaze_family_t fam);
  /* Prints the maze {M} to {wr}, with the encoding of family {fam}. */

/* POSTSCRIPT PLOTTING */

void tmaze_gen_plot_maze
  ( char *dir,
    char *name,
    tmaze_t *M,
    tmaze_family_t fam,
    bool_t cells,  /* TRUE paints the background of each tile. */
    bool_t grids,  /* TRUE paints the cell grid. */
    double tsize,  /* Cell size (mm) */
    bool_t curved  /* If TRUE uses curved roads when appropriate, FALSE uses polygonal ones. */
  );
  /* Writes a picture of maze {M} of family {fam} as an Encapsulated
    Postscript file, called "{dir}/{name}.eps". */

/* GRAPHS OF MAZES */

dgraph_t tmaze_gen_graph_make(tmaze_t *M, tmaze_family_t fam);
  /* Creates the undirected graph for the maze {M}, interrpeted as
    being of family {fam}.  The graph's size and topology depends 
    on {M} and {fam}. */

tmaze_size_t tmaze_gen_tot_size(tmaze_t *M, tmaze_family_t fam);
  /* Returs the total size of any maze of family {fam}
    with the dimensions and topology of {M}. */

tmaze_size_t tmaze_gen_max_comp_size(tmaze_t *M, tmaze_family_t fam);
  /* Returs the maximum possible size of any component of
    any maze of family {fam} with the dimensions and topology of {M}. */

void tmaze_gen_graph_count_components_by_size
  ( dgraph_t *G, 
    tmaze_t *M, 
    tmaze_family_t fam,
    tmaze_size_t *msP,
    tmaze_comp_count_t **ctP
  );
  /* Takes the graph {G} of maze {M}, assumed to be of family {fam}.
    Returns in {*msP} the maximum size {ms} of any
    connected component of {G}, and in {*ctP} the address of a newly
    allocated vector {ct[0..ms]} such that {ct[k]} is the number of
    connected components of {G} with exactly {k} vertices, for each
    {k}. The definition of `size' depends on the family. */

double tmaze_gen_predicted_comp_count(tmaze_t *M, tmaze_family_t fam, tmaze_size_t size);
  /* Statistically expected count of components with size {size}
    in a maze of family {fam} with {M}'s size and topology
    but with random tiles.  The tiles in {M} are ignored.
    May return NAN when there is no prediction. */
    
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
  );
  /* Gathers statistics on component sizes for tile mazes of family {fam}
    with {M}'s size and topology.  
    
    Specifically, let {ms} be the maximum possible component count for
    any tiling of {M}. The procedure generates {trials} random tile
    assignments for {M}, based on the given {seed}. For each
    assignment, counts the number {ct[s]} of components with each size
    {s} in {0..ms}. Then computes the average and deviation of
    {ct[s]}, over all trials.
    
    The largest component in every trial is tallied separately from
    the other components. Thus the procedure computes two means,
    {avg_large[s]}, referring to the largest components in each
    instance, and {avg_small[s]} referring to all components except
    the largest ones. Ditto for the deviations {dev_large[s]} and
    {dev_small[s]}.
    
    Returns in {*msP} the value of {ms}, and in {*ct_avg_smallP},
    {*ct_dev_smallP}, {*ct_avg_largeP}, {*ct_dev_largeP} the addresses
    of four newly allocated vectors with the above statistics. */
 
void tmaze_gen_print_comp_size_counts
  ( FILE *wr, 
    tmaze_t *M, 
    tmaze_family_t fam,
    tmaze_size_t ms,
    tmaze_comp_count_t ct[]
  );
  /* Given a maze {M} of family {fam}, prints to {wr} the observed
    and expected counts of connected maze components of each size.
    
    Assumes that {ms} is an upper bound to the sze of any component,
    and that the vector {ct} has {ms+1} elements each; and, for {s} in
    {0..ms}, {ct[s]} is the count of graph components that have size
    {s}.
    
    The output has one line for each component size, in the format
    "{SIZE} {CT} {TSZ} {CT_CUM} {TSZ_CUM} {ECT}" where
    
      {SIZE} is the size of a component. 
       
      {CT} is the count of components with that size.
      
      {TSZ = SIZE*CT} is the total size of those components.
      
      {CT_CUM} is the sum of {CT} up to this line.
      
      {TSZ_CUM} is the same for {TSZ}.
      
      {ECT} is the theoretical prediction for {CT}, or "nan",
        as computed by {tmaze_gen_predicted_comp_count}..
        
    The lines are sorted by increasing {SIZE}. LInes with zero {CT} and 
    NAN or zero {ECT} are omitted. Also prints #-comments
    with some global maze attributes. */

void tmaze_gen_print_comp_size_distr
  ( FILE *wr, 
    tmaze_t *M, 
    tmaze_family_t fam,
    tmaze_size_t ms,
    double ct_avg[],
    double ct_dev[]
  );
  /* Prints to {wr} average and deviation of the counts of connected
    maze components of each size, averaged over several random mazes
    of family {fam} with the size and topology of {M}.
    
    Assumes that {ms} is an upper bound to all component sizes seen,
    and that the vectors {ct_avg,ct_dev} have {ms+1} elements each;
    and, for {s} in {0..ms},
    
      {ct_avg[s]} is the count of graph components that
        have size {s}, averaged over several random tilings of {M};
      
      {ct_dev[s]} is the standard deviation of those counts.
    
    The output has one line for each component size, in the format
    "{SIZE} {ACT} {DCT} {TSZ} {ACT_CUM} {TSZ_CUM} {ECT}" where
    
      {SIZE} is the size of a component. 
       
      {ACT} is the count of components with that size, averaged over all trials.
      
      {DCT} is the observed deviation of those counts among all trials.
      
      {TSZ = SIZE*ACT} is the total size those components, averaged over all trials.
      
      {ACT_CUM} is the sum of {ACT} up to this line.
      
      {TSZ_CUM} is the same for {TSZ}.
      
      {ECT} is the theoretical prediction for {ACT}, or "nan",
        as computed by {tmaze_gen_predicted_comp_count}.
        
    The lines are sorted by increasing {SIZE}.  Also prints #-comments
    with some global maze attributes. */
 
#endif
